"""Setup MRF do conversor de torque CT_V1 e varredura de razões de velocidade (PyFluent, Fluent 2024 R1).

Uso (rodar com PYTHONIOENCODING=utf-8 e o .venv241):
    python tools/conversor_mrf.py setup               monta o caso, testa imagens e grava CT_V1_MRF_setup.cas.h5
    python tools/conversor_mrf.py sentido [N]         SR = 0 com os dois sentidos de rotação, N iterações cada
    python tools/conversor_mrf.py varredura S [SRs] [--ini X.dat.h5]
                                                      S = +1 ou -1; SRs opcionais (padrão SRS); --ini: partida do 1º ponto
    opções globais (qualquer comando): --malha <arquivo .msh.h5> --saida <pasta em results/conversor>
    python tools/conversor_mrf.py sensibilidade [sst_sem_cc ke_realizable fou meio_passo]
                                                      variantes do setup no estol (DIR/sens/<variante>)
    python tools/conversor_mrf.py refazer_figuras [SRs]
                                                      relê o .cas/.dat de cada ponto e refaz as imagens do Fluent
                                                      (na varredura elas saem em branco, só com a legenda)

Física (confirmada pelo Murilo em 2026-11-20; justificativas no C_TG.md e no TG):
    bomba a 2000 rpm; SR = 0; 0,2; 0,4; 0,6; 0,7; 0,8; 0,85; 0,9; óleo ρ = 830 kg/m³, μ = 0,007 Pa·s;
    k-ω SST com correção de curvatura; PRESTO!; SIMPLE (o Coupled não cabe na RAM); formulação de velocidade relativa;
    inicialização standard no referencial absoluto.

Saídas em results/conversor/CT_V1_MRF/: fluent/ (diretório de trabalho do Fluent e caso de setup),
SRxxx/{fluent,dados,graficos} (um ponto de operação por pasta), dados/ (monitores contínuos) e logs/.
"""
import json, math, re, shutil, sys, time, traceback
from pathlib import Path
import ansys.fluent.core as pf

sys.path.insert(0, str(Path(__file__).resolve().parent))
from casos import ROOT                     # noqa: E402
from run_caso import encerrar, ler_rep     # noqa: E402

DIR = ROOT / "results" / "conversor" / "CT_V1_MRF"
MALHA = ROOT / "mesh" / "CT_V1_3zonas.msh.h5"
CASO_SETUP = "CT_V1_MRF_setup"
PROCS = 6
RPM_B = 2000.0
OMEGA_B = RPM_B * 2 * math.pi / 60
SRS = [0.0, 0.2, 0.4, 0.6, 0.7, 0.8, 0.85, 0.9]
RHO, MU = 830.0, 0.007
K0, W0 = 1.0, 500.0            # k [m²/s²] e ω [1/s] iniciais: μt/μ ≈ 240
REP = time.strftime("CT_V1_MRF_rep_%Y%m%d-%H%M%S.out")  # monitores desta sessão (todas as SR), em DIR/dados
ITER_FOU = 200                  # iterações iniciais em primeira ordem (só no 1º ponto)
BLOCO = 100                     # iterações entre checagens de convergência
JANELA = 200                    # janela da estabilidade dos torques
TOL = 0.005                     # variação relativa máxima na janela (0,5 %)
ITER_MIN = {"primeiro": 1000, "demais": 600}
ITER_MAX = {"primeiro": 3000, "demais": 2000}   # ~4,5 s/iteração com SIMPLE: 3,75 h e 2,5 h no pior caso

ZONAS = ("fluido_bomba", "fluido_turbina", "fluido_estator")
PECA = {"fluido_bomba": "bomba", "fluido_turbina": "turbina", "fluido_estator": "estator"}
PARES = {  # interface: (lado da 1ª zona, lado da 2ª zona)
    "if_bt": ("if_bomba_turbina_lado_bomba", "if_bomba_turbina_lado_turbina"),
    "if_te": ("if_turbina_estator_lado_turbina", "if_turbina_estator_lado_estator"),
    "if_eb": ("if_estator_bomba_lado_estator", "if_estator_bomba_lado_bomba"),
}
# Paredes que pertencem a uma peça mas tocam o fluido de outra zona (gira com a peça dona, referencial absoluto).
PAREDES_ABS = {"parede_bomba_lado_turbina": "bomba", "parede_bomba_lado_estator": "bomba",
               "parede_estator_lado_turbina": "estator", "parede_estator_lado_bomba": "estator"}


def log(*a):
    """Imprime com horário e esvazia o buffer.

    Args:
        *a: Itens a imprimir.
    """
    print(time.strftime("[%H:%M:%S]"), *a, flush=True)


def pasta_sr(sr):
    """Devolve (e cria) a pasta de um ponto de operação.

    Args:
        sr (float): Razão de velocidade ω_t/ω_b.

    Returns:
        Path: ``DIR/SRxxx`` com as subpastas ``fluent``, ``dados`` e ``graficos``.
    """
    p = DIR / f"SR{round(sr * 100):03d}"
    for sub in ("fluent", "dados", "graficos"):
        (p / sub).mkdir(parents=True, exist_ok=True)
    return p


def abrir():
    """Abre o Fluent 3D, dupla precisão, ``PROCS`` processos, sem GUI, com diretório de trabalho ``DIR/fluent``.

    Returns:
        ansys.fluent.core.session_solver.Solver: Sessão conectada.
    """
    for sub in ("fluent", "dados", "graficos", "logs"):
        (DIR / sub).mkdir(parents=True, exist_ok=True)
    return pf.launch_fluent(product_version=pf.FluentVersion.v241, dimension=pf.Dimension.THREE,
                            precision=pf.Precision.DOUBLE, processor_count=PROCS, ui_mode=pf.UIMode.NO_GUI,
                            mode="solver", cwd=str(DIR / "fluent"), cleanup_on_exit=False)


def adjacencia(st, nomes):
    """Descobre a zona de células adjacente a cada parede (renomeia por adjacência e desfaz).

    Args:
        st: Raiz ``session.settings``.
        nomes (list[str]): Paredes a consultar.

    Returns:
        dict[str, str]: Parede -> zona de células adjacente.
    """
    bc = st.setup.boundary_conditions
    out = {}
    for n in nomes:                       # uma por vez: duas "wall-N" juntas dariam nomes ambíguos
        antes = set(bc.wall.keys())
        bc.rename_by_adjacency(zone_name=[n], abbreviate_types=False, exclude=False)
        m = (set(bc.wall.keys()) - antes) or {n}
        m = m.pop()
        out[n] = next(z for z in ZONAS if z in m)
        if m != n:
            bc.set_zone_name(zonename=m, newname=n)
    return out


def criar_interfaces(st):
    """Converte os pares ``if_*`` em interface e cria as três mesh interfaces não conformes.

    O pareamento one-to-one fica ativo; as faces sem correspondente viram paredes, renomeadas para
    ``sobra_<if>_lado_<peça>``. Cada sobra pertence à peça do *outro* lado da interface (a casca sólida
    que o fluido toca ali) e gira com ela.

    Args:
        st: Raiz ``session.settings``.

    Returns:
        dict: ``sobras`` (parede -> peça dona) e ``interfaces`` (estado lido da sessão).
    """
    bc = st.setup.boundary_conditions
    bc.set_zone_type(zone_list=list(bc.symmetry.keys()), new_type="interface")
    it = st.setup.mesh_interfaces.interface
    sobras = {}
    for nome, (a, b) in PARES.items():
        antes = set(bc.wall.keys())
        it.create(nome)
        it[nome].zone_names = [a, b]     # zone1/zone2 são só leitura na 24.1; o TUI create entra em laço de prompt
        novas = sorted(set(bc.wall.keys()) - antes)
        adj = adjacencia(st, novas) if novas else {}
        for w, zona in adj.items():
            lado = PECA[zona]
            dono = PECA[[z for z in ZONAS if PECA[z] in nome_lado(nome) and z != zona][0]]
            novo = f"sobra_{nome}_lado_{lado}"
            bc.set_zone_name(zonename=w, newname=novo)
            sobras[novo] = dono
    return {"sobras": sobras, "interfaces": {k: it[k]() for k in it.keys()}}


def nome_lado(nome):
    """Peças ligadas por uma interface.

    Args:
        nome (str): ``if_bt``, ``if_te`` ou ``if_eb``.

    Returns:
        str: Nomes das duas peças concatenados (para busca por substring).
    """
    return {"if_bt": "bomba turbina", "if_te": "turbina estator", "if_eb": "estator bomba"}[nome]


def girar_parede(st, zona, omega, sinal):
    """Parede móvel rotacional em referencial absoluto, em torno do eixo x pela origem.

    Args:
        st: Raiz ``session.settings``.
        zona (str): Parede.
        omega (float): Velocidade angular absoluta [rad/s] no sentido de ``sinal``·x.
        sinal (int): +1 ou −1, sentido do eixo de rotação.
    """
    m = st.setup.boundary_conditions.wall[zona].momentum
    m.motion_bc = "Moving Wall"
    m.relative = False
    m.rotating = True
    m.rotation_axis_origin = [0.0, 0.0, 0.0]
    m.rotation_axis_direction = [float(sinal), 0.0, 0.0]
    try:
        m.omega = omega
    except Exception:
        m.omega = {"option": "value", "value": omega}


def frame(st, zona, omega, sinal):
    """Liga o Frame Motion de uma zona de células em torno de ``sinal``·x pela origem.

    Args:
        st: Raiz ``session.settings``.
        zona (str): Zona de células.
        omega (float): Velocidade angular [rad/s].
        sinal (int): +1 ou −1.
    """
    rf = st.setup.cell_zone_conditions.fluid[zona].reference_frame
    rf.frame_motion = True
    rf.mrf_omega = {"option": "value", "value": omega}
    rf.reference_frame_axis_origin = [0.0, 0.0, 0.0]
    rf.reference_frame_axis_direction = [float(sinal), 0.0, 0.0]


def aplicar_rotacao(st, sr, sinal, sobras):
    """Aplica ω_b à bomba e SR·ω_b à turbina (zonas e paredes absolutas); o estator fica parado.

    Args:
        st: Raiz ``session.settings``.
        sr (float): Razão de velocidade.
        sinal (int): +1 ou −1.
        sobras (dict[str, str]): Paredes de sobra das interfaces e a peça dona de cada uma.

    Returns:
        dict: Estado de Frame Motion e das paredes absolutas lido da sessão.
    """
    w = {"bomba": OMEGA_B, "turbina": sr * OMEGA_B, "estator": 0.0}
    frame(st, "fluido_bomba", w["bomba"], sinal)
    frame(st, "fluido_turbina", w["turbina"], sinal)
    for zona, dono in {**PAREDES_ABS, **sobras}.items():
        girar_parede(st, zona, w[dono], sinal)
    cz = st.setup.cell_zone_conditions.fluid
    bc = st.setup.boundary_conditions.wall
    return {"frames": {z: cz[z].reference_frame() for z in ZONAS},
            "paredes_abs": {z: bc[z].momentum() for z in {**PAREDES_ABS, **sobras}}}


def setup_do_caso(st):
    """Reconstrói ``sobras`` e ``grupos_torque`` a partir das paredes de um caso já montado.

    Serve quando ``dados/setup.json`` não existe (a varredura grava o arquivo só no fim). Cada
    ``sobra_<if>_lado_<peça>`` pertence à outra peça da interface, como em ``criar_interfaces``.

    Args:
        st: Raiz ``session.settings`` com o caso lido.

    Returns:
        dict: ``sobras`` (parede -> peça dona) e ``grupos_torque`` (peça -> paredes).

    Raises:
        AssertionError: Se alguma parede ficar sem peça (de ``paredes_por_peca``).
    """
    sobras = {}
    for w in st.setup.boundary_conditions.wall.keys():
        m = re.fullmatch(r"sobra_(if_[a-z]{2})_lado_([a-z]+)", w)
        if m:
            nome, lado = m.groups()
            sobras[w] = next(x for x in nome_lado(nome).split() if x != lado)
    return {"sobras": sobras, "grupos_torque": paredes_por_peca(st, sobras)}


def paredes_por_peca(st, sobras):
    """Agrupa as paredes por peça para os torques.

    Args:
        st: Raiz ``session.settings``.
        sobras (dict[str, str]): Paredes de sobra e a peça dona.

    Returns:
        dict[str, list[str]]: Peça -> paredes.

    Raises:
        AssertionError: Se alguma parede ficar sem peça.
    """
    g = {"bomba": ["parede_bomba"], "turbina": ["parede_turbina"], "estator": ["parede_estator"]}
    for z, dono in {**PAREDES_ABS, **sobras}.items():
        g[dono].append(z)
    todas = set(st.setup.boundary_conditions.wall.keys())
    assert todas == {z for v in g.values() for z in v}, f"paredes sem peça: {todas - {z for v in g.values() for z in v}}"
    return g


def expressao_lx(s, st):
    """Cria o campo do momento angular específico em torno de x, r·c_θ = y·w − z·v (velocidade absoluta).

    Tenta named expression; se o Fluent não aceitar, usa custom field function.

    Args:
        s (ansys.fluent.core.session_solver.Solver): Sessão.
        st: Raiz ``session.settings``.

    Returns:
        str | None: Nome do campo criado, ou None se nenhuma forma funcionou.
    """
    ne = st.setup.named_expressions
    for defin in ("y*Velocity.z - z*Velocity.y", "y*VelocityZ - z*VelocityY"):
        try:
            if "lx" not in ne.keys():
                ne.create("lx")
            ne["lx"].definition = defin
            ne["lx"].output_parameter = False
            st.solution.report_definitions.surface.create("teste_lx")
            st.solution.report_definitions.surface["teste_lx"] = {"report_type": "surface-flowrate", "field": "expr:lx",
                                                                  "surface_names": ["if_bomba_turbina_lado_bomba"]}
            # compute() fica inativo antes da inicialização: basta o campo ser aceito pelo report
            st.solution.report_definitions.surface.delete(name_list=["teste_lx"])
            log("lx =", defin)
            return "expr:lx"           # named expressions entram como campo com o prefixo "expr:"
        except Exception as e:
            log("lx com", repr(defin), "falhou:", repr(e)[:300])
            try:
                st.solution.report_definitions.surface.delete(name_list=["teste_lx"])
            except Exception:
                pass
    return None


def momento_em_nm(obj):
    """Faz um report ``moment`` sair em N·m (sem isso o Fluent entrega o coeficiente Cm).

    A opção é descoberta em ``report_output_type.allowed_values()``: a que não contém "Coefficient".
    Na 2024 R1 o report de momento oferece ``['Drag Coefficient', 'Drag Force']`` (rótulos herdados do
    report de força); ``Drag Force`` dá o momento em N·m, exatamente 0,6125 × o coeficiente com os
    valores de referência padrão (sonda no SR000 em 2026-11-21).

    Args:
        obj: Objeto ``report_definitions.moment[nome]``.

    Returns:
        str: Opção aplicada.

    Raises:
        RuntimeError: Se o Fluent não oferecer uma opção de momento dimensional.
    """
    opcoes = list(obj.report_output_type.allowed_values())
    escolha = next((o for o in opcoes if "coef" not in o.lower()), None)
    if escolha is None:
        raise RuntimeError(f"sem opção de momento dimensional em {opcoes}")
    obj.report_output_type = escolha
    return escolha


def corrigir_momentos(st):
    """Aplica ``momento_em_nm`` a todos os reports ``moment`` (casos gravados antes da correção).

    Args:
        st: Raiz ``session.settings``.

    Returns:
        dict[str, str]: Report -> opção aplicada.
    """
    rd = st.solution.report_definitions.moment
    return {n: momento_em_nm(rd[n]) for n in rd.keys()}


def reports(s, st, grupos, sinal):
    """Cria os torques por peça, o total, as vazões e os fluxos de momento angular nas interfaces.

    Args:
        s (ansys.fluent.core.session_solver.Solver): Sessão.
        st: Raiz ``session.settings``.
        grupos (dict[str, list[str]]): Peça -> paredes.
        sinal (int): Sentido do eixo; os torques são projetados em ``sinal``·x.

    Returns:
        list[str]: Nomes dos reports, na ordem do Report File.
    """
    rd = st.solution.report_definitions
    nomes = []
    for p, zs in [*grupos.items(), ("total", [z for v in grupos.values() for z in v])]:
        n = f"t_{p}"
        rd.moment.create(n)
        rd.moment[n] = {"mom_center": [0.0, 0.0, 0.0], "mom_axis": [float(sinal), 0.0, 0.0], "zones": zs}
        momento_em_nm(rd.moment[n])
        nomes.append(n)
    log("MOMENT exemplo:", rd.moment["t_bomba"]())
    lx = expressao_lx(s, st)
    for lados in PARES.values():
        for z in lados:
            n = "m_" + z.replace("if_", "").replace("_lado_", "_")
            rd.surface.create(n)
            rd.surface[n] = {"report_type": "surface-massflowrate", "surface_names": [z]}
            nomes.append(n)
            if lx:
                n2 = "h_" + n[2:]
                rd.surface.create(n2)
                rd.surface[n2] = {"report_type": "surface-flowrate", "field": lx, "surface_names": [z]}
                nomes.append(n2)
    rf = st.solution.monitor.report_files
    rf.create("mon")
    rf["mon"] = {"file_name": f"../dados/{REP}", "report_defs": nomes}
    rp = st.solution.monitor.report_plots
    rp.create("torques")
    rp["torques"] = {"report_defs": [n for n in nomes if n.startswith("t_")]}
    return nomes


def montar(s, sinal):
    """Monta o setup MRF completo e inicializa.

    Args:
        s (ansys.fluent.core.session_solver.Solver): Sessão recém-aberta.
        sinal (int): +1 ou −1, sentido do eixo de rotação.

    Returns:
        dict: Setup efetivo lido da sessão.
    """
    global PAREDES_ABS
    st = s.settings
    st.file.read_mesh(file_name=str(MALHA))
    st.mesh.check()
    ef = {"malha": MALHA.name, "sinal": sinal, "rpm_b": RPM_B, "omega_b": OMEGA_B}
    existentes = set(st.setup.boundary_conditions.wall.keys())
    ausentes = [z for z in PAREDES_ABS if z not in existentes]
    if ausentes:                    # malhas grossas do GCI podem perder zonas pequenas (62 faces na G0)
        log("ATENÇÃO: paredes ausentes nesta malha:", ausentes)
        PAREDES_ABS = {z: d for z, d in PAREDES_ABS.items() if z in existentes}
    ef["paredes_abs_ausentes"] = ausentes
    ef.update(criar_interfaces(st))
    adj = adjacencia(st, list(PAREDES_ABS))
    ef["adjacencia_paredes_lado"] = adj
    for z, viz in adj.items():
        esperado = "fluido_" + z.split("_lado_")[1]
        if viz != esperado:
            log(f"ATENÇÃO: {z} é adjacente a {viz}, não a {esperado}")

    g = st.setup.general
    g.solver.velocity_formulation = "relative"
    g.operating_conditions.reference_pressure_location = [-0.03, 0.30, 0.0]
    v = st.setup.models.viscous
    v.model = "k-omega"
    v.k_omega_model = "sst"
    v.options.curvature_correction = True
    st.setup.models.energy.enabled = False

    fl = st.setup.materials.fluid
    fl.make_a_copy(from_="air", to="oleo-atf")
    fl["oleo-atf"].density = {"option": "constant", "value": RHO}
    fl["oleo-atf"].viscosity = {"option": "constant", "value": MU}
    for z in ZONAS:
        st.setup.cell_zone_conditions.fluid[z].material = "oleo-atf"

    ef["rotacao"] = aplicar_rotacao(st, 0.0, sinal, ef["sobras"])
    grupos = paredes_por_peca(st, ef["sobras"])
    ef["grupos_torque"] = grupos
    ef["reports"] = reports(s, st, grupos, sinal)

    m = st.solution.methods
    # Coupled precisa de ~17,6 GB nesta malha (16 GB de RAM: swap, 15 s/iteração). TG 4.4: o Coupled vale
    # "quando a memória disponível o permite"; a família SIMPLE é a alternativa mais leve. URFs padrão (TG 2.8).
    m.p_v_coupling.flow_scheme = "SIMPLE"
    m.gradient_scheme = "least-square-cell-based"
    m.discretization_scheme["pressure"] = "presto!"
    for k in m.discretization_scheme.keys():
        if k != "pressure":
            m.discretization_scheme[k] = "second-order-upwind"
    res = st.solution.monitor.residual
    for eq in res.equations.keys():
        res.equations[eq] = {"monitor": True, "check_convergence": False}  # parada pelos torques (TG 2.8, 4.4)

    ini = st.solution.initialization
    ini.initialization_type = "standard"
    ini.reference_frame = "absolute"            # TG 2.8: recomendação do User's Guide para MRF
    ini.defaults = {"pressure": 0, "x-velocity": 0, "y-velocity": 0, "z-velocity": 0, "k": K0, "omega": W0}
    ini.initialize()
    ef.update({"general": g(), "viscous": v(), "methods": m(), "init": ini(),
               "urf": st.solution.controls.under_relaxation(),
               "rho": RHO, "mu": MU, "pyfluent": pf.__version__})
    return ef


# ----------------------------------------------------------------------------------------------- imagens
def camera(st, direita, cima):
    """Posiciona a câmera a partir da vista frontal do Fluent, orbitando em torno do alvo.

    Na 2024 R1 sem GUI, ``camera.position``/``target`` fora do eixo z e ``views.auto_scale`` deixam a
    imagem só com a legenda (sondas no SR000 em 2026-11-21); ``restore_view("front")`` enquadra a
    geometria e ``camera.orbit`` gira a vista sem perder o enquadramento.

    Args:
        st: Raiz ``session.settings``.
        direita (float): Órbita para a direita, em graus.
        cima (float): Órbita para cima, em graus.
    """
    views = st.results.graphics.views
    views.restore_view(view_name="front")      # câmera em +z olhando o plano meridional (z = 0)
    if direita or cima:
        views.camera.orbit(right=direita, up=cima)


VISTAS = {  # nome: (órbita para a direita, órbita para cima) em graus, a partir da vista frontal (+z)
    "meridional": (0.0, 0.0),
    "de_menos_x": (-90.0, 0.0),
    "de_mais_x": (90.0, 0.0),
    "iso": (35.0, 25.0),
    "de_mais_y": (0.0, 90.0),
}


def preparar_picture(st):
    """Formato PNG, 1600 × 1200, fundo branco.

    Args:
        st: Raiz ``session.settings``.

    Returns:
        dict: Estado de ``picture``.
    """
    pic = st.results.graphics.picture
    for nome, valor in (("driver_options.hardcopy_format", "png"), ("use_window_resolution", False),
                        ("x_resolution", 1600), ("y_resolution", 1200), ("invert_background", True)):
        obj = pic
        try:
            for parte in nome.split(".")[:-1]:
                obj = getattr(obj, parte)
            setattr(obj, nome.split(".")[-1], valor)
        except Exception as e:     # sem GUI alguns itens ficam inativos
            log(f"picture.{nome}:", repr(e)[:150])
    try:
        return pic()
    except Exception:
        return {}


def salvar(st, arquivo):
    """Grava a janela gráfica atual em PNG.

    Args:
        st: Raiz ``session.settings``.
        arquivo (Path): Destino.

    Returns:
        bool: True se o arquivo foi criado.
    """
    arquivo.unlink(missing_ok=True)
    st.results.graphics.picture.save_picture(file_name=str(arquivo))
    return arquivo.exists() and arquivo.stat().st_size > 0


def contorno(st, nome, campo, superficies, vista, arquivo, faixa=None):
    """Cria (ou atualiza) um contorno, exibe na vista pedida e salva.

    Args:
        st: Raiz ``session.settings``.
        nome (str): Nome do objeto gráfico.
        campo (str): Campo do Fluent.
        superficies (list[str]): Superfícies.
        vista (str): Chave de ``VISTAS``.
        arquivo (Path): PNG de saída.
        faixa (tuple[float, float], optional): Faixa fixa; senão automática.

    Returns:
        bool: True se a imagem foi gravada.
    """
    c = st.results.graphics.contour
    if nome not in c.keys():
        c.create(nome)
    c[nome].field = campo
    c[nome].surfaces_list = superficies
    c[nome].filled = True
    c[nome].node_values = True
    if faixa:
        c[nome].range_option = {"option": "auto-range-off", "auto_range_off": {"minimum": faixa[0], "maximum": faixa[1]}}
    else:
        c[nome].range_option = {"option": "auto-range-on", "auto_range_on": {"global_range": False}}
    c[nome].display()
    camera(st, *VISTAS[vista])
    return salvar(st, arquivo)


def vetores(st, nome, superficies, vista, arquivo, campo_vetor="relative-velocity", escala=1.0, pulo=2):
    """Vetores de velocidade (relativa por padrão) coloridos pelo módulo.

    Args:
        st: Raiz ``session.settings``.
        nome (str): Nome do objeto gráfico.
        superficies (list[str]): Superfícies.
        vista (str): Chave de ``VISTAS``.
        arquivo (Path): PNG de saída.
        campo_vetor (str, optional): Campo vetorial. Padrão ``"relative-velocity"``.
        escala (float, optional): Escala das setas.
        pulo (int, optional): Pula vetores para não poluir.

    Returns:
        bool: True se a imagem foi gravada.
    """
    v = st.results.graphics.vector
    if nome not in v.keys():
        v.create(nome)
    v[nome].vector_field = campo_vetor
    v[nome].field = "rel-velocity-magnitude" if campo_vetor.startswith("rel") else "velocity-magnitude"
    v[nome].surfaces_list = superficies
    v[nome].scale = {"auto_scale": True, "scale_f": escala}
    v[nome].skip = pulo
    v[nome].display()
    camera(st, *VISTAS[vista])
    return salvar(st, arquivo)


def superficies_pos(st):
    """Cria as superfícies de pós-processamento: plano meridional, cortes axiais e raios constantes.

    Args:
        st: Raiz ``session.settings``.

    Returns:
        dict[str, list[str]]: Grupos de superfícies.
    """
    sf = st.results.surfaces
    out = {}
    if "plano_meridional" not in sf.plane_surface.keys():
        sf.plane_surface.create("plano_meridional")
        sf.plane_surface["plano_meridional"] = {"method": "xy-plane", "z": 0.0}
    out["meridional"] = ["plano_meridional"]
    for n, x in (("corte_x_bomba", -0.09), ("corte_x_turbina", 0.08), ("corte_x_estator", 0.0)):
        if n not in sf.plane_surface.keys():
            sf.plane_surface.create(n)
            sf.plane_surface[n] = {"method": "yz-plane", "x": x}
    out["cortes"] = ["corte_x_bomba", "corte_x_turbina", "corte_x_estator"]
    for n, r, zs in (("raio_bomba", 0.29, ["fluido_bomba"]), ("raio_turbina", 0.29, ["fluido_turbina"]),
                     ("raio_estator", 0.17, ["fluido_estator"])):
        if n not in sf.iso_surface.keys():
            sf.iso_surface.create(n)
            sf.iso_surface[n] = {"field": "radial-coordinate", "zones": zs, "iso_values": [r]}
    out["raios"] = ["raio_bomba", "raio_turbina", "raio_estator"]
    return out


def teste_imagem(st, destino):
    """Testa se ``save_picture`` funciona sem GUI (malha das paredes da bomba).

    Args:
        st: Raiz ``session.settings``.
        destino (Path): Pasta do PNG de teste.

    Returns:
        bool: True se gravou uma imagem não vazia.
    """
    try:
        preparar_picture(st)
        gm = st.results.graphics.mesh
        if "teste_malha" not in gm.keys():
            gm.create("teste_malha")
        gm["teste_malha"].surfaces_list = ["parede_bomba", "parede_turbina", "parede_estator"]
        gm["teste_malha"].display()
        camera(st, *VISTAS["iso"])
        ok = salvar(st, destino / "00_teste_save_picture.png")
    except Exception as e:
        log("save_picture falhou:", repr(e))
        ok = False
    log("TESTE save_picture sem GUI:", ok)
    return ok


# ----------------------------------------------------------------------------------------------- solução
def convergiu(linhas, n0, colunas, n_min):
    """Testa a estabilidade dos torques pela média móvel.

    Os torques oscilam ~1–2 % de iteração para iteração (SIMPLE, MRF, interfaces não conformes) sem
    deriva; o critério é a média das últimas ``JANELA`` iterações mudar menos que ``TOL`` em relação à
    média das ``JANELA`` anteriores (T_e relativo a |T_b|, pois cruza zero perto do acoplamento).

    Args:
        linhas (list[list[float]]): Linhas numéricas do Report File (iteração, valores...).
        n0 (int): Iteração em que o ponto começou.
        colunas (dict[str, int]): Report -> índice da coluna.
        n_min (int): Iterações mínimas no ponto.

    Returns:
        tuple[bool, dict]: (convergiu, deriva relativa da média por torque).
    """
    d = [l for l in linhas if l[0] > n0]
    if len(d) < max(n_min, 2 * JANELA):
        return False, {}
    w1, w0 = d[-JANELA:], d[-2 * JANELA:-JANELA]
    med = lambda w, n: sum(l[colunas[n]] for l in w) / len(w)
    ref = max(abs(med(w1, "t_bomba")), 1e-12)
    var = {}
    for n in ("t_bomba", "t_turbina", "t_estator"):
        m1, m0 = med(w1, n), med(w0, n)
        var[n] = abs(m1 - m0) / (max(abs(m1), 1e-12) if n != "t_estator" else ref)
    return all(x < TOL for x in var.values()), var


def media_janela(arq, n0):
    """Média e desvio-padrão de todos os reports nas últimas ``JANELA`` iterações do ponto.

    Args:
        arq (Path): Report File.
        n0 (int): Iteração em que o ponto começou.

    Returns:
        tuple[dict[str, float], dict[str, float]]: (médias, desvios-padrão) por report.
    """
    nomes, linhas = ler_tabela(arq)
    w = [l for l in linhas if l[0] > n0][-JANELA:]
    med, dp = {}, {}
    for i, n in enumerate(nomes[1:], 1):
        v = [l[i] for l in w]
        m = sum(v) / len(v)
        med[n] = m
        dp[n] = (sum((x - m) ** 2 for x in v) / max(len(v) - 1, 1)) ** 0.5
    return med, dp


def ler_tabela(arquivo):
    """Lê o Report File com todas as colunas.

    Args:
        arquivo (Path): Arquivo ``.out``.

    Returns:
        tuple[list[str], list[list[float]]]: Nomes das colunas e linhas numéricas.
    """
    nomes, linhas = [], []
    if not arquivo.exists():
        return nomes, linhas
    for ln in arquivo.read_text(errors="replace").splitlines():
        if '"Iteration"' in ln:
            nomes = re.findall(r'"([^"]+)"', ln)
        p = ln.split()
        if p and p[0].isdigit():
            try:
                linhas.append([float(x) for x in p])
            except ValueError:
                pass
    return nomes, linhas


def iterar_ponto(s, sr, primeiro, rotulo, n_min=None):
    """Itera um ponto de operação até a estabilidade dos torques ou o teto.

    Args:
        s (ansys.fluent.core.session_solver.Solver): Sessão com o ponto aplicado.
        sr (float): Razão de velocidade (só para o log).
        primeiro (bool): Primeiro ponto da varredura (começa em primeira ordem e tem teto maior).
        rotulo (str): Texto do log.
        n_min (int, optional): Iterações mínimas; padrão ``ITER_MIN`` do tipo de ponto.

    Returns:
        dict: Iteração inicial e final, convergência, variações, tempo de parede.

    Raises:
        RuntimeError: Se um bloco de iterações não acrescentar linhas ao report file.
    """
    st = s.settings
    arq = DIR / "dados" / REP
    nomes, linhas = ler_tabela(arq)
    n0 = int(linhas[-1][0]) if linhas else 0
    chave = "primeiro" if primeiro else "demais"
    t0 = time.time()
    m = st.solution.methods
    if primeiro:
        for k in m.discretization_scheme.keys():
            if k != "pressure":
                m.discretization_scheme[k] = "first-order-upwind"
        log(f"{rotulo}: {ITER_FOU} iterações em primeira ordem")
        st.solution.run_calculation.iterate(iter_count=ITER_FOU)
        for k in m.discretization_scheme.keys():
            if k != "pressure":
                m.discretization_scheme[k] = "second-order-upwind"
    ok, var, ultimo = False, {}, -1
    while True:
        st.solution.run_calculation.iterate(iter_count=BLOCO)
        nomes, linhas = ler_tabela(arq)
        col = {n: i for i, n in enumerate(nomes)}
        feito = int(linhas[-1][0]) - n0
        if feito <= ultimo:
            raise RuntimeError(f"{rotulo}: o report file {arq.name} não ganhou linhas em {BLOCO} iterações "
                               "(o Fluent pode estar gravando em outro arquivo)")
        ultimo = feito
        ok, var = convergiu(linhas, n0, col, n_min if n_min is not None else ITER_MIN[chave])
        u = linhas[-1]
        log(f"{rotulo} it+{feito}: T_b={u[col['t_bomba']]:.4g} T_t={u[col['t_turbina']]:.4g} "
            f"T_e={u[col['t_estator']]:.4g} soma={u[col['t_total']]:.3g} var={ {k: f'{x:.2e}' for k, x in var.items()} }")
        if ok or feito >= ITER_MAX[chave]:
            break
    med, dp = media_janela(arq, n0)
    return {"sr": sr, "it_ini": n0, "it_fim": int(linhas[-1][0]), "convergiu": ok, "deriva_media": var,
            "media_janela": med, "desvio_janela": dp, "janela": JANELA, "tempo_s": round(time.time() - t0, 1)}


def finais(st, nomes):
    """Calcula o valor atual de todos os reports.

    Args:
        st: Raiz ``session.settings``.
        nomes (list[str]): Reports.

    Returns:
        dict[str, float]: Report -> valor.
    """
    out = {}
    for d in st.solution.report_definitions.compute(report_defs=nomes):
        out.update({k: v[0] for k, v in d.items()})
    return out


def aceitos(cmd, campos):
    """Filtra a lista de campos pelos valores aceitos pelo Fluent (quando ele informa).

    Args:
        cmd: Objeto de settings com ``allowed_values()``.
        campos (list[str]): Campos desejados.

    Returns:
        list[str]: Campos aceitos (todos, se o Fluent não informar).
    """
    try:
        ok = set(cmd.allowed_values())
    except Exception:
        return campos
    fora = [c for c in campos if c not in ok]
    if fora:
        log("campos não aceitos:", fora)
    return [c for c in campos if c in ok]


def torques_por_zona(st, sobras, sinal):
    """Torque do fluido nas paredes adjacentes a cada zona de fluido (para comparar com Euler na zona).

    Args:
        st: Raiz ``session.settings``.
        sobras (dict[str, str]): Paredes de sobra (o sufixo ``_lado_<peça>`` diz a zona adjacente).
        sinal (int): Sentido do eixo.

    Returns:
        dict[str, float]: Zona -> torque [N·m] no eixo ``sinal``·x.
    """
    rd = st.solution.report_definitions.moment
    grupos = {z: [f"parede_{PECA[z]}"] for z in ZONAS}
    for w in [*PAREDES_ABS, *sobras]:
        grupos["fluido_" + w.split("_lado_")[1]].append(w)
    nomes = []
    for z, ws in grupos.items():
        n = "z_" + z
        if n not in rd.keys():
            rd.create(n)
        rd[n] = {"mom_center": [0.0, 0.0, 0.0], "mom_axis": [float(sinal), 0.0, 0.0], "zones": ws}
        momento_em_nm(rd[n])
        nomes.append(n)
    f = finais(st, nomes)
    return {z: f["z_" + z] for z in ZONAS}


def exportar(st, p, sufixo):
    """Exporta interfaces e paredes (ASCII; valores da célula adjacente a cada face) para a análise em Python.

    Args:
        st: Raiz ``session.settings``.
        p (Path): Pasta do ponto.
        sufixo (str): Rótulo do ponto (``SRxxx``).
    """
    campos_if = ["pressure", "x-velocity", "y-velocity", "z-velocity", "axial-velocity", "radial-velocity",
                 "tangential-velocity", "rel-tangential-velocity", "density", "x-face-area", "y-face-area",
                 "z-face-area", "face-area-magnitude", "radial-coordinate", "expr:lx"]
    campos_if = aceitos(st.file.export.ascii.cell_func_domain, campos_if)
    campos_par = aceitos(st.file.export.ascii.cell_func_domain, ["pressure", "wall-shear", "y-plus", "radial-coordinate"])
    ifs = [z for v in PARES.values() for z in v]
    for z in ifs:
        try:
            st.file.export.ascii(file_name=str(p / "dados" / f"{sufixo}_{z}.csv"), surface_name_list=[z],
                                 delimiter="comma", location="cell-center", cell_func_domain=campos_if)
        except Exception as e:
            log("export interface", z, e)
    for peca in ("bomba", "turbina", "estator"):
        try:
            st.file.export.ascii(file_name=str(p / "dados" / f"{sufixo}_parede_{peca}.csv"),
                                 surface_name_list=[f"parede_{peca}"], delimiter="comma", location="cell-center",
                                 cell_func_domain=campos_par)
        except Exception as e:
            log("export parede", peca, e)


def figuras(st, p, sufixo, grupos):
    """Gera o conjunto de imagens do Fluent de um ponto; cada figura falha sozinha, sem parar o resto.

    Args:
        st: Raiz ``session.settings``.
        p (Path): Pasta do ponto.
        sufixo (str): Rótulo do ponto.
        grupos (dict[str, list[str]]): Peça -> paredes.

    Returns:
        dict[str, bool]: Figura -> gravada.
    """
    g = p / "graficos"
    sup = superficies_pos(st)
    res = {}

    def tenta(nome, f, *a, **k):
        """Executa uma figura e registra o resultado.

        Args:
            nome (str): Nome do arquivo (sem extensão).
            f (callable): Função de figura.
            *a: Argumentos posicionais.
            **k: Argumentos nomeados.
        """
        try:
            res[nome] = f(*a, **k)
        except Exception as e:
            res[nome] = False
            log("figura", nome, "falhou:", repr(e)[:300])

    preparar_picture(st)
    mer = sup["meridional"]
    for i, (campo, rot) in enumerate([("pressure", "pressao"), ("total-pressure", "pressao_total"),
                                      ("velocity-magnitude", "vabs"), ("rel-velocity-magnitude", "vrel"),
                                      ("axial-velocity", "v_axial"), ("radial-velocity", "v_radial"),
                                      ("tangential-velocity", "v_tangencial"), ("turb-kinetic-energy", "k"),
                                      ("viscosity-ratio", "mut_mu")], 1):
        n = f"{i:02d}_meridional_{rot}"
        tenta(n, contorno, st, "c_mer", campo, mer, "meridional", g / f"{sufixo}_{n}.png")
    tenta("10_meridional_vetores_vrel", vetores, st, "v_mer", mer, "meridional", g / f"{sufixo}_10_meridional_vetores_vrel.png")
    for j, s_ in enumerate(sup["raios"]):
        tenta(f"2{j}_raio_{s_}_pressao", contorno, st, "c_raio", "pressure", [s_], "de_mais_y",
              g / f"{sufixo}_2{j}_{s_}_pressao.png")
        tenta(f"2{j}_raio_{s_}_vetores", vetores, st, "v_raio", [s_], "de_mais_y",
              g / f"{sufixo}_2{j}_{s_}_vetores_vrel.png", pulo=1)
    for j, s_ in enumerate(sup["cortes"]):
        for campo, rot in (("pressure", "pressao"), ("rel-velocity-magnitude", "vrel")):
            tenta(f"3{j}_{s_}_{rot}", contorno, st, "c_corte", campo, [s_], "de_mais_x", g / f"{sufixo}_3{j}_{s_}_{rot}.png")
    for peca, zs in grupos.items():
        base = [f"parede_{peca}"]
        for campo, rot in (("pressure", "pressao"), ("wall-shear", "tau_w"), ("y-plus", "yplus")):
            for vista in ("de_menos_x", "de_mais_x"):
                tenta(f"4_{peca}_{rot}_{vista}", contorno, st, "c_par", campo, base, vista,
                      g / f"{sufixo}_4_{peca}_{rot}_{vista}.png")
    for z in [z for v in PARES.values() for z in v]:
        for campo, rot in (("axial-velocity", "v_axial"), ("tangential-velocity", "v_tangencial"),
                           ("radial-velocity", "v_radial")):
            vista = "de_mais_x" if "bomba_turbina" in z else "iso"
            tenta(f"5_{z}_{rot}", contorno, st, "c_if", campo, [z], vista, g / f"{sufixo}_5_{z}_{rot}.png")
    # linhas de corrente a partir da saída da bomba
    try:
        pl = st.results.graphics.pathline
        if "pl_bt" not in pl.keys():
            pl.create("pl_bt")
        pl["pl_bt"].field = "velocity-magnitude"
        pl["pl_bt"].release_from_surfaces = ["if_bomba_turbina_lado_bomba"]
        pl["pl_bt"].skip = 40
        pl["pl_bt"].step = 4000
        pl["pl_bt"].display()
        for vista in ("iso", "de_mais_x", "meridional"):
            camera(st, *VISTAS[vista])
            res[f"6_pathlines_{vista}"] = salvar(st, g / f"{sufixo}_6_pathlines_{vista}.png")
    except Exception as e:
        log("pathlines falharam:", repr(e)[:300])
    # vórtices: iso-superfície de critério Q
    try:
        sf = st.results.surfaces.iso_surface
        if "q_iso" not in sf.keys():
            sf.create("q_iso")
        q = OMEGA_B ** 2 * 4
        sf["q_iso"] = {"field": "q-criterion", "zones": list(ZONAS), "iso_values": [q]}
        tenta("7_vortices_q", contorno, st, "c_q", "velocity-magnitude", ["q_iso"], "iso", g / f"{sufixo}_7_vortices_q.png")
    except Exception as e:
        log("critério Q falhou:", repr(e)[:300])
    return res


# ----------------------------------------------------------------------------------------------- comandos
def cmd_setup():
    """Monta o setup com sinal +1, testa imagens, mede o custo por iteração e grava o caso de setup."""
    s = abrir()
    try:
        ef = montar(s, +1)
        st = s.settings
        ef["teste_imagem"] = teste_imagem(st, DIR / "graficos")
        st.file.write_case(file_name=CASO_SETUP)
        (DIR / "dados" / "setup.json").write_text(json.dumps(ef, indent=1, default=str))
        log("SETUP EFETIVO:", json.dumps(ef, default=str)[:4000])
    except Exception:
        traceback.print_exc()
        raise
    finally:
        encerrar(s)


def preparar_sessao(sinal):
    """Abre o Fluent, monta o setup com o sinal pedido e devolve a sessão e o setup efetivo.

    Args:
        sinal (int): +1 ou −1.

    Returns:
        tuple: (sessão, setup efetivo).
    """
    s = abrir()
    ef = montar(s, sinal)
    return s, ef


def cmd_sentido(n_it=300):
    """SR = 0 com os dois sentidos, ``n_it`` iterações cada, a partir da mesma inicialização.

    Critério (registrado no C_TG.md): no sentido certo o conversor multiplica torque no estol, ou seja,
    TR = −T_t/T_b > 1, com T_b resistente (T_b < 0 no eixo de ω_b) e T_t motor (T_t > 0).

    Args:
        n_it (int, optional): Iterações por sentido. Padrão 300.
    """
    out = {}
    for sinal in (+1, -1):
        (DIR / "dados" / REP).unlink(missing_ok=True)
        s, ef = preparar_sessao(sinal)
        try:
            st = s.settings
            if sinal == +1:
                ef["teste_imagem"] = teste_imagem(st, DIR / "graficos")
            st.file.write_case(file_name=CASO_SETUP + ("" if sinal > 0 else "_menos"))
            t0 = time.time()
            m = st.solution.methods
            for k in m.discretization_scheme.keys():
                if k != "pressure":
                    m.discretization_scheme[k] = "first-order-upwind"
            st.solution.run_calculation.iterate(iter_count=ITER_FOU)
            for k in m.discretization_scheme.keys():
                if k != "pressure":
                    m.discretization_scheme[k] = "second-order-upwind"
            st.solution.run_calculation.iterate(iter_count=max(n_it - ITER_FOU, 0))
            dt = time.time() - t0
            f = finais(st, ef["reports"])
            tr = -f["t_turbina"] / f["t_bomba"] if f["t_bomba"] else float("nan")
            out[sinal] = {"reports": f, "TR": tr, "s_por_it": dt / n_it, "teste_imagem": ef.get("teste_imagem")}
            log(f"SENTIDO {sinal:+d}: TR={tr:.4f} T_b={f['t_bomba']:.2f} T_t={f['t_turbina']:.2f} "
                f"T_e={f['t_estator']:.2f} soma={f['t_total']:.3g} ({dt / n_it:.2f} s/it)")
            if sinal > 0:                         # testa todas as figuras antes da varredura
                pt = DIR / "teste_figuras"
                (pt / "graficos").mkdir(parents=True, exist_ok=True)
                (pt / "dados").mkdir(parents=True, exist_ok=True)
                try:
                    out[sinal]["figuras"] = figuras(st, pt, "teste", ef["grupos_torque"])
                    exportar(st, pt, "teste")
                except Exception as e:
                    log("teste de figuras falhou:", repr(e)[:300])
            st.file.write_case_data(file_name=f"CT_V1_MRF_sentido_{'mais' if sinal > 0 else 'menos'}")
            try:                                  # o Fluent mantém o .out aberto: copiar, não mover
                shutil.copy2(DIR / "dados" / REP, DIR / "dados" / f"sentido_{'mais' if sinal > 0 else 'menos'}_rep.out")
            except OSError as e:
                log("cópia do report file:", e)
            (DIR / "dados" / "setup.json").write_text(json.dumps(ef, indent=1, default=str)) if sinal > 0 else None
        except Exception:
            traceback.print_exc()
            out[sinal] = {"erro": traceback.format_exc()[-1500:]}
        finally:
            encerrar(s)
    (DIR / "dados" / "sentido.json").write_text(json.dumps(out, indent=1, default=str))
    log("SENTIDO RESUMO:", json.dumps(out, default=str))


def cmd_varredura(sinal, srs, ini=None):
    """Varre as razões de velocidade em sequência, cada ponto partindo do anterior.

    Retoma: um ponto com ``SRxxx/dados/resultado.json`` é pulado e o próximo parte do .dat dele.
    O primeiro ponto parte do caso do teste de sentido (``CT_V1_MRF_sentido_mais|menos``) quando existe.

    Args:
        sinal (int): +1 ou −1.
        srs (list[float]): Razões de velocidade.
        ini (str, optional): ``.dat.h5`` de partida do primeiro ponto (por exemplo um checkpoint); sem FOU.
    """
    s, ef = preparar_sessao(sinal)
    st = s.settings
    try:
        ef["teste_imagem"] = teste_imagem(st, DIR / "graficos")
        anterior = None
        for i, sr in enumerate(srs):
            p = pasta_sr(sr)
            rot = p.name
            if (p / "dados" / "resultado.json").exists():
                anterior = p
                log(f"{rot} já feito; pulando")
                continue
            if anterior is not None:
                st.file.read_data(file_name=str(anterior / "fluent" / f"CT_V1_MRF_{anterior.name}.dat.h5"))
                log(f"{rot}: parte de {anterior.name}")
            elif ini:
                st.file.read_data(file_name=str(ini))
                log(f"{rot}: parte de {ini}")
            ef_rot = aplicar_rotacao(st, sr, sinal, ef["sobras"])
            primeiro = anterior is None and not ini
            info = iterar_ponto(s, sr, primeiro, rot, n_min=2 * JANELA if (anterior is None and ini) else None)
            st.file.write_case_data(file_name=str(p / "fluent" / f"CT_V1_MRF_{rot}"))
            f = finais(st, ef["reports"])                # valores instantâneos da última iteração
            try:
                f["paredes_zona"] = torques_por_zona(st, ef["sobras"], sinal)
            except Exception as e:
                log("torques por zona falharam:", repr(e)[:300])
            info["reports_inst"] = f
            f = info["media_janela"]                      # valores oficiais: média das últimas JANELA iterações
            info.update({"reports": f, "rotacao": ef_rot, "omega_b": OMEGA_B, "omega_t": sr * OMEGA_B,
                         "rpm_b": RPM_B, "sinal": sinal, "unidade_torque": "N·m"})
            log(f"{rot} FINAL: TR={-f['t_turbina'] / f['t_bomba']:.4f} T_b={f['t_bomba']:.2f} "
                f"T_t={f['t_turbina']:.2f} T_e={f['t_estator']:.2f} soma={f['t_total']:.3g} conv={info['convergiu']}")
            exportar(st, p, rot)
            (p / "dados" / "resultado.json").write_text(json.dumps(info, indent=1, default=str))
            (DIR / "dados" / "setup.json").write_text(json.dumps(ef, indent=1, default=str))
            try:                                          # depois dos números: um travamento aqui não perde o ponto
                info["figuras"] = figuras(st, p, rot, ef["grupos_torque"])
            except Exception as e:
                log("figuras falharam:", repr(e)[:300])
            (p / "dados" / "resultado.json").write_text(json.dumps(info, indent=1, default=str))
            anterior = p
        (DIR / "dados" / "setup.json").write_text(json.dumps(ef, indent=1, default=str))
    except Exception:
        traceback.print_exc()
        raise
    finally:
        encerrar(s)


def cmd_figuras_teste(nome_caso):
    """Relê um caso com dados e gera todas as figuras e exportações em ``DIR/teste_figuras`` (depuração).

    Args:
        nome_caso (str): Caso em ``DIR/fluent`` (sem extensão), por exemplo ``CT_V1_MRF_sentido_mais``.
    """
    s = abrir()
    try:
        st = s.settings
        st.file.read_case_data(file_name=nome_caso)
        ef = json.loads((DIR / "dados" / "setup.json").read_text())
        pt = DIR / "teste_figuras"
        for sub in ("graficos", "dados"):
            (pt / sub).mkdir(parents=True, exist_ok=True)
        log("teste_imagem:", teste_imagem(st, pt / "graficos"))
        res = figuras(st, pt, "teste", ef["grupos_torque"])
        log("FIGURAS:", res)
        exportar(st, pt, "teste")
        log("EXPORT:", sorted(f.name for f in (pt / "dados").glob("*.csv")))
    except Exception:
        traceback.print_exc()
        raise
    finally:
        encerrar(s)


PAS = {"bomba": 20, "turbina": 18, "estator": 18}   # contadas pela harmônica angular das faces de parede
VARIANTES = {  # nome: descrição (todas no SR = 0; as de modelo e esquema partem da solução do SR000)
    "sst_sem_cc": "k-ω SST sem correção de curvatura",
    "ke_realizable": "k-ε Realizable com funções de parede",
    "fou": "primeira ordem (upwind) em momento e turbulência; PRESTO! na pressão",
    "meio_passo": "turbina girada de meio passo (10°) em relação à bomba e ao estator; setup refeito da malha",
}


def girar_turbina(s, graus):
    """Gira a zona de célula da turbina (com suas faces) em torno do eixo x, antes de criar as interfaces.

    Tenta o objeto de settings ``mesh.modify_zones.rotate_zone``; sem ele, usa o TUI
    ``/mesh/modify-zones/rotate-zone`` com todos os argumentos na chamada (sem prompts).

    Args:
        s (ansys.fluent.core.session_solver.Solver): Sessão com a malha lida.
        graus (float): Ângulo de rotação [°], positivo no sentido de +x.

    Returns:
        str: Caminho usado (``settings`` ou ``tui``).

    Raises:
        RuntimeError: Se nenhum dos caminhos existir.
    """
    st = s.settings
    mz = getattr(getattr(st, "mesh", None), "modify_zones", None)
    if mz is not None and hasattr(mz, "rotate_zone"):
        log("rotate_zone (settings) argumentos:", mz.rotate_zone.argument_names)
        mz.rotate_zone(zone_names=["fluido_turbina"], rotation_angle=graus, origin=[0.0, 0.0, 0.0],
                       axis=[1.0, 0.0, 0.0])
        return "settings"
    tui = s.tui.mesh.modify_zones
    if hasattr(tui, "rotate_zone"):
        tui.rotate_zone("fluido_turbina", "()", graus, 0, 0, 0, 1, 0, 0)
        return "tui"
    raise RuntimeError("rotate_zone não encontrado")


def cmd_sensibilidade(variante):
    """Roda uma variante do setup no SR = 0 (estol) e grava em ``DIR/sens/<variante>``.

    As variantes de modelo e esquema relêem o ``.cas/.dat`` do SR000, mudam só o item em estudo e iteram
    pelo mesmo critério da varredura (mínimo de 2 janelas). ``meio_passo`` refaz o setup da malha com a
    turbina girada e parte do repouso, como o SR000.

    Args:
        variante (str): Chave de ``VARIANTES``.

    Raises:
        KeyError: Variante desconhecida.
    """
    global MALHA, REP
    desc = VARIANTES[variante]
    # Report file próprio: se o arquivo já tem linhas além da iteração atual (variante anterior), o Fluent
    # não o sobrescreve e cria "<nome>_<it+1>.out", e o laço ficaria lendo o arquivo antigo para sempre.
    REP = time.strftime(f"CT_V1_MRF_sens_{variante}_%Y%m%d-%H%M%S.out")
    p = DIR / "sens" / variante
    for sub in ("fluent", "dados", "graficos"):
        (p / sub).mkdir(parents=True, exist_ok=True)
    s = abrir()
    try:
        st = s.settings
        ef = {"variante": variante, "descricao": desc}
        if variante == "meio_passo":
            st.file.read_mesh(file_name=str(MALHA))
            ef["rotacao_turbina"] = {"graus": 180.0 / PAS["turbina"], "via": girar_turbina(s, 180.0 / PAS["turbina"])}
            st.file.write_mesh(file_name=str(p / "fluent" / "CT_V1_3zonas_meio_passo.msh.h5"))
            malha_orig, MALHA = MALHA, p / "fluent" / "CT_V1_3zonas_meio_passo.msh.h5"
            try:
                ef.update(montar(s, +1))       # relê a malha girada e monta tudo de novo
            finally:
                MALHA = malha_orig
            primeiro, n_min = True, None
        else:
            st.file.read_case_data(file_name=str(DIR / "SR000" / "fluent" / "CT_V1_MRF_SR000.cas.h5"))
            st.solution.monitor.report_files["mon"].file_name = f"../dados/{REP}"
            ef["momentos"] = corrigir_momentos(st)
            v = st.setup.models.viscous
            m = st.solution.methods
            if variante == "sst_sem_cc":
                v.options.curvature_correction = False
            elif variante == "ke_realizable":
                v.model = "k-epsilon"
                v.k_epsilon_model = "realizable"
                wt = v.near_wall_treatment.wall_treatment
                opcoes = list(wt.allowed_values())
                ef["wall_treatment_opcoes"] = opcoes
                wt.set_state(next((o for o in opcoes if "scalable" in o), "standard-wall-fn"))
            elif variante == "fou":
                for k in m.discretization_scheme.keys():
                    if k != "pressure":
                        m.discretization_scheme[k] = "first-order-upwind"
            if variante != "fou":
                # Ao trocar o modelo, o Fluent põe as equações novas (k, ε) em 1ª ordem sozinho
                # ("Changing Discretization Scheme ... to First Order Upwind"); a variante só muda o modelo.
                for k in m.discretization_scheme.keys():
                    if k != "pressure":
                        m.discretization_scheme[k] = "second-order-upwind"
            ef.update({"viscous": v(), "methods": m()})
            primeiro, n_min = False, 2 * JANELA
        log(f"SENS {variante}: {desc}; setup:", json.dumps(ef, default=str)[:3000])
        info = iterar_ponto(s, 0.0, primeiro, f"sens_{variante}", n_min=n_min)
        st.file.write_case_data(file_name=str(p / "fluent" / f"CT_V1_MRF_sens_{variante}"))
        f = info["media_janela"]
        info.update({"reports": f, "setup": ef, "unidade_torque": "N·m"})
        try:
            sobras = ef.get("sobras") or setup_do_caso(st)["sobras"]
            info["paredes_zona"] = torques_por_zona(st, sobras, +1)
        except Exception as e:
            log("torques por zona falharam:", repr(e)[:300])
        log(f"SENS {variante} FINAL: TR={-f['t_turbina'] / f['t_bomba']:.4f} T_b={f['t_bomba']:.2f} "
            f"T_t={f['t_turbina']:.2f} T_e={f['t_estator']:.2f} soma={f['t_total']:.3g} conv={info['convergiu']}")
        exportar(st, p, f"sens_{variante}")
        (p / "dados" / "resultado.json").write_text(json.dumps(info, indent=1, default=str))
        try:
            shutil.copy2(DIR / "dados" / REP, p / "dados" / f"sens_{variante}_rep.out")
        except OSError as e:
            log("cópia do report file:", e)
    except Exception:
        traceback.print_exc()
        raise
    finally:
        encerrar(s)


def verificar_ponto(st, p, rot, sobras, euler=False):
    """Confere a conversão Cm → N·m, mede o torque das paredes de sobra e exporta as paredes com área.

    Lê os valores de referência, calcula os torques como estão no caso (Cm, se o caso é anterior à
    correção), troca para N·m e recalcula. Exporta as paredes com ``face-area-magnitude`` (para o y⁺
    ponderado por área). Com ``euler``, exporta os dois lados de cada interface com velocidades, áreas,
    viscosidades e (se aceitos) os gradientes de velocidade, para o balanço de Euler com o termo difusivo.

    Args:
        st: Raiz ``session.settings``.
        p (Path): Pasta do ponto.
        rot (str): Rótulo (``SRxxx``).
        sobras (list[str]): Paredes de sobra.
        euler (bool, optional): Também faz a parte das zonas interior.

    Returns:
        dict: Valores de referência, torques antes/depois, razão, sobras e (se pedido) as zonas interior.
    """
    out = {}
    try:
        out["reference_values"] = st.setup.reference_values()
    except Exception as e:
        out["reference_values"] = repr(e)[:200]
    nomes = ["t_bomba", "t_turbina", "t_estator", "t_total"]
    rd = st.solution.report_definitions.moment
    out["tipo_antes"] = {n: rd[n].report_output_type() for n in nomes}
    antes = finais(st, nomes)
    out["opcoes"] = corrigir_momentos(st)
    depois = finais(st, nomes)
    out.update({"antes": antes, "depois_nm": depois,
                "razao_nm_sobre_antes": {n: depois[n] / antes[n] for n in nomes if antes[n]}})
    log(f"{rot} VERIFICA Cm→N·m:", json.dumps(out, default=str)[:1500])
    try:
        for w in sobras:
            n = "s_" + w[:40]
            if n not in rd.keys():
                rd.create(n)
            rd[n] = {"mom_center": [0.0, 0.0, 0.0], "mom_axis": [1.0, 0.0, 0.0], "zones": [w]}
            momento_em_nm(rd[n])
        out["sobras_nm"] = finais(st, ["s_" + w[:40] for w in sobras])
        log(f"{rot} SOBRAS N·m:", out["sobras_nm"])
    except Exception as e:
        log("sobras falharam:", repr(e)[:300])
    campos = aceitos(st.file.export.ascii.cell_func_domain,
                     ["pressure", "wall-shear", "y-plus", "radial-coordinate", "face-area-magnitude",
                      "x-face-area", "y-face-area", "z-face-area"])
    for peca in ("bomba", "turbina", "estator"):
        try:
            st.file.export.ascii(file_name=str(p / "dados" / f"{rot}_parede_{peca}_area.csv"),
                                 surface_name_list=[f"parede_{peca}"], delimiter="comma", location="cell-center",
                                 cell_func_domain=campos)
        except Exception as e:
            log("export parede com área", peca, repr(e)[:200])
    if euler:
        # As zonas interior criadas pelas interfaces (if_bt, if_te, if_eb) não são aceitas como superfície
        # de report nem de exportação na 24.1 ("Values contain disallowed entries"); usam-se os dois
        # lados de cada interface, com viscosidade efetiva e gradientes para o termo difusivo.
        try:
            out["interior"] = list(st.setup.boundary_conditions.interior.keys())
            grad = [f"d{a}-velocity-d{b}" for a in "xyz" for b in "xyz"]
            ce = aceitos(st.file.export.ascii.cell_func_domain,
                         ["x-velocity", "y-velocity", "z-velocity", "density", "x-face-area", "y-face-area",
                          "z-face-area", "viscosity-eff", "viscosity-lam", "viscosity-ratio", "radial-coordinate",
                          "tangential-velocity"] + grad)
            out["euler_campos"] = ce
            out["euler_gradientes"] = [g for g in grad if g in ce]
            for z in [z for v in PARES.values() for z in v]:
                st.file.export.ascii(file_name=str(p / "dados" / f"{rot}_euler_{z}.csv"), surface_name_list=[z],
                                     delimiter="comma", location="cell-center", cell_func_domain=ce)
            log(f"{rot} Euler: lados das interfaces exportados; gradientes aceitos:", out["euler_gradientes"])
        except Exception as e:
            log("parte Euler/interfaces falhou:", repr(e)[:300])
    return out


def cmd_refazer_figuras(srs):
    """Refaz as figuras do Fluent de pontos prontos, relendo o ``.cas/.dat`` de cada um numa sessão nova.

    Na varredura as imagens saem em branco (só a legenda): a janela gráfica vem dos plots de monitor
    usados durante as iterações. Relendo o caso com dados, como em ``cmd_figuras_teste``, elas saem certas.
    As imagens antigas são sobrescritas e ``resultado.json`` ganha ``figuras_refeitas`` e ``verificacao``
    (conversão Cm → N·m, torque das sobras, paredes com área e zonas interior das interfaces).
    Sem ``dados/setup.json``, as sobras e os grupos de torque saem do primeiro caso lido
    (``setup_do_caso``) e o arquivo é gravado.

    Args:
        srs (list[float]): Razões de velocidade cujos pontos têm ``resultado.json``.
    """
    arq_setup = DIR / "dados" / "setup.json"
    ef = json.loads(arq_setup.read_text()) if arq_setup.exists() else None
    s = abrir()
    try:
        st = s.settings
        for sr in srs:
            p = pasta_sr(sr)
            rot = p.name
            caso = p / "fluent" / f"CT_V1_MRF_{rot}.cas.h5"
            if not (p / "dados" / "resultado.json").exists() or not caso.exists():
                log(f"{rot}: sem resultado ou caso; pulando")
                continue
            try:
                st.file.read_case_data(file_name=str(caso))
                if ef is None:                     # sem setup.json: tira sobras e grupos do próprio caso
                    ef = setup_do_caso(st)
                    arq_setup.write_text(json.dumps(ef, indent=1, default=str))
                    log("setup.json reconstruído do caso:", json.dumps(ef["sobras"]))
                try:
                    ver = verificar_ponto(st, p, rot, list(ef["sobras"]), euler=True)
                    info = json.loads((p / "dados" / "resultado.json").read_text())
                    info["verificacao"] = ver
                    (p / "dados" / "resultado.json").write_text(json.dumps(info, indent=1, default=str))
                except Exception:
                    traceback.print_exc()
                    log(f"{rot}: verificação falhou")
                log(f"{rot}: caso relido; teste_imagem:", teste_imagem(st, p / "graficos"))
                res = figuras(st, p, rot, ef["grupos_torque"])
                ok = sum(1 for v in res.values() if v)
                log(f"{rot} FIGURAS REFEITAS: {ok}/{len(res)}")
                info = json.loads((p / "dados" / "resultado.json").read_text())
                info["figuras_refeitas"] = res
                (p / "dados" / "resultado.json").write_text(json.dumps(info, indent=1, default=str))
            except Exception:
                traceback.print_exc()
                log(f"{rot}: falhou ao refazer figuras")
    finally:
        encerrar(s)


def opcoes_globais(argv):
    """Retira ``--malha <arquivo>`` e ``--saida <pasta>`` da linha de comando e ajusta ``MALHA`` e ``DIR``.

    Serve às malhas do GCI: ``--malha mesh/CT_V1_3zonas_G1.msh.h5 --saida CT_V1_GCI_G1`` grava em
    ``results/conversor/CT_V1_GCI_G1`` sem tocar nos resultados da malha base.

    Args:
        argv (list[str]): ``sys.argv``.

    Returns:
        list[str]: Argumentos restantes.
    """
    global MALHA, DIR
    argv = list(argv)
    for opc in ("--malha", "--saida"):
        if opc in argv:
            k = argv.index(opc)
            valor = argv[k + 1]
            del argv[k:k + 2]
            if opc == "--malha":
                MALHA = (ROOT / valor) if not Path(valor).is_absolute() else Path(valor)
            else:
                DIR = ROOT / "results" / "conversor" / valor
    log("MALHA =", MALHA, "| DIR =", DIR)
    return argv


if __name__ == "__main__":
    sys.argv = opcoes_globais(sys.argv)
    c = sys.argv[1]
    if c == "setup":
        cmd_setup()
    elif c == "sentido":
        cmd_sentido(int(sys.argv[2]) if len(sys.argv) > 2 else 300)
    elif c == "figuras_teste":
        cmd_figuras_teste(sys.argv[2])
    elif c == "refazer_figuras":
        cmd_refazer_figuras([float(x) for x in sys.argv[2:]] or SRS)
    elif c == "sensibilidade":
        for var in sys.argv[2:] or list(VARIANTES):
            try:
                cmd_sensibilidade(var)
            except Exception:
                log(f"SENS {var}: falhou; seguindo para a próxima")
    elif c == "varredura":
        args = sys.argv[3:]
        ini = None
        if "--ini" in args:
            k = args.index("--ini")
            ini = args[k + 1]
            del args[k:k + 2]
        cmd_varredura(int(sys.argv[2]), [float(x) for x in args] or SRS, ini)
    else:
        raise SystemExit(__doc__)
