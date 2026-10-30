"""Runner PyFluent: monta, resolve e exporta uma execução de ``casos.CASES``.

Uso: python tools/run_caso.py R02             (resolve e exporta)
     python tools/run_caso.py --exportar R02  (relê .cas/.dat e refaz só a exportação)
     python tools/run_caso.py --continuar R04 7000 [1e-6]
         (relê .cas/.dat e itera mais; com o 3º número, até o monitor variar menos que ele em 50 it)

Requer .venv241 (ansys-fluent-core 0.37.2, último com suporte à 2024 R1). Saídas em
results/<bloco>/<nome>/{fluent,dados,logs}.
"""
import json, sys, time, threading, traceback, subprocess
from pathlib import Path
import ansys.fluent.core as pf
from casos import ROOT, caso

ITER_SINAL = 50


def ler_rep(path):
    """Lê as linhas numéricas de um Report File do Fluent.

    Args:
        path (str | Path): Caminho do arquivo ``.out`` gravado pelo Report File.

    Returns:
        list[tuple[int, float]]: Pares (iteração, valor), na ordem do arquivo.
    """
    out = []
    for ln in Path(path).read_text().splitlines():
        p = ln.split()
        if len(p) >= 2 and p[0].isdigit():
            out.append((int(p[0]), float(p[1])))
    return out


def parede_movel(st, zona, u, relativa=False):
    """Define uma parede como móvel, translação em (1, 0).

    Args:
        st: Raiz ``session.settings``.
        zona (str): Nome da zona de parede.
        u (float): Velocidade da parede [m/s].
        relativa (bool, optional): Se True, a velocidade é relativa à zona adjacente
            (referencial rotativo); senão, absoluta. Padrão False.
    """
    w = st.setup.boundary_conditions.wall[zona].momentum
    w.motion_bc = "Moving Wall"
    w.relative = relativa
    w.vmag = {"option": "value", "value": u}
    w.wall_translation = [1, 0]


def linha_monitor(st, nome, p0, p1, campo="x-velocity"):
    """Cria uma linha, um report de média por área nela, Report File e Report Plot.

    Args:
        st: Raiz ``session.settings``.
        nome (str): Nome-base do monitor (a linha e o report usam o mesmo nome).
        p0 (list[float]): Ponto inicial (x, y, z).
        p1 (list[float]): Ponto final (x, y, z).
        campo (str, optional): Campo monitorado. Padrão ``"x-velocity"``.
    """
    st.results.surfaces.line_surface.create(nome)
    st.results.surfaces.line_surface[nome] = {"p0": p0, "p1": p1}
    rd = st.solution.report_definitions.surface
    rd.create("rd_" + nome)
    rd["rd_" + nome] = {"report_type": "surface-areaavg", "field": campo, "surface_names": [nome]}


def setup(s, c):
    """Monta o caso na sessão a partir dos parâmetros de ``casos.caso``.

    Comum a todos: malha, laminar, energia off, material ``fluido-tg``, métodos, URFs,
    resíduos absolutos de 1e-6, reports e inicialização standard com zeros. Depois aplica
    o contorno do tipo (``cp_per``, ``cp_io``, ``ldc``, ``cavr`` ou ``bfs``). Os valores são lidos
    de volta da sessão para o registro.

    Args:
        s (ansys.fluent.core.session_solver.Solver): Sessão do solver já conectada.
        c (dict): Parâmetros completos do caso.

    Returns:
        dict: Setup efetivo lido da sessão.

    Raises:
        AssertionError: Se a periodicidade ou o gradiente periódico não forem aplicados.
    """
    st = s.settings
    st.file.read_mesh(file_name=str(ROOT / "mesh" / f"{c['malha']}.msh"))
    bc = st.setup.boundary_conditions
    t = c["tipo"]
    ef = {}

    if t == "cp_per":
        # make-periodic refeito a cada malha (guia, passo Fluent 5)
        s.tui.mesh.modify_zones.make_periodic("periodic_a", "periodic_b", "no", "yes", "yes")
        assert sorted(bc.periodic.keys()) == ["periodic_a", "periodic_b"]

    st.setup.models.viscous.model = "laminar"           # padrão da sessão nova é k-omega
    st.setup.models.energy.enabled = False
    fl = st.setup.materials.fluid
    fl.make_a_copy(from_="air", to="fluido-tg")
    fl["fluido-tg"].density = {"option": "constant", "value": c["rho"]}
    fl["fluido-tg"].viscosity = {"option": "constant", "value": c["mu"]}
    zona = list(st.setup.cell_zone_conditions.fluid.keys())[0]
    st.setup.cell_zone_conditions.fluid[zona].material = "fluido-tg"

    if t == "cp_per":
        parede_movel(st, "wall_top", c["u_parede"])
        # Sem equivalente na API de settings da 241: prompts = gradiente, dir x, dir y.
        s.execute_tui(f"/define/periodic-conditions/pressure-gradient-specification {c['dpdx']} -1 0")
        ef["periodic_dpdx"] = s.scheme.eval("(rpgetvar 'periodic/pressure-derivative)")
        assert abs(ef["periodic_dpdx"] - c["dpdx"]) < 1e-12, f"gradiente não aplicado: {ef['periodic_dpdx']}"
        linha_monitor(st, "line_x2p5", [2.5, -1, 0], [2.5, 0, 0])
    elif t == "cp_io":
        parede_movel(st, "wall_top", c["u_parede"])
        vel = bc.velocity_inlet["inlet"].momentum.velocity
        if "perfil" in c:
            st.file.read_profile(file_name=str(ROOT / "refs" / c["perfil"]))
            vel.set_state({"option": "profile", "profile_name": "cp-inlet", "field_name": "u"})
        else:
            vel.set_state({"option": "value", "value": c["uin"]})
        bc.pressure_outlet["outlet"].momentum.gauge_pressure = {"option": "value", "value": 0}
        linha_monitor(st, "line_x2p5", [2.5, -1, 0], [2.5, 0, 0])
        for nome, sup in (("p_in", "inlet"), ("p_out", "outlet")):
            st.solution.report_definitions.surface.create(nome)
            st.solution.report_definitions.surface[nome] = {"report_type": "surface-areaavg",
                                                           "field": "pressure", "surface_names": [sup]}
        ef["inlet"] = bc.velocity_inlet["inlet"].momentum()
        ef["outlet"] = bc.pressure_outlet["outlet"].momentum.gauge_pressure()
    elif t == "ldc":
        parede_movel(st, "lid", c["u_parede"])
        linha_monitor(st, "line_x0p5", [0.5, 0, 0], [0.5, 1, 0])
    elif t == "cavr":
        # Frame Motion em 2D planar (eixo z implícito); nomes lidos pela API em tools/teste_frame_motion.py.
        rf = st.setup.cell_zone_conditions.fluid[zona].reference_frame
        rf.frame_motion = True
        rf.mrf_omega = {"option": "value", "value": c["omega"]}
        rf.reference_frame_axis_origin = c["eixo"]
        # Demais paredes ficam Stationary Wall: no Fluent, paradas em relação à zona adjacente (que gira).
        parede_movel(st, "lid", c["u_parede"], relativa=True)
        linha_monitor(st, "line_x0p5", [0.5, 0, 0], [0.5, 1, 0], campo="relative-x-velocity")
        ef["reference_frame"] = rf()
    elif t == "bfs":
        st.file.read_profile(file_name=str(ROOT / "refs" / c["perfil"]))
        bc.velocity_inlet["inlet"].momentum.velocity.set_state(
            {"option": "profile", "profile_name": "bfs-inlet", "field_name": "u"})
        bc.pressure_outlet["outlet"].momentum.gauge_pressure = {"option": "value", "value": 0}
        linha_monitor(st, "line_x2", [2, 0, 0], [2, 1.5, 0])
        ef["inlet"] = bc.velocity_inlet["inlet"].momentum()

    if t.startswith("cp"):
        fo = st.solution.report_definitions.force
        for nome, z in (("fx_top", "wall_top"), ("fx_bottom", "wall_bottom")):
            fo.create(nome)
            fo[nome].set_state({"force_vector": [1, 0, 0], "zones": [z]})

    rf = st.solution.monitor.report_files
    rf.create("monitor-rfile")
    mon = [k for k in st.solution.report_definitions.surface.keys() if k.startswith("rd_")]
    rf["monitor-rfile"] = {"file_name": f"../dados/{c['nome']}_rep.out", "report_defs": mon}
    st.solution.monitor.report_plots.create("monitor-rplot")
    st.solution.monitor.report_plots["monitor-rplot"] = {"report_defs": mon}

    m = st.solution.methods
    m.p_v_coupling.flow_scheme = c["alg"]                # padrão da sessão nova é Coupled
    m.gradient_scheme = "least-square-cell-based"
    m.discretization_scheme["pressure"] = c["pres"]
    m.discretization_scheme["mom"] = c["mom"]
    m.high_order_term_relaxation.enable = False
    if c["alg"] != "Coupled":                            # no Coupled o objeto de URF fica inativo
        urf = st.solution.controls.under_relaxation
        urf["pressure"], urf["mom"] = c["urf_p"], c["urf_mom"]
        ef["urf"] = urf()

    res = st.solution.monitor.residual
    res.options.criterion_type = "absolute"
    for eq in ("continuity", "x-velocity", "y-velocity"):
        res.equations[eq] = {"monitor": True, "check_convergence": True, "absolute_criteria": 1e-6}

    ini = st.solution.initialization
    ini.initialization_type = "standard"
    ini.defaults = {"pressure": 0, "x-velocity": 0, "y-velocity": 0}
    ini.initialize()

    ef.update({
        "malha": c["malha"], "zona_fluido": zona, "general": st.setup.general.solver(),
        "viscous": st.setup.models.viscous.model(), "energy": st.setup.models.energy.enabled(),
        "rho": fl["fluido-tg"].density.value(), "mu": fl["fluido-tg"].viscosity.value(),
        "paredes": {k: v.momentum() for k, v in bc.wall.items()},
        "methods": {"p_v": m.p_v_coupling.flow_scheme(), "grad": m.gradient_scheme(),
                    "disc": m.discretization_scheme()},
        "residual": {k: v["absolute_criteria"] for k, v in res.equations().items()},
        "monitores": mon, "init": ini.defaults(), "iter_max": c["iter_max"],
    })
    print("SETUP EFETIVO:", ef, flush=True)
    return ef


def exportar(s, c):
    """Exporta o campo (Cell Center) e os reports finais.

    O histórico de resíduos não sai daqui: a análise o remonta dos logs, que são completos
    (o serviço de monitores do PyFluent devolve o histórico com lacunas após uma continuação).

    Args:
        s (ansys.fluent.core.session_solver.Solver): Sessão com a solução carregada.
        c (dict): Parâmetros completos do caso.

    Returns:
        dict: Valores finais dos reports de força e pressão (vazio se não houver).
    """
    st = s.settings
    exportar_campo(st, c, f"../dados/{c['nome']}_campo.csv")
    rds = [k for k in ("fx_top", "fx_bottom", "p_in", "p_out")
           if k in list(st.solution.report_definitions.force.keys()) + list(st.solution.report_definitions.surface.keys())]
    finais = {}
    if rds:
        for d in st.solution.report_definitions.compute(report_defs=rds):
            finais.update({k: v[0] for k, v in d.items()})
    print("REPORTS FINAIS:", finais, flush=True)
    return finais


def exportar_campo(st, c, arquivo):
    """Exporta o campo de todas as células (Cell Center, vírgula) para um CSV.

    Args:
        st: Raiz ``session.settings``.
        c (dict): Parâmetros completos do caso.
        arquivo (str): Caminho relativo ao diretório de trabalho do Fluent (``<dir>/fluent``).
    """
    # Guia, passo 14: Surfaces vazio exporta todas as células do domínio.
    st.file.export.ascii(file_name=arquivo, surface_name_list=[], delimiter="comma", location="cell-center",
                         cell_func_domain=["pressure", "total-pressure", "velocity-magnitude", "x-velocity",
                                           "y-velocity", "stream-function"]
                         + (["relative-x-velocity", "relative-y-velocity", "rel-velocity-magnitude"]
                            if c["tipo"] == "cavr" else []))


def iterar_com_quadros(s, c, rep, it):
    """Itera em blocos de ``c["anim_cada"]`` e exporta o campo ao fim de cada bloco (quadros da animação).

    Para quando o Fluent converge (o bloco termina antes do previsto) ou ao atingir ``it`` iterações.

    Args:
        s (ansys.fluent.core.session_solver.Solver): Sessão com o caso inicializado.
        c (dict): Parâmetros completos do caso.
        rep (Path): Report File do monitor, usado para saber a iteração atual.
        it (int): Teto de iterações.
    """
    anim = c["dir"] / "dados" / "anim"
    anim.mkdir(exist_ok=True)
    for q in anim.glob("*_it*.csv"):  # quadros de uma tentativa anterior não podem entrar no GIF
        q.unlink()
    feito = 0
    while feito < it:
        k = min(c["anim_cada"], it - feito)
        s.settings.solution.run_calculation.iterate(iter_count=k)
        n = ler_rep(rep)[-1][0]
        exportar_campo(s.settings, c, f"../dados/anim/{c['nome']}_it{n:05d}.csv")
        if n - feito < k:  # convergiu dentro do bloco
            break
        feito = n


def encerrar(s):
    """Encerra o Fluent lançado por este script, sem travar o cliente.

    Na 24.1 com a 0.37.2, nem ``exit()`` nem ``/exit`` terminam os processos:
    ``/exit`` só derruba o servidor gRPC (e nunca responde à chamada, por isso vai
    numa thread). Se os processos não terminarem em 30 s, a árvore do Cortex
    desta sessão é encerrada pelo PID.

    Args:
        s (ansys.fluent.core.session_solver.Solver): Sessão aberta com ``cleanup_on_exit=False``.

    Returns:
        bool: True se, ao final, não restar processo desta sessão.
    """
    cp = s._fluent_connection.connection_properties
    pids = (cp.cortex_pid, cp.fluent_host_pid)
    t = threading.Thread(target=lambda: s.execute_tui("/exit yes"), daemon=True)
    t.start()
    t.join(15)
    s.exit()
    if not s.wait_process_finished(wait=30):
        print(f"Fluent não terminou; encerrando a árvore desta sessão (PIDs {pids})", flush=True)
        for pid in pids:
            subprocess.run(["taskkill", "/T", "/F", "/PID", str(pid)], capture_output=True)
    ok = s.wait_process_finished(wait=15)
    print("Fluent encerrado:", ok, flush=True)
    return ok


def lancar(c):
    """Cria as pastas da execução e abre o Fluent 2D, dupla precisão, serial, sem GUI.

    O diretório de trabalho do Fluent é ``<dir>/fluent``: .cas/.dat e transcrições caem lá.

    Args:
        c (dict): Parâmetros completos do caso.

    Returns:
        ansys.fluent.core.session_solver.Solver: Sessão conectada.
    """
    for sub in ("fluent", "dados", "graficos", "logs"):
        (c["dir"] / sub).mkdir(parents=True, exist_ok=True)
    return pf.launch_fluent(product_version=pf.FluentVersion.v241, dimension=pf.Dimension.TWO,
                            precision=pf.Precision.DOUBLE, processor_count=1, ui_mode=pf.UIMode.NO_GUI,
                            mode="solver", cwd=str(c["dir"] / "fluent"), cleanup_on_exit=False)


def arrumar(c):
    """Move transcrições para ``logs/`` e apaga scripts de limpeza que o Fluent deixa.

    Args:
        c (dict): Parâmetros completos do caso.
    """
    for f in (c["dir"] / "fluent").glob("*"):
        if f.suffix == ".trn" or f.name.endswith("error.log"):
            f.replace(c["dir"] / "logs" / f.name)
        elif f.name.startswith("cleanup-fluent"):
            f.unlink()


def run(cid):
    """Executa um caso: setup, solução, gravação e exportação.

    No CP periódico roda ``ITER_SINAL`` iterações para checar o sinal de ū antes de
    seguir até ``iter_max`` (o Fluent para sozinho ao convergir).

    Args:
        cid (str): ID do registro, por exemplo ``"R02"``.

    Raises:
        RuntimeError: Se ū não for positivo após ``ITER_SINAL`` iterações no CP periódico.
    """
    c = caso(cid)
    s = lancar(c)
    rep = c["dir"] / "dados" / f"{c['nome']}_rep.out"
    rep.unlink(missing_ok=True)
    try:
        ef = setup(s, c)
        t0 = time.time()
        it = c["iter_max"]
        if c["tipo"] == "cp_per":
            s.settings.solution.run_calculation.iterate(iter_count=ITER_SINAL)
            ub = ler_rep(rep)[-1][1]
            print(f"ubar após {ITER_SINAL} it = {ub:.6f}", flush=True)
            if ub <= 0:
                raise RuntimeError(f"ubar negativo após {ITER_SINAL} iterações: sinal do gradiente/direção errado")
            it -= ITER_SINAL
        if c.get("anim_cada"):
            iterar_com_quadros(s, c, rep, it)
        else:
            s.settings.solution.run_calculation.iterate(iter_count=it)
        ef["tempo_parede_s"] = round(time.time() - t0, 1)
        ef["pyfluent"] = pf.__version__
        s.settings.file.write_case_data(file_name=c["nome"])
        ef["reports_finais"] = exportar(s, c)
        (c["dir"] / "dados" / f"{c['nome']}_setup.json").write_text(json.dumps(ef, indent=1, default=str))
        print(f"tempo de parede (iterações) = {ef['tempo_parede_s']} s", flush=True)
    except Exception:
        traceback.print_exc()
        sys.stdout.flush()
        raise
    finally:
        encerrar(s)
        arrumar(c)


def reexportar(cid):
    """Relê o caso e os dados gravados e refaz só a exportação, sem resolver de novo.

    Args:
        cid (str): ID do registro, por exemplo ``"R02"``.
    """
    c = caso(cid)
    s = lancar(c)
    try:
        s.settings.file.read_case_data(file_name=c["nome"])
        exportar(s, c)
    except Exception:
        traceback.print_exc()
        sys.stdout.flush()
        raise
    finally:
        encerrar(s)
        arrumar(c)


def continuar(cid, n_iter, tol_mon=None):
    """Relê o caso e os dados de uma execução e segue iterando.

    Sem ``tol_mon``: itera até ``n_iter`` a mais (o Fluent para antes se os resíduos
    convergirem). Com ``tol_mon``: aplica a 2ª metade do critério do guia ("monitor estável
    por 50 iterações"): desliga o corte por resíduo e itera em blocos de 250 até a variação
    do monitor nas últimas 50 iterações ficar abaixo de ``tol_mon`` ou somar ``n_iter``.
    Os desvios do protocolo ficam em ``_setup.json`` (``continuacoes`` e tempo somado).
    Regrava .cas/.dat e refaz a exportação.

    Args:
        cid (str): ID do registro.
        n_iter (int): Máximo de iterações adicionais.
        tol_mon (float, optional): Variação máxima do monitor em 50 iterações [m/s].
    """
    c = caso(cid)
    sj = c["dir"] / "dados" / f"{c['nome']}_setup.json"
    ef = json.loads(sj.read_text())
    rep = c["dir"] / "dados" / f"{c['nome']}_rep.out"
    s = lancar(c)
    try:
        s.settings.file.read_case_data(file_name=c["nome"])
        t0 = time.time()
        it0 = ler_rep(rep)[-1][0]
        if tol_mon is None:
            s.settings.solution.run_calculation.iterate(iter_count=n_iter)
        else:
            res = s.settings.solution.monitor.residual
            for eq in ("continuity", "x-velocity", "y-velocity"):
                res.equations[eq].absolute_criteria = 1e-12
            while True:
                s.settings.solution.run_calculation.iterate(iter_count=250)
                v = [x for _, x in ler_rep(rep)[-50:]]
                feito = ler_rep(rep)[-1][0] - it0
                print(f"CONTINUACAO it+{feito}: variação do monitor em 50 it = {max(v) - min(v):.3e}", flush=True)
                if max(v) - min(v) < tol_mon or feito >= n_iter:
                    break
        extra = ler_rep(rep)[-1][0] - it0
        ef.setdefault("continuacoes", []).append(
            {"iter_extra": extra, "tol_monitor": tol_mon,
             "residuo_abs": 1e-12 if tol_mon else ef["residual"]["continuity"]})
        ef["iter_extra"] = ef.get("iter_extra", 0) + extra
        ef["tempo_parede_s"] = round(ef["tempo_parede_s"] + time.time() - t0, 1)
        s.settings.file.write_case_data(file_name=c["nome"])
        ef["reports_finais"] = exportar(s, c)
        sj.write_text(json.dumps(ef, indent=1, default=str))
    except Exception:
        traceback.print_exc()
        sys.stdout.flush()
        raise
    finally:
        encerrar(s)
        arrumar(c)


if __name__ == "__main__":
    if sys.argv[1] == "--exportar":
        reexportar(sys.argv[2])
    elif sys.argv[1] == "--continuar":
        continuar(sys.argv[2], int(sys.argv[3]), float(sys.argv[4]) if len(sys.argv) > 4 else None)
    else:
        run(sys.argv[1])
