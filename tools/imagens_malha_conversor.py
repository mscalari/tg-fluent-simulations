"""Imagens das malhas G0, G1 e G2 do conversor CT_V1 para o relatório, com enquadramento idêntico (sem GUI).

Uso (rodar com PYTHONIOENCODING=utf-8 e o .venv241):
    python tools/imagens_malha_conversor.py G2 [G1 G0] [--processos N]
    python tools/imagens_malha_conversor.py --so-montagem

Abre uma sessão do solver por malha (uma depois da outra, ``UIMode.NO_GUI``), lê só a malha, normaliza os nomes
das zonas (``conversor_mrf.normalizar_zonas``) e grava as vistas de ``VISTAS_MALHA``. Antes de cada sessão confere
a RAM livre (mínimo ``RAM_MIN_MB``). No fim monta, para cada vista, uma figura matplotlib lado a lado com as malhas
que existirem (G2 | G1 | G0, da grossa para a fina), rotuladas com o número de células lido do solver
(``/mesh/size-info``).

Enquadramento: todas as vistas partem de ``views.restore_view("front")`` (câmera em +z), giram com
``camera.orbit`` e fixam a janela com ``camera.dolly`` (deslocamento em metros, no plano da tela) e
``camera.field`` (largura × altura visíveis, em metros). Como os valores são fixos e a geometria é a mesma, o
enquadramento é idêntico entre as malhas. O plano z = 0 e o corte x = −10 mm são iso-superfícies (a
``plane_surface`` sai em branco na exibição de malha sem GUI). Nessas vistas de corte a imagem é composta em
Python: com a câmera parada, o Fluent grava a máscara de cada zona (contorno de ``cell-zone`` só daquela zona) e
as arestas da malha; ``compor_zonas`` pinta as zonas com ``CORES`` e multiplica pelas arestas (o ``draw_mesh``
do contorno e os overlays não funcionam sem GUI na 24.1).

Saídas em results/conversor/CT_V1_GCI/:
    graficos/malha_<nivel>_<nn>_<nome>.png     uma por vista e por malha
    graficos/malha_comparacao_<nome>.png       montagem lado a lado
    dados/malha_<nivel>_info.json              células, faces, nós, tempo e RAM da sessão
    fluent/imagens_malha/                      diretório de trabalho do Fluent (transcripts)
"""
import argparse, json, re, subprocess, sys, time, traceback
from pathlib import Path
import ansys.fluent.core as pf

sys.path.insert(0, str(Path(__file__).resolve().parent))
from casos import ROOT                     # noqa: E402
from run_caso import encerrar              # noqa: E402
import conversor_mrf as cm                 # noqa: E402

SAIDA = ROOT / "results" / "conversor" / "CT_V1_GCI"
GRAF, DADOS, TRAB = SAIDA / "graficos", SAIDA / "dados", SAIDA / "fluent" / "imagens_malha"
MALHAS = {"G0": ROOT / "mesh" / "CT_V1_3zonas.msh.h5",
          "G1": ROOT / "mesh" / "CT_V1_3zonas_G1.msh.h5",
          "G2": ROOT / "mesh" / "CT_V1_3zonas_G2.msh.h5"}
ORDEM = ("G2", "G1", "G0")      # da mais grossa para a mais fina (colunas da montagem)
RAM_MIN_MB = 3000
PROCESSOS = 2
# Alvo da vista frontal = centro da caixa da geometria exibida: x de −0,179 a 0,159 m, y de −0,36 a 0,36 m.
CENTRO_FRONTAL = (-0.010, 0.0)

CORES = {"bomba": (140, 180, 230), "turbina": (245, 170, 120), "estator": (150, 210, 150)}   # RGB
# Superfícies de cada grupo (só as que existirem na malha entram).
PAREDES = {
    "bomba": ["parede_bomba", "parede_bomba_lado_turbina", "parede_bomba_lado_estator"],
    "turbina": ["parede_turbina", "parede_turbina_lado_estator", "parede_turbina_lado_bomba"],
    "estator": ["parede_estator", "parede_estator_lado_turbina", "parede_estator_lado_bomba"],
}

# nn, nome, título, conteúdo, vista (chave de cm.VISTAS), centro, campo (largura, altura) em m.
# centro: na vista "meridional", (x, y) do mundo; nas outras, (direita, cima) da tela a partir do alvo.
VISTAS_MALHA = [
    ("01", "meridional", "Plano z = 0 (meridional), cor por zona", "z0", "meridional", (-0.010, 0.190), (0.48, 0.36)),
    ("02", "zoom_if_bomba_turbina", "Interface bomba | turbina (não conforme), plano z = 0", "z0", "meridional",
     (0.000, 0.290), (0.10, 0.075)),
    ("03", "zoom_estator", "Estator e interfaces com bomba e turbina, plano z = 0", "z0", "meridional",
     (0.050, 0.105), (0.24, 0.18)),
    ("04", "zoom_if_estator_bomba", "Interface estator | bomba (não conforme), plano z = 0", "z0", "meridional",
     (-0.050, 0.045), (0.06, 0.045)),
    ("05", "corte_x_bomba", "Corte axial x = −10 mm (canais da bomba), de +x", "x_bomba", "de_mais_x",
     (0.0, 0.0), (0.96, 0.72)),
    ("06", "zoom_corte_x_bomba", "Detalhe do corte x = −10 mm (canal entre pás da bomba), de +x", "x_bomba",
     "de_mais_x", (0.0, 0.27), (0.16, 0.12)),
    ("07", "superficie_turbina", "Malha de superfície da turbina, de −x", "turbina", "de_menos_x",
     (0.0, 0.17), (0.50, 0.375)),
    ("08", "superficie_estator", "Malha de superfície do estator, de +x", "estator", "de_mais_x",
     (0.0, 0.0), (0.48, 0.36)),
    ("09", "superficie_bomba", "Malha de superfície da bomba, de +x", "bomba", "de_mais_x",
     (0.0, 0.17), (0.50, 0.375)),
    ("10", "zoom_superficie_bomba", "Detalhe da superfície da bomba (pás), de +x", "bomba", "de_mais_x",
     (0.0, 0.24), (0.20, 0.15)),
]


def log(*a):
    """Imprime com horário e esvazia o buffer.

    Args:
        *a: Itens a imprimir.
    """
    print(time.strftime("[%H:%M:%S]"), *a, flush=True)


def ram_livre_mb():
    """Lê a memória física livre do Windows.

    Returns:
        int: RAM livre em MB.
    """
    r = subprocess.run(["powershell", "-NoProfile", "-Command",
                        "[int]((Get-CimInstance Win32_OperatingSystem).FreePhysicalMemory/1024)"],
                       capture_output=True, text=True)
    return int(r.stdout.strip())


def ram_sessao_mb(pids):
    """Soma o working set da árvore de processos de uma sessão do Fluent.

    Args:
        pids (tuple[int, ...]): PIDs raiz (Cortex e host do Fluent).

    Returns:
        int: Memória em MB (0 se a consulta falhar).
    """
    ps = ("$p = Get-CimInstance Win32_Process | Select-Object ProcessId,ParentProcessId,WorkingSetSize; "
          "$p | ForEach-Object { \"$($_.ProcessId) $($_.ParentProcessId) $($_.WorkingSetSize)\" }")
    try:
        r = subprocess.run(["powershell", "-NoProfile", "-Command", ps], capture_output=True, text=True, timeout=60)
        filhos, ws = {}, {}
        for lin in r.stdout.splitlines():
            pid, pai, mem = (int(x) for x in lin.split())
            filhos.setdefault(pai, []).append(pid)
            ws[pid] = mem
        vistos, pilha = set(), list(pids)
        while pilha:
            p = pilha.pop()
            if p not in vistos:
                vistos.add(p)
                pilha.extend(filhos.get(p, []))
        return int(sum(ws.get(p, 0) for p in vistos) / 2 ** 20)
    except Exception:              # noqa: BLE001
        return 0


def tamanho_malha(s, nivel):
    """Lê células, faces e nós da malha pelo ``/mesh/size-info`` gravado num transcript.

    Args:
        s (ansys.fluent.core.session_solver.Solver): Sessão com a malha lida.
        nivel (str): ``"G0"``, ``"G1"`` ou ``"G2"`` (nome do transcript).

    Returns:
        dict: ``celulas``, ``faces``, ``nos`` (``None`` se a leitura falhar).
    """
    trn = TRAB / f"size_info_{nivel}.trn"
    trn.unlink(missing_ok=True)
    s.tui.file.start_transcript(f'"{trn.as_posix()}"')
    try:
        s.tui.mesh.size_info()
    finally:
        s.tui.file.stop_transcript()
    txt = trn.read_text(encoding="utf-8", errors="replace") if trn.exists() else ""
    m = re.search(r"Level\s+Cells\s+Faces\s+Nodes.*\n\s*\d+\s+(\d+)\s+(\d+)\s+(\d+)", txt)
    if not m:
        log("AVISO: /mesh/size-info não lido; ver", trn)
        return {"celulas": None, "faces": None, "nos": None}
    return {"celulas": int(m.group(1)), "faces": int(m.group(2)), "nos": int(m.group(3))}


def criar_superficies(st):
    """Cria as iso-superfícies z = 0 e x = −10 mm (uma por zona) e lista as paredes de cada peça.

    A ``plane_surface`` (``xy-plane``) sai em branco na exibição de malha sem GUI (sonda de 2026-11-21 na G2);
    a iso-superfície de ``z-coordinate`` funciona e permite restringir a zona, o que dá a máscara de cada zona.

    Args:
        st: Raiz ``session.settings`` com a malha lida e as zonas normalizadas.

    Returns:
        dict: ``z0`` e ``x_bomba`` -> {peça: superfície}; ``bomba``, ``turbina``, ``estator`` -> lista de paredes.
    """
    iso = st.results.surfaces.iso_surface
    out = {}
    for chave, campo, valor in (("z0", "z-coordinate", 0.0), ("x_bomba", "x-coordinate", -0.010)):
        out[chave] = {}
        for z in cm.ZONAS:
            n = f"malha_{chave}_{cm.PECA[z]}"
            try:
                if n not in iso.keys():
                    iso.create(n)
                iso[n] = {"field": campo, "zones": [z], "iso_values": [valor]}
                out[chave][cm.PECA[z]] = n
            except Exception as e:          # noqa: BLE001  (o corte pode não cruzar a zona)
                log(f"iso {n}:", repr(e)[:150])
    paredes = set(st.setup.boundary_conditions.wall.keys())
    for peca, nomes in PAREDES.items():
        out[peca] = [n for n in nomes if n in paredes]
    return out


def objetos_graficos(st):
    """Cria os objetos gráficos: arestas (malha sem faces), paredes (faces + arestas) e contorno de ``cell-zone``.

    O contorno só serve de máscara de cada zona: o ``draw_mesh`` do contorno fica inativo sem GUI e o
    ``/display/set/overlays yes`` não sobrepõe os objetos de settings (sonda de 2026-11-21 na G2). Por isso a
    legenda é desligada e a imagem final é composta em ``compor_zonas``.

    Args:
        st: Raiz ``session.settings``.

    Returns:
        tuple: (objeto de malha das arestas, objeto de malha das paredes, contorno de zonas).
    """
    g = st.results.graphics
    for n in ("malha_arestas", "malha_paredes"):
        if n not in g.mesh.keys():
            g.mesh.create(n)
    arestas, paredes = g.mesh["malha_arestas"], g.mesh["malha_paredes"]
    # Só arestas (faces=False) sai em branco sem GUI; com faces, as da iso-superfície ficam brancas.
    arestas.options = {"nodes": False, "edges": True, "faces": True, "partitions": False}
    # Com a cor por tipo, as arestas da iso-superfície saem brancas (somem no fundo branco); por id saem pretas.
    arestas.coloring = {"option": "automatic", "automatic": {"option": "id"}}
    arestas.edge_type = {"option": "all"}
    paredes.options = {"nodes": False, "edges": True, "faces": True, "partitions": False}
    paredes.edge_type = {"option": "all"}
    if "zonas" not in g.contour.keys():
        g.contour.create("zonas")
    c = g.contour["zonas"]
    c.field = "cell-zone"
    c.filled = True
    c.node_values = False
    c.boundary_values = False
    try:
        c.color_map.visible = False
    except Exception as e:                    # noqa: BLE001
        log("color_map.visible:", repr(e)[:150])
    return arestas, paredes, c


def enquadrar(st, vista, centro, campo):
    """Aplica a vista padrão do projeto e fixa o centro e o tamanho da janela.

    Args:
        st: Raiz ``session.settings``.
        vista (str): Chave de ``conversor_mrf.VISTAS``.
        centro (tuple[float, float]): Na vista ``meridional``, (x, y) do mundo; nas demais, deslocamento
            (direita, cima) da tela a partir do alvo, em m.
        campo (tuple[float, float]): Largura e altura visíveis, em m.
    """
    cm.camera(st, *cm.VISTAS[vista])
    cam = st.results.graphics.views.camera
    if vista == "meridional":
        direita, cima = centro[0] - CENTRO_FRONTAL[0], centro[1] - CENTRO_FRONTAL[1]
    else:
        direita, cima = centro
    if direita or cima:
        cam.dolly(right=-direita, up=cima, in_=0.0)   # na 24.1 right > 0 leva a vista para −x (conferido na G2)
    cam.field(width=campo[0], height=campo[1])


def compor_zonas(mascaras, arestas, destino):
    """Pinta cada zona com a sua cor e multiplica pelas arestas da malha (mesma câmera em todas as imagens).

    Args:
        mascaras (dict[str, Path]): Peça -> PNG do contorno só daquela zona (fundo branco).
        arestas (Path): PNG só com as arestas da malha (traço escuro em fundo branco).
        destino (Path): PNG composto.

    Returns:
        bool: True se gravou.
    """
    import numpy as np
    from PIL import Image, ImageDraw, ImageFont

    base = np.asarray(Image.open(arestas).convert("L"), dtype=float) / 255.0
    rgb = np.ones(base.shape + (3,))
    presentes = []
    for peca, arq in mascaras.items():
        m = np.asarray(Image.open(arq).convert("RGB"), dtype=float) / 255.0
        dentro = m.min(axis=2) < 0.9                     # qualquer cor não branca é a zona
        if dentro.sum() > 50:
            presentes.append(peca)
        rgb[dentro] = np.array(CORES[peca]) / 255.0
    img = Image.fromarray((rgb * base[..., None] * 255).astype("uint8"))
    d = ImageDraw.Draw(img)
    try:
        fonte = ImageFont.truetype("arial.ttf", 30)
    except OSError:
        fonte = ImageFont.load_default()
    if presentes:
        d.rectangle([10, 10, 190, 20 + 44 * len(presentes)], fill=(255, 255, 255), outline=(150, 150, 150))
    y = 20
    for peca in presentes:
        d.rectangle([20, y, 60, y + 30], fill=CORES[peca], outline=(0, 0, 0))
        d.text((72, y - 2), peca, fill=(0, 0, 0), font=fonte)
        y += 44
    destino.unlink(missing_ok=True)
    img.save(destino)
    return destino.exists()


def vista_por_zona(st, objs, sup, vista, centro, campo, destino):
    """Grava uma vista de corte colorida por zona com a malha por cima.

    Exibe o contorno de todas as zonas e enquadra; depois, sem mexer na câmera, grava o contorno de cada zona
    sozinha (máscara) e as arestas de todas; ``compor_zonas`` junta as imagens.

    Args:
        st: Raiz ``session.settings``.
        objs (tuple): Saída de ``objetos_graficos``.
        sup (dict[str, str]): Peça -> iso-superfície daquela zona.
        vista (str): Chave de ``conversor_mrf.VISTAS``.
        centro (tuple[float, float]): Ver ``enquadrar``.
        campo (tuple[float, float]): Ver ``enquadrar``.
        destino (Path): PNG final.

    Returns:
        bool: True se a imagem composta foi gravada.
    """
    arestas, _paredes, contorno = objs
    todas = list(sup.values())
    contorno.surfaces_list = todas
    contorno.display()
    enquadrar(st, vista, centro, campo)
    tmp = TRAB / "tmp"
    tmp.mkdir(exist_ok=True)
    mascaras = {}
    for peca, n in sup.items():
        contorno.surfaces_list = [n]
        contorno.display()
        if cm.salvar(st, tmp / f"zona_{peca}.png"):
            mascaras[peca] = tmp / f"zona_{peca}.png"
    arestas.surfaces_list = todas
    arestas.display()
    if not cm.salvar(st, tmp / "arestas.png"):
        return False
    return compor_zonas(mascaras, tmp / "arestas.png", destino)


def gerar(nivel, processos):
    """Abre o solver, lê a malha de um nível e grava todas as vistas de ``VISTAS_MALHA``.

    Args:
        nivel (str): ``"G0"``, ``"G1"`` ou ``"G2"``.
        processos (int): Processos MPI do solver.

    Returns:
        dict: Resumo gravado em ``DADOS/malha_<nivel>_info.json``.

    Raises:
        FileNotFoundError: Se a malha do nível não existir.
        RuntimeError: Se a RAM livre estiver abaixo de ``RAM_MIN_MB``.
    """
    malha = MALHAS[nivel]
    if not malha.exists():
        raise FileNotFoundError(malha)
    livre = ram_livre_mb()
    log(f"{nivel}: RAM livre {livre} MB")
    if livre < RAM_MIN_MB:
        raise RuntimeError(f"RAM livre {livre} MB < {RAM_MIN_MB} MB; não abro o Fluent para {nivel}")
    for d in (GRAF, DADOS, TRAB):
        d.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    s = pf.launch_fluent(product_version=pf.FluentVersion.v241, dimension=pf.Dimension.THREE,
                         precision=pf.Precision.DOUBLE, processor_count=processos, ui_mode=pf.UIMode.NO_GUI,
                         mode="solver", cwd=str(TRAB), cleanup_on_exit=False)
    cp = s._fluent_connection.connection_properties
    pids = (cp.cortex_pid, cp.fluent_host_pid)
    info = {"nivel": nivel, "malha": malha.name, "processos": processos, "pids": pids, "imagens": {}}
    ram = []
    try:
        st = s.settings
        st.file.read_mesh(file_name=str(malha))
        cm.normalizar_zonas(st)
        cm.renomear_paredes_stl(st)
        info.update(tamanho_malha(s, nivel))
        log(f"{nivel}: {info['celulas']} células, {info['faces']} faces, {info['nos']} nós "
            f"(lida em {time.time() - t0:.0f} s)")
        ram.append(ram_sessao_mb(pids))
        cm.preparar_picture(st)
        sup = criar_superficies(st)
        objs = objetos_graficos(st)
        for nn, nome, _titulo, conteudo, vista, centro, campo in VISTAS_MALHA:
            arq = GRAF / f"malha_{nivel}_{nn}_{nome}.png"
            try:
                if conteudo in ("z0", "x_bomba"):
                    ok = vista_por_zona(st, objs, sup[conteudo], vista, centro, campo, arq)
                else:
                    objs[1].surfaces_list = sup[conteudo]
                    objs[1].display()
                    enquadrar(st, vista, centro, campo)
                    ok = cm.salvar(st, arq)
            except Exception as e:        # noqa: BLE001
                log(f"{nivel} {nn}_{nome}: falhou", repr(e)[:300])
                ok = False
            info["imagens"][f"{nn}_{nome}"] = ok
            log(f"{nivel} {nn}_{nome}: {'ok' if ok else 'FALHOU'}")
        ram.append(ram_sessao_mb(pids))
    except Exception:
        traceback.print_exc()
        info["erro"] = traceback.format_exc()[-2000:]
    finally:
        info["tempo_s"] = round(time.time() - t0, 1)
        info["ram_sessao_mb_max"] = max(ram) if ram else None
        (DADOS / f"malha_{nivel}_info.json").write_text(json.dumps(info, indent=2, ensure_ascii=False),
                                                        encoding="utf-8")
        encerrar(s)
    log(f"{nivel}: {info['tempo_s']} s, RAM da sessão até {info['ram_sessao_mb_max']} MB")
    return info


def montagem():
    """Monta, para cada vista, as imagens das malhas existentes lado a lado (G2 | G1 | G0).

    Returns:
        list[Path]: Figuras gravadas em ``GRAF/malha_comparacao_<nome>.png``.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.image as mpimg

    def rotulo(nivel):
        """Rótulo da coluna com o número de células lido do solver.

        Args:
            nivel (str): ``"G0"``, ``"G1"`` ou ``"G2"``.

        Returns:
            str: ``"G2 — 459 592 células"`` (ou só o nível, se o JSON não existir).
        """
        arq = DADOS / f"malha_{nivel}_info.json"
        n = json.loads(arq.read_text(encoding="utf-8")).get("celulas") if arq.exists() else None
        return f"{nivel} — {n:,} células".replace(",", " ") if n else nivel

    feitas = []
    for nn, nome, titulo, *_ in VISTAS_MALHA:
        imgs = [(nv, GRAF / f"malha_{nv}_{nn}_{nome}.png") for nv in ORDEM]
        imgs = [(nv, a) for nv, a in imgs if a.exists()]
        if not imgs:
            continue
        fig, axs = plt.subplots(1, len(imgs), figsize=(6.0 * len(imgs), 5.0), squeeze=False)
        for ax, (nv, a) in zip(axs[0], imgs):
            ax.imshow(mpimg.imread(a))
            ax.set_title(rotulo(nv), fontsize=13)
            ax.set_xticks([])
            ax.set_yticks([])
            for sp in ax.spines.values():
                sp.set_color("0.7")
        fig.suptitle(titulo, fontsize=14)
        fig.tight_layout(rect=(0, 0, 1, 0.93))   # reserva o topo para o suptitle
        dest = GRAF / f"malha_comparacao_{nome}.png"
        fig.savefig(dest, dpi=150)
        plt.close(fig)
        feitas.append(dest)
        log("montagem:", dest.name, "com", [nv for nv, _ in imgs])
    return feitas


def main():
    """Lê os argumentos, gera as imagens de cada malha pedida (uma sessão por vez) e monta as comparações."""
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("niveis", nargs="*", choices=list(MALHAS), help="malhas a processar, na ordem dada")
    ap.add_argument("--processos", type=int, default=PROCESSOS, help=f"processos MPI (padrão {PROCESSOS})")
    ap.add_argument("--so-montagem", action="store_true", help="só refaz as montagens com os PNGs existentes")
    a = ap.parse_args()
    if not a.niveis and not a.so_montagem:
        ap.error("informe ao menos um nível (G0, G1, G2) ou --so-montagem")
    if not a.so_montagem:
        for nivel in a.niveis:
            gerar(nivel, a.processos)
    montagem()


if __name__ == "__main__":
    main()
