"""Gera as malhas grossas G1 e G2 do conversor CT_V1 (estudo de independência de malha pelo GCI) -- RASCUNHO.

Reproduz, pelo PyFluent (Fluent Meshing 2024 R1, ``ansys-fluent-core==0.37.2``), o mesmo Fault-tolerant Meshing
usado na malha de referência ``mesh/CT_V1_3zonas.msh.h5`` (G0, 1 900 736 células), trocando só os tamanhos.
Fontes dos nomes e valores: ``mesh/CT_V1_3zonas_workflow_files/CT_V1_3zonas.wft`` (tarefas e argumentos),
seção "Volume travado na GUI, refeito sem GUI e sem prismas" do ``C_TG.md`` (zonas, TUI de junção) e
``tools/conversor_mrf.py`` (nomes de zona exigidos pelo setup MRF).

Uso (NÃO foi executado; rodar com ``PYTHONIOENCODING=utf-8`` e o ``.venv241``)::

    python tools/malha_gci_conversor.py G1 G2                  gera as duas malhas, uma depois da outra
    python tools/malha_gci_conversor.py G1 --so-superficie     para depois da malha de superfície (calibração)
    python tools/malha_gci_conversor.py G1 --partir-da-superficie
                                                               copia CT_V1_3zonas_superficie.msh.h5 (G0) para a pasta de
                                                               trabalho, relê (restaura o workflow até Update Boundaries) e
                                                               refaz a partir do controle global e do dimensionamento
    python tools/malha_gci_conversor.py G1 --retomar-superficie
                                                               parte de CT_V1_3zonas_G1_superficie.msh.h5 e só gera o volume
    python tools/malha_gci_conversor.py G1 --fator 1.4         substitui o fator (só com um nível; calibração)

Receita (f = fator de escala do nível em relação à G0; G1: f = 1,35; G2: f = 1,35² = 1,8225)::

    comprimentos (MinSize, MaxSize, GlobalMax, GlobalMin)               x f
    taxa de crescimento                                                 1 + 0,2 f   (limite 1,34)
    ângulo normal de curvatura                                          18° x f
    células por folga (proximidade)                                     3 / f       (limite inferior 1,775)
    Wrap*: 2/3 dos valores acima para tamanhos; demais iguais aos da superfície (como no .wft)

Saídas: ``mesh/CT_V1_3zonas_G1.msh.h5`` e ``mesh/CT_V1_3zonas_G2.msh.h5`` (mais ``*_workflow_files``) e logs em
``results/conversor/CT_V1_GCI/logs``. Nomes finais idênticos aos da G0 (3 cell zones, 13 face zones de contorno
+ 3 interiores), para o ``conversor_mrf.py`` rodar sem mudança (só trocar ``MALHA``).
"""
import argparse
import fnmatch
import json
import math
import re
import shutil
import subprocess
import sys
import threading
import time
from pathlib import Path

import ansys.fluent.core as pf

sys.path.insert(0, str(Path(__file__).resolve().parent))

RAIZ = Path(__file__).resolve().parents[1]
MESH = RAIZ / "mesh"
GEOM = RAIZ.parent / "Conversor de torque V1"            # STLs e interfaces
LOGS = RAIZ / "results" / "conversor" / "CT_V1_GCI" / "logs"
TRAB = RAIZ / "results" / "conversor" / "CT_V1_GCI" / "fluent"
PROCESSOS = 2                                            # igual à G0; 12 núcleos lógicos, mas a RAM é o limite

# ---------------------------------------------------------------------------------------------------------------
# G0 (referência), lida do .wft e do C_TG.md
# ---------------------------------------------------------------------------------------------------------------
BASE = {"MinSize": 1.5, "MaxSize": 20.0, "GrowthRate": 1.2, "CurvatureNormalAngle": 18.0, "CellsPerGap": 3.0,
        "GlobalMax": 18.39594116210937, "GlobalMin": 0.071859145164489741,   # GlobalMin = GlobalMax / 256
        "WrapRatio": 2.0 / 3.0}
CELULAS_G0 = 1_900_736
FACES_SUPERFICIE_G0 = 1_440_968                          # soma das face zones da malha de superfície (sem origin-*)
VOLUME_G0_M3 = 8.576847e-2                               # "total volume" da G0 no solver (só referência da G0)
FATORES = {"G1": 1.35, "G2": 1.35 ** 2}                  # f nominal; r efetivo é medido com V e N reais (ver estatisticas)
GROWTH_MAX = 1.34                                        # a semelhança pede 1 + 0,2 f; G2 precisaria de 1,365 com f = 1,8225
CPG_MIN = 1.775                                          # CONFERIR: o campo aceita valor não inteiro (1,775)?

STLS = ["Bomba_com_carcaca.stl", "estator.stl", "interface_bomba_turbina.stl", "interface_estator_bomba.stl",
        "interface_turbina_estator.stl", "turbine.stl"]  # mesma ordem do .wft (define .1 = bomba, .2 = estator, .3 = turbina)
OBJETOS = ["bomba_com_carcaca", "estator", "interface_bomba_turbina", "interface_estator_bomba",
           "interface_turbina_estator", "turbine"]
OBJETOS_SOLIDOS = ["bomba_com_carcaca", "estator", "turbine"]
BASE_MPT = {"bomba": ((-130.0, 164.1, 162.6), 20), "turbina": ((87.0, 282.0, 1.2), 18),
            "estator": ((80.0, 59.0, -19.2), 18)}        # ponto-base (mm) e número de pás/canais (pontos_material.py)
S = "stlb_atf_15.15.0.0_color_____"                        # nome da zona de parede vindo dos STLs (+ .1 bomba, .2 estator, .3 turbina)
GRUPOS_FACE = {                                           # nome final: padrão das face zones a juntar (passo_faces.py da G0)
    "parede_bomba": f"{S}.1:fluido_bomba*", "parede_bomba_lado_turbina": f"{S}.1:fluido_turbina*",
    "parede_bomba_lado_estator": f"{S}.1:fluido_estator*",
    "parede_estator": f"{S}.2:fluido_estator*", "parede_estator_lado_turbina": f"{S}.2:fluido_turbina*",
    "parede_estator_lado_bomba": f"{S}.2:fluido_bomba*",
    "parede_turbina": f"{S}.3:fluido_turbina*",
    "if_bomba_turbina_lado_bomba": "interface_bomba_turbina.1:fluido_bomba*",
    "if_bomba_turbina_lado_turbina": "interface_bomba_turbina.1:fluido_turbina*",
    "if_turbina_estator_lado_turbina": "interface_turbina_estator.1:fluido_turbina*",
    "if_turbina_estator_lado_estator": "interface_turbina_estator.1:fluido_estator*",
    "if_estator_bomba_lado_estator": "interface_estator_bomba.1:fluido_estator*",
    "if_estator_bomba_lado_bomba": "interface_estator_bomba.1:fluido_bomba*",
}
CELL_ZONES = ("fluido_bomba", "fluido_turbina", "fluido_estator")
FACE_ZONES_ESPERADAS = tuple(GRUPOS_FACE)                 # 13 de contorno; os 3 ``interior--fluido_*`` vêm da junção das células
# Zonas pequenas que o conversor_mrf.py lista em PAREDES_ABS: podem sumir numa malha grossa (G0 tinha 62 faces).
ZONAS_PEQUENAS = ("parede_estator_lado_bomba", "parede_estator_lado_turbina", "parede_bomba_lado_estator",
                  "parede_bomba_lado_turbina")     # faces na G0, nesta ordem: 62 / 973 / 830 / 1744 (opcionais: só aviso)
N_REGIOES, N_FACE_ZONES, N_SYMMETRY = 56, 16, 112


def log(*a):
    """Imprime com horário e esvazia o buffer.

    Args:
        *a (object): Itens a imprimir.
    """
    print(time.strftime("[%H:%M:%S]"), *a, flush=True)


def parametros_nivel(fator):
    """Calcula os parâmetros de tamanho de um nível a partir do fator de escala em relação à G0.

    Comprimentos multiplicam por ``fator``; a taxa de crescimento cresce com ``fator`` (para a lei h = h0 + (GR-1)·d
    continuar semelhante); células por folga dividem por ``fator``; o ângulo de curvatura multiplica por ``fator``.

    Args:
        fator (float): Razão de tamanho em relação à G0 (G1 = 1,3; G2 = 1,69).

    Returns:
        dict: ``MinSize``, ``MaxSize``, ``GrowthRate``, ``CurvatureNormalAngle``, ``CellsPerGap`` e os respectivos
        ``Wrap*``, mais ``GlobalMax``, ``GlobalMin`` e ``GlobalGrowthRate``; tamanhos em mm.
    """
    gr = min(1.0 + (BASE["GrowthRate"] - 1.0) * fator, GROWTH_MAX)
    cpg = max(BASE["CellsPerGap"] / fator, CPG_MIN)
    ang = BASE["CurvatureNormalAngle"] * fator
    mn, mx = BASE["MinSize"] * fator, BASE["MaxSize"] * fator
    w = BASE["WrapRatio"]
    return {"fator": fator, "MinSize": mn, "MaxSize": mx, "GrowthRate": gr, "CurvatureNormalAngle": ang,
            "CellsPerGap": cpg, "WrapMin": mn * w, "WrapMax": mx * w, "WrapGrowthRate": gr,
            "WrapCurvatureNormalAngle": ang, "WrapCellsPerGap": cpg,
            "GlobalMax": BASE["GlobalMax"] * fator, "GlobalMin": BASE["GlobalMin"] * fator, "GlobalGrowthRate": gr}


def pontos_material():
    """Gera os 56 pontos de material (um por canal entre pás), girando o ponto-base em torno de X.

    Returns:
        list[tuple[str, float, float, float]]: ``(nome, x, y, z)`` em mm; ``fluido_<peça>`` para k = 0 e
        ``fluido_<peça>_02`` ... para os demais (mesmos nomes e coordenadas de ``pontos_material_CT_V1.csv``).
    """
    out = []
    for peca, ((x, y, z), n) in BASE_MPT.items():
        for k in range(n):
            a = 2 * math.pi * k / n
            c, s = math.cos(a), math.sin(a)
            nome = f"fluido_{peca}" if k == 0 else f"fluido_{peca}_{k + 1:02d}"
            out.append((nome, x, round(c * y - s * z, 2), round(s * y + c * z, 2)))
    return out


def lancar_meshing(processos):
    """Abre o Fluent Meshing 3D, dupla precisão, sem GUI (mesmos argumentos do ``ms.py`` da G0).

    Args:
        processos (int): Número de processos MPI.

    Returns:
        tuple: ``(sessão, pids)``; ``sessão`` é ``ansys.fluent.core.session_meshing.Meshing`` (aberta com
        ``cleanup_on_exit=False``) e ``pids`` é ``(cortex_pid, fluent_host_pid)``, guardados já aqui porque
        ``exit()`` zera ``_fluent_connection``.
    """
    TRAB.mkdir(parents=True, exist_ok=True)
    m = pf.launch_fluent(product_version=pf.FluentVersion.v241, precision=pf.Precision.DOUBLE,
                         processor_count=processos, ui_mode=pf.UIMode.NO_GUI, mode="meshing", cwd=str(TRAB),
                         cleanup_on_exit=False, start_watchdog=False)
    cp = m._fluent_connection.connection_properties
    return m, (cp.cortex_pid, cp.fluent_host_pid)


def pid_vivo(pid):
    """Diz se um PID ainda aparece no ``tasklist``.

    Args:
        pid (int): PID do processo.

    Returns:
        bool: True se o processo existe.
    """
    r = subprocess.run(["tasklist", "/FI", f"PID eq {pid}", "/NH"], capture_output=True, text=True)
    return re.search(rf"\b{pid}\b", r.stdout) is not None


def encerrar_meshing(m, pids):
    """Encerra o Fluent Meshing lançado por este script (o ``exit()`` não termina os processos na 24.1).

    Tenta ``exit()`` numa thread e, se os PIDs guardados na abertura continuarem vivos, mata as árvores deles.
    A conferência é por PID (e não por nome de imagem), porque outra simulação pode estar com ``fluent.exe``,
    ``cx2410.exe`` e ``fl_mpi2410.exe`` abertos; os processos que restam no sistema só são registrados no log.

    Args:
        m (ansys.fluent.core.session_meshing.Meshing): Sessão aberta com ``cleanup_on_exit=False``.
        pids (tuple[int, int]): ``(cortex_pid, fluent_host_pid)`` guardados logo após o launch.

    Returns:
        bool: True só se nenhum dos PIDs desta sessão estiver vivo ao final.
    """
    t = threading.Thread(target=lambda: m.exit(), daemon=True)
    t.start()
    t.join(30)
    for pid in pids:
        if pid_vivo(pid):
            subprocess.run(["taskkill", "/T", "/F", "/PID", str(pid)], capture_output=True)
    time.sleep(3)
    vivos = [p for p in pids if pid_vivo(p)]
    r = subprocess.run(["tasklist"], capture_output=True, text=True).stdout.lower()
    resto = {n: r.count(n) for n in ("fluent.exe", "cx2410.exe", "fl_mpi2410.exe")}
    log(f"PIDs desta sessão vivos: {vivos}; processos Fluent no sistema (podem ser de outra simulação): {resto}")
    return not vivos


def tarefa(m, nome):
    """Devolve a tarefa do workflow pelo nome de exibição.

    Args:
        m (ansys.fluent.core.session_meshing.Meshing): Sessão do Meshing.
        nome (str): Nome da tarefa, p. ex. ``"Add Boundary Layers"`` (ou o nome de um filho, p. ex. ``"curvature-1"``).

    Returns:
        Tarefa do workflow clássico (``m.workflow.TaskObject[nome]``).
    """
    return m.workflow.TaskObject[nome]


def executar(m, nome, args=None, filho=False):
    """Define os argumentos de uma tarefa e a executa (ou cria um filho com ``AddChildAndUpdate``).

    Args:
        m (ansys.fluent.core.session_meshing.Meshing): Sessão do Meshing.
        nome (str): Nome da tarefa.
        args (dict | None): Argumentos (nomes do ``.wft``); ``None`` mantém os atuais.
        filho (bool): True para tarefas compostas (Identify Regions, Add Local Sizing, Extract Edge Features):
            cria um filho em vez de executar a tarefa-mãe.

    Returns:
        bool: Resultado do comando do Fluent.

    Raises:
        RuntimeError: Se o comando devolver falso ou a tarefa ficar com erros.
    """
    t = tarefa(m, nome)
    if args is not None:
        t.Arguments.set_state(args)
    t0 = time.time()
    ok = t.AddChildAndUpdate() if filho else t.Execute()
    log(f"{nome}: ok={ok} estado={t.State()} ({(time.time() - t0) / 60:.1f} min)")
    erros = t.Errors()
    if not ok or erros:
        raise RuntimeError(f"Tarefa '{nome}' falhou: ok={ok}, erros={erros}, avisos={t.Warnings()}")
    return ok


def transcrever(m, arquivo, comando):
    """Roda um comando de listagem da TUI gravando a saída num transcript da TUI e devolve o texto.

    O transcript do cliente (``m.transcript``) é outro mecanismo e fica separado, em ``gerar_malha``.

    Args:
        m (ansys.fluent.core.session_meshing.Meshing): Sessão do Meshing.
        arquivo (Path): Arquivo do transcript.
        comando (callable): Função da TUI sem argumentos, p. ex. ``m.tui.boundary.manage.list``.

    Returns:
        str: Conteúdo do transcript.

    Raises:
        FileNotFoundError: Se o Fluent não gravou o arquivo do transcript.
    """
    arquivo.parent.mkdir(parents=True, exist_ok=True)
    arquivo.unlink(missing_ok=True)
    m.tui.file.start_transcript(f'"{arquivo.as_posix()}"')
    try:
        comando()
    finally:
        m.tui.file.stop_transcript()
    if not arquivo.exists():
        raise FileNotFoundError(f"Transcript não gravado: {arquivo}")
    return arquivo.read_text(encoding="utf-8", errors="replace")


def contar_faces(texto):
    """Lê as face zones listadas por ``/boundary/manage/list``, ignorando as ``origin-*``.

    O regex não exige fim de linha, porque a listagem do volume poly pode ganhar colunas.

    Args:
        texto (str): Transcript da listagem (colunas: id, nome, tipo, count, tri, quad, ...).

    Returns:
        tuple[int, list[tuple[str, str, int]]]: Total de faces e a lista ``(nome, tipo, faces)``.
    """
    total, zonas = 0, []
    for lin in texto.splitlines():
        mt = re.match(r"^\s*(\d+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\d+)", lin)
        if mt and not mt.group(2).startswith("origin"):
            total += int(mt.group(4))
            zonas.append((mt.group(2), mt.group(3), int(mt.group(4))))
    return total, zonas


def importar_geometria(m, nivel):
    """Inicia o Fault-tolerant Meshing, importa os 6 STLs (um objeto por arquivo) e confere os objetos criados.

    Args:
        m (ansys.fluent.core.session_meshing.Meshing): Sessão do Meshing.
        nivel (str): ``"G1"`` ou ``"G2"`` (nome dos logs).

    Raises:
        FileNotFoundError: Se algum STL não existir.
        RuntimeError: Se a importação falhar ou faltar algum dos objetos esperados.
    """
    arqs = [GEOM / n for n in STLS]
    faltam = [str(a) for a in arqs if not a.exists()]
    if faltam:
        raise FileNotFoundError(f"STLs ausentes: {faltam}")
    m.workflow.InitializeWorkflow(WorkflowType="Fault-tolerant Meshing")
    lista = ";".join(a.as_posix() for a in arqs)
    # CONFERIR: InputFileChanged/LoadFiles existem no datamodel 24.1, mas não foram rodados com vários STL separados por ';'
    m.PartManagement.InputFileChanged(FilePath=lista, IgnoreSolidNames=False, PartPerBody=False)
    m.PMFileManagement.FileManager.LoadFiles()
    executar(m, "Import CAD and Part Management", {      # só as chaves gravadas no .wft; o mm vem de GlobalSettings
        "Append": False, "Context": 0, "FMDFileName": lista, "FileLoaded": "yes",
        "ObjectSetting": "DefaultObjectSetting"})
    txt = transcrever(m, LOGS / f"{nivel}_objetos.trn", m.tui.objects.list)  # CONFERIR: nome do comando (/objects/list)
    ausentes = [o for o in OBJETOS if not re.search(rf"\b{o}\b", txt)]
    if ausentes:
        raise RuntimeError(f"Objetos não criados na importação: {ausentes}")
    try:                      # CONFERIR: comando da caixa delimitadora (esperado x de -187,5 a 377,6 mm)
        caixa = transcrever(m, LOGS / f"{nivel}_caixa.trn", m.tui.report.bounding_box)
        log("caixa delimitadora (esperado x de -187,5 a 377,6 mm):", " | ".join(caixa.split("\n")[-6:]))
    except Exception as e:                                 # noqa: BLE001
        log("AVISO: caixa delimitadora não lida:", repr(e)[:150])


def descrever_e_extrair_arestas(m):
    """Executa Describe Geometry and Flow (escoamento interno) e extrai as arestas das 3 peças sólidas.

    Args:
        m (ansys.fluent.core.session_meshing.Meshing): Sessão do Meshing.

    Raises:
        RuntimeError: Se alguma das tarefas falhar.
    """
    executar(m, "Describe Geometry and Flow", {
        "AddEnclosure": "No", "CloseCaps": "No", "FlowType": "Internal flow through the object",
        "DescribeGeometryAndFlowOptions": {"AdvancedOptions": True, "ExtractEdgeFeatures": "Yes"}})
    executar(m, "Extract Edge Features", {
        "ExtractEdgesName": "edge-group-1", "ObjectSelectionList": OBJETOS_SOLIDOS,
        "CompleteObjectSelectionList": OBJETOS_SOLIDOS}, filho=True)


def identificar_regioes(m):
    """Cria os 56 pontos de material (Numerical Inputs), um por canal entre pás.

    Args:
        m (ansys.fluent.core.session_meshing.Meshing): Sessão do Meshing.

    Raises:
        RuntimeError: Se algum ponto não puder ser criado.
    """
    for nome, x, y, z in pontos_material():
        executar(m, "Identify Regions", {
            "GraphicalSelection": False, "MaterialPointsName": nome, "MptMethodType": "Numerical Inputs",
            "ShowCoordinates": True, "X": x, "Y": y, "Z": z}, filho=True)


def listas_regioes(nomes_estado):
    """Monta as listas de Update Region Settings: 6 objetos void e 56 regiões fluid/wrap/poly-hexcore.

    Args:
        nomes_estado (list[str] | None): ``AllRegionNameList`` lida da sessão (preserva a ordem do Fluent); se
            vazia ou ``None``, usa os 6 objetos seguidos dos 56 pontos de material.

    Returns:
        dict: Argumentos de ``UpdateRegionSettings`` com os mesmos campos gravados no ``.wft``.
    """
    flu = [p[0] for p in pontos_material()]
    nomes = list(nomes_estado) if nomes_estado else OBJETOS + flu
    ehflu = [n.startswith("fluido_") for n in nomes]
    nf = [n for n in nomes if n.startswith("fluido_")]
    return {
        "AllRegionNameList": nomes,
        "AllRegionTypeList": ["fluid" if f else "void" for f in ehflu],
        "AllRegionMeshMethodList": ["wrap" if f else "none" for f in ehflu],
        "AllRegionVolumeFillList": ["poly-hexcore" if f else "none" for f in ehflu],
        "AllRegionLeakageSizeList": ["none"] * len(nomes),
        "AllRegionOversetComponenList": ["no"] * len(nomes),
        "AllRegionLinkedConstructionSurfaceList": ["no" if f else "n/a" for f in ehflu],
        "AllRegionSourceList": ["mpt" if f else "object" for f in ehflu],
        "AllRegionFilterCategories": ["1" if f else "2" for f in ehflu],
        "RegionNameList": nf, "RegionTypeList": ["fluid"] * len(nf), "RegionMeshMethodList": ["wrap"] * len(nf),
        "RegionVolumeFillList": ["poly-hexcore"] * len(nf), "RegionLeakageSizeList": [""] * len(nf),
        "RegionOversetComponenList": ["no"] * len(nf),
        "OldRegionNameList": nf, "OldRegionTypeList": ["fluid"] * len(nf), "OldRegionMeshMethodList": ["wrap"] * len(nf),
        "OldRegionVolumeFillList": ["hexcore"] * len(nf), "OldRegionLeakageSizeList": [""] * len(nf),
        "OldRegionOversetComponenList": ["no"] * len(nf)}


def atualizar_regioes(m):
    """Executa Update Region Settings: objetos void; regiões fluid, wrap e poly-hexcore.

    Se o ``Execute`` recusar (p. ex. por causa de Define Leakage Threshold, que foi pulada), tenta
    ``ExecuteUpstreamNonExecutedAndThisTask``.

    Args:
        m (ansys.fluent.core.session_meshing.Meshing): Sessão do Meshing.

    Raises:
        RuntimeError: Se as duas formas de execução falharem.
    """
    t = tarefa(m, "Update Region Settings")
    nomes = None
    try:      # CONFERIR: caminho de dados das regiões (FTMRegionData); se falhar, usa o estado da tarefa e depois o fallback
        nomes = list(m.meshing.GlobalSettings.FTMRegionData.AllRegionNameList())
    except Exception:                                              # noqa: BLE001
        try:
            nomes = t.Arguments.get_state().get("AllRegionNameList")
        except Exception:                                          # noqa: BLE001
            nomes = None
    args = listas_regioes(nomes)
    try:
        executar(m, "Update Region Settings", args)
    except RuntimeError as e:
        log("Execute recusado, tentando ExecuteUpstreamNonExecutedAndThisTask:", str(e)[:200])
        t.Arguments.set_state(args)
        ok = t.ExecuteUpstreamNonExecutedAndThisTask()             # CONFERIR: nome no TaskObject 24.1
        if not ok or t.Errors():
            raise RuntimeError(f"Update Region Settings falhou: ok={ok}, erros={t.Errors()}") from e
    est = t.Arguments.get_state()
    tipos = est.get("AllRegionTypeList", [])
    log(f"regiões: {tipos.count('fluid')} fluid, {tipos.count('void')} void, "
        f"preenchimento {sorted(set(est.get('AllRegionVolumeFillList', [])))}")


def controle_global(m, p):
    """Executa Choose Mesh Control Options (Custom) com crescimento e limites do nível.

    Args:
        m (ansys.fluent.core.session_meshing.Meshing): Sessão do Meshing.
        p (dict): Parâmetros do nível (``parametros_nivel``).

    Raises:
        RuntimeError: Se a tarefa falhar.
    """
    t = tarefa(m, "Choose Mesh Control Options")
    t.Arguments.update_dict({"CreationMethod": "Custom", "GlobalGrowthRate": p["GlobalGrowthRate"],
                             "GlobalMax": p["GlobalMax"], "GlobalMin": p["GlobalMin"]})  # CONFERIR: update_dict em tarefa condicional
    executar(m, "Choose Mesh Control Options")


def configurar_dimensionamento(m, p):
    """Define o controle global (Custom) e os dois dimensionamentos locais (curvatura e proximidade) do nível.

    Args:
        m (ansys.fluent.core.session_meshing.Meshing): Sessão do Meshing.
        p (dict): Parâmetros do nível (``parametros_nivel``).

    Raises:
        RuntimeError: Se alguma tarefa falhar.
    """
    controle_global(m, p)
    for tipo in ("curvature", "proximity"):
        par = {"AdvancedOptions": False, "CellsPerGap": p["CellsPerGap"],
               "CurvatureNormalAngle": p["CurvatureNormalAngle"], "GrowthRate": p["GrowthRate"],
               "IgnoreSelf": True, "InitialSizeControl": False, "MaxSize": p["MaxSize"], "MinSize": p["MinSize"],
               "ScopeProximityTo": "faces-and-edges", "SizingType": tipo, "TargetSizeControl": False,
               "WrapCellsPerGap": p["WrapCellsPerGap"], "WrapCurvatureNormalAngle": p["WrapCurvatureNormalAngle"],
               "WrapGrowthRate": p["WrapGrowthRate"], "WrapMax": p["WrapMax"], "WrapMin": p["WrapMin"]}
        executar(m, "Add Local Sizing", {
            "LocalSettingsName": f"{tipo}-1", "SelectionType": "object", "ObjectSelectionList": OBJETOS,
            "CompleteObjectSelectionList": OBJETOS,
            "ValueChanged": "NotDefault", "LocalSizeControlParameters": par}, filho=True)


def refazer_dimensionamento(m, p):
    """Reaplica controle global e dimensionamento num workflow restaurado de uma malha de superfície.

    Args:
        m (ansys.fluent.core.session_meshing.Meshing): Sessão do Meshing com o workflow restaurado.
        p (dict): Parâmetros do nível (``parametros_nivel``).

    Raises:
        RuntimeError: Se um filho de dimensionamento ou o controle global não puder ser reexecutado.
    """
    controle_global(m, p)                             # reescala o global antes dos filhos
    for tipo in ("curvature", "proximity"):
        nome = f"{tipo}-1"
        t = tarefa(m, nome)
        est = t.Arguments.get_state()
        par = dict(est.get("LocalSizeControlParameters", {}))
        par.update({"CellsPerGap": p["CellsPerGap"], "CurvatureNormalAngle": p["CurvatureNormalAngle"],
                    "GrowthRate": p["GrowthRate"], "MaxSize": p["MaxSize"], "MinSize": p["MinSize"],
                    "WrapCellsPerGap": p["WrapCellsPerGap"], "WrapCurvatureNormalAngle": p["WrapCurvatureNormalAngle"],
                    "WrapGrowthRate": p["WrapGrowthRate"], "WrapMax": p["WrapMax"], "WrapMin": p["WrapMin"]})
        t.Arguments.update_dict({"LocalSizeControlParameters": par})   # update_dict preserva ObjectSelectionList etc.
        ok = t.Execute()                          # CONFERIR: o GUI usa "Update" no nó; Execute() do filho é equivalente?
        log(f"{nome} reexecutado: ok={ok} estado={t.State()}")
        if not ok:
            raise RuntimeError(f"Não consegui reexecutar {nome}")


def interfaces_como_symmetry(m):
    """Troca as zonas ``interface*`` para symmetry (o ``conversor_mrf.py`` lê ``bc.symmetry`` e converte em interface).

    Args:
        m (ansys.fluent.core.session_meshing.Meshing): Sessão do Meshing.
    """
    # CONFERIR: na G0 foi digitado no console do GUI como "/boundary/manage/type (interface*) symmetry" (112 zonas)
    m.tui.boundary.manage.type("(interface*)", "symmetry")


def gerar_superficie(m, nivel, p):
    """Gera a malha de superfície, troca as interfaces para symmetry, grava o backup e confere as zonas.

    O backup ``mesh/CT_V1_3zonas_<nivel>_superficie.msh.h5`` é gravado antes de contar faces, já com as 112 zonas
    de interface como symmetry (a retomada com ``--retomar-superficie`` parte dele).

    Args:
        m (ansys.fluent.core.session_meshing.Meshing): Sessão do Meshing.
        nivel (str): ``"G1"`` ou ``"G2"``.
        p (dict): Parâmetros do nível.

    Returns:
        int: Total de faces das zonas de superfície (G0: 1 440 968).

    Raises:
        RuntimeError: Se uma tarefa falhar, não houver 112 zonas symmetry ou as regiões não forem 56 (dois pontos de
            material na mesma região).
    """
    padroes = tarefa(m, "Generate the Surface Mesh").Arguments.get_state()
    log("padrões do Generate the Surface Mesh (G0 usou os padrões; skewness alvo 0,8):",
        json.dumps(padroes, default=str)[:600])
    executar(m, "Generate the Surface Mesh")
    executar(m, "Update Boundaries")              # CONFERIR: Execute() sem argumentos; na G0 foi pela GUI (estado Up-to-date)
    interfaces_como_symmetry(m)
    arq = MESH / f"CT_V1_3zonas_{nivel}_superficie.msh.h5"
    m.tui.file.write_mesh(f'"{arq.as_posix()}"')
    total, zonas = contar_faces(transcrever(m, LOGS / f"{nivel}_zonas_superficie.trn", m.tui.boundary.manage.list))
    sym = [z for z in zonas if z[1] == "symmetry"]
    regioes = {z[0].split(":", 1)[1] for z in zonas if ":" in z[0] and z[0].split(":", 1)[1].startswith("fluido_")}
    prev = FACES_SUPERFICIE_G0 / p["fator"] ** 2
    log(f"superfície {nivel}: {total} faces em {len(zonas)} zonas ({len(sym)} symmetry, {len(regioes)} regiões); "
        f"esperado ~{prev:,.0f} (r_superfície efetivo = {math.sqrt(FACES_SUPERFICIE_G0 / max(total, 1)):.3f}, "
        f"alvo {p['fator']:.3f})")
    if len(sym) != N_SYMMETRY:
        raise RuntimeError(f"Esperava {N_SYMMETRY} zonas symmetry, li {len(sym)}")
    if len(regioes) != N_REGIOES:
        raise RuntimeError(f"Esperava {N_REGIOES} regiões, li {len(regioes)}: dois pontos de material na mesma região?")
    return total


def gerar_volume(m):
    """Desliga os prismas e gera o volume poly-hexcore (sem mexer nos argumentos do volume, como na G0).

    Args:
        m (ansys.fluent.core.session_meshing.Meshing): Sessão do Meshing.

    Raises:
        RuntimeError: Se alguma tarefa falhar.
    """
    executar(m, "Add Boundary Layers", {"AddChild": "no"})   # sem prismas (as folgas não comportam a 1ª camada)
    est = tarefa(m, "Generate the Volume Mesh").Arguments.get_state()
    log("padrões do volume (OrthogonalQuality etc.):", json.dumps({k: v for k, v in est.items()
        if not k.startswith(("All", "Old", "Region"))}, default=str)[:400])
    executar(m, "Generate the Volume Mesh")


def gravar_malha(m, nome):
    """Grava uma malha em ``mesh/`` (junto com a pasta ``_workflow_files``).

    Args:
        m (ansys.fluent.core.session_meshing.Meshing): Sessão do Meshing.
        nome (str): Nome do arquivo, p. ex. ``"CT_V1_3zonas_G1.msh.h5"``.

    Returns:
        Path: Caminho gravado.
    """
    arq = MESH / nome
    m.tui.file.write_mesh(f'"{arq.as_posix()}"')
    log("gravado", arq.name)
    return arq


def juntar_zonas(m, nivel):
    """Remove os objetos, junta as 56 cell zones em 3 e as face zones em 13 (sequência validada na G0).

    Grupos de ``GRUPOS_FACE`` sem zona correspondente na listagem são pulados com AVISO (zonas pequenas podem
    não existir numa malha grossa).

    Args:
        m (ansys.fluent.core.session_meshing.Meshing): Sessão do Meshing.
        nivel (str): ``"G1"`` ou ``"G2"`` (nome dos logs).

    Returns:
        list[str]: Nomes finais dos grupos que foram pulados.
    """
    m.tui.objects.delete_all()                     # remove só os objetos; a malha fica (senão o merge de cell zones recusa)
    for peca in ("bomba", "turbina", "estator"):
        m.tui.mesh.manage.merge(f"fluido_{peca}*", "()")
    bm = m.tui.boundary.manage
    _, zonas = contar_faces(transcrever(m, LOGS / f"{nivel}_face_zones_antes_merge.trn", bm.list))
    nomes = [z[0] for z in zonas]
    pulados = []
    for novo, padrao in GRUPOS_FACE.items():
        if not any(fnmatch.fnmatchcase(n, padrao) for n in nomes):
            log(f"AVISO: nenhuma face zone casa com '{padrao}'; '{novo}' não será criada")
            pulados.append(novo)
            continue
        bm.merge(padrao, "()")
        bm.name(padrao, novo)                      # CONFERIR: a G0 renomeou com o próprio padrão como zona atual
        log("face zone:", novo)
    bm.delete("(origin-*)", "yes")                 # forma validada na G0
    return pulados


def validar_zonas(m, nivel):
    """Confere que os nomes finais são os mesmos da G0 e registra a listagem das zonas.

    As quatro zonas pequenas (``ZONAS_PEQUENAS``) são opcionais: só geram aviso.

    Args:
        m (ansys.fluent.core.session_meshing.Meshing): Sessão do Meshing.
        nivel (str): ``"G1"`` ou ``"G2"``.

    Returns:
        dict: ``n_face_zones``, ``faces_por_zona`` e ``texto_cell_zones`` (listagem crua das cell zones).

    Raises:
        RuntimeError: Se faltar cell zone ou face zone obrigatória, ou sobrar zona com nome de origem.
    """
    cel = transcrever(m, LOGS / f"{nivel}_cell_zones.trn", m.tui.mesh.manage.list)
    fac = transcrever(m, LOGS / f"{nivel}_face_zones.trn", m.tui.boundary.manage.list)
    _, zonas = contar_faces(fac)
    faces = {z[0]: z[2] for z in zonas}
    obrig = [n for n in FACE_ZONES_ESPERADAS if n not in ZONAS_PEQUENAS]
    faltam = [n for n in CELL_ZONES if not re.search(rf"\b{n}\b", cel)]
    faltam += [n for n in obrig if n not in faces]
    sobras = sorted(set(re.findall(r"\b(?:stlb_\S+|interface_\S+|fluido_\S+_\d\d\S*|origin-\S+)", cel + fac)))
    if faltam or sobras:
        raise RuntimeError(f"Zonas de {nivel} diferentes da G0: faltam {faltam}; sobram {sobras[:6]}")
    for n in ZONAS_PEQUENAS:
        if n not in faces:
            log(f"AVISO: zona pequena '{n}' ausente em {nivel}: o conversor_mrf.py (PAREDES_ABS) vai precisar tratá-la")
    log(f"{nivel}: {len(zonas)} face zones (esperado {N_FACE_ZONES}); zonas pequenas "
        f"(G0: 62/973/830/1744): {[faces.get(n, 0) for n in ZONAS_PEQUENAS]}")
    return {"n_face_zones": len(zonas), "faces_por_zona": faces, "texto_cell_zones": cel}


def estatisticas(m, nivel, p, val):
    """Mede células, volume e tamanho médio reais da malha final e compara com a G0.

    Args:
        m (ansys.fluent.core.session_meshing.Meshing): Sessão do Meshing.
        nivel (str): ``"G1"`` ou ``"G2"``.
        p (dict): Parâmetros do nível.
        val (dict): Saída de ``validar_zonas``.

    Returns:
        dict: ``celulas_por_zona``, ``celulas``, ``volume_m3``, ``h_mm``, ``r_efetivo`` (valores ``None`` se a
        leitura falhar; a listagem crua fica nos transcripts de ``LOGS``).
    """
    cel = {}
    for lin in val["texto_cell_zones"].splitlines():
        mt = re.match(r"^\s*(\d+)\s+(\S+)\s+(.*)", lin)
        if mt and mt.group(2) in CELL_ZONES:
            ints = re.findall(r"\d+", mt.group(3))
            cel[mt.group(2)] = int(ints[0]) if ints else None   # CONFERIR: coluna de contagem da /mesh/manage/list
    n = sum(v for v in cel.values() if v) if len(cel) == 3 and all(cel.values()) else None
    vol = None
    try:                                           # CONFERIR: comando que imprime "total volume" (/mesh/check-mesh?)
        txt = transcrever(m, LOGS / f"{nivel}_check.trn", m.tui.mesh.check_mesh)
        mt = re.search(r"total volume \((\w+)\):\s*([\d.eE+-]+)", txt)
        if mt:
            vol = float(mt.group(2)) * (1e-9 if mt.group(1) == "mm3" else 1.0)
    except Exception as e:                         # noqa: BLE001
        log("AVISO: volume total não lido:", repr(e)[:150])
    h = 1e3 * (vol / n) ** (1 / 3) if vol and n else None
    h0 = 1e3 * (VOLUME_G0_M3 / CELULAS_G0) ** (1 / 3)
    r = h / h0 if h else None
    log(f"{nivel}: células {n} {cel}; volume {vol} m3 (G0 {VOLUME_G0_M3}); h = {h} mm (G0 {h0:.2f}); "
        f"r efetivo = {r} (nominal {p['fator']:.3f})")
    if vol and abs(vol / VOLUME_G0_M3 - 1) > 0.02:
        log("AVISO: volume difere mais de 2% da G0 (wrap pode ter engolido folgas)")
    return {"celulas_por_zona": cel, "celulas": n, "volume_m3": vol, "h_mm": h, "r_efetivo": r}


def gerar_malha(nivel, fator, so_superficie=False, origem="stl", processos=PROCESSOS):
    """Gera uma malha completa (ou só a de superfície) para um nível do GCI.

    Args:
        nivel (str): ``"G1"`` ou ``"G2"``.
        fator (float): Fator de escala em relação à G0.
        so_superficie (bool): Para depois da malha de superfície (calibração do fator pela contagem de faces).
        origem (str): ``"stl"`` (importa os STLs), ``"g0"`` (copia e relê a superfície da G0 para ``TRAB``) ou
            ``"retomar"`` (relê ``CT_V1_3zonas_<nivel>_superficie.msh.h5`` e só gera o volume).
        processos (int): Processos MPI do Meshing.

    Returns:
        dict: Resumo gravado em ``LOGS/<nivel>_parametros.json`` (parâmetros, regiões, face zones, células, faces, volume).

    Raises:
        FileNotFoundError: Se a malha de superfície a reler não existir.
        RuntimeError: Se alguma tarefa ou conferência falhar.
    """
    LOGS.mkdir(parents=True, exist_ok=True)
    p = parametros_nivel(fator)
    log(f"== {nivel}: fator {fator:.3f} origem {origem}", json.dumps(p))
    t0 = time.time()
    m, pids = lancar_meshing(processos)
    resumo = {"nivel": nivel, "parametros": p, "origem": origem, "n_regioes": N_REGIOES}
    try:
        m.transcript.stop()                        # CONFERIR: transcript do cliente; a TUI usa o dela em transcrever()
        m.transcript.start(file_name=str(LOGS / f"{nivel}_sessao.log"), write_to_stdout=True)
        if origem == "stl":
            importar_geometria(m, nivel)
            descrever_e_extrair_arestas(m)
            identificar_regioes(m)
            atualizar_regioes(m)
            configurar_dimensionamento(m, p)
        else:
            nome = "CT_V1_3zonas_superficie" if origem == "g0" else f"CT_V1_3zonas_{nivel}_superficie"
            if origem == "g0":                     # copia para TRAB para não sobrescrever a G0
                TRAB.mkdir(parents=True, exist_ok=True)
                shutil.copy2(MESH / f"{nome}.msh.h5", TRAB / f"{nome}.msh.h5")
                shutil.copytree(MESH / f"{nome}_workflow_files", TRAB / f"{nome}_workflow_files", dirs_exist_ok=True)
                base = TRAB / f"{nome}.msh.h5"
            else:
                base = MESH / f"{nome}.msh.h5"
            if not base.exists():
                raise FileNotFoundError(base)
            m.tui.file.read_mesh(f'"{base.as_posix()}"')   # restaura o workflow (verificado na G0); _workflow_files ao lado
            if origem == "g0":
                refazer_dimensionamento(m, p)
        if origem != "retomar":
            resumo["faces_superficie"] = gerar_superficie(m, nivel, p)
        if not so_superficie:
            gerar_volume(m)
            gravar_malha(m, f"CT_V1_3zonas_{nivel}_56regioes.msh.h5")
            resumo["face_zones_puladas"] = juntar_zonas(m, nivel)
            m.tui.report.cell_quality_limits()
            gravar_malha(m, f"CT_V1_3zonas_{nivel}.msh.h5")
            val = validar_zonas(m, nivel)
            resumo["n_face_zones"] = val["n_face_zones"]
            resumo["faces_por_zona"] = val["faces_por_zona"]
            resumo["estatisticas"] = estatisticas(m, nivel, p, val)
            log(f"alvo ~{CELULAS_G0 / fator ** 3:,.0f} células (G0: {CELULAS_G0:,}). Conferir a contagem real no solver.")
    finally:
        resumo["minutos"] = round((time.time() - t0) / 60, 1)
        (LOGS / f"{nivel}_parametros.json").write_text(json.dumps(resumo, indent=1, default=str), encoding="utf-8")
        resumo["encerrado"] = encerrar_meshing(m, pids)
        log("Fluent encerrado:", resumo["encerrado"])
    return resumo


def main(argv=None):
    """Interpreta a linha de comando e gera os níveis pedidos, um de cada vez (um Fluent por vez, para limitar a RAM).

    Args:
        argv (list[str] | None): Argumentos; ``None`` usa ``sys.argv``.

    Returns:
        int: Código de saída (0 se todos os níveis terminarem).
    """
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("niveis", nargs="+", choices=sorted(FATORES))
    ap.add_argument("--so-superficie", action="store_true", help="para após a malha de superfície")
    g = ap.add_mutually_exclusive_group()
    g.add_argument("--partir-da-superficie", action="store_true",
                   help="copia e relê CT_V1_3zonas_superficie.msh.h5 (G0) e refaz a partir do dimensionamento")
    g.add_argument("--retomar-superficie", action="store_true",
                   help="relê CT_V1_3zonas_<nivel>_superficie.msh.h5 e só gera o volume")
    ap.add_argument("--fator", type=float, help="substitui o fator (só com um nível)")
    ap.add_argument("--processos", type=int, default=PROCESSOS)
    a = ap.parse_args(argv)
    if a.fator and len(a.niveis) != 1:
        ap.error("--fator só vale com um nível")
    origem = "g0" if a.partir_da_superficie else "retomar" if a.retomar_superficie else "stl"
    for nv in a.niveis:
        gerar_malha(nv, a.fator or FATORES[nv], a.so_superficie, origem, a.processos)
    return 0


if __name__ == "__main__":
    sys.exit(main())
