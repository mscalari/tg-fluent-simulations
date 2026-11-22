"""Mapas de campo das variantes de sensibilidade do conversor CT_V1 (MRF), com a mesma escala de cores por campo.

Relê, numa única sessão do Fluent, o ``.cas/.dat`` da configuração base (G0, SR 0), das variantes de
sensibilidade (SST sem correção de curvatura, k-ε Realizable, FOU e meio passo) e do par SIMPLE × Coupled da G2,
todos guardados em ``D:/TG_backup_fluent``. Nada é iterado: só pós-processamento.

Em duas passadas:

1. lê cada caso, obtém o mínimo e o máximo de cada campo nas superfícies (``field_data``) e grava os mapas com
   escala automática (``*_auto.png``, só para conferência);
2. com a faixa comum a todos os casos (``faixas.json``), relê cada caso e grava os mapas com a escala fixa
   (``<caso>_<campo>.png``), que são os usados no TG.

Campos: velocidade tangencial e razão de viscosidade turbulenta no plano meridional (z = 0) e velocidade axial no
lado da bomba da interface bomba-turbina, nas mesmas vistas das figuras da varredura.

Uso (rodar com ``PYTHONIOENCODING=utf-8`` e o ``.venv241``)::

    python tools/campos_sensibilidade_conversor.py                 as duas passadas
    python tools/campos_sensibilidade_conversor.py --so-segunda    só a 2ª, com os extremos de ``faixas.json``

Saída: ``results/conversor/CT_V1_MRF/sens_campos/`` (PNGs e ``faixas.json``). As figuras do TG são montadas por
``tools/figuras_tg_campos_conversor.py`` (nomes ``sens_*``).
"""
import json
import sys
import traceback
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import conversor_mrf as cm                 # noqa: E402
from run_caso import encerrar              # noqa: E402

BACKUP = Path("D:/TG_backup_fluent")
CASOS = {  # rótulo: arquivo .cas.h5 (o .dat.h5 de mesmo nome é lido junto)
    "base": BACKUP / "CT_V1_MRF/SR000/fluent/CT_V1_MRF_SR000.cas.h5",
    "sst_sem_cc": BACKUP / "CT_V1_MRF/sens/sst_sem_cc/fluent/CT_V1_MRF_sens_sst_sem_cc.cas.h5",
    "ke_realizable": BACKUP / "CT_V1_MRF/sens/ke_realizable/fluent/CT_V1_MRF_sens_ke_realizable.cas.h5",
    "fou": BACKUP / "CT_V1_MRF/sens/fou/fluent/CT_V1_MRF_sens_fou.cas.h5",
    "meio_passo": BACKUP / "CT_V1_MRF/sens/meio_passo/fluent/CT_V1_MRF_sens_meio_passo.cas.h5",
    "g2_simple": BACKUP / "CT_V1_GCI_G2/SR000/fluent/CT_V1_MRF_SR000.cas.h5",
    "g2_coupled": BACKUP / "CT_V1_GCI_G2/sens/coupled/fluent/CT_V1_MRF_sens_coupled.cas.h5",
}
CAMPOS = {  # rótulo: (campo do Fluent, superfície, vista, faixa de reserva se o field_data falhar)
    "vtheta": ("tangential-velocity", "plano_meridional", "meridional", (-45.0, 100.0)),
    "mut": ("viscosity-ratio", "plano_meridional", "meridional", (0.0, 3000.0)),
    "ifbt_vx": ("axial-velocity", "if_bomba_turbina_lado_bomba", "de_menos_x", (-60.0, 60.0)),
}
SAIDA = cm.DIR / "sens_campos"
# Mapas usados no TG (2ª passada): campo -> casos que dividem a mesma escala. A faixa de cada campo sai só dos
# casos do seu grupo, para que um extremo de outra variante (o k-ε chega a 160 m/s na if_bt) não achate a escala.
GRUPOS = {"vtheta": ("base", "fou", "g2_simple", "g2_coupled"),
          "mut": ("base", "sst_sem_cc", "ke_realizable", "fou"),
          "ifbt_vx": ("base", "meio_passo")}
# Faixa imposta: no k-ε o máximo de μt/μ (~1,2e4) é quatro vezes o do SST e apagaria o campo do SST; a escala
# fica na do SST e o k-ε satura (declarado na legenda da figura).
FAIXA_FIXA = {"mut": [0.0, 3000.0]}


def extremos(s, campo, superficie):
    """Mínimo e máximo de um campo (valores nos nós) numa superfície, pelo ``field_data`` da sessão.

    Args:
        s: Sessão do Fluent (``Solver``).
        campo (str): Nome do campo.
        superficie (str): Nome da superfície.

    Returns:
        tuple[float, float] | None: (mínimo, máximo), ou None se a leitura falhar.
    """
    try:
        dados = s.fields.field_data.get_scalar_field_data(field_name=campo, surfaces=[superficie],
                                                          node_value=True)
    except Exception as e:
        cm.log(f"field_data {campo}/{superficie} falhou:", repr(e)[:200])
        return None

    def valores(obj):
        """Junta recursivamente os números de um retorno do ``field_data`` (dict, lista ou objeto com ``data``).

        Args:
            obj: Retorno (ou parte dele).

        Returns:
            list[np.ndarray]: Vetores de valores encontrados.
        """
        if isinstance(obj, np.ndarray):
            return [obj.ravel().astype(float)]
        if isinstance(obj, dict):
            return [v for x in obj.values() for v in valores(x)]
        if isinstance(obj, (list, tuple)):
            if obj and all(isinstance(x, (int, float)) for x in obj):
                return [np.asarray(obj, dtype=float)]
            return [v for x in obj for v in valores(x)]
        for atributo in ("data", "scalar_data"):
            if hasattr(obj, atributo):
                return valores(getattr(obj, atributo))
        return []

    v = [x for x in valores(dados) if x.size]
    if not v:
        cm.log(f"field_data {campo}/{superficie}: retorno sem valores ({type(dados)})")
        return None
    v = np.concatenate(v)
    return float(np.nanmin(v)), float(np.nanmax(v))


def ler_caso(s, rotulo):
    """Lê o ``.cas/.dat`` de um caso, cria as superfícies de pós-processamento e prepara a janela gráfica.

    Args:
        s: Sessão do Fluent.
        rotulo (str): Chave de ``CASOS``.

    Returns:
        bool: True se o caso foi lido e a imagem de teste gravada.
    """
    st = s.settings
    st.file.read_case_data(file_name=str(CASOS[rotulo]))
    cm.sincronizar_paredes_abs(st)
    cm.superficies_pos(st)
    return cm.teste_imagem(st, SAIDA)


def passada(s, faixas=None):
    """Lê todos os casos e grava os mapas; sem ``faixas``, mede os extremos e usa escala automática.

    Args:
        s: Sessão do Fluent.
        faixas (dict[str, list[float]] | None): Faixa fixa por campo; None na 1ª passada.

    Returns:
        dict[str, dict[str, list[float]]]: Extremos medidos por caso e campo (vazio na 2ª passada).
    """
    st = s.settings
    medidos = {}
    for rot in CASOS:
        if faixas is not None and not any(rot in g for g in GRUPOS.values()):
            continue
        try:
            cm.log(f"lendo {rot}: {CASOS[rot]}")
            cm.log(f"{rot}: teste de imagem", ler_caso(s, rot))
            for nome, (campo, sup, vista, _) in CAMPOS.items():
                if faixas is not None and rot not in GRUPOS[nome]:
                    continue
                if faixas is None:
                    e = extremos(s, campo, sup)
                    if e:
                        medidos.setdefault(rot, {})[nome] = list(e)
                    arq, faixa = SAIDA / f"{rot}_{nome}_auto.png", None
                else:
                    arq, faixa = SAIDA / f"{rot}_{nome}.png", tuple(faixas[nome])
                ok = cm.contorno(st, f"c_{nome}", campo, [sup], vista, arq, faixa)
                cm.log(f"{rot} {nome}: imagem {ok}", medidos.get(rot, {}).get(nome, ""))
        except Exception:
            traceback.print_exc()
            cm.log(f"{rot}: falhou")
    return medidos


def faixa_comum(medidos):
    """Faixa comum por campo: do menor mínimo ao maior máximo entre os casos do grupo, arredondada para fora.

    Args:
        medidos (dict[str, dict[str, list[float]]]): Extremos por caso e campo.

    Returns:
        dict[str, list[float]]: Faixa por campo (a de reserva de ``CAMPOS`` se nenhum caso foi medido).
    """
    out = {}
    for nome, (_, _, _, reserva) in CAMPOS.items():
        if nome in FAIXA_FIXA:
            out[nome] = list(FAIXA_FIXA[nome])
            continue
        ex = [m[nome] for c, m in medidos.items() if nome in m and c in GRUPOS[nome]]
        if not ex:
            out[nome] = list(reserva)
            continue
        lo, hi = min(e[0] for e in ex), max(e[1] for e in ex)
        passo = 10.0 ** np.floor(np.log10(max(abs(lo), abs(hi)))) / 2
        out[nome] = [float(np.floor(lo / passo) * passo), float(np.ceil(hi / passo) * passo)]
        if nome == "mut":
            out[nome][0] = 0.0
    return out


def main(so_segunda=False):
    """Executa as duas passadas e grava ``faixas.json`` com os extremos medidos e as faixas usadas.

    Args:
        so_segunda (bool): Pula a 1ª passada e usa os extremos já medidos em ``faixas.json``.
    """
    SAIDA.mkdir(parents=True, exist_ok=True)
    s = cm.abrir()
    try:
        cm.preparar_picture(s.settings)
        if so_segunda:
            medidos = json.loads((SAIDA / "faixas.json").read_text())["medidos"]
        else:
            medidos = passada(s)
        faixas = faixa_comum(medidos)
        (SAIDA / "faixas.json").write_text(json.dumps({"medidos": medidos, "faixas": faixas}, indent=1))
        cm.log("faixas comuns:", json.dumps(faixas))
        passada(s, faixas)
        cm.log("FIM campos_sensibilidade")
    finally:
        encerrar(s)


if __name__ == "__main__":
    main("--so-segunda" in sys.argv[1:])
