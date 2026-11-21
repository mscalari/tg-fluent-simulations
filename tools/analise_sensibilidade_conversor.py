"""Sensibilidade do ponto de estol (SR = 0) do conversor CT_V1 ao setup do MRF.

Compara a base (``results/conversor/CT_V1_MRF/SR000``) com as variantes gravadas em
``<dir-sens>/<variante>/dados/resultado.json``: ``sst_sem_cc`` (SST sem correção de curvatura),
``ke_realizable`` (k-ε Realizable), ``fou`` (esquemas de 1ª ordem) e ``meio_passo`` (turbina girada
de meio passo de pá, 10°). Variante que ainda não existe é avisada e pulada; o script roda com o que houver.

Uso:
    python tools/analise_sensibilidade_conversor.py
    python tools/analise_sensibilidade_conversor.py --dir-sens <pasta> --saida <pasta>

Opções:
    --dir-sens: pasta com as subpastas de variante (padrão: results/conversor/CT_V1_MRF/sens).
    --saida: pasta onde gravar ``dados/sensibilidade.csv`` e ``graficos/g6_sensibilidade.png``
        (padrão: results/conversor/CT_V1_MRF).

Grandezas: TR = -T_t/T_b; K = n_b/sqrt(|T_b|) com n_b em rpm (``rpm_b``, padrão 2000);
balanço = (T_b+T_t+T_e)/|T_b|·100. Torques gravados em Cm são convertidos para N·m com
``fator_torque`` de ``analise_conversor``.

Incerteza (por ponto, para T_b, T_t, TR e K):
    * SEM (erro-padrão da média) pela série por iteração da janela, como em ``gci_conversor``:
      ``serie_janela`` acha o Report File ``.out`` cuja janela bate com ``media_janela`` e ``grandezas_ponto``
      calcula TR e K iteração a iteração e o SEM por batch means (5 lotes). O τ integral sai junto
      (colunas ``tau_*``, como conferência). Se o ``.out`` não for achado, cai no σ por iteração
      (``desvio_janela``, método antigo, que superestima a incerteza) com aviso; ``metodo_incerteza`` marca isso.
    * Diferença entre variante e base: U = √(SEM_variante² + SEM_base²) e ``abs_d_sobre_U_*`` = |Δ|/U.
      Valor ≥ 2 indica diferença acima do ruído iterativo (cerca de 2 desvios-padrão).
      As barras de erro do gráfico são ±2U.
    * As colunas ``sigma_*``, ``abs_dT_b_sigma`` e ``abs_dTR_sigma`` (σ por iteração, sem correlação) são
      mantidas só para referência.

Os ``.out`` das variantes são achados registrando a pasta de cada variante em ``PASTAS`` de ``gci_conversor``
(chave ``VAR_<nome>``), pois ``serie_janela`` procura os arquivos pelo nome do nível. O import de
``gci_conversor`` altera ``rcParams`` do matplotlib; essa mudança é desfeita logo após o import.
"""
import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from analise_conversor import fator_torque  # noqa: E402

_RC_ANTES = matplotlib.rcParams.copy()
from gci_conversor import N_LOTES, PASTAS, grandezas_ponto, serie_janela  # noqa: E402
matplotlib.rcParams.update(_RC_ANTES)

ROOT = Path(__file__).resolve().parents[1]
DIR = ROOT / "results" / "conversor" / "CT_V1_MRF"
BASE = DIR / "SR000" / "dados" / "resultado.json"
VARIANTES = ["sst_sem_cc", "ke_realizable", "fou", "meio_passo"]
RPM_B_PADRAO = 2000.0
NOME_BASE = "base (SR000)"
DESCRICAO_BASE = "Base: SST com correção de curvatura, PRESTO!, SIMPLE, SOU"
CHAVE_BASE = "G0"  # chave de gci_conversor.PASTAS da base (malha fina, SR000 em CT_V1_MRF)
QUANTIDADES = ("T_b", "T_t", "TR", "K")

# Colunas do console: (chave na linha, título curto, casas decimais). Títulos em ASCII para o console cp1252.
COLUNAS_CONSOLE = [("variante", "variante", 0), ("convergiu", "conv", 0), ("iteracoes", "it", 0),
                   ("T_b [N·m]", "T_b", 1), ("T_t [N·m]", "T_t", 1), ("TR", "TR", 4), ("K [rpm/√(N·m)]", "K", 2),
                   ("dT_b_pct", "dTb%", 2), ("dT_t_pct", "dTt%", 2), ("dTR_pct", "dTR%", 2), ("dK_pct", "dK%", 2),
                   ("abs_d_sobre_U_T_b", "|dTb|/U", 2), ("abs_d_sobre_U_T_t", "|dTt|/U", 2),
                   ("abs_d_sobre_U_TR", "|dTR|/U", 2), ("abs_d_sobre_U_K", "|dK|/U", 2)]
COLUNAS_CONSOLE_SEM = [("variante", "variante", 0), ("sem_T_b", "SEM_Tb", 2), ("sem_T_t", "SEM_Tt", 2),
                       ("sem_TR", "SEM_TR", 5), ("sem_K", "SEM_K", 3), ("tau_T_b", "tau_Tb", 1),
                       ("tau_T_t", "tau_Tt", 1), ("tau_TR", "tau_TR", 1), ("tau_K", "tau_K", 1)]


def ler_resultado(arq):
    """Lê o ``resultado.json`` de um ponto.

    Args:
        arq (Path): Arquivo JSON do ponto.

    Returns:
        dict: Conteúdo do arquivo.
    """
    return json.loads(arq.read_text(encoding="utf-8"))


def registrar_variante(nome, pasta):
    """Registra a pasta de uma variante em ``PASTAS`` (de ``gci_conversor``) para ``serie_janela`` achar o ``.out``.

    Efeito colateral: acrescenta a chave ``VAR_<nome>`` ao dicionário ``PASTAS`` do módulo importado.

    Args:
        nome (str): Nome da variante.
        pasta (Path): Pasta da variante (contém ``dados/resultado.json`` e ``dados/*.out``).

    Returns:
        str: Chave registrada em ``PASTAS`` (``VAR_<nome>``).
    """
    chave = f"VAR_{nome}"
    PASTAS[chave] = pasta
    return chave


def incerteza_ponto(chave, r):
    """SEM e τ de T_b, T_t, TR e K de um ponto, pela série por iteração da janela, com fallback.

    Usa ``serie_janela`` (acha o ``.out`` da janela em ``PASTAS[chave]``) e ``grandezas_ponto`` (TR e K iteração
    a iteração; SEM por batch means). Sem o ``.out``, ``grandezas_ponto`` devolve o σ por iteração (método antigo).

    Args:
        chave (str): Chave de ``PASTAS`` onde procurar o ``.out`` (``CHAVE_BASE`` ou ``VAR_<nome>``).
        r (dict): ``resultado.json`` do ponto. ``rpm_b`` ausente ou nulo é preenchido com ``RPM_B_PADRAO``, como em
            ``metricas``.

    Returns:
        dict: ``sem`` e ``tau`` (dicionários por grandeza de ``QUANTIDADES``; τ em iterações, NaN no fallback),
        ``metodo`` (texto para ``metodo_incerteza``), ``achou`` (bool) e ``info`` (nome do ``.out`` usado ou
        motivo de não ter achado).

    Raises:
        KeyError: Se faltar ``media_janela`` ou ``sr`` no ``resultado.json``.
    """
    r2 = dict(r)
    r2["rpm_b"] = r.get("rpm_b") or RPM_B_PADRAO
    serie, info = serie_janela(chave, r2, float(r2["sr"]))
    g = grandezas_ponto(r2, serie)
    achou = serie is not None
    metodo = f"batch means ({N_LOTES} lotes, serie .out)" if achou else "sigma por iteracao (fallback)"
    return {"sem": {q: g[q]["sem"] for q in QUANTIDADES}, "tau": {q: g[q]["tau"] for q in QUANTIDADES},
            "metodo": metodo, "achou": achou, "info": info}


def metricas(nome, r, descricao_padrao, inc):
    """Grandezas de desempenho de um ponto, já convertidas para N·m, com as incertezas.

    Args:
        nome (str): Nome da variante (ou da base).
        r (dict): ``resultado.json`` do ponto.
        descricao_padrao (str): Descrição usada quando o JSON não traz ``setup.descricao``.
        inc (dict): Saída de ``incerteza_ponto`` para o mesmo ponto.

    Returns:
        dict: Variante, SR, convergiu, iterações, T_b, T_t, T_e [N·m], TR, K, balanço [%], σ absolutos por iteração
        (``s_T_b``, ``s_T_t``, ``s_TR``, ``s_K``; NaN quando ``desvio_janela`` faltar), e ``sem``, ``tau``,
        ``metodo``, ``achou`` e ``info`` vindos de ``inc``.
    """
    c = fator_torque(r)
    m = r["media_janela"]
    dp = r.get("desvio_janela", {})
    tb, tt, te = m["t_bomba"] * c, m["t_turbina"] * c, m["t_estator"] * c
    sb = abs(dp.get("t_bomba", float("nan")) * c)
    st = abs(dp.get("t_turbina", float("nan")) * c)
    tr = -tt / tb
    k = float(r.get("rpm_b") or RPM_B_PADRAO) / np.sqrt(abs(tb))
    setup = r.get("setup") or {}
    return {
        "variante": nome,
        "SR": r["sr"],
        "descricao": setup.get("descricao", descricao_padrao),
        "convergiu": bool(r.get("convergiu")),
        "iteracoes": r["it_fim"] - r["it_ini"],
        "T_b": tb, "T_t": tt, "T_e": te, "TR": tr, "K": k,
        "balanco_pct": 100 * (tb + tt + te) / abs(tb),
        "s_T_b": sb, "s_T_t": st,
        "s_TR": abs(tr) * np.hypot(st / tt, sb / tb),
        "s_K": 0.5 * k * sb / abs(tb),
        "sem": inc["sem"], "tau": inc["tau"], "metodo": inc["metodo"], "achou": inc["achou"], "info": inc["info"],
    }


def delta_pct(x, ref):
    """Variação percentual de um valor em relação à referência.

    Args:
        x (float): Valor da variante.
        ref (float): Valor da base.

    Returns:
        float: 100·(x − ref)/|ref|.
    """
    return 100 * (x - ref) / abs(ref)


def sigma_pct(s, x):
    """Desvio-padrão da janela em porcentagem do próprio valor.

    Args:
        s (float): Desvio-padrão absoluto.
        x (float): Valor médio.

    Returns:
        float: 100·s/|x|.
    """
    return 100 * s / abs(x)


def razao_sigma(delta, s_a, s_b):
    """Diferença entre duas médias em número de desvios combinados, |Δ|/sqrt(s_a² + s_b²).

    Usada com σ por iteração (colunas antigas) e com SEM (``abs_d_sobre_U_*``, onde é |Δ|/U).

    Args:
        delta (float): Diferença (variante − base).
        s_a (float): Incerteza (SEM ou σ) da variante.
        s_b (float): Incerteza (SEM ou σ) da base.

    Returns:
        float: |Δ|/combinado, ou NaN se o combinado for zero ou indefinido.
    """
    s = float(np.hypot(s_a, s_b))
    return abs(delta) / s if np.isfinite(s) and s > 0 else float("nan")


def linha(m, base):
    """Monta a linha da tabela de sensibilidade de uma variante (valores, Δ%, σ antigos, SEM, τ e |Δ|/U).

    Args:
        m (dict): Métricas da variante (saída de ``metricas``).
        base (dict): Métricas da base.

    Returns:
        dict: Linha com valores absolutos, Δ% em relação à base, σ antigos em %, |Δ|/σ antigos, ``sem_*``, ``tau_*``,
        ``metodo_incerteza`` e ``abs_d_sobre_U_*`` (|Δ|/U com U = √(SEM_variante² + SEM_base²)).
    """
    lin = {
        "variante": m["variante"],
        "descricao": m["descricao"],
        "convergiu": "sim" if m["convergiu"] else "não",
        "iteracoes": m["iteracoes"],
        "T_b [N·m]": m["T_b"],
        "T_t [N·m]": m["T_t"],
        "T_e [N·m]": m["T_e"],
        "TR": m["TR"],
        "K [rpm/√(N·m)]": m["K"],
        "balanco_pct": m["balanco_pct"],
        "dT_b_pct": delta_pct(m["T_b"], base["T_b"]),
        "dT_t_pct": delta_pct(m["T_t"], base["T_t"]),
        "dTR_pct": delta_pct(m["TR"], base["TR"]),
        "dK_pct": delta_pct(m["K"], base["K"]),
        "sigma_T_b_pct": sigma_pct(m["s_T_b"], m["T_b"]),
        "sigma_T_t_pct": sigma_pct(m["s_T_t"], m["T_t"]),
        "sigma_TR_pct": sigma_pct(m["s_TR"], m["TR"]),
        "abs_dT_b_sigma": razao_sigma(m["T_b"] - base["T_b"], m["s_T_b"], base["s_T_b"]),
        "abs_dTR_sigma": razao_sigma(m["TR"] - base["TR"], m["s_TR"], base["s_TR"]),
        "metodo_incerteza": f"base: {base['metodo']}; variante: {m['metodo']}",
    }
    for q in QUANTIDADES:
        lin[f"sem_{q}"] = m["sem"][q]
    for q in QUANTIDADES:
        lin[f"tau_{q}"] = m["tau"][q]
    for q in QUANTIDADES:
        lin[f"abs_d_sobre_U_{q}"] = razao_sigma(m[q] - base[q], m["sem"][q], base["sem"][q])
    return lin


def arredonda(v):
    """Arredonda números para 4 casas e troca NaN por célula vazia, para o CSV.

    Args:
        v (object): Valor da célula.

    Returns:
        object: Número arredondado, texto ou "" (para NaN).
    """
    if isinstance(v, float):
        return round(v, 4) if np.isfinite(v) else ""
    return v


def grafico(base, variantes, arq):
    """Gráfico de barras de Δ% de TR, T_b, T_t e K por variante, com barras de erro ±2U e linha zero.

    Args:
        base (dict): Métricas da base.
        variantes (list[dict]): Métricas de cada variante, na ordem do gráfico.
        arq (Path): PNG de saída.
    """
    quantidades = [("TR", "Razão de torque TR"), ("T_b", "Torque bomba T_b"),
                   ("T_t", "Torque turbina T_t"), ("K", "Fator de capacidade K")]
    nomes = [v["variante"] for v in variantes]
    x = np.arange(len(variantes))
    w = 0.2
    fig, ax = plt.subplots(figsize=(max(7.0, 2.2 * len(variantes) + 3), 4.2))
    for i, (q, rotulo) in enumerate(quantidades):
        vals = [delta_pct(v[q], base[q]) for v in variantes]
        erros = [200 * np.hypot(v["sem"][q], base["sem"][q]) / abs(base[q]) for v in variantes]
        ax.bar(x + (i - 1.5) * w, vals, w, yerr=np.nan_to_num(erros), capsize=2, label=rotulo, color=f"C{i}")
    ax.axhline(0, color="gray", lw=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(nomes)
    ax.set(ylabel="Δ% em relação à base (SR = 0)",
           title="Sensibilidade do estol ao setup; barras de erro = ±2U (SEM por batch means)")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(arq, dpi=130)
    plt.close(fig)


def descobrir_variantes(pasta, avisos):
    """Lê as variantes disponíveis em ``pasta``, calcula a incerteza de cada uma e registra as que faltam.

    Args:
        pasta (Path): Pasta com ``<variante>/dados/resultado.json``.
        avisos (list[str]): Lista onde são acrescentados os avisos (preenchida no lugar).

    Returns:
        tuple[list[dict], list[str]]: (métricas das variantes encontradas, nomes das variantes ausentes).
    """
    arqs = {p.parents[1].name: p for p in sorted(pasta.glob("*/dados/resultado.json"))} if pasta.exists() else {}
    if not pasta.exists():
        avisos.append(f"pasta de variantes não existe: {pasta}")
    faltam = [n for n in VARIANTES if n not in arqs]
    ordem = VARIANTES + sorted(n for n in arqs if n not in VARIANTES)
    achadas = []
    for nome in ordem:
        if nome not in arqs:
            continue
        r = ler_resultado(arqs[nome])
        if abs(r["sr"]) > 1e-9:
            avisos.append(f"variante {nome} não está em SR = 0 (sr = {r['sr']}); ignorada")
            continue
        inc = incerteza_ponto(registrar_variante(nome, arqs[nome].parents[1]), r)
        if not inc["achou"]:
            avisos.append(f"variante {nome}: .out da janela não achado ({inc['info']}); "
                          "SEM = σ por iteração (método antigo, conservador)")
        achadas.append(metricas(nome, r, "", inc))
    return achadas, faltam


def formata(v, casas):
    """Número formatado para o console, ou "-" quando indefinido.

    Args:
        v (object): Valor.
        casas (int): Casas decimais.

    Returns:
        str: Texto formatado.
    """
    if isinstance(v, (float, int, np.floating)) and not isinstance(v, bool):
        return f"{v:.{casas}f}" if np.isfinite(v) else "-"
    return str(v)


def imprimir_tabela(rows, colunas):
    """Imprime uma tabela do console com as colunas dadas (uma linha por ponto).

    Args:
        rows (list[dict]): Linhas da tabela (base primeiro).
        colunas (list[tuple[str, str, int]]): (chave, título, casas decimais) de cada coluna.
    """
    def largura(titulo):
        """Largura de uma coluna do console.

        Args:
            titulo (str): Título da coluna.

        Returns:
            int: Largura mínima de 10 caracteres ou o título mais um espaço.
        """
        return max(10, len(titulo) + 1)

    print("".join(f"{t:>{largura(t)}}" for _, t, _ in colunas))
    for r in rows:
        print("".join(f"{formata(r[k], c):>{largura(t)}}" for k, t, c in colunas))
    print()


def imprimir(rows, avisos):
    """Imprime as tabelas de sensibilidade (valores e razões; SEM e τ), os métodos e os avisos no console.

    Args:
        rows (list[dict]): Linhas da tabela (base primeiro).
        avisos (list[str]): Avisos a mostrar ao final.
    """
    imprimir_tabela(rows, COLUNAS_CONSOLE)
    imprimir_tabela(rows, COLUNAS_CONSOLE_SEM)
    for r in rows[1:]:
        print(f"metodo {r['variante']}: {r['metodo_incerteza']}")
    for a in avisos:
        print("AVISO:", a)


def main(argv=None):
    """Lê a base e as variantes disponíveis, grava a tabela e o gráfico e imprime a tabela.

    Args:
        argv (list[str] | None): Argumentos de linha de comando; None usa ``sys.argv``.

    Raises:
        SystemExit: Se o ``resultado.json`` da base não existir.
    """
    ap = argparse.ArgumentParser(description="Tabela e gráfico de sensibilidade do ponto de estol (SR = 0).")
    ap.add_argument("--dir-sens", type=Path, default=DIR / "sens",
                    help="pasta com as variantes (<variante>/dados/resultado.json)")
    ap.add_argument("--saida", type=Path, default=DIR, help="pasta de saída (dados/ e graficos/)")
    args = ap.parse_args(argv)
    if not BASE.exists():
        sys.exit(f"Base não encontrada: {BASE}")
    avisos = []
    r_base = ler_resultado(BASE)
    inc_base = incerteza_ponto(CHAVE_BASE, r_base)
    if not inc_base["achou"]:
        avisos.append(f"base: .out da janela não achado ({inc_base['info']}); "
                      "SEM = σ por iteração (método antigo, conservador)")
    base = metricas(NOME_BASE, r_base, DESCRICAO_BASE, inc_base)
    variantes, faltam = descobrir_variantes(args.dir_sens, avisos)
    rows = [linha(base, base)] + [linha(v, base) for v in variantes]

    arq_csv = args.saida / "dados" / "sensibilidade.csv"
    arq_csv.parent.mkdir(parents=True, exist_ok=True)
    with open(arq_csv, "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter=";")
        w.writeheader()
        w.writerows({k: arredonda(v) for k, v in r.items()} for r in rows)

    if variantes:
        arq_png = args.saida / "graficos" / "g6_sensibilidade.png"
        arq_png.parent.mkdir(parents=True, exist_ok=True)
        grafico(base, variantes, arq_png)
        print("gráfico:", arq_png)
    else:
        avisos.append("gráfico não gerado: nenhuma variante encontrada")
    if faltam:
        avisos.append("variantes ausentes (puladas): " + ", ".join(faltam))
    print("tabela:", arq_csv)
    imprimir(rows, avisos)


if __name__ == "__main__":
    main()
