"""Distribuição de y+ na parede do conversor CT_V1 (MRF): faixas, percentis, comparação e gráfico.

Uso: PYTHONIOENCODING=utf-8 python tools/analise_yplus_conversor.py

Entradas, por SR e peça (bomba, turbina, estator), em results/conversor/CT_V1_MRF/SRxxx/dados/:
- SRxxx_parede_<peca>_area.csv (quando existe): exportação ASCII do Fluent com as colunas
  y-plus e face-area-magnitude. Cada face pesa a sua área.
- SRxxx_parede_<peca>.csv (sempre; usado quando não há _area.csv): uma linha por face de parede,
  sem área. Cada face pesa 1.

Ponderação (coluna 'ponderacao' do CSV de saída):
- 'area': frações por faixa são fração da ÁREA da parede; percentis e média são ponderados por área.
- 'face': frações são fração do número de FACES; percentis e média não ponderados.

Percentil ponderado: ordena os y+ e calcula a distribuição acumulada com pontos médios,
F_i = (soma dos pesos até a face i - w_i / 2) / soma total dos pesos. Depois interpola
linearmente o y+ em F = p / 100 (np.interp). Com pesos iguais, isso aproxima o np.percentile
(diferença desprezível para o número de faces desta malha).

Saídas:
- CT_V1_MRF/dados/yplus_faixas.csv (separador ';'), uma linha por SR e peça.
- CT_V1_MRF/graficos/g4_yplus_faixas.png, com o SR marcado com '*' quando ponderado por face.
- Impressão no console: tabela resumida e comparação face x área para os pontos com os dois arquivos.
"""
import csv, sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
DIR = ROOT / "results" / "conversor" / "CT_V1_MRF"
PECAS = ["bomba", "turbina", "estator"]
# Faixas [inicio, fim) de y+. A última é aberta.
FAIXAS = [("<5", 0.0, 5.0), ("5-30", 5.0, 30.0), ("30-300", 30.0, 300.0),
          ("300-1000", 300.0, 1000.0), (">1000", 1000.0, np.inf)]
CORES_FAIXAS = ["#9ecae1", "#fdae6b", "#a1d99b", "#fc9272", "#de2d26"]
NOMES_FAIXAS_PT = {"<5": "< 5 (subcamada viscosa)", "5-30": "5-30 (amortecimento)",
                   "30-300": "30-300 (faixa log)", "300-1000": "300-1000", ">1000": "> 1000"}
plt.rcParams.update({"figure.dpi": 130, "axes.grid": True, "grid.alpha": 0.3, "font.size": 9})


def ler_parede(arquivo):
    """Lê y+ e, quando existir, a área de cada face de um CSV de parede exportado pelo Fluent.

    Args:
        arquivo (Path): Caminho do CSV com cabeçalho e separador vírgula (cabeçalho com espaços).

    Returns:
        tuple[np.ndarray, np.ndarray | None]: y+ de cada face e área de cada face. A área é None se
        a coluna 'face-area-magnitude' não existir. Faces com y+ não finito, ou com área não finita
        ou não positiva, são descartadas (o y+ e a área ficam alinhados).

    Raises:
        ValueError: Se o cabeçalho não tiver a coluna 'y-plus'.
    """
    with arquivo.open(newline="", errors="replace") as fh:
        leitor = csv.reader(fh, skipinitialspace=True)
        cabecalho = [c.strip() for c in next(leitor)]
        if "y-plus" not in cabecalho:
            raise ValueError(f"Coluna 'y-plus' ausente em {arquivo.name}: {cabecalho}")
        idx_y = cabecalho.index("y-plus")
        idx_a = cabecalho.index("face-area-magnitude") if "face-area-magnitude" in cabecalho else None
        yplus, area = [], []
        for linha in leitor:
            if len(linha) <= idx_y or (idx_a is not None and len(linha) <= idx_a):
                continue
            try:
                y = float(linha[idx_y])
                a = float(linha[idx_a]) if idx_a is not None else 1.0
            except ValueError:
                continue
            if np.isfinite(y) and np.isfinite(a) and a > 0:
                yplus.append(y)
                area.append(a)
    y_arr = np.array(yplus, dtype=float)
    a_arr = np.array(area, dtype=float) if idx_a is not None else None
    return y_arr, a_arr


def percentil_ponderado(valores, pesos, p):
    """Calcula o percentil p de uma amostra com pesos (ponderação por área ou por face).

    Método: ordena os valores; a distribuição acumulada usa pontos médios
    F_i = (soma dos pesos até i - w_i / 2) / soma total; o resultado é interpolado linearmente
    em F = p / 100. Com pesos iguais, reduz-se a uma aproximação do np.percentile.

    Args:
        valores (np.ndarray): Valores da amostra (y+), não vazios.
        pesos (np.ndarray): Peso positivo de cada valor (área da face ou 1), mesmo tamanho.
        p (float): Percentil desejado, de 0 a 100.

    Returns:
        float: Valor do percentil p na escala de 'valores'.
    """
    ordem = np.argsort(valores, kind="mergesort")
    v = valores[ordem]
    w = pesos[ordem]
    acumulado = np.cumsum(w)
    frac_acum = (acumulado - w / 2.0) / acumulado[-1]
    return float(np.interp(p / 100.0, frac_acum, v))


def estatisticas(yplus, pesos):
    """Calcula frações por faixa de y+, média e percentis ponderados pelos pesos.

    Args:
        yplus (np.ndarray): Valores de y+ das faces da parede.
        pesos (np.ndarray): Peso de cada face (área, no método por área, ou 1, no método por face).
            Mesmo tamanho de yplus.

    Returns:
        dict: Chaves 'n' (número de faces), 'peso_total' (soma dos pesos), 'media' (média ponderada),
        'p5', 'p50', 'p95' (percentis ponderados) e 'frac' (dicionário faixa -> fração 0 a 1 do peso
        total). Se não houver faces, os valores numéricos são NaN e as frações são zero.
    """
    n = yplus.size
    if n == 0:
        return {"n": 0, "peso_total": 0.0, "media": float("nan"), "p5": float("nan"),
                "p50": float("nan"), "p95": float("nan"),
                "frac": {nome: 0.0 for nome, _, _ in FAIXAS}}
    total = float(pesos.sum())
    frac = {}
    for nome, ini, fim in FAIXAS:
        mascara = (yplus >= ini) & (yplus < fim)
        frac[nome] = float(pesos[mascara].sum() / total)
    return {
        "n": n,
        "peso_total": total,
        "media": float(np.average(yplus, weights=pesos)),
        "p5": percentil_ponderado(yplus, pesos, 5),
        "p50": percentil_ponderado(yplus, pesos, 50),
        "p95": percentil_ponderado(yplus, pesos, 95),
        "frac": frac,
    }


def listar_sr():
    """Lista os pontos SR que têm as exportações de parede (face) das três peças.

    Returns:
        list[str]: Nomes dos SR (ex.: 'SR000'), em ordem crescente.

    Raises:
        SystemExit: Se nenhum ponto SR com dados de parede for encontrado.
    """
    srs = []
    for pasta in sorted(DIR.glob("SR*")):
        dados = pasta / "dados"
        if all((dados / f"{pasta.name}_parede_{p}.csv").exists() for p in PECAS):
            srs.append(pasta.name)
    if not srs:
        sys.exit("Nenhum SR com exportações de parede encontrado.")
    return srs


def carregar_peca(sr, peca):
    """Carrega y+ e pesos de uma peça de um SR, usando área quando o _area.csv existir.

    Args:
        sr (str): Nome do ponto SR (ex.: 'SR000').
        peca (str): Nome da peça ('bomba', 'turbina' ou 'estator').

    Returns:
        dict: Chaves 'ponderacao' ('area' ou 'face'), 'yplus' (np.ndarray) e 'pesos' (np.ndarray:
        área das faces, ou 1 por face).
    """
    dados = DIR / sr / "dados"
    arq_area = dados / f"{sr}_parede_{peca}_area.csv"
    if arq_area.exists():
        yplus, area = ler_parede(arq_area)
        return {"ponderacao": "area", "yplus": yplus, "pesos": area}
    yplus, _ = ler_parede(dados / f"{sr}_parede_{peca}.csv")
    return {"ponderacao": "face", "yplus": yplus, "pesos": np.ones(yplus.size)}


def coletar(srs):
    """Percorre os SR e as peças e calcula as estatísticas de y+ pelo método disponível.

    Args:
        srs (list[str]): Nomes dos SR a processar.

    Returns:
        list[dict]: Uma linha por (SR, peça), com as estatísticas de estatisticas() mais as chaves
        'sr', 'peca' e 'ponderacao'.
    """
    linhas = []
    for sr in srs:
        for peca in PECAS:
            dados = carregar_peca(sr, peca)
            est = estatisticas(dados["yplus"], dados["pesos"])
            est.update({"sr": sr, "peca": peca, "ponderacao": dados["ponderacao"]})
            linhas.append(est)
            print(f"{sr} {peca:8s} [{dados['ponderacao']:4s}] faces={est['n']:7d} "
                  f"y+ med={est['media']:8.2f} p50={est['p50']:8.2f}")
    return linhas


def comparar_face_area(srs):
    """Compara o método por face com o método por área nos pontos que têm o _area.csv.

    Para cada SR e peça com _area.csv, usa o mesmo y+ nos dois métodos: por face (peso 1) e por
    área (peso face-area-magnitude). Assim a diferença vem só da ponderação.

    Args:
        srs (list[str]): Nomes dos SR a comparar.

    Returns:
        list[dict]: Uma linha por (SR, peça) com 'sr', 'peca', 'f3030_face', 'f3030_area' (fração
        da faixa 30-300 em cada método), 'p50_face' e 'p50_area'.
    """
    comparacao = []
    for sr in srs:
        for peca in PECAS:
            arq_area = DIR / sr / "dados" / f"{sr}_parede_{peca}_area.csv"
            if not arq_area.exists():
                continue
            yplus, area = ler_parede(arq_area)
            por_face = estatisticas(yplus, np.ones(yplus.size))
            por_area = estatisticas(yplus, area)
            comparacao.append({"sr": sr, "peca": peca,
                               "f3030_face": por_face["frac"]["30-300"],
                               "f3030_area": por_area["frac"]["30-300"],
                               "p50_face": por_face["p50"], "p50_area": por_area["p50"]})
    return comparacao


def gravar_csv(linhas):
    """Grava a tabela de faixas, percentis e ponderação em CT_V1_MRF/dados/yplus_faixas.csv.

    Args:
        linhas (list[dict]): Resultado de coletar().

    Returns:
        Path: Caminho do CSV gravado.
    """
    saida = DIR / "dados" / "yplus_faixas.csv"
    cab = ["SR", "peca", "ponderacao", "n_faces", "area_total_m2", "yplus_medio",
           "yplus_p5", "yplus_p50", "yplus_p95"]
    cab += [f"frac_{nome}" for nome, _, _ in FAIXAS]
    with saida.open("w", newline="", encoding="utf-8") as fh:
        escritor = csv.writer(fh, delimiter=";")
        escritor.writerow(cab)
        for l in linhas:
            area_total = f"{l['peso_total']:.6e}" if l["ponderacao"] == "area" else ""
            escritor.writerow([l["sr"], l["peca"], l["ponderacao"], l["n"], area_total,
                               f"{l['media']:.4f}", f"{l['p5']:.4f}",
                               f"{l['p50']:.4f}", f"{l['p95']:.4f}"]
                              + [f"{l['frac'][nome]:.6f}" for nome, _, _ in FAIXAS])
    return saida


def imprimir_tabela(linhas):
    """Imprime a tabela resumida: frações por faixa (%), P50 e ponderação por SR e peça.

    Args:
        linhas (list[dict]): Resultado de coletar().

    Returns:
        None
    """
    cab = f"{'SR':7s}{'peca':9s}{'pond':6s}" + "".join(f"{nome:>10s}" for nome, _, _ in FAIXAS)
    print(cab + f"{'P50':>9s}")
    for l in linhas:
        valores = "".join(f"{l['frac'][nome] * 100:10.2f}" for nome, _, _ in FAIXAS)
        print(f"{l['sr']:7s}{l['peca']:9s}{l['ponderacao']:6s}{valores}{l['p50']:9.2f}")


def imprimir_comparacao(comparacao):
    """Imprime a comparação face x área (fração 30-300 e P50) e a diferença área - face.

    Args:
        comparacao (list[dict]): Resultado de comparar_face_area().

    Returns:
        None
    """
    if not comparacao:
        print("Nenhum ponto com _area.csv para comparar.")
        return
    print(f"{'SR':7s}{'peca':9s}{'30-300 face%':>13s}{'30-300 area%':>13s}{'delta pp':>10s}"
          f"{'P50 face':>10s}{'P50 area':>10s}{'delta P50':>10s}")
    for c in comparacao:
        f_face = c["f3030_face"] * 100
        f_area = c["f3030_area"] * 100
        print(f"{c['sr']:7s}{c['peca']:9s}{f_face:13.2f}{f_area:13.2f}{f_area - f_face:10.2f}"
              f"{c['p50_face']:10.2f}{c['p50_area']:10.2f}{c['p50_area'] - c['p50_face']:10.2f}")


def grafico(linhas, srs):
    """Gera o gráfico de barras empilhadas das frações por faixa de y+.

    Uma barra por peça, agrupadas por SR. O eixo y e o título dizem se a fração é de área ou de faces.
    Se um SR for ponderado por face, o rótulo do SR recebe '*' e uma nota de rodapé explica.

    Args:
        linhas (list[dict]): Resultado de coletar().
        srs (list[str]): Nomes dos SR, na ordem do eixo x.

    Returns:
        Path: Caminho do PNG gravado.
    """
    fig, ax = plt.subplots(figsize=(max(9, 1.1 * len(srs) * len(PECAS)), 5.5))
    largura = 0.26
    espaco_grupo = 1.0
    xticks, xlabels = [], []
    rotulos_usados = set()
    ponderacoes = {l["ponderacao"] for l in linhas}

    def desenhar_barra(x, base_frac, rotulo_faixa):
        """Desenha uma barra empilhada para uma peça de um SR.

        Args:
            x (float): Posição horizontal da barra.
            base_frac (dict): Frações por faixa (faixa -> 0 a 1) da peça.
            rotulo_faixa (bool): Se True, inclui a faixa na legenda (só na primeira barra).

        Returns:
            None
        """
        base = 0.0
        for (nome, _, _), cor in zip(FAIXAS, CORES_FAIXAS):
            valor = base_frac[nome] * 100
            ax.bar(x, valor, width=largura, bottom=base, color=cor,
                   edgecolor="white", linewidth=0.4,
                   label=NOMES_FAIXAS_PT[nome] if rotulo_faixa and nome not in rotulos_usados else None)
            if rotulo_faixa:
                rotulos_usados.add(nome)
            if valor >= 6:
                ax.text(x, base + valor / 2, f"{valor:.0f}", ha="center", va="center", fontsize=7)
            base += valor

    for i, sr in enumerate(srs):
        centro = i * espaco_grupo * (len(PECAS) * largura + 0.5)
        por_face = any(l["ponderacao"] == "face" for l in linhas if l["sr"] == sr)
        for j, peca in enumerate(PECAS):
            x = centro + (j - 1) * largura
            linha = next(l for l in linhas if l["sr"] == sr and l["peca"] == peca)
            desenhar_barra(x, linha["frac"], rotulo_faixa=True)
            xticks.append(x)
            xlabels.append(peca.capitalize())
        rotulo_sr = sr + ("*" if por_face else "")
        ax.text(centro, -14, rotulo_sr, ha="center", va="top", fontsize=9, fontweight="bold")

    if ponderacoes == {"area"}:
        eixo_y, sufixo = "Fração da ÁREA da parede (%)", "fração de área"
    elif ponderacoes == {"face"}:
        eixo_y, sufixo = "Fração das FACES da parede (%)", "fração de faces"
    else:
        eixo_y, sufixo = "Fração (%): área, exceto * (faces)", "área; * = fração de faces"

    ax.set_xticks(xticks)
    ax.set_xticklabels(xlabels, fontsize=8)
    ax.set_ylabel(eixo_y)
    ax.set_ylim(0, 100)
    ax.set_title(f"CT_V1 MRF: distribuição de y+ na parede por peça e ponto SR ({sufixo})")
    ax.legend(title="Faixa de y+", bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=8)
    ax.grid(axis="x", visible=False)
    if "face" in ponderacoes:
        fig.text(0.01, 0.01, "* SR ponderado por face (sem arquivo _area.csv para essa peça/ponto).",
                 fontsize=7, style="italic")
    fig.subplots_adjust(bottom=0.2, right=0.78)
    saida = DIR / "graficos" / "g4_yplus_faixas.png"
    fig.savefig(saida)
    plt.close(fig)
    return saida


def main():
    """Executa a análise completa: estatísticas, CSV, comparação face x área e gráfico.

    Returns:
        None
    """
    srs = listar_sr()
    linhas = coletar(srs)
    print()
    imprimir_tabela(linhas)
    print()
    print("Comparação face x área (mesmo y+; só muda o peso):")
    imprimir_comparacao(comparar_face_area(srs))
    csv_saida = gravar_csv(linhas)
    png_saida = grafico(linhas, srs)
    print(f"\nCSV: {csv_saida}")
    print(f"Gráfico: {png_saida}")


if __name__ == "__main__":
    main()
