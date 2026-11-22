"""Gráficos do conversor CT_V1 (MRF) para a Seção 4.5 do TG, a partir dos CSVs já gerados.

Uso (rodar com ``PYTHONIOENCODING=utf-8`` e o ``.venv241``)::

    python tools/figuras_tg_conversor.py            gera todos os gráficos
    python tools/figuras_tg_conversor.py gci lit    gera só os escolhidos (nomes em ``FIGURAS``)

Entradas (nada é recalculado no Fluent)::

    results/conversor/CT_V1_MRF/dados/desempenho.csv, sensibilidade.csv, euler_interfaces.csv, yplus_faixas.csv
    results/conversor/CT_V1_MRF/dados/CT_V1_MRF_rep_20261121-022556.out   (histórico da varredura da G0, em Cm)
    results/conversor/CT_V1_GCI/dados/gci.csv
    results/conversor/CT_V1_GCI_G2/sens/coupled/dados/resultado.json      (variante Coupled, G2)
    results/conversor/CT_V1_GCI_G2/SR000/dados/resultado.json            (base SIMPLE da G2)
    refs/conversor_literatura.csv                                         (via comparacao_literatura_conversor)
    results/conversor/CT_V1_MRF/dados/retorno_if_bt.csv                   (gerado aqui na 1ª vez, a partir dos
                                                                           exports da if_bt em D:/TG_backup_drive)

Saídas: PNG em ``tg-fem-unicamp/Figuras/conversor_res_*.png`` (rótulos em português, vírgula decimal).

Convenções: torques em N·m (o report file da varredura da G0 está em Cm e é convertido por 0,6125); TR = −T_t/T_b;
η = SR·TR; K = n_b/√|T_b|. O ponto de acoplamento (T_e = 0) é interpolado linearmente em T_e.
"""
import csv
import json
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

sys.path.insert(0, str(Path(__file__).resolve().parent))
from comparacao_literatura_conversor import ler_literatura  # noqa: E402

RAIZ = Path(__file__).resolve().parents[1]
CONV = RAIZ / "results" / "conversor"
DADOS = CONV / "CT_V1_MRF" / "dados"
SAIDA = RAIZ.parent / "tg-fem-unicamp" / "Figuras"
BACKUP_CSV = Path("D:/TG_backup_drive/03_dados_csv/results/conversor/CT_V1_MRF")
CM_PARA_NM = 0.6125
SRS_CAMPOS = ("SR000", "SR020", "SR040", "SR060", "SR070", "SR080", "SR085", "SR090", "SR095", "SR100")
COR = {"bomba": "#1f77b4", "turbina": "#d62728", "estator": "#2ca02c"}

plt.rcParams.update({"font.size": 9, "axes.grid": True, "grid.alpha": 0.3, "figure.dpi": 100,
                     "savefig.dpi": 220, "axes.titlesize": 9, "legend.fontsize": 7.5, "mathtext.fontset": "cm"})


def virgula(ax):
    """Troca o ponto decimal por vírgula nos rótulos dos eixos x e y.

    Args:
        ax (matplotlib.axes.Axes): Eixo a formatar.
    """
    def fmt(v, _):
        """Formata um valor de marcação com vírgula decimal.

        Args:
            v (float): Valor.
            _ (int): Posição (não usada).

        Returns:
            str: Texto do rótulo.
        """
        s = f"{v:.10g}"
        if "e" in s:
            s = f"{v:.3g}"
        return s.replace(".", ",").replace("-", "−")
    ax.xaxis.set_major_formatter(FuncFormatter(fmt))
    ax.yaxis.set_major_formatter(FuncFormatter(fmt))


def ler_csv(caminho):
    """Lê um CSV separado por ``;``.

    Args:
        caminho (Path): Arquivo.

    Returns:
        list[dict]: Linhas.
    """
    with open(caminho, encoding="utf-8-sig") as f:
        return list(csv.DictReader(f, delimiter=";"))


def desempenho():
    """Tabela de desempenho da G0 como arrays ordenados por SR.

    Returns:
        dict[str, np.ndarray]: Colunas numéricas de ``desempenho.csv``.
    """
    ls = sorted(ler_csv(DADOS / "desempenho.csv"), key=lambda l: float(l["SR"]))
    cols = ("SR", "T_b", "T_t", "T_e", "TR", "eta", "K", "balanco_pct", "mdot_bt", "dp_TR", "dp_T_b_pct",
            "dp_T_t_pct")
    return {c: np.array([float(l[c]) for l in ls]) for c in cols}


def sr_acoplamento(d):
    """SR em que T_e cruza zero, por interpolação linear.

    Args:
        d (dict): Saída de ``desempenho``.

    Returns:
        float: SR de acoplamento.
    """
    sr, te = d["SR"], d["T_e"]
    for i in range(len(sr) - 1):
        if te[i] * te[i + 1] <= 0:
            return sr[i] - te[i] * (sr[i + 1] - sr[i]) / (te[i + 1] - te[i])
    return float("nan")


def sombrear(ax, sr_ac):
    """Sombreia a faixa acima do acoplamento, em que o estator travado não representa o conversor real.

    Args:
        ax (matplotlib.axes.Axes): Eixo.
        sr_ac (float): SR de acoplamento.
    """
    ax.axvspan(sr_ac, 1.03, color="0.88", alpha=0.6, lw=0, zorder=0)
    ax.axvline(sr_ac, color="0.4", ls=":", lw=0.9)


def salvar(fig, nome):
    """Grava a figura em ``SAIDA``.

    Args:
        fig (matplotlib.figure.Figure): Figura.
        nome (str): Nome do arquivo sem extensão.
    """
    SAIDA.mkdir(parents=True, exist_ok=True)
    fig.savefig(SAIDA / f"{nome}.png", bbox_inches="tight")
    plt.close(fig)
    print("gravado", SAIDA / f"{nome}.png")


# ---------------------------------------------------------------------------------------------------------------
def fig_convergencia():
    """Histórico dos torques e do balanço no SR 0 da G0, com as duas janelas do critério de parada."""
    arq = DADOS / "CT_V1_MRF_rep_20261121-022556.out"
    linhas = [l.split() for l in arq.read_text().splitlines()[3:] if l.strip()]
    a = np.array([[float(x) for x in l] for l in linhas])
    a = a[a[:, 0] <= 600]
    it, tb, tt, te = a[:, 0], a[:, 1] * CM_PARA_NM, a[:, 2] * CM_PARA_NM, a[:, 3] * CM_PARA_NM
    fig, axs = plt.subplots(1, 2, figsize=(7.2, 2.7))
    ax = axs[0]
    ax.plot(it, tt / 1e3, color=COR["turbina"], lw=0.8, label="$T_t$")
    ax.plot(it, -tb / 1e3, color=COR["bomba"], lw=0.8, label="$-T_b$")
    for w0, w1, c in ((200, 400, "0.85"), (400, 600, "0.7")):
        ax.axvspan(w0, w1, color=c, alpha=0.35, lw=0)
        for y, cor in ((tt, COR["turbina"]), (-tb, COR["bomba"])):
            m = y[(it > w0) & (it <= w1)].mean() / 1e3
            ax.hlines(m, w0, w1, color=cor, lw=1.6)
    ax.set_ylim(15, 35)
    ax.set_xlabel("Iteração")
    ax.set_ylabel("Torque [kN·m]")
    ax.set_title("(a) Torques na bomba e na turbina")
    ax.legend(loc="lower right")
    ax = axs[1]
    bal = (tb + tt + te) / np.abs(tb) * 100
    ax.plot(it, bal, color="k", lw=0.7)
    ax.axhline(0, color="0.4", lw=0.8)
    ax.set_ylim(-6, 6)
    ax.set_xlabel("Iteração")
    ax.set_ylabel(r"$(T_b+T_t+T_e)/|T_b|$ [%]")
    ax.set_title("(b) Balanço de torque instantâneo")
    for ax in axs:
        ax.set_xlim(0, 600)
        virgula(ax)
    fig.tight_layout()
    salvar(fig, "conversor_res_convergencia")


def fig_desempenho():
    """TR, η e K em função de SR na G0, com o desvio-padrão na janela final e a faixa acima do acoplamento."""
    d = desempenho()
    sr_ac = sr_acoplamento(d)
    fig, axs = plt.subplots(1, 3, figsize=(7.4, 2.5))
    ax = axs[0]
    ax.errorbar(d["SR"], d["TR"], yerr=d["dp_TR"], fmt="o-", color="k", ms=3.5, lw=1, capsize=2)
    ax.axhline(1, color="0.4", ls="--", lw=0.8)
    ax.set_ylabel("TR")
    ax.set_title("(a) Razão de torque")
    ax = axs[1]
    ax.plot(d["SR"], d["eta"], "o-", color="k", ms=3.5, lw=1)
    ax.set_ylabel(r"$\eta$")
    ax.set_title("(b) Eficiência")
    ax = axs[2]
    ax.plot(d["SR"], d["K"], "o-", color="k", ms=3.5, lw=1)
    ax.set_ylabel(r"$K$ [rpm/(N·m)$^{1/2}$]")
    ax.set_title("(c) Fator de capacidade")
    for ax in axs:
        sombrear(ax, sr_ac)
        ax.set_xlim(-0.03, 1.03)
        ax.set_xlabel("SR")
        virgula(ax)
    fig.tight_layout()
    salvar(fig, "conversor_res_desempenho")


def fig_torques():
    """Torques nas três rodas, balanço de torque e vazão no circuito em função de SR."""
    d = desempenho()
    sr_ac = sr_acoplamento(d)
    fig, axs = plt.subplots(1, 3, figsize=(7.4, 2.5))
    ax = axs[0]
    for k, n in (("T_b", "bomba"), ("T_t", "turbina"), ("T_e", "estator")):
        ax.plot(d["SR"], d[k] / 1e3, "o-", ms=3, lw=1, color=COR[n], label=f"$T_{n[0]}$ ({n})")
    ax.axhline(0, color="0.4", lw=0.8)
    ax.set_ylabel("Torque [kN·m]")
    ax.set_title("(a) Torques")
    ax.legend(loc="upper right", fontsize=7)
    ax = axs[1]
    ax.plot(d["SR"], d["balanco_pct"], "o-", ms=3, lw=1, color="k")
    ax.axhspan(-1, 1, color="tab:green", alpha=0.12, lw=0)
    ax.axhline(0, color="0.4", lw=0.8)
    ax.set_ylabel(r"$(T_b+T_t+T_e)/|T_b|$ [%]")
    ax.set_title("(b) Balanço de torque")
    ax = axs[2]
    ax.plot(d["SR"], d["mdot_bt"], "o-", ms=3, lw=1, color="k")
    ax.set_ylabel(r"$\dot{m}$ [kg/s]")
    ax.set_title("(c) Vazão no circuito")
    for ax in axs:
        sombrear(ax, sr_ac)
        ax.set_xlim(-0.03, 1.03)
        ax.set_xlabel("SR")
        virgula(ax)
    fig.tight_layout()
    salvar(fig, "conversor_res_torques")


def fig_euler():
    """Razão entre o torque de Euler (fluxo convectivo de momento angular nas interfaces) e o torque nas paredes."""
    ls = ler_csv(DADOS / "euler_interfaces.csv")
    fig, ax = plt.subplots(figsize=(4.6, 2.7))
    for zona, n in (("fluido_bomba", "bomba"), ("fluido_turbina", "turbina"), ("fluido_estator", "estator")):
        sel = sorted((float(l["SR"]), float(l["razao_conv"])) for l in ls if l["zona"] == zona)
        ax.plot([s for s, _ in sel], [r for _, r in sel], "o-", ms=3.5, lw=1, color=COR[n], label=n.capitalize())
    ax.axhline(1, color="0.4", lw=0.8)
    ax.axhspan(0.97, 1.03, color="0.85", alpha=0.5, lw=0)
    ax.set_xlim(-0.03, 1.03)
    ax.set_ylim(0.86, 1.04)
    ax.set_xlabel("SR")
    ax.set_ylabel(r"$T_{Euler}/T_{parede}$")
    ax.legend(loc="lower left", ncol=3)
    virgula(ax)
    fig.tight_layout()
    salvar(fig, "conversor_res_euler")


def fig_yplus():
    """Fração da área de parede de cada roda por faixa de y⁺, em função de SR."""
    ls = ler_csv(DADOS / "yplus_faixas.csv")
    faixas = [("frac_<5", "$y^+<5$", "#08306b"), ("frac_5-30", "$5$–$30$", "#6baed6"),
              ("frac_30-300", "$30$–$300$", "#c6dbef"), ("frac_300-1000", "$300$–$1000$", "#fd8d3c"),
              ("frac_>1000", "$>1000$", "#a63603")]
    fig, axs = plt.subplots(1, 3, figsize=(7.4, 2.6), sharey=True)
    for ax, peca in zip(axs, ("bomba", "turbina", "estator")):
        sel = sorted((l for l in ls if l["peca"] == peca), key=lambda l: l["SR"])
        x = np.arange(len(sel))
        base = np.zeros(len(sel))
        for k, rot, cor in faixas:
            v = np.array([float(l[k]) for l in sel]) * 100
            ax.bar(x, v, bottom=base, color=cor, width=0.8, label=rot, edgecolor="0.3", lw=0.3)
            base += v
        ax.set_xticks(x, [f"{int(l['SR'][2:]) / 100:g}".replace(".", ",") for l in sel], rotation=90, fontsize=7)
        ax.set_title(peca.capitalize())
        ax.set_xlabel("SR")
        ax.set_ylim(0, 100)
        ax.grid(False)
    axs[0].set_ylabel("Fração da área [%]")
    axs[2].legend(loc="center left", bbox_to_anchor=(1.02, 0.5), fontsize=7)
    fig.tight_layout()
    salvar(fig, "conversor_res_yplus")


def fig_literatura():
    """Tendências adimensionais do CT_V1 e de conversores da literatura (ensaios e RANS sem cavitação)."""
    d = desempenho()
    sr_ac = sr_acoplamento(d)
    series = [s for s in ler_literatura()
              if s["origem"] == "exp" or (s["fonte"] == "guo2021" and "sem cav" in s["rotulo"])]
    nomes = {"zhang2023": "Zhang et al. (2023)", "guo2021": "Guo et al. (2021)", "ran2021": "Ran et al. (2021)",
             "nadeem2020": "Nadeem (2020)"}
    cores = {"zhang2023": "#9467bd", "guo2021": "#8c564b", "ran2021": "#e377c2", "nadeem2020": "#7f7f7f"}
    fig, axs = plt.subplots(1, 3, figsize=(7.4, 2.9))
    for s in series:
        estilo = dict(marker="o", ls="-", ms=2.8, lw=0.9) if s["origem"] == "exp" else dict(marker="^", ls="--", ms=2.8,
                                                                                          lw=0.9)
        rot = nomes[s["fonte"]] + (", ensaio" if s["origem"] == "exp" else ", RANS sem cavitação")
        tr0 = s["TR"][np.argmin(s["SR"])]
        axs[0].plot(s["SR"], s["TR"] / tr0, color=cores[s["fonte"]], label=rot, **estilo)
        axs[1].plot(s["SR"], s["eta"], color=cores[s["fonte"]], **estilo)
        ok = ~np.isnan(s["KK"])
        axs[2].plot(s["SR"][ok], s["KK"][ok], color=cores[s["fonte"]], **estilo)
    kw = dict(color="k", lw=1.8, marker="D", ms=3.5, label="CT_V1, este trabalho")
    axs[0].plot(d["SR"], d["TR"] / d["TR"][0], **kw)
    axs[1].plot(d["SR"], d["eta"], **kw)
    axs[2].plot(d["SR"], d["K"] / d["K"][0], **kw)
    axs[0].set_ylabel(r"TR/TR$_{estol}$")
    axs[1].set_ylabel(r"$\eta$")
    axs[2].set_ylabel(r"$K/K_{estol}$")
    axs[2].set_ylim(0.8, 2.0)
    for ax, t in zip(axs, ("(a) Razão de torque normalizada", "(b) Eficiência", "(c) Fator de capacidade normalizado")):
        sombrear(ax, sr_ac)
        ax.set_xlim(-0.03, 1.03)
        ax.set_xlabel("SR")
        ax.set_title(t)
        virgula(ax)
    h, lab = axs[0].get_legend_handles_labels()
    fig.legend(h, lab, loc="lower center", ncol=3, fontsize=7, bbox_to_anchor=(0.5, -0.1))
    fig.tight_layout()
    salvar(fig, "conversor_res_literatura")


def fig_sensibilidade():
    """Variação de TR, |T_b| e T_t em relação à configuração base, por variante, com barras de ±2U."""
    ls = ler_csv(DADOS / "sensibilidade.csv")
    nomes = {"sst_sem_cc": "SST sem\ncorreção de\ncurvatura", "ke_realizable": r"$k$-$\varepsilon$" + "\nRealizable",
             "fou": "Primeira\nordem", "meio_passo": "Meio passo\n(10°)"}
    var, dados = [], []
    for l in ls:
        v = l["\ufeffvariante"] if "\ufeffvariante" in l else l["variante"]
        if v not in nomes:
            continue
        lin = []
        for k, ku in (("dTR_pct", "abs_d_sobre_U_TR"), ("dT_b_pct", "abs_d_sobre_U_T_b"), ("dT_t_pct", "abs_d_sobre_U_T_t")):
            dv = float(l[k])
            if k == "dT_b_pct":
                dv = -dv                                   # variação do módulo de T_b (T_b < 0)
            u = abs(float(l[k])) / float(l[ku]) if float(l[ku]) > 0 else 0.0
            lin.append((dv, 2 * u))
        var.append(nomes[v])
        dados.append(lin)
    base = json.loads((CONV / "CT_V1_GCI_G2" / "SR000" / "dados" / "resultado.json").read_text(encoding="utf-8"))
    cp = json.loads((CONV / "CT_V1_GCI_G2" / "sens" / "coupled" / "dados" / "resultado.json").read_text(encoding="utf-8"))
    def tor(r):
        """Torques médios (N·m) de um resultado.json.

        Args:
            r (dict): Conteúdo do resultado.json.

        Returns:
            tuple[float, float]: T_b e T_t.
        """
        m = r["media_janela"]
        f = 1.0 if r.get("unidade_torque") == "N·m" else CM_PARA_NM
        return m["t_bomba"] * f, m["t_turbina"] * f
    b0, t0 = tor(base)
    b1, t1 = tor(cp)
    var.append("Coupled\n(malha G2)")
    dados.append([((t1 / -b1) / (t0 / -b0) * 100 - 100, 0), ((abs(b1) / abs(b0) - 1) * 100, 0),
                  ((t1 / t0 - 1) * 100, 0)])
    fig, ax = plt.subplots(figsize=(6.6, 2.8))
    x = np.arange(len(var))
    w = 0.26
    for j, (rot, cor) in enumerate((("TR", "k"), (r"$|T_b|$", COR["bomba"]), (r"$T_t$", COR["turbina"]))):
        v = [dd[j][0] for dd in dados]
        e = [dd[j][1] for dd in dados]
        ax.bar(x + (j - 1) * w, v, w, color=cor, alpha=0.85, label=rot, yerr=e, capsize=2, error_kw=dict(lw=0.8))
    ax.axhline(0, color="0.3", lw=0.8)
    ax.set_xticks(x, var, fontsize=7.5)
    ax.set_ylabel("Variação em relação à base [%]")
    ax.legend(ncol=3, loc="upper left")
    ax.grid(axis="x", visible=False)
    virgula(ax)
    ax.set_xticks(x, var, fontsize=7.5)
    fig.tight_layout()
    salvar(fig, "conversor_res_sensibilidade")


def fig_gci():
    """TR, T_b, T_t e K nas malhas G2, G1 e G0, em SR 0 e 0,4, com ±U_iter e a banda do GCI21 na G0."""
    ls = ler_csv(CONV / "CT_V1_GCI" / "dados" / "gci.csv")
    chave = "\ufeffSR" if "\ufeffSR" in ls[0] else "SR"
    grand = (("TR", "TR", 1), ("T_b", r"$T_b$ [kN·m]", 1e-3), ("T_t", r"$T_t$ [kN·m]", 1e-3),
             ("K", r"$K$ [rpm/(N·m)$^{1/2}$]", 1))
    fig, axs = plt.subplots(2, 4, figsize=(7.6, 4.0))
    for i, sr in enumerate(("0", "0.4")):
        for j, (g, rot, esc) in enumerate(grand):
            l = next(r for r in ls if r[chave] == sr and r["grandeza"] == g)
            h = np.array([float(l[f"h{k}_mm"]) for k in (1, 2, 3)])
            phi = np.array([float(l[f"phi{k}"]) for k in (1, 2, 3)]) * esc
            sem = np.array([float(l[f"sem{k}"]) for k in (1, 2, 3)]) * abs(esc)
            ax = axs[i, j]
            ax.errorbar(h, phi, yerr=2 * sem, fmt="o", color="k", ms=4, capsize=2.5, lw=1)
            ax.plot(h, phi, color="0.6", lw=0.8, zorder=0)
            gci = float(l["GCI21_pct"]) / 100 * abs(phi[0])
            ax.add_patch(plt.Rectangle((h[0] - 0.12, phi[0] - gci), 0.24, 2 * gci, color="tab:orange", alpha=0.35,
                                       lw=0))
            for hk, pk, n in zip(h, phi, ("G0", "G1", "G2")):
                ax.annotate(n, (hk, pk), textcoords="offset points", xytext=(5, 3), fontsize=7)
            ax.set_xlim(3.0, 6.3)
            ax.set_title(f"{rot.split(' [')[0]}, SR = {sr.replace('.', ',')}")
            ax.set_ylabel(rot)
            if i == 1:
                ax.set_xlabel(r"$h$ [mm]")
            virgula(ax)
    fig.tight_layout()
    salvar(fig, "conversor_res_gci")


# ---------------------------------------------------------------------------------------------------------------
def calcular_retorno():
    """Distribuição radial da vazão líquida na interface bomba-turbina, a partir dos exports por célula.

    A vazão de cada face é ρ(u·A), com A apontando para fora da zona que exporta; o sinal é ajustado para que
    positivo seja no sentido bomba -> turbina. Grava ``retorno_if_bt.csv`` (SR, lado, faixa de raio, vazão).
    Os valores por centro de célula não conservam a massa entre os lados (desvio de até ~15 % na vazão líquida),
    por isso o gráfico usa a média dos dois lados.

    Returns:
        Path: CSV gravado.
    """
    bordas = np.concatenate([np.linspace(0.213, 0.270, 9), np.linspace(0.270, 0.360, 14)[1:]])  # aresta em 270 mm
    saida = DADOS / "retorno_if_bt.csv"
    with open(saida, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter=";")
        w.writerow(["SR", "lado", "r_min_m", "r_max_m", "mdot_kg_s"])
        for s in SRS_CAMPOS:
            for lado in ("bomba", "turbina"):
                arq = BACKUP_CSV / s / "dados" / f"{s}_euler_if_bomba_turbina_lado_{lado}.csv"
                a = np.genfromtxt(arq, delimiter=",", skip_header=1)
                m = a[:, 7] * np.einsum("ij,ij->i", a[:, 4:7], a[:, 8:11])
                if lado == "turbina":
                    m = -m
                r = a[:, 14]
                for k in range(len(bordas) - 1):
                    sel = (r >= bordas[k]) & (r < bordas[k + 1]) if k < len(bordas) - 2 else (r >= bordas[k])
                    w.writerow([int(s[2:]) / 100, lado, f"{bordas[k]:.5f}", f"{bordas[k + 1]:.5f}", f"{m[sel].sum():.4f}"])
    print("gravado", saida)
    return saida


def fig_retorno():
    """(a) Vazão líquida por faixa de raio na if_bt para alguns SR; (b) vazão de retorno junto ao núcleo e vazão no
    circuito em função de SR."""
    arq = DADOS / "retorno_if_bt.csv"
    if not arq.exists():
        calcular_retorno()
    ls = ler_csv(arq)
    d = desempenho()
    srs = sorted({float(l["SR"]) for l in ls})
    tab = {}
    for l in ls:
        tab.setdefault(float(l["SR"]), {}).setdefault((float(l["r_min_m"]), float(l["r_max_m"])), []).append(
            float(l["mdot_kg_s"]))
    fig, axs = plt.subplots(1, 2, figsize=(7.2, 2.8))
    ax = axs[0]
    cores = {0.0: "#08306b", 0.4: "#2171b5", 0.9: "#6baed6", 1.0: "#fd8d3c"}
    for sr, cor in cores.items():
        fx = sorted(tab[sr])
        rc = np.array([(a + b) / 2 for a, b in fx])
        dr = np.array([b - a for a, b in fx])
        m = np.array([np.mean(tab[sr][k]) for k in fx])
        ax.step(rc * 1e3, m / dr / 1e3, where="mid", color=cor, lw=1.2, label=f"SR = {sr:g}".replace(".", ","))
    ax.axhline(0, color="0.3", lw=0.8)
    ax.axvline(270, color="0.5", ls=":", lw=0.9)
    ax.set_xlabel(r"$r$ [mm]")
    ax.set_ylabel(r"$d\dot{m}/dr$ [t/(s·m)]")
    ax.set_title("(a) Vazão líquida na interface bomba–turbina")
    ax.legend(fontsize=7)
    ax = axs[1]
    ret, ext = [], []
    for sr in srs:
        fx = tab[sr]
        ret.append(sum(np.mean(v) for (a, b), v in fx.items() if b <= 0.2701))
        ext.append(sum(np.mean(v) for (a, b), v in fx.items() if b > 0.2701))
    ax.plot(srs, ext, "o-", ms=3, lw=1, color="#d62728", label=r"$r > 270$ mm (bomba $\rightarrow$ turbina)")
    ax.plot(srs, -np.array(ret), "s-", ms=3, lw=1, color="#1f77b4", label=r"$r < 270$ mm (turbina $\rightarrow$ bomba)")
    ax.plot(d["SR"], d["mdot_bt"], "^-", ms=3, lw=1, color="k", label="vazão no circuito (solver)")
    ax.axhline(0, color="0.3", lw=0.8)
    ax.set_xlabel("SR")
    ax.set_ylabel(r"$\dot{m}$ [kg/s]")
    ax.set_title("(b) Vazões em função de SR")
    ax.legend(fontsize=6.5)
    ax.set_xlim(-0.03, 1.03)
    for ax in axs:
        virgula(ax)
    fig.tight_layout()
    salvar(fig, "conversor_res_retorno")


def ler_historico(caminho, fator=1.0, it_max=None):
    """Lê um report file do Fluent e devolve iteração, T_b, T_t e T_e.

    Args:
        caminho (Path): Arquivo ``.out`` (três linhas de cabeçalho, colunas Iteration, t_bomba, t_turbina,
            t_estator, ...).
        fator (float): Fator de conversão dos torques para N·m (0,6125 para os arquivos em Cm).
        it_max (int | None): Última iteração mantida; ``None`` mantém todas.

    Returns:
        tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]: Iteração, T_b, T_t e T_e [N·m].
    """
    a = np.array([[float(x) for x in l.split()[:4]] for l in caminho.read_text().splitlines()[3:] if l.strip()])
    if it_max is not None:
        a = a[a[:, 0] <= it_max]
    return a[:, 0], a[:, 1] * fator, a[:, 2] * fator, a[:, 3] * fator


def fig_sens_historico():
    """Históricos de TR e |T_b| das variantes de sensibilidade, cada uma ao lado da sua configuração base.

    As variantes que releram o caso base (SST sem correção de curvatura, k-ε, FOU, Coupled) continuam a contagem
    de iterações a partir do fim da base; a de meio passo partiu do repouso. A linha horizontal de cada série é a
    média na janela final de 200 iterações.
    """
    sens = CONV / "CT_V1_MRF" / "sens"
    base_g0 = ler_historico(DADOS / "CT_V1_MRF_rep_20261121-022556.out", CM_PARA_NM, 600)
    base_g2 = ler_historico(CONV / "CT_V1_GCI_G2" / "dados" / "CT_V1_MRF_rep_20261122-000054.out", 1.0, 1000)
    var = {v: ler_historico(sens / v / "dados" / f"sens_{v}_rep.out")
           for v in ("fou", "ke_realizable", "sst_sem_cc", "meio_passo")}
    var["coupled"] = ler_historico(CONV / "CT_V1_GCI_G2" / "sens" / "coupled" / "dados" / "sens_coupled_rep.out")
    # (rótulo, histórico, cor, espessura, primeira iteração mostrada); a ordem define quem fica por cima
    grupos = [("(a) Acoplamento (G2)", [("SIMPLE", base_g2, "k", 0.7, 300), ("Coupled", var["coupled"], "#d62728",
                                                                             0.7, 1010)]),
              ("(b) Esquema convectivo", [("SOU (base)", base_g0, "k", 0.7, 150), ("FOU", var["fou"], "#d62728",
                                                                                   0.7, 610)]),
              ("(c) Modelo de turbulência", [(r"$k$-$\varepsilon$ Realizable", var["ke_realizable"], "#d62728", 0.5,
                                              610),
                                             ("SST (base)", base_g0, "k", 0.7, 150),
                                             ("SST sem corr. curv.", var["sst_sem_cc"], "#1f77b4", 0.9, 610)]),
              ("(d) Posição das rodas", [("Original (base)", base_g0, "k", 0.7, 150), ("Meio passo", var["meio_passo"],
                                                                                      "#d62728", 0.7, 150)])]
    fig, axs = plt.subplots(2, 4, figsize=(7.6, 4.6), sharex="col")
    for j, (titulo, series) in enumerate(grupos):
        for rot, (it, tb, tt, _), cor, lw, it0 in series:
            ok = it >= it0
            tr = -tt / tb
            jan = it > it[-1] - 200
            for i, y in enumerate((tr, -tb / 1e3)):
                axs[i, j].plot(it[ok], y[ok], color=cor, lw=lw, label=rot)
                axs[i, j].hlines(y[jan].mean(), it[jan][0], it[-1], color=cor, lw=1.8)
        axs[0, j].set_title(titulo)
        axs[1, j].legend(loc="upper right" if j == 0 else "lower left", fontsize=6.3, handlelength=1.2)
        axs[1, j].set_xlabel("Iteração")
    for i in range(2):
        for j, ax in enumerate(axs[i]):
            ax.set_ylim(*((0.9, 1.3) if (i, j) == (0, 2) else (0.98, 1.16) if i == 0 else (18, 32)))
            virgula(ax)
    axs[0, 0].set_ylabel("TR")
    axs[1, 0].set_ylabel(r"$|T_b|$ [kN·m]")
    fig.tight_layout()
    salvar(fig, "conversor_res_sens_historico")


FIGURAS = {"conv": fig_convergencia, "des": fig_desempenho, "tor": fig_torques, "euler": fig_euler,
           "yplus": fig_yplus, "lit": fig_literatura, "sens": fig_sensibilidade, "gci": fig_gci, "ret": fig_retorno,
           "senshist": fig_sens_historico}

if __name__ == "__main__":
    for nome in (sys.argv[1:] or FIGURAS):
        FIGURAS[nome]()
