"""Compara as curvas adimensionais do CT_V1 (MRF) com conversores da literatura (tendência, não valor).

Uso: python tools/comparacao_literatura_conversor.py

Lê ``results/conversor/CT_V1_MRF/dados/desempenho.csv`` (gerado por ``analise_conversor.py``) e
``refs/conversor_literatura.csv`` (curvas digitalizadas; fontes e limitações em ``refs/conversor_literatura.md``).
Grava ``results/conversor/CT_V1_MRF/graficos/g5_comparacao_literatura.png`` com quatro painéis:
TR × SR, TR/TR_estol × SR, η × SR e K/K_estol × SR.

As geometrias são outras (o CT_V1 é dado), então o gráfico mostra se o setup reproduz as tendências
esperadas: TR decrescente, η quase linear até o acoplamento, K crescente.
"""
import csv, re
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
DES = ROOT / "results" / "conversor" / "CT_V1_MRF" / "dados" / "desempenho.csv"
LIT = ROOT / "refs" / "conversor_literatura.csv"
SAIDA = ROOT / "results" / "conversor" / "CT_V1_MRF" / "graficos" / "g5_comparacao_literatura.png"
ESTILO = {"exp": dict(marker="o", ls="-", lw=1.0, ms=3.5), "cfd_rans": dict(marker="^", ls="--", lw=0.9, ms=3.5),
          "cfd_les": dict(marker="s", ls=":", lw=0.9, ms=3.5), "modelo": dict(marker=None, ls="-.", lw=1.0)}
CORES = {"zhang2023": "#9467bd", "guo2021": "#8c564b", "ran2021": "#e377c2", "nadeem2020": "#7f7f7f",
         "kotwicki_forma_nadeem2020": "#bcbd22"}


def ler_ct():
    """Lê a tabela de desempenho do CT_V1.

    Returns:
        dict[str, np.ndarray]: ``SR``, ``TR``, ``eta``, ``K`` e ``dp_TR`` ordenados por SR.

    Raises:
        FileNotFoundError: Se a tabela ainda não existir.
    """
    with open(DES, encoding="utf-8") as f:
        linhas = list(csv.DictReader(f, delimiter=";"))
    linhas.sort(key=lambda l: float(l["SR"]))
    return {k: np.array([float(l[k]) for l in linhas]) for k in ("SR", "TR", "eta", "K", "dp_TR")}


def ler_literatura():
    """Lê as curvas digitalizadas e agrupa por série (fonte, origem, figura).

    O K/K_estol vem do campo ``obs`` (``K/K_estol(derivado)=x``) quando o artigo dá o fator de
    capacidade adimensional; senão, é K dividido pelo K de estol da série (ou do experimento da mesma
    fonte, no caso do modelo, que não vale no estol).

    Returns:
        list[dict]: Séries com ``fonte``, ``origem``, ``rotulo`` e arrays ``SR``, ``TR``, ``eta``, ``KK``.
    """
    with open(LIT, encoding="utf-8") as f:
        linhas = list(csv.DictReader(f, delimiter=";"))

    def num(s):
        """Converte texto em float.

        Args:
            s (str): Valor do CSV.

        Returns:
            float: Valor, ou NaN se vazio ou inválido.
        """
        try:
            return float(s)
        except (TypeError, ValueError):
            return float("nan")

    series = {}
    for l in linhas:
        chave = (l["fonte"], l["origem"], l["figura_ou_tabela"])
        series.setdefault(chave, []).append(l)
    k_estol_exp = {}
    for (fonte, origem, _), ls in series.items():
        if origem == "exp":
            ls0 = min(ls, key=lambda l: num(l["SR"]))
            k_estol_exp[fonte] = num(ls0["K"])
    out = []
    for (fonte, origem, fig), ls in series.items():
        ls.sort(key=lambda l: num(l["SR"]))
        sr = np.array([num(l["SR"]) for l in ls])
        tr = np.array([num(l["TR"]) for l in ls])
        eta = np.array([num(l["eta"]) for l in ls])
        kk = []
        for l in ls:
            m = re.search(r"K/K_estol\(derivado\)=([0-9.]+)", l["obs"] or "")
            if m:
                kk.append(float(m.group(1)))
            elif "n_bomba" in (l["unidade_K"] or ""):
                base = k_estol_exp.get(fonte.replace("kotwicki_forma_", "")) if origem == "modelo" else num(ls[0]["K"])
                kk.append(num(l["K"]) / base if base else float("nan"))
            else:
                kk.append(float("nan"))
        eta = np.where(np.isnan(eta), sr * tr, eta) if origem != "cfd_les" else eta
        sub = re.search(r"\(([^)]*)\)", fig)
        rot = f"{fonte} {origem}" + (f" ({sub.group(1)[:22]})" if sub and origem != "exp" else "")
        out.append({"fonte": fonte, "origem": origem, "rotulo": rot, "SR": sr, "TR": tr, "eta": eta, "KK": np.array(kk)})
    return out


def acoplamento(sr, tr):
    """SR em que TR cruza 1 por interpolação linear.

    Args:
        sr (np.ndarray): Razões de velocidade crescentes.
        tr (np.ndarray): Razões de torque.

    Returns:
        float: SR do cruzamento, ou NaN se TR não cruzar 1 no intervalo.
    """
    for i in range(len(sr) - 1):
        if (tr[i] - 1) * (tr[i + 1] - 1) <= 0 and tr[i] != tr[i + 1]:
            return sr[i] + (1 - tr[i]) * (sr[i + 1] - sr[i]) / (tr[i + 1] - tr[i])
    return float("nan")


def main():
    """Gera a figura de comparação e imprime o resumo (TR de estol e SR de acoplamento por série)."""
    ct = ler_ct()
    lit = ler_literatura()
    fig, ax = plt.subplots(2, 2, figsize=(12, 8.5))
    a_tr, a_trn, a_eta, a_k = ax[0, 0], ax[0, 1], ax[1, 0], ax[1, 1]
    for s in lit:
        st = dict(ESTILO.get(s["origem"], {}), color=CORES.get(s["fonte"], "k"), label=s["rotulo"], alpha=0.85)
        ok = ~np.isnan(s["TR"])
        a_tr.plot(s["SR"][ok], s["TR"][ok], **st)
        if s["origem"] != "modelo" and ok.any():
            a_trn.plot(s["SR"][ok], s["TR"][ok] / s["TR"][ok][0], **st)
        ok = ~np.isnan(s["eta"])
        a_eta.plot(s["SR"][ok], s["eta"][ok], **st)
        ok = ~np.isnan(s["KK"])
        if ok.any():
            a_k.plot(s["SR"][ok], s["KK"][ok], **st)
        print(f"{s['rotulo']:<55} TR(SR mín)={s['TR'][0]:.3f}  acoplamento SR≈{acoplamento(s['SR'], s['TR']):.3f}")
    est = dict(color="k", marker="D", ms=5, lw=2.2, label="CT_V1 MRF (este trabalho)", zorder=5)
    a_tr.errorbar(ct["SR"], ct["TR"], yerr=ct["dp_TR"], capsize=3, **est)
    a_trn.plot(ct["SR"], ct["TR"] / ct["TR"][0], **est)
    a_eta.plot(ct["SR"], ct["eta"], **est)
    a_k.plot(ct["SR"], ct["K"] / ct["K"][0], **est)
    print(f"{'CT_V1 MRF':<55} TR(0)={ct['TR'][0]:.3f}  acoplamento SR≈{acoplamento(ct['SR'], ct['TR']):.3f}")
    a_tr.axhline(1, color="gray", lw=0.8, ls=":")
    a_tr.set(xlabel="SR = ω_t/ω_b", ylabel="TR = −T_t/T_b", title="Razão de torque")
    a_trn.set(xlabel="SR", ylabel="TR / TR_estol", title="Razão de torque normalizada pelo estol (forma)")
    a_eta.set(xlabel="SR", ylabel="η = SR·TR", title="Eficiência")
    a_k.set(xlabel="SR", ylabel="K / K_estol", title="Fator de capacidade normalizado (K = n_b/√T_b)")
    for a in ax.flat:
        a.set_xlim(-0.02, 1.0)
    h, l = a_tr.get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=3, fontsize=7, frameon=False)
    fig.suptitle("CT_V1 (MRF, k-ω SST) × conversores da literatura — comparação de tendência; geometrias diferentes",
                 fontsize=10)
    fig.tight_layout(rect=(0, 0.12, 1, 0.97))
    SAIDA.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(SAIDA, dpi=140)
    plt.close(fig)
    print("gravado:", SAIDA)


if __name__ == "__main__":
    main()
