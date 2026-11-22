"""Derivada das curvas de desempenho do CT_V1 (MRF) contra a literatura, para validar a TENDÊNCIA.

Uso: python tools/tendencia_derivada_conversor.py

Calcula dTR/dSR, dη/dSR e d(K/K_estol)/dSR do CT_V1 e de cada série da literatura (``refs/conversor_literatura.csv``,
lida por import de ``comparacao_literatura_conversor``) e grava:

* ``results/conversor/CT_V1_MRF/graficos/g8_derivada_tendencia.png`` (3 painéis, barras de erro no CT_V1);
* ``results/conversor/CT_V1_MRF/dados/derivada_tendencia.csv`` (separador ``;``, vírgula decimal).

Método: diferença finita pelo PONTO MÉDIO de cada intervalo, d = (y_{i+1} - y_i)/(SR_{i+1} - SR_i), plotada em
SR_{i+1/2}. Foi escolhida em vez da diferença centrada nos pontos internos porque a malha de SR é não uniforme
(passos de 0,2 a 0,05): a fórmula centrada de 3 pontos em malha irregular perde precisão e não dá derivada nas
pontas, enquanto o quociente do intervalo usa todos os intervalos e tem propagação de incerteza simples.

Incerteza (CT_V1): U_i = 2·SEM_i (k = 2, ~95 %), com o SEM de TR, η e K por batch means das séries do report file
(``gci_conversor.serie_janela`` + ``grandezas_ponto``). Para o intervalo: U_d = sqrt(U_i² + U_{i+1}²)/ΔSR (pontos
independentes). Isso cobre só o ruído iterativo da janela, não o erro de discretização nem o de modelo (sistemáticos
e em boa parte comuns aos pontos vizinhos). Para K/K_estol divide-se por K(SR=0), tratado como fator exato (escala
comum). A literatura é digitalizada e não tem barras (a incerteza de leitura está em ``refs/conversor_literatura.csv``).
"""
import csv
import json
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
import comparacao_literatura_conversor as cl  # noqa: E402
import gci_conversor as gci  # noqa: E402

ROOT = cl.ROOT
BASE = ROOT / "results" / "conversor" / "CT_V1_MRF"
SAIDA_PNG = BASE / "graficos" / "g8_derivada_tendencia.png"
SAIDA_CSV = BASE / "dados" / "derivada_tendencia.csv"
K_SIGMA = 2.0  # fator de abrangência (~95 %)
ROTULO = {"TR": "dTR/dSR", "eta": "dη/dSR", "KK": "d(K/K_estol)/dSR"}


def derivada(sr, y, u=None):
    """Derivada por diferença finita no ponto médio de cada intervalo, com a incerteza propagada.

    Args:
        sr (np.ndarray): Razões de velocidade (crescentes).
        y (np.ndarray): Valores da grandeza (pontos com NaN são descartados).
        u (np.ndarray | None): Incerteza expandida de cada ``y`` (mesma unidade) ou ``None``.

    Returns:
        tuple[np.ndarray, np.ndarray, np.ndarray]: ``(sr_medio, dy_dsr, u_deriv)``; ``u_deriv`` é NaN sem ``u``.
    """
    sr, y = np.asarray(sr, float), np.asarray(y, float)
    uu = np.full_like(y, np.nan) if u is None else np.asarray(u, float)
    ok = ~np.isnan(y) & ~np.isnan(sr)
    sr, y, uu = sr[ok], y[ok], uu[ok]
    ds = np.diff(sr)
    return 0.5 * (sr[:-1] + sr[1:]), np.diff(y) / ds, np.hypot(uu[:-1], uu[1:]) / ds


def incertezas_ct():
    """Lê o SEM de TR, η e K de cada ponto do CT_V1 pelas séries dos report files.

    Returns:
        tuple[dict[float, dict], list[str]]: ``{SR: {"TR","eta","K": sem}}`` e avisos dos pontos sem série
        (nesses o SEM cai no σ por iteração, conservador).
    """
    sem, avisos = {}, []
    for pasta in sorted(BASE.glob("SR???")):
        r = json.loads((pasta / "dados" / "resultado.json").read_text(encoding="utf-8"))
        sr = float(r.get("sr", int(pasta.name[2:]) / 100))
        serie, info = gci.serie_janela("G0", r, sr)
        g = gci.grandezas_ponto(r, serie)
        if serie is None:
            avisos.append(f"SR={sr:g}: sem série ({info}); SEM = σ por iteração (conservador)")
        sem[round(sr, 4)] = {q: g[q]["sem"] for q in ("TR", "eta", "K")}
    return sem, avisos


def fmt(x, casas=4):
    """Formata número com vírgula decimal.

    Args:
        x (float): Valor.
        casas (int): Casas decimais.

    Returns:
        str: Texto (vazio para NaN).
    """
    return "" if x is None or np.isnan(x) else f"{x:.{casas}f}".replace(".", ",")


def gravar_csv(der_ct, dk, uk, lit):
    """Grava a tabela de derivadas (CT_V1 e literatura).

    Args:
        der_ct (dict): ``{"TR"|"eta"|"KK": (sr_medio, d, u)}`` do CT_V1.
        dk (np.ndarray): dK/dSR do CT_V1 (K em rpm/√(N·m)).
        uk (np.ndarray): Incerteza expandida de ``dk``.
        lit (list[dict]): Séries da literatura (``ler_literatura``).

    Returns:
        None
    """
    cab = ["fonte", "origem", "SR_medio", "dTR_dSR", "U_dTR_dSR", "deta_dSR", "U_deta_dSR",
           "d(K/Kestol)_dSR", "U_d(K/Kestol)_dSR", "dK_dSR", "U_dK_dSR"]
    with open(SAIDA_CSV, "w", encoding="utf-8-sig", newline="") as f:
        w = csv.writer(f, delimiter=";")
        w.writerow(cab)
        for i, sm in enumerate(der_ct["TR"][0]):
            w.writerow(["CT_V1", "mrf", fmt(sm), fmt(der_ct["TR"][1][i]), fmt(der_ct["TR"][2][i]),
                        fmt(der_ct["eta"][1][i]), fmt(der_ct["eta"][2][i]), fmt(der_ct["KK"][1][i]),
                        fmt(der_ct["KK"][2][i]), fmt(dk[i]), fmt(uk[i])])
        for s in lit:
            o = {q: derivada(s["SR"], s[q]) for q in ("TR", "eta", "KK")}
            for sm in o["TR"][0]:
                vals = []
                for q in ("TR", "eta", "KK"):
                    k = np.where(np.isclose(o[q][0], sm))[0]
                    vals += [fmt(o[q][1][k[0]]) if len(k) else "", ""]
                w.writerow([s["fonte"], s["origem"], fmt(sm)] + vals + ["", ""])


def figura(der_ct, lit):
    """Desenha a figura g8 (3 painéis: dTR/dSR, dη/dSR, d(K/K_estol)/dSR).

    Args:
        der_ct (dict): Derivadas do CT_V1 com incerteza.
        lit (list[dict]): Séries da literatura.

    Returns:
        None
    """
    fig, axs = plt.subplots(1, 3, figsize=(15, 5.6), sharex=True)
    for ax, q in zip(axs, ("TR", "eta", "KK")):
        for s in lit:
            m, d, _ = derivada(s["SR"], s[q])
            if len(m):
                ax.plot(m, d, color=cl.CORES.get(s["fonte"], "gray"), alpha=0.8, label=s["rotulo"],
                        **{**cl.ESTILO[s["origem"]], "lw": 0.8, "ms": 2.5})
        m, d, ud = der_ct[q]
        ax.errorbar(m, d, yerr=ud, color="k", marker="D", ms=5, lw=1.8, capsize=3, zorder=5,
                    label="CT_V1 (MRF), barras = 2·SEM")
        ax.axhline(0, color="gray", lw=1)
        ax.set_xlabel("SR (ponto médio do intervalo)")
        ax.set_ylabel(ROTULO[q])
        ax.set_title(ROTULO[q])
        ax.grid(alpha=0.3)
    axs[2].set_ylim(-1, 4)  # corta o disparo do nadeem2020 perto de SR = 0,95 (K/K_estol ~ 13, ruído de digitalização)
    axs[2].text(0.98, 0.97, "eixo cortado em 4 (nadeem2020 sobe a 12,8 perto de SR = 0,95)", transform=axs[2].transAxes,
                ha="right", va="top", fontsize=7)
    h, l = axs[0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=4, fontsize=7.5, frameon=False, bbox_to_anchor=(0.5, 0.05))
    fig.suptitle("Derivada das curvas de desempenho: CT_V1 (MRF) × literatura", fontsize=12)
    fig.text(0.5, 0.012, "A derivada confirma o sinal da tendência, não a forma: a forma depende da geometria (a do CT_V1 é "
             "dada e difere da literatura).\nDiferenças finitas nos pontos médios; barras do CT_V1: ruído iterativo "
             "(batch means, k = 2).", ha="center", fontsize=8)
    fig.tight_layout(rect=(0, 0.15, 1, 0.96))
    fig.savefig(SAIDA_PNG, dpi=150)
    plt.close(fig)


def main():
    """Calcula as derivadas, grava tabela e figura e imprime o resumo.

    Returns:
        None
    """
    ct = cl.ler_ct()
    sem, avisos = incertezas_ct()
    sr = ct["SR"]
    u = {q: np.array([K_SIGMA * sem[round(s, 4)][q] for s in sr]) for q in ("TR", "eta", "K")}
    k0 = ct["K"][0]
    y = {"TR": ct["TR"], "eta": ct["eta"], "KK": ct["K"] / k0}
    uy = {"TR": u["TR"], "eta": u["eta"], "KK": u["K"] / k0}
    der = {q: derivada(sr, y[q], uy[q]) for q in y}
    _, dk, uk = derivada(sr, ct["K"], u["K"])
    lit = cl.ler_literatura()
    gravar_csv(der, dk, uk, lit)
    figura(der, lit)

    print("Avisos:", avisos or "nenhum (SEM por batch means em todos os pontos)")
    print("U(TR) por ponto:", np.round(u["TR"], 4), "\nU(eta):", np.round(u["eta"], 4), "\nU(K):", np.round(u["K"], 3))
    for q in ("TR", "eta", "KK"):
        m, d, ud = der[q]
        incerto = [f"{a:g}-{b:g}" for a, b, di, ui in zip(sr[:-1], sr[1:], d, ud) if abs(di) <= ui]
        print(f"CT_V1 {ROTULO[q]}: min {d.min():.3f}, max {d.max():.3f}; sinal incerto em: {incerto or 'nenhum'}")
        print("   por intervalo:", np.round(d, 3))
        for s in lit:
            dd = derivada(s["SR"], s[q])[1]
            if len(dd):
                print(f"   {s['rotulo']}: [{dd.min():.3f}; {dd.max():.3f}]")
    print("tabela:", SAIDA_CSV, "\nfigura:", SAIDA_PNG)


if __name__ == "__main__":
    main()
