"""Figuras-resumo de comparação com a referência, uma por bloco (A–H), para as rodadas R01–R43.

Cada figura sobrepõe todas as rodadas do bloco e as compara com a referência disponível em ``refs/``,
no código (``analise_caso``/``casos``) ou nos docs. Nada é inventado: sem referência, a figura mostra só as rodadas.
Não roda o Fluent: lê os ``_campo.csv``, ``_parede_*.csv`` e ``_metricas.json`` já gravados.

Uso:
    python tools/comparacao_blocos.py            (todos os blocos)
    python tools/comparacao_blocos.py A C H      (só os blocos indicados; H gera duto e degrau)

Saída em ``results/comparacao_blocos/``:
    bloco_A.png, bloco_B.png, bloco_C.png, bloco_D.png, bloco_E.png, bloco_F.png, bloco_G.png,
    bloco_H_duto.png, bloco_H_degrau.png e resumo.csv (``;``, uma linha por rodada).

Rodadas com ``Convergiu_S_N = N`` em ``00_admin/registro_pyfluent.csv`` aparecem com ``*`` no rótulo,
hachura nas barras e linha tracejada nas curvas.
"""
import json
import math
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from casos import ROOT, CASES, caso  # noqa: E402
from analise_caso import (ler_campo, arquivos, coluna, perto, linha_central, comparar_ghia,  # noqa: E402
                          paredes_bft, u_analitica, GHIA, DRIVER, XR_H_EXP, XR_H_INC)

SAIDA = ROOT / "results" / "comparacao_blocos"
DPI = 140
DP_REF = 20.0           # Pa, G·L = 4 Pa/m · 5 m
UB_DNS = 14.7           # m/s, DNS de Gavrilakis (1992), Re_b ≈ 4410 (casos.py, tipo "duto")
CF_PIROZZOLI = 9.38e-3  # Pirozzoli et al. (2018), Tabela 1 caso A (refs/duto_secundario.csv)
VS_UB_FAIXA = (1.0, 2.0)  # % de u_b, Pirozzoli et al. (2018), resumo (refs/duto_secundario.csv)
H_DEGRAU = 0.5          # m, altura do degrau laminar (bloco E)
RESUMO = []             # linhas do resumo.csv, preenchidas pelos blocos
NAO_CONV = set()        # IDs com Convergiu_S_N = N


# ---------------------------------------------------------------- utilidades

def avisar(msg):
    """Imprime um aviso no console (item pulado ou dado ausente).

    Args:
        msg (str): Texto do aviso.
    """
    print(f"AVISO: {msg}")


def virg(v, nd=2):
    """Formata número com vírgula decimal (pt-BR).

    Args:
        v (float): Valor.
        nd (int, optional): Casas decimais. Padrão 2.

    Returns:
        str: Número formatado, por exemplo ``"4,35"``.
    """
    return f"{v:.{nd}f}".replace(".", ",")


def carregar_registro():
    """Lê ``00_admin/registro_pyfluent.csv`` e preenche ``NAO_CONV`` com as rodadas não convergidas.

    Returns:
        pd.DataFrame: Registro indexado por ID.
    """
    reg = pd.read_csv(ROOT / "00_admin" / "registro_pyfluent.csv", sep=";", encoding="utf-8-sig")
    reg = reg.set_index("ID")
    NAO_CONV.clear()
    NAO_CONV.update(reg.index[reg["Convergiu_S_N"].astype(str).str.strip() == "N"])
    return reg


def rotulo(cid, texto=None):
    """Rótulo de uma rodada, com ``*`` se não convergiu.

    Args:
        cid (str): ID da rodada.
        texto (str, optional): Texto após o ID (por exemplo, o modelo).

    Returns:
        str: Por exemplo ``"R37* BSL"``.
    """
    s = cid + ("*" if cid in NAO_CONV else "")
    return f"{s} {texto}" if texto else s


def metricas(cid):
    """Lê o ``_metricas.json`` de uma rodada.

    Args:
        cid (str): ID da rodada.

    Returns:
        dict | None: Métricas, ou None (com aviso) se o arquivo não existir.
    """
    p = arquivos(caso(cid))["met"]
    if not p.exists():
        avisar(f"{cid}: sem {p.name}; rodada pulada")
        return None
    return json.loads(p.read_text(encoding="utf-8"))


def campo(cid):
    """Lê o campo de uma rodada.

    Args:
        cid (str): ID da rodada.

    Returns:
        pd.DataFrame | None: Campo, ou None (com aviso) se não existir.
    """
    p = arquivos(caso(cid))["campo"]
    if not p.exists():
        avisar(f"{cid}: sem campo CSV; rodada pulada")
        return None
    return ler_campo(p)


def linha(bloco, cid, metrica, valor, ref, desvio, nota=""):
    """Acrescenta uma linha ao resumo.

    Args:
        bloco (str): Letra do bloco.
        cid (str): ID da rodada.
        metrica (str): Nome e unidade da métrica principal.
        valor (float | None): Valor obtido.
        ref (str): Referência (valor e fonte).
        desvio (str): Desvio em relação à referência (com unidade).
        nota (str, optional): Observação.
    """
    RESUMO.append(dict(bloco=bloco, rodada=cid, metrica_principal=metrica,
                       valor="" if valor is None else f"{valor:.6g}", referencia=ref, desvio=desvio,
                       convergiu="N" if cid in NAO_CONV else "S", nota=nota))


def salvar(fig, nome):
    """Grava a figura em ``SAIDA`` e a fecha.

    Args:
        fig (matplotlib.figure.Figure): Figura.
        nome (str): Nome do arquivo PNG.
    """
    SAIDA.mkdir(parents=True, exist_ok=True)
    fig.savefig(SAIDA / nome, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"gerado: {SAIDA / nome}")


def barras(ax, ids, valores, rotulos, cores=None):
    """Barras por rodada, com hachura nas não convergidas e o valor escrito em cima.

    Args:
        ax (matplotlib.axes.Axes): Eixo.
        ids (list[str]): IDs das rodadas (definem a hachura).
        valores (list[float]): Alturas.
        rotulos (list[str]): Rótulos do eixo x.
        cores (list | None, optional): Cor por barra. Padrão: cores do ciclo.
    """
    cores = cores or [f"C{i}" for i in range(len(ids))]
    xs = np.arange(len(ids))
    for x, cid, v, c in zip(xs, ids, valores, cores):
        ax.bar(x, v, color=c, edgecolor="k", hatch="//" if cid in NAO_CONV else None, alpha=0.85)
    ax.set_xticks(xs)
    ax.set_xticklabels(rotulos, rotation=30, ha="right", fontsize=8)
    ax.set_ylim(0, 1.28 * max(max(valores), 1e-30))
    ax.grid(axis="y", alpha=0.3)


def rodape(fig, texto):
    """Escreve uma nota de rodapé na figura.

    Args:
        fig (matplotlib.figure.Figure): Figura.
        texto (str): Nota.
    """
    fig.text(0.01, -0.09, texto, fontsize=8, va="top", ha="left", style="italic")


def perfil_cp(cid):
    """Perfil u(d) do canal na coluna de células indicada por ``x0`` nas métricas.

    Args:
        cid (str): ID da rodada (blocos A, B, F).

    Returns:
        tuple[np.ndarray, np.ndarray] | None: d = −y [m] e u [m/s], ou None se faltar dado.
    """
    m, df = metricas(cid), campo(cid)
    if m is None or df is None:
        return None
    _, col = coluna(df, "x-coordinate", m["x0"])
    return -col["y-coordinate"].to_numpy(), col["x-velocity"].to_numpy()


# ---------------------------------------------------------------- bloco A

def figura_a():
    """Bloco A (R01–R05): perfis u(d) por malha contra a analítica e e_L2 × h em log-log."""
    ids = ["R01", "R02", "R03", "R04", "R05"]
    malhas = {"R01": "M1 25×10", "R02": "M2 50×20", "R03": "M3 100×40", "R04": "M4 200×80", "R05": "M2, FOU"}
    fig, ax = plt.subplots(1, 3, figsize=(17, 5.2), gridspec_kw=dict(width_ratios=[1.1, 1.1, 1]))
    dd = np.linspace(0, 1, 200)
    ax[0].plot(u_analitica(dd), dd, "k-", lw=1.5, label="Analítico: u = 2d + 20d(1−d)")
    ax[1].plot(dd, np.zeros_like(dd), "k-", lw=1.0, label="Analítico (erro = 0)")
    mar = dict(R01="o", R02="s", R03="^", R04="v")
    h, e = [], []
    for cid in ids:
        m, p = metricas(cid), perfil_cp(cid)
        if m is None or p is None:
            continue
        d, u = p
        ny = math.sqrt(m["n_celulas"] / 2.5)
        if cid == "R05":
            ax[0].plot(u, d, "*", color="tab:red", ms=13, mfc="none", mew=1.6, label=rotulo(cid, malhas[cid]) + " (destaque)")
            ax[1].plot(d, u - u_analitica(d), "*--", color="tab:red", ms=11, mfc="none", lw=1.2,
                       label=rotulo(cid, malhas[cid]) + " (destaque)")
            ax[2].plot(1 / 20, m["eL2"], "*", color="tab:red", ms=16, mfc="none", mew=1.8, label=rotulo(cid, "FOU (M2)"))
            ref = metricas("R02")
            linha("A", cid, "e_L2 [%]", m["eL2"], "0 (analítico u = 2d+20d(1−d))", f"{virg(m['eL2'], 4)} % (= R02; max|Δu| vs R02 = {m.get('max_du_R02', float('nan')):.1e} m/s)",
                  "FOU igual ao SOU em escoamento desenvolvido")
            continue
        ax[0].plot(u, d, mar[cid], ms=5, mfc="none", label=rotulo(cid, malhas[cid]))
        ax[1].plot(d, u - u_analitica(d), mar[cid] + "-", ms=4, lw=0.8, mfc="none", label=rotulo(cid, malhas[cid]))
        h.append(1 / ny)
        e.append(m["eL2"])
        ax[2].plot(1 / ny, m["eL2"], mar[cid], ms=9, mfc="none", label=rotulo(cid, malhas[cid]))
        linha("A", cid, "e_L2 [%]", m["eL2"], "0 (analítico u = 2d+20d(1−d))", f"{virg(m['eL2'], 4)} %", f"h = 1/{ny:.0f} m")
    ax[0].set(xlabel="u [m/s]", ylabel="d = distância à parede parada [m]", title="Perfil de velocidade em x = 2,5 m")
    ax[1].set(xlabel="d [m]", ylabel="u(Fluent) − u(analítico) [m/s]", title="Erro do perfil em relação ao analítico")
    if len(h) >= 2:
        h, e = np.array(h), np.array(e)
        p_obs = np.polyfit(np.log(h), np.log(e), 1)[0]
        par = [math.log(e[i] / e[i + 1]) / math.log(h[i] / h[i + 1]) for i in range(len(h) - 1)]
        ax[2].plot(h, e, "-", color="gray", lw=0.8)
        hh = np.array([h.min() * 0.8, h.max() * 1.2])
        ax[2].plot(hh, e[0] * (hh / h[0]) ** 2, "k--", label="Inclinação 2 (ideal)")
        ax[2].text(0.04, 0.95, f"Ordem observada p = {virg(p_obs)} (ajuste)\n" +
                   "\n".join(f"R0{i + 1}→R0{i + 2}: p = {virg(q)}" for i, q in enumerate(par)),
                   transform=ax[2].transAxes, va="top", fontsize=9, bbox=dict(fc="w", ec="gray", alpha=0.9))
    ax[2].set(xscale="log", yscale="log", xlabel="h = 1/N_y [m]", ylabel="e_L2 relativo [%]",
              title="Convergência de malha (R05 em destaque)")
    for a_ in ax:
        a_.grid(alpha=0.3, which="both")
        a_.legend(fontsize=8)
    fig.suptitle("Bloco A (R01–R05): Couette-Poiseuille periódico contra a solução analítica", fontsize=13)
    salvar(fig, "bloco_A.png")


# ---------------------------------------------------------------- blocos B e F

def referencia_dp_entrada_uniforme(m):
    """Referências de Δp do R17 (entrada uniforme), que não tem 20 Pa como valor esperado.

    Os 20 Pa (G·L) valem só para o escoamento plenamente desenvolvido. Com entrada uniforme há duas
    parcelas a mais, que o Roteiro (ficha F) não quantifica: o fluxo de quantidade de movimento sai
    do uniforme para o perfil CP (ρ(β−1)ū², exato) e o atrito é maior na região de entrada.

    Args:
        m (dict): Métricas do R17 (``ubar``, ``fx_top``, ``fx_bottom``).

    Returns:
        dict: ``mom`` (ρ(β−1)ū² [Pa]), ``minimo`` (20 + mom: atrito desenvolvido, sem excesso de entrada)
        e ``balanco`` (atrito do Fluent + mom: balanço global de quantidade de movimento).
    """
    dd = np.linspace(0, 1, 20001)
    u = u_analitica(dd)
    mom = float(np.trapezoid(u ** 2, dd) - m["ubar"] ** 2)
    return dict(mom=mom, minimo=DP_REF + mom, balanco=m["fx_top"] + m["fx_bottom"] + mom)


def figura_io(cid, bloco, titulo, nome, nota, entrada_uniforme=False):
    """Bloco B (R06) ou F (R17): perfil de saída contra o analítico e Δp contra a referência.

    Em B a referência é G·L = 20 Pa. Em F (``entrada_uniforme``) 20 Pa é só o limite desenvolvido; a
    referência é o balanço de quantidade de movimento do domínio (ver ``referencia_dp_entrada_uniforme``).

    Args:
        cid (str): ID da rodada.
        bloco (str): Letra do bloco.
        titulo (str): Título da figura.
        nome (str): Nome do PNG.
        nota (str): Nota de rodapé (interpretação).
        entrada_uniforme (bool, optional): True para o R17. Padrão False.
    """
    m, p = metricas(cid), perfil_cp(cid)
    if m is None or p is None:
        return
    d, u = p
    dd = np.linspace(0, 1, 200)
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.8))
    ax[0].plot(u_analitica(dd), dd, "k-", lw=1.5, label="Analítico")
    ax[0].plot(u, d, "o", mfc="none", color="tab:blue", label=rotulo(cid))
    ax[0].set(xlabel="u [m/s]", ylabel="d = distância à parede parada [m]", title=f"Perfil em x = {virg(m['x0'])} m")
    ax[1].axhline(0, color="k", lw=1.2, label="Analítico (erro = 0)")
    ax[1].plot(d, u - u_analitica(d), "o-", mfc="none", color="tab:blue", label=rotulo(cid))
    ax[1].set(xlabel="d [m]", ylabel="u(Fluent) − u(analítico) [m/s]",
              title=f"Erro do perfil: e_L2 = {virg(m['eL2'], 3)} %, max = {virg(m['einf'], 3)} m/s")
    dp = m["dp"]
    if entrada_uniforme:
        r = referencia_dp_entrada_uniforme(m)
        ref, rot_ref = r["balanco"], f"Balanço de quantidade de movimento: {virg(r['balanco'], 2)} Pa"
        txt_ref = (f"{virg(ref, 2)} Pa (balanço: atrito do Fluent {virg(m['fx_top'] + m['fx_bottom'], 2)} Pa "
                   f"+ ρ(β−1)ū² = {virg(r['mom'], 2)} Pa)")
    else:
        ref, rot_ref, txt_ref = DP_REF, "Esperado: G·L = 4 × 5 = 20 Pa", "20 Pa (G·L = 4 Pa/m · 5 m)"
    desv = 100 * (dp - ref) / ref
    ax[2].bar([0], [dp], color="tab:blue", edgecolor="k", width=0.5, label=f"{rotulo(cid)}: {virg(dp)} Pa")
    ax[2].axhline(ref, color="k", ls="--", label=rot_ref)
    topo = max(dp, ref)
    if entrada_uniforme:
        ax[2].axhline(DP_REF, color="gray", ls=":", label="Limite desenvolvido (G·L): 20 Pa, não é a referência do F")
        ax[2].axhline(r["minimo"], color="tab:green", ls="-.", label=f"Mínimo físico 20 + ρ(β−1)ū² = {virg(r['minimo'], 2)} Pa")
    ax[2].text(0, dp * 0.5, f"{virg(dp)} Pa\n({'+' if desv >= 0 else '−'}{virg(abs(desv), 1)} %)", ha="center",
               color="w", fontsize=12, fontweight="bold")
    ax[2].set(xticks=[0], xticklabels=[cid], ylabel="Δp = p_entrada − p_saída [Pa]", ylim=(0, topo * 1.2),
              title=f"Queda de pressão (desvio {'+' if desv >= 0 else '−'}{virg(abs(desv), 1)} % da referência)")
    for a_ in ax:
        a_.grid(alpha=0.3)
        a_.legend(fontsize=8 if entrada_uniforme else 9, loc="best")
    fig.suptitle(titulo, fontsize=13)
    rodape(fig, nota)
    salvar(fig, nome)
    linha(bloco, cid, "Δp [Pa]", dp, txt_ref, f"{'+' if desv >= 0 else '−'}{virg(abs(desv), 1)} %",
          f"e_L2 do perfil = {m['eL2']:.4g} %" + (f"; mínimo físico {virg(r['minimo'], 2)} Pa; +28,6 % contra G·L não é desvio" if entrada_uniforme else ""))


def figura_b():
    """Bloco B (R06): entrada com perfil analítico e saída a 0 Pa."""
    figura_io("R06", "B", "Bloco B (R06): Couette-Poiseuille com entrada (perfil) e saída (0 Pa) contra o analítico",
              "bloco_B.png", "O desvio de Δp (−1,7 %) é aceito como limitação do caso: a saída a 0 Pa e a entrada por perfil "
              "não impõem o gradiente exato de 4 Pa/m; o perfil em x ≈ 2,5 m segue o analítico (e_L2 ≈ 0,5 %).")


def figura_f():
    """Bloco F (R17): entrada uniforme ū = 4,3333 m/s (região de desenvolvimento)."""
    figura_io("R17", "F", "Bloco F (R17): Couette-Poiseuille com entrada uniforme contra o analítico",
              "bloco_F.png", "Os 20 Pa (G·L) valem só para o escoamento desenvolvido; o Roteiro (ficha F) não dá Δp para a entrada uniforme. "
              "A referência usada é o balanço de quantidade de movimento no domínio (atrito do Fluent + ρ(β−1)ū², β do perfil CP). "
              "O perfil em x ≈ 2,5 m já é praticamente o analítico.", entrada_uniforme=True)


# ---------------------------------------------------------------- blocos C e D

def figura_ldc(ids, rot, bloco, titulo, nome, nota, extra_ref=None):
    """Blocos C e D: u na vertical central e v na horizontal central contra Ghia (1982), com max|Δ| por rodada.

    Args:
        ids (list[str]): Rodadas do bloco.
        rot (dict[str, str]): Rótulo de cada rodada (esquema ou algoritmo).
        bloco (str): Letra do bloco.
        titulo (str): Título da figura.
        nome (str): Nome do PNG.
        nota (str): Nota de rodapé (interpretação).
        extra_ref (str | None, optional): Rodada-base desenhada junto como referência de comparação (D usa R08).
    """
    fig, ax = plt.subplots(2, 2, figsize=(14, 10))
    gu = GHIA[GHIA["tabela"] == "I"]
    gv = GHIA[GHIA["tabela"] == "II"]
    ax[0, 0].plot(gu["vel"], gu["pos"], "ks", mfc="none", ms=6, label="Ghia et al. (1982)", zorder=5)
    ax[0, 1].plot(gv["pos"], gv["vel"], "ks", mfc="none", ms=6, label="Ghia et al. (1982)", zorder=5)
    lista = ([extra_ref] if extra_ref else []) + ids
    estilos = ["-", "--", "-.", ":", (0, (5, 1)), (0, (3, 1, 1, 1))]
    dados = {}
    for i, cid in enumerate(lista):
        df, m = campo(cid), metricas(cid)
        if df is None or m is None:
            continue
        s, u = linha_central(df, "x-coordinate", "x-velocity")
        s2, v = linha_central(df, "y-coordinate", "y-velocity")
        txt = rot.get(cid, "SIMPLE (base)")
        kw = dict(color=f"C{i}", ls=estilos[i % len(estilos)], lw=1.6, label=rotulo(cid, txt))
        ax[0, 0].plot(u, s, **kw)
        ax[0, 1].plot(s2, v, **kw)
        g = comparar_ghia(df)
        ax[1, 0].plot(g["I"]["pos"], g["I"]["erro"], marker="o", ms=3, **kw)
        dados[cid] = (m, i, txt)
    ax[0, 0].set(xlabel="u [m/s]", ylabel="y [m]", title="u na vertical central (x = 0,5 m)")
    ax[0, 1].set(xlabel="x [m]", ylabel="v [m/s]", title="v na horizontal central (y = 0,5 m)")
    ax[1, 0].axhline(0, color="k", lw=0.8)
    ax[1, 0].set(xlabel="y [m]", ylabel="u(Fluent) − u(Ghia) [m/s]", title="Erro de u nos pontos de Ghia (interior)")
    ids_ok = [c for c in lista if c in dados]
    x = np.arange(len(ids_ok))
    w = 0.38
    for k, (chave, nomev, cor) in enumerate((("ghia_du_max", "max|Δu|", "tab:blue"), ("ghia_dv_max", "max|Δv|", "tab:orange"))):
        vals = [dados[c][0][chave] for c in ids_ok]
        for xi, c, v in zip(x, ids_ok, vals):
            ax[1, 1].bar(xi + (k - 0.5) * w, v, w, color=cor, edgecolor="k", hatch="//" if c in NAO_CONV else None,
                         label=nomev if xi == 0 else None)
            ax[1, 1].text(xi + (k - 0.5) * w, v, f"{v * 1000:.2f}".replace(".", ","), ha="center", va="bottom", fontsize=7)
    ax[1, 1].set_xticks(x)
    ax[1, 1].set_xticklabels([rotulo(c, dados[c][2].split(" (")[0]) for c in ids_ok], rotation=25, ha="right", fontsize=8)
    ax[1, 1].set_ylim(0, 1.18 * max(dados[c][0][k] for c in ids_ok for k in ("ghia_du_max", "ghia_dv_max")))
    ax[1, 1].set(ylabel="Desvio máximo em relação a Ghia [m/s]", title="max|Δ| nos pontos de Ghia (valores em 10⁻³ m/s)")
    for a_ in ax.ravel():
        a_.grid(alpha=0.3)
        a_.legend(fontsize=8)
    fig.suptitle(titulo, fontsize=13)
    rodape(fig, "Referência: Ghia, Ghia e Shin (1982), Re = 100, malha 129 × 129. " + nota)
    salvar(fig, nome)
    for cid in ids:
        if cid not in dados:
            continue
        m = dados[cid][0]
        linha(bloco, cid, "max|Δu| vs Ghia [m/s]", m["ghia_du_max"], "Ghia et al. (1982), Tabela I (Re = 100)",
              f"{virg(m['ghia_du_max'] * 1000, 2)}e-3 m/s", f"max|Δv| = {m['ghia_dv_max']:.4g} m/s; "
              f"rms Δu = {m['ghia_du_rms']:.4g}; centro = ({m['xc']:.4f}; {m['yc']:.4f}) vs Ghia (0,6172; 0,7344)")


def figura_c():
    """Bloco C (R07–R11): esquemas de momento e malha fina contra Ghia (1982)."""
    rot = {"R07": "FOU N040", "R08": "SOU N040", "R09": "QUICK N040", "R10": "MUSCL N040", "R11": "SOU N128"}
    figura_ldc(["R07", "R08", "R09", "R10", "R11"], rot, "C",
               "Bloco C (R07–R11): cavidade Re = 100, esquemas de momento contra Ghia (1982)", "bloco_C.png",
               "Só o FOU (1ª ordem) se afasta de Ghia (máx. ≈ 18e-3 m/s); SOU, QUICK e MUSCL em N040 ficam em ≈ 5e-3 m/s e a N128 não zera o desvio.")


def figura_d():
    """Bloco D (R12–R14): acoplamento pressão-velocidade contra Ghia (1982), com R08 (SIMPLE) como base."""
    rot = {"R12": "SIMPLEC", "R13": "PISO", "R14": "Coupled"}
    figura_ldc(["R12", "R13", "R14"], rot, "D",
               "Bloco D (R12–R14): cavidade Re = 100, algoritmos pressão-velocidade contra Ghia (1982)", "bloco_D.png",
               "SIMPLEC, PISO e Coupled convergem para a mesma solução do SIMPLE (R08): as curvas e os desvios coincidem.", extra_ref="R08")


# ---------------------------------------------------------------- bloco E

def figura_e():
    """Bloco E (R15–R16): u junto ao fundo e comprimento de recolamento x_r/h por comprimento de saída."""
    fig, ax = plt.subplots(1, 2, figsize=(14, 5), gridspec_kw=dict(width_ratios=[1.5, 1]))
    leg = {"R15": "saída x = 2,5 m", "R16": "saída x = 11 m"}
    xr_h, ids = [], []
    for i, cid in enumerate(("R15", "R16")):
        df, m = campo(cid), metricas(cid)
        if df is None or m is None:
            continue
        y0 = df["y-coordinate"][df["y-coordinate"] > 0].min()
        lin = df[perto(df["y-coordinate"], y0) & (df["x-coordinate"] > 1)].sort_values("x-coordinate")
        ax[0].plot(lin["x-coordinate"], lin["x-velocity"], color=f"C{i}", label=rotulo(cid, leg[cid]))
        if m.get("xr") is not None:
            ax[0].axvline(1 + m["xr"], color=f"C{i}", ls="--", label=f"{cid}: recolamento x = {virg(1 + m['xr'], 3)} m")
            xr_h.append(m["xr_h"])
            linha("E", cid, "x_r/h", m["xr_h"], "sem referência em refs/ (Roteiro não traz x_r)", "n/d",
                  f"x_r = {m['xr']:.4f} m a partir do degrau; h = 0,5 m")
        else:
            xr_h.append(m["xr_maior_que"] / H_DEGRAU)
            ax[0].text(0.98, 0.1, "R15: u < 0 até a saída (bolha atravessa a saída)", transform=ax[0].transAxes,
                       ha="right", color=f"C{i}", fontsize=9)
            linha("E", cid, "x_r/h (cota inferior)", m["xr_maior_que"] / H_DEGRAU, "sem referência em refs/", "n/d",
                  f"x_r > {m['xr_maior_que']:.3f} m: a bolha atravessa a saída, x_r indefinido")
        ids.append(cid)
    ax[0].axhline(0, color="k", lw=0.6)
    ax[0].set(xlabel="x [m]", ylabel="u [m/s]", title="u na 1ª fileira de células acima do fundo a jusante do degrau")
    xs = np.arange(len(ids))
    for x, cid, v in zip(xs, ids, xr_h):
        cruza = metricas(cid).get("xr") is None
        ax[1].bar(x, v, color=f"C{int(cid == 'R16')}", edgecolor="k", hatch="//" if cruza else None, alpha=0.85)
        ax[1].text(x, v * (1.2 if cruza else 1.0), (">" if cruza else "") + virg(v, 2), ha="center", va="bottom", fontsize=11)
        if cruza:
            ax[1].annotate("", xy=(x, v * 1.18), xytext=(x, v), arrowprops=dict(arrowstyle="->", lw=2))
            ax[1].text(x, v * 0.5, "bolha\natravessa\na saída", ha="center", color="w", fontweight="bold")
    ax[1].set_xticks(xs)
    ax[1].set_xticklabels([rotulo(c, leg[c]) for c in ids])
    ax[1].set(ylabel="x_r / h  (h = 0,5 m)", title="Comprimento de recolamento", ylim=(0, max(xr_h) * 1.45 if xr_h else 1))
    for a_ in ax:
        a_.grid(alpha=0.3)
        a_.legend(fontsize=8) if a_ is ax[0] else None
    fig.suptitle("Bloco E (R15–R16): degrau laminar, influência do comprimento da saída", fontsize=13)
    rodape(fig, "Não há x_r de referência em refs/; a comparação é entre as duas saídas. Com a saída curta (R15) o recolamento "
                "não ocorre dentro do domínio, então a saída de 2,5 m é insuficiente.")
    salvar(fig, "bloco_E.png")


# ---------------------------------------------------------------- bloco G

def figura_g():
    """Bloco G (R18–R30): interpolação de pressão; controle Ω = 0, cavidade rotativa Ω = 1 rad/s e duto 3D."""
    fig, ax = plt.subplots(2, 3, figsize=(18, 10.5), constrained_layout=True)
    # (a) controle Omega = 0
    ids0, rot0 = ["R18", "R19", "R20"], {"R18": "PRESTO!", "R19": "Second Order", "R20": "Standard"}
    ms = {c: metricas(c) for c in ids0}
    if all(ms.values()):
        umax = ms["R19"]["u_max"]
        w = 0.38
        for k, (chave, nome_, cor) in enumerate((("max_sem_tampa", "max|Δu| (sem a fileira da tampa)", "tab:blue"),
                                                 ("rms", "rms Δu", "tab:orange"))):
            for i, c in enumerate(ids0):
                if c == "R19":
                    v = 0.0
                else:
                    v = 100 * ms[c]["du_R19"][chave] / umax
                ax[0, 0].bar(i + (k - 0.5) * w, v, w, color=cor, edgecolor="k", label=nome_ if i == 0 else None)
                ax[0, 0].text(i + (k - 0.5) * w, v, virg(v, 2), ha="center", va="bottom", fontsize=8)
        ax[0, 0].axhline(1.0, color="r", ls="--", label="Limite de G1 (1 %)")
        ax[0, 0].set_xticks(range(3))
        ax[0, 0].set_xticklabels([rotulo(c, rot0[c]) for c in ids0], fontsize=8)
        ax[0, 0].set(ylabel="Diferença para R19 [% de u_máx]", title="Ω = 0 (controle, N040): diferença para R19 (Second Order)")
        ax[0, 0].legend(fontsize=8)
        for c in ids0:
            if c == "R19":
                linha("G", c, "max|Δu| vs R19 sem tampa [% u_max]", 0.0, "R19 (Second Order, Ω = 0)", "0 (própria referência)",
                      "controle Ω = 0")
            else:
                v = 100 * ms[c]["du_R19"]["max_sem_tampa"] / umax
                linha("G", c, "max|Δu| vs R19 sem tampa [% u_max]", v, "R19 (Second Order, Ω = 0); limite G1 = 1 %",
                      f"{'+' if v >= 1 else ''}{virg(v - 1, 2)} p.p. em relação ao limite de 1 %",
                      f"rms = {100 * ms[c]['du_R19']['rms'] / umax:.3g} %; Ω = 0 (controle)")
    else:
        avisar("bloco G: faltam R18–R20; painel Ω = 0 pulado")
    # (b, c, d) cavidade rotativa N040, referência R26
    ids1 = ["R21", "R22", "R23", "R24", "R25"]
    nomes = {"R21": "Linear", "R22": "Standard", "R23": "Second\nOrder", "R24": "Body Force\nWeighted", "R25": "PRESTO!"}
    m1 = {c: metricas(c) for c in ids1}
    m26, m27 = metricas("R26"), metricas("R27")
    ok = [c for c in ids1 if m1[c]]
    for a_, chave, yl, tit in ((ax[0, 1], "eL2_u", "e_L2 de u contra R26 [%]", "Ω = 1 rad/s, N040: erro e_L2 (u) contra R26"),
                               (ax[0, 2], "cb", "CB (oscilação de p em y = 0,5) [–]", "Ω = 1 rad/s: oscilação de pressão (cb)"),
                               (ax[1, 0], "overshoot", "Overshoot de u contra R26 [m/s]", "Ω = 1 rad/s: overshoot de u em y = 0,5")):
        vals = [m1[c][chave] for c in ok]
        rots = [nomes[c] for c in ok]
        ids = list(ok)
        if m27 and chave in m27:
            vals.append(m27[chave])
            rots.append("N128\nSecond Order")
            ids.append("R27")
        barras(a_, ids, vals, [rotulo(c, r.replace("\n", " ")) for c, r in zip(ids, rots)])
        for i, v in enumerate(vals):
            a_.text(i, v, f"{v:.3g}".replace(".", ","), ha="center", va="bottom", fontsize=8)
        a_.set(ylabel=yl, title=tit)
        if chave == "cb" and m26:
            a_.axhline(m26["cb"], color="k", ls="--", label=f"R26 (N128, PRESTO!): {virg(m26['cb'], 4)}")
            a_.legend(fontsize=8)
        if chave == "eL2_u":
            a_.axhline(0, color="k", ls="--", label="R26 (referência: e_L2 = 0 por definição)")
            a_.legend(fontsize=8)
        if chave == "overshoot":
            a_.axhline(0, color="k", ls="--", label="R26 (referência: 0)")
            a_.legend(fontsize=8)
    for c in ok:
        m = m1[c]
        linha("G", c, "e_L2 de u contra R26 [%]", m["eL2_u"], "R26 (N128, PRESTO!, Ω = 1, referência interna)",
              f"+{virg(m['eL2_u'], 2)} %", f"cb = {m['cb']:.4g} (R26: {m26['cb']:.4g} se existir); overshoot = {m['overshoot']:.4g} m/s" if m26 else "")
    if m27:
        linha("G", "R27", "e_L2 de u contra R26 [%]", m27["eL2_u"], "R26 (N128, PRESTO!)", f"+{virg(m27['eL2_u'], 2)} %",
              f"checagem da referência; max|Δu|/|V| sem fileira da tampa = {m27.get('ref_dif_sem_tampa', float('nan')):.3g} %; cb = {m27['cb']:.4g}")
    if m26:
        linha("G", "R26", "cb [–]", m26["cb"], "própria referência interna (N128, PRESTO!)", "0 (referência)", "Ω = 1 rad/s")
    # (e, f) duto 3D
    idsd = ["R28", "R29", "R30"]
    nd = {"R28": "Second Order", "R29": "Standard", "R30": "PRESTO!"}
    md = {c: metricas(c) for c in idsd}
    okd = [c for c in idsd if md[c]]
    if okd:
        barras(ax[1, 1], okd, [md[c]["cf"] for c in okd], [rotulo(c, nd[c]) for c in okd])
        ax[1, 1].axhline(CF_PIROZZOLI, color="k", ls="--", label="Pirozzoli et al. (2018): Cf = 9,38e-3 (Re_b = 4410)")
        ax[1, 1].set(ylabel="Cf = 2/(ρ U_b²) [–]", title="Duto 3D + SST: Cf por pressão", ylim=(0, 0.0125))
        ax[1, 1].legend(fontsize=8)
        barras(ax[1, 2], okd, [md[c]["ub"] for c in okd], [rotulo(c, nd[c]) for c in okd])
        ax[1, 2].axhline(UB_DNS, color="k", ls="--", label="DNS Gavrilakis (1992): U_b ≈ 14,7 m/s")
        ax[1, 2].set(ylabel="U_b [m/s]", title="Duto 3D + SST: velocidade média", ylim=(0, 17))
        ax[1, 2].legend(fontsize=8)
        for i, c in enumerate(okd):
            ax[1, 1].text(i, md[c]["cf"], f"{md[c]['cf'] * 1000:.2f}e-3".replace(".", ","), ha="center", va="bottom", fontsize=8)
            ax[1, 2].text(i, md[c]["ub"], virg(md[c]["ub"], 3), ha="center", va="bottom", fontsize=8)
            desv = 100 * (md[c]["cf"] - CF_PIROZZOLI) / CF_PIROZZOLI
            linha("G", c, "Cf [–]", md[c]["cf"], "Pirozzoli et al. (2018), Tabela 1 caso A: 9,38e-3 (ordem de grandeza)",
                  f"{'+' if desv >= 0 else '−'}{virg(abs(desv), 1)} %",
                  f"U_b = {md[c]['ub']:.4f} m/s (DNS ≈ 14,7); vs/U_b = {md[c]['vs_ub']:.2g} %; duto 3D, pressão {nd[c]}")
    for a_ in ax.ravel():
        a_.grid(alpha=0.3)
    fig.suptitle("Bloco G (R18–R30): esquemas de interpolação de pressão (Ω = 0, Ω = 1 rad/s e duto 3D)", fontsize=13)
    rodape(fig, "Ω = 0: R19 é a referência de comparação (R18 e R20 diferem sobretudo nos cantos singulares da tampa). "
                "Ω = 1: referência interna R26 (N128, PRESTO!); R27 (N128, SO) checa a referência. "
                "Duto: as três pressões dão resultados iguais (SST sem escoamento secundário), Cf e U_b só como ordem de grandeza. "
                "* = não convergiu.")
    salvar(fig, "bloco_G.png")


# ---------------------------------------------------------------- bloco H (duto)

def texto_valor(chave, v):
    """Texto do valor escrito sobre a barra, conforme a grandeza.

    Args:
        chave (str): ``"ub"``, ``"vs_ub"`` ou ``"cf"``.
        v (float): Valor.

    Returns:
        str: Número com vírgula decimal.
    """
    if chave == "ub":
        return virg(v, 2)
    if chave == "vs_ub":
        return f"{v:.2g}".replace(".", ",")
    return f"{v * 1000:.2f}e-3".replace(".", ",")


def figura_h_duto():
    """Bloco H, duto 3D (R31–R34): U_b, intensidade do escoamento secundário e Cf por modelo."""
    ids = ["R31", "R32", "R33", "R34"]
    nome = {"R31": "k-ε Standard", "R32": "k-ω SST", "R33": "RSM linear", "R34": "SST + Corner\nFlow Corr."}
    ms = {c: metricas(c) for c in ids}
    ok = [c for c in ids if ms[c]]
    if not ok:
        avisar("duto: sem rodadas")
        return
    rots = [rotulo(c, nome[c].replace("\n", " ")) for c in ok]
    fig, ax = plt.subplots(1, 3, figsize=(18, 6), constrained_layout=True)
    barras(ax[0], ok, [ms[c]["ub"] for c in ok], rots)
    ax[0].axhline(UB_DNS, color="k", ls="--", label="DNS Gavrilakis (1992): U_b ≈ 14,7 m/s")
    ax[0].set(ylabel="Velocidade média U_b [m/s]", title="Velocidade média do duto", ylim=(0, 17))
    barras(ax[1], ok, [ms[c]["vs_ub"] for c in ok], rots)
    ax[1].axhspan(*VS_UB_FAIXA, color="gold", alpha=0.35, label="Pirozzoli et al. (2018): 1–2 % de U_b")
    ax[1].set(ylabel="|V_s|máx / U_b [%]", title="Intensidade do escoamento secundário (plano xz)")
    barras(ax[2], ok, [ms[c]["cf"] for c in ok], rots)
    ax[2].axhline(CF_PIROZZOLI, color="k", ls="--", label="Pirozzoli et al. (2018): Cf = 9,38e-3 (Re_b = 4410)")
    ax[2].set(ylabel="Cf = 2/(ρ U_b²) [–]", title="Coeficiente de atrito", ylim=(0, 0.0145))
    for a_, chave in ((ax[0], "ub"), (ax[1], "vs_ub"), (ax[2], "cf")):
        for i, c in enumerate(ok):
            a_.text(i, ms[c][chave], texto_valor(chave, ms[c][chave]), ha="center", va="bottom", fontsize=9)
        a_.legend(fontsize=8)
        a_.grid(alpha=0.3)
    fig.suptitle("Bloco H, duto quadrado 3D (R31–R34): modelos de turbulência contra a ordem de grandeza da DNS", fontsize=13)
    rodape(fig, "Modelos de Boussinesq (k-ε, SST) não geram escoamento secundário (≈ 0 %); o RSM (2,3 %, pouco acima da faixa) e o SST com correção de canto "
                "(1,5 %, dentro) são os únicos com escoamento secundário na ordem de 1–2 % de U_b. Cf e U_b só como ordem de grandeza (Re_b do Fluent difere de 4410). * = não convergiu.")
    salvar(fig, "bloco_H_duto.png")
    for c in ok:
        m = ms[c]
        v = m["vs_ub"]
        desv = 0.0 if VS_UB_FAIXA[0] <= v <= VS_UB_FAIXA[1] else (v - VS_UB_FAIXA[0] if v < VS_UB_FAIXA[0] else v - VS_UB_FAIXA[1])
        linha("H-duto", c, "|V_s|máx/U_b [%]", v, "Pirozzoli et al. (2018): 1–2 % de U_b",
              f"{'dentro da faixa' if desv == 0 else ('+' if desv > 0 else '−') + virg(abs(desv), 2) + ' p.p. fora da faixa'}",
              f"U_b = {m['ub']:.3f} m/s (DNS ≈ 14,7); Cf = {m['cf']:.3e} (ref. 9,38e-3)")


# ---------------------------------------------------------------- bloco H (degrau)

def figura_h_degrau():
    """Bloco H, degrau turbulento (R35–R43): x_r/H, Cf(x/H) e Cp(x/H) contra Driver e Seegmiller (1985)."""
    nome = {"R35": "k-ε Std", "R36": "k-ε Realiz.", "R37": "k-ω BSL", "R38": "k-ω SST", "R39": "SST + Curv.", "R40": "RSM",
            "R41": "SST (y+≈30)", "R42": "k-ε RNG", "R43": "GEKO C_SEP=2,5"}
    ordem = [f"R{n}" for n in range(35, 44)]
    fig, ax = plt.subplots(2, 2, figsize=(17, 12), constrained_layout=True)
    ex = DRIVER.set_index("serie")
    ax[0, 1].errorbar(ex.loc["cf", "pos"], ex.loc["cf", "val"], yerr=abs(ex.loc["cf", "val"] * ex.loc["cf", "dcf_rel"]),
                      fmt="ks", mfc="none", ms=5, label="Driver e Seegmiller (1985)", zorder=10)
    ax[1, 0].plot(ex.loc["cp_step", "pos"], ex.loc["cp_step", "val"], "ks", mfc="none", ms=5, label="Driver e Seegmiller (1985)", zorder=10)
    ax[1, 1].plot(ex.loc["cp_top", "pos"], ex.loc["cp_top", "val"], "k^", mfc="none", ms=5, label="Driver e Seegmiller (1985)", zorder=10)
    xr_ok, xr_val = [], []
    for i, cid in enumerate(ordem):
        m = metricas(cid)
        if m is None:
            continue
        try:
            lado, top = paredes_bft(caso(cid))
        except (FileNotFoundError, KeyError) as exc:
            avisar(f"{cid}: paredes ausentes ({exc}); curvas de Cf/Cp puladas")
            lado = top = None
        cor = plt.get_cmap("tab10")(i)
        if cid in NAO_CONV:
            ls = "--"
        elif i % 3 == 1:
            ls = "-."
        else:
            ls = "-"
        r = rotulo(cid, nome[cid])
        if lado is not None:
            lado, top = lado.sort_values("x"), top.sort_values("x")
            ax[0, 1].plot(lado["x"], lado["cf"], color=cor, ls=ls, lw=1.4, label=r)
            ax[1, 0].plot(lado["x"], lado["cp"], color=cor, ls=ls, lw=1.4, label=r)
            ax[1, 1].plot(top["x"], top["cp"], color=cor, ls=ls, lw=1.4, label=r)
        if m.get("xr_h") is not None:
            xr_ok.append(cid)
            xr_val.append(m["xr_h"])
            desv = 100 * (m["xr_h"] - XR_H_EXP) / XR_H_EXP
            linha("H-degrau", cid, "x_r/H", m["xr_h"], f"{virg(XR_H_EXP)} ± {virg(XR_H_INC)} (Driver e Seegmiller, 1985)",
                  f"{'+' if desv >= 0 else '−'}{virg(abs(desv), 1)} %",
                  f"Cf_rms = {m['cf_rms']:.2e}; Cp_rms degrau/oposta = {m['cp_step_rms']:.3f}/{m['cp_top_rms']:.3f}; "
                  f"{m['iteracoes']} iterações" + ("; NÃO convergiu" if cid in NAO_CONV else ""))
        else:
            avisar(f"{cid}: sem x_r (sem recolamento detectado); fora da barra")
    # x_r/H em barras
    ax0 = ax[0, 0]
    barras(ax0, xr_ok, xr_val, [rotulo(c, nome[c]) for c in xr_ok], cores=[plt.get_cmap("tab10")(ordem.index(c)) for c in xr_ok])
    ax0.axhspan(XR_H_EXP - XR_H_INC, XR_H_EXP + XR_H_INC, color="k", alpha=0.35, zorder=5)
    ax0.axhline(XR_H_EXP, color="k", ls="--", label=f"Exp.: {virg(XR_H_EXP)} ± {virg(XR_H_INC)} (Driver e Seegmiller, 1985)")
    for i, (c, v) in enumerate(zip(xr_ok, xr_val)):
        ax0.text(i, v, virg(v, 2) + f"\n({'+' if v >= XR_H_EXP else '−'}{virg(abs(100 * (v - XR_H_EXP) / XR_H_EXP), 1)} %)",
                 ha="center", va="center", fontsize=9, bbox=dict(fc="w", ec="none", alpha=0.8), zorder=8)
    ax0.set(ylabel="x_r / H [–]", title="Comprimento de recolamento (faixa preta = incerteza do experimento)",
            ylim=(0, max(xr_val + [XR_H_EXP]) * 1.3 if xr_val else 8))
    ax0.legend(fontsize=8, loc="upper left")
    ax[0, 1].axhline(0, color="k", lw=0.5)
    ax[0, 1].set(xlim=(-8, 30), ylim=(-0.004, 0.006), xlabel="x/H [–]", ylabel="Cf [–]", title="Cf na parede do lado do degrau")
    ax[1, 0].set(xlim=(-8, 30), xlabel="x/H [–]", ylabel="Cp (ref. x/H = −4) [–]", title="Cp na parede do lado do degrau")
    ax[1, 1].set(xlim=(-8, 30), xlabel="x/H [–]", ylabel="Cp (ref. x/H = −4) [–]", title="Cp na parede oposta (superior)")
    for a_ in ax.ravel():
        a_.grid(alpha=0.3)
    for a_ in (ax[0, 1], ax[1, 0], ax[1, 1]):
        a_.legend(fontsize=7, ncol=2)
    fig.suptitle("Bloco H, degrau turbulento (R35–R43): modelos de turbulência contra Driver e Seegmiller (1985), Re_H = 36 000",
                 fontsize=13)
    rodape(fig, "* (hachura e traço tracejado) = rodada com Convergiu_S_N = N (R37–R40 e R43: não atingiram o critério de resíduos 1e-6 "
                "em 20 000–22 000 iterações); seus valores são indicativos. R41 usa a malha de parede grosseira BFYW (y+ ≈ 30). "
                "O perfil de entrada não é o medido em x/H = −4; vale mais a comparação entre modelos que o valor absoluto.")
    salvar(fig, "bloco_H_degrau.png")


# ---------------------------------------------------------------- principal

def gravar_resumo():
    """Grava ``resumo.csv`` (separador ``;``) com uma linha por rodada."""
    SAIDA.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(RESUMO)
    df["_n"] = df["rodada"].str[1:].astype(int)
    df = df.sort_values("_n").drop(columns="_n")
    df.to_csv(SAIDA / "resumo.csv", sep=";", index=False, encoding="utf-8-sig")
    print(f"gerado: {SAIDA / 'resumo.csv'} ({len(df)} linhas)")


def principal(argv):
    """Gera as figuras dos blocos pedidos e o resumo.

    Args:
        argv (list[str]): Letras dos blocos (A–H); vazio = todos.
    """
    carregar_registro()
    plt.rcParams.update({"font.size": 10, "axes.titlesize": 10.5, "figure.dpi": 100})
    pedidos = [a.upper() for a in argv] or list("ABCDEFGH")
    tabela = {"A": [figura_a], "B": [figura_b], "C": [figura_c], "D": [figura_d], "E": [figura_e], "F": [figura_f],
              "G": [figura_g], "H": [figura_h_duto, figura_h_degrau]}
    for b in pedidos:
        for f in tabela[b]:
            try:
                f()
            except Exception as exc:  # um bloco com dado faltando não impede os demais
                avisar(f"bloco {b} ({f.__name__}) falhou: {type(exc).__name__}: {exc}")
    gravar_resumo()


if __name__ == "__main__":
    principal(sys.argv[1:])
