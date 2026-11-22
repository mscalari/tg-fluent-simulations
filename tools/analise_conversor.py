"""Análise da varredura MRF do conversor CT_V1: gráficos por ponto, curvas globais, Euler e nota do Obsidian.

Uso: python tools/analise_conversor.py            (analisa todos os SRxxx com resultado.json)
     python tools/analise_conversor.py --sem-nota (só gráficos e tabela)

Lê results/conversor/CT_V1_MRF/{dados,logs,SRxxx}. Grava gráficos em SRxxx/graficos e
CT_V1_MRF/graficos, a tabela em CT_V1_MRF/dados/desempenho.csv e a nota em
<vault>/C_TG/05_Conversor/CT_V1_MRF.md (imagens copiadas para <vault>/C_TG/img).
"""
import csv, json, re, shutil, sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
DIR = ROOT / "results" / "conversor" / "CT_V1_MRF"
VAULT = Path(r"C:\Users\murilo\Repositorios\obsidian-vault-murilo-scalari\CLAUDE\C_TG")
CORES = {"bomba": "#1f77b4", "turbina": "#d62728", "estator": "#2ca02c", "total": "#555555"}
# Interfaces de cada zona (lado que pertence à zona). As normais das faces apontam para fora da zona.
LADOS = {"fluido_bomba": ["bomba_turbina_bomba", "estator_bomba_bomba"],
         "fluido_turbina": ["bomba_turbina_turbina", "turbina_estator_turbina"],
         "fluido_estator": ["turbina_estator_estator", "estator_bomba_estator"]}
# Os reports "moment" da varredura foram criados sem ``report_output_type`` e saíram como coeficiente Cm,
# com os valores de referência padrão do Fluent (ρ = 1,225 kg/m³, v = 1 m/s, A = 1 m², L = 1 m):
# M = Cm · ½·ρ·v²·A·L = 0,6125·Cm. Resultados novos gravam ``unidade_torque = "N·m"`` e não são convertidos.
CM_PARA_NM = 0.5 * 1.225 * 1.0 ** 2 * 1.0 * 1.0


def fator_torque(r):
    """Fator que converte os torques gravados num ``resultado.json`` para N·m.

    Args:
        r (dict): ``resultado.json`` de um ponto.

    Returns:
        float: 1 se o ponto já gravou em N·m; ``CM_PARA_NM`` se gravou o coeficiente Cm (varredura original).
    """
    return 1.0 if r.get("unidade_torque") == "N·m" else CM_PARA_NM


plt.rcParams.update({"figure.dpi": 130, "axes.grid": True, "grid.alpha": 0.3, "font.size": 9})


def ler_tabela(arquivo):
    """Lê um Report File do Fluent.

    Args:
        arquivo (Path): Arquivo ``.out``.

    Returns:
        tuple[list[str], np.ndarray]: Nomes das colunas e matriz (linhas × colunas).
    """
    nomes, linhas = [], []
    for ln in arquivo.read_text(errors="replace").splitlines():
        if '"Iteration"' in ln:
            nomes = re.findall(r'"([^"]+)"', ln)
        p = ln.split()
        if p and p[0].isdigit():
            try:
                linhas.append([float(x) for x in p])
            except ValueError:
                pass
    return nomes, np.array(linhas) if linhas else np.zeros((0, len(nomes)))


def monitores():
    """Junta todos os Report Files da varredura (uma por sessão), sem iterações repetidas.

    Returns:
        tuple[list[str], np.ndarray]: Colunas e matriz ordenada por iteração.
    """
    nomes, tudo = None, []
    for f in sorted((DIR / "dados").glob("CT_V1_MRF_rep_*.out")) + sorted((DIR / "dados").glob("sentido_mais_rep.out")):
        n, m = ler_tabela(f)
        if len(m) and (nomes is None or n == nomes):
            nomes = n
            tudo.append(m)
    if not tudo:
        return [], np.zeros((0, 0))
    m = np.vstack(tudo)
    _, idx = np.unique(m[:, 0], return_index=True)
    return nomes, m[idx]


def residuos(log):
    """Extrai o histórico de resíduos do log do Fluent (linhas "iter continuity ...").

    Args:
        log (Path): Log da sessão.

    Returns:
        tuple[list[str], np.ndarray]: Equações e matriz (iteração, resíduos...).
    """
    eqs, linhas = [], []
    for ln in log.read_text(errors="replace").splitlines():
        p = ln.split()
        if p and p[0] == "iter" and "continuity" in p:
            eqs = [x for x in p[1:] if x not in ("time/iter",)]
            eqs = eqs[:eqs.index("time/iter")] if "time/iter" in eqs else eqs
            eqs = [e for e in eqs if not e.startswith(("t_", "m_", "h_"))]
            continue
        if eqs and p and p[0].isdigit() and len(p) > len(eqs):
            try:
                linhas.append([float(p[0])] + [float(x) for x in p[1:1 + len(eqs)]])
            except ValueError:
                pass
    if not linhas:
        return eqs, np.zeros((0, 0))
    m = np.array(linhas)
    _, idx = np.unique(m[:, 0], return_index=True)
    return eqs, m[idx]


def pontos():
    """Lista os pontos concluídos.

    Returns:
        list[tuple[Path, dict]]: (pasta, resultado) por SR crescente.
    """
    out = []
    for p in sorted(DIR.glob("SR[0-9][0-9][0-9]")):
        f = p / "dados" / "resultado.json"
        if f.exists():
            out.append((p, json.loads(f.read_text())))
    return sorted(out, key=lambda x: x[1]["sr"])


def desempenho(r):
    """Grandezas de desempenho de um ponto.

    Args:
        r (dict): ``resultado.json``.

    Returns:
        dict: SR, torques, TR, η, K, vazão, balanço e Euler por zona.
    """
    c = fator_torque(r)
    def conv(dd):
        """Converte para N·m os valores ``t_*``/``z_*`` de um dicionário de reports.

        Args:
            dd (dict): Reports (chave -> valor).

        Returns:
            dict: Cópia com os torques multiplicados por ``c`` e o resto inalterado.
        """
        return {k: (v * c if isinstance(v, (int, float)) and k.startswith(("t_", "z_")) else v)
                for k, v in dd.items()}
    f = conv(r["reports"])                 # médias das últimas JANELA iterações, em N·m
    fi = conv(r.get("reports_inst", r["reports"]))  # última iteração (Euler e torques por zona são instantâneos)
    if isinstance(fi.get("paredes_zona"), dict):
        fi["paredes_zona"] = {k: v * c for k, v in fi["paredes_zona"].items()}
    dp = conv(r.get("desvio_janela", {}))
    tb, tt, te, tot = f["t_bomba"], f["t_turbina"], f["t_estator"], f["t_total"]
    sr = r["sr"]
    tr = -tt / tb
    d = {"SR": sr, "T_b": tb, "T_t": tt, "T_e": te, "soma": tot, "TR": tr, "eta": sr * tr,
         "K": r["rpm_b"] / np.sqrt(abs(tb)), "balanco_pct": 100 * tot / abs(tb),
         "mdot_bt": abs(f.get("m_bomba_turbina_bomba") or float("nan")), "iter": r["it_fim"] - r["it_ini"], "conv": r["convergiu"],
         "dp_T_b_pct": 100 * dp.get("t_bomba", float("nan")) / abs(tb), "dp_T_t_pct": 100 * dp.get("t_turbina", float("nan")) / abs(tt),
         "dp_TR": tr * ((dp.get("t_bomba", 0) / abs(tb)) ** 2 + (dp.get("t_turbina", 0) / abs(tt)) ** 2) ** 0.5}
    sinal = r.get("sinal", 1)
    for zona, lados in LADOS.items():
        hs = [fi.get("h_" + l) for l in lados]
        if all(h is not None for h in hs) and all(abs(h) > 0 for h in hs):
            # Reports do Fluent: fluxo positivo ENTRANDO na zona (as vazões confirmam: −298,7/+298,7 kg/s em if_bt).
            # Torque do fluido nas paredes da zona = −(fluxo líquido que sai) = +Σ fluxo que entra.
            d["euler_" + zona] = sinal * sum(hs)
        ws = fi.get("paredes_zona", {}).get(zona) if isinstance(fi.get("paredes_zona"), dict) else None
        if ws is not None:
            d["parede_" + zona] = ws
    return d


def fluxos_lado(f):
    """Fluxos de massa e de momento angular (eixo +x) num lado de interface exportado.

    Args:
        f (Path): CSV de um lado (valores da célula adjacente, área orientada para fora da zona).

    Returns:
        tuple[float, float]: (ṁ, H) = (Σ ρ V·A, Σ ρ V·A (y·w − z·v)), positivos saindo da zona.
    """
    d = ler_csv(f)
    vA = d["x-velocity"] * d["x-face-area"] + d["y-velocity"] * d["y-face-area"] + d["z-velocity"] * d["z-face-area"]
    lx = d["y-coordinate"] * d["z-velocity"] - d["z-coordinate"] * d["y-velocity"]
    return float(np.sum(d["density"] * vA)), float(np.sum(d["density"] * vA * lx))


def euler_csv(p, sinal):
    """Torque nas paredes de cada zona por Euler, a partir das interfaces exportadas.

    O export ASCII das interfaces não conformes traz o valor da célula adjacente de cada lado, não o da
    face; a média dos dois lados (diferença central) recupera a vazão exata dos reports do Fluent
    (SR000: 298,9 contra 298,7 kg/s). O fluxo de A para B é F = (H_A − H_B)/2, com H de cada lado
    orientado para fora da própria zona; T_zona = −Σ fluxo que sai = −Σ F.

    Args:
        p (Path): Pasta do ponto.
        sinal (int): Sentido do eixo de ω_b (+1 ou −1).

    Returns:
        dict[str, float]: Zona -> torque por Euler [N·m] (projeção em sinal·x), mais ``mdot_<if>``; vazio se faltar arquivo.
    """
    pares = {"bomba_turbina": ("fluido_bomba", "fluido_turbina"),
             "turbina_estator": ("fluido_turbina", "fluido_estator"),
             "estator_bomba": ("fluido_estator", "fluido_bomba")}
    sai = {z: 0.0 for z in LADOS}
    out = {}
    for par, (za, zb) in pares.items():
        fa = p / "dados" / f"{p.name}_if_{par}_lado_{za.split('_')[1]}.csv"
        fb = p / "dados" / f"{p.name}_if_{par}_lado_{zb.split('_')[1]}.csv"
        if not (fa.exists() and fb.exists()):
            return {}
        (ma, ha), (mb, hb) = fluxos_lado(fa), fluxos_lado(fb)
        F = (ha - hb) / 2
        sai[za] += F
        sai[zb] -= F
        out["mdot_" + par] = (ma - mb) / 2
    for z in LADOS:
        out[z] = -sinal * sai[z]
    return out


def graf_monitores(nomes, m, p, r):
    """Torques, soma e vazões do ponto ao longo das iterações.

    Args:
        nomes (list[str]): Colunas.
        m (np.ndarray): Monitores.
        p (Path): Pasta do ponto.
        r (dict): Resultado do ponto.
    """
    ini = r["it_ini"] + (50 if r["sr"] == 0 else 0)     # o 1º ponto parte do repouso: o transitório esmaga a escala
    sel = (m[:, 0] > ini) & (m[:, 0] <= r["it_fim"])
    if not sel.any():
        return
    c = {n: i for i, n in enumerate(nomes)}
    it = m[sel, 0]
    k = fator_torque(r)
    fig, ax = plt.subplots(1, 3, figsize=(13, 3.6))
    for n in ("bomba", "turbina", "estator"):
        ax[0].plot(it, k * m[sel, c["t_" + n]], label=f"T_{n[0]}", color=CORES[n])
    ax[0].set(xlabel="iteração", ylabel="torque [N·m]", title="Torques (eixo de ω_b)")
    ax[0].legend()
    ax[1].plot(it, 100 * m[sel, c["t_total"]] / np.abs(m[sel, c["t_bomba"]]), color="k")
    ax[1].set(xlabel="iteração", ylabel="(T_b+T_t+T_e)/|T_b| [%]", title="Balanço de torque")
    for n in [n for n in nomes if n.startswith("m_")]:
        ax[2].plot(it, m[sel, c[n]], label=n[2:], lw=1)
    ax[2].set(xlabel="iteração", ylabel="vazão [kg/s]", title="Vazão nas interfaces (normal para fora da zona)")
    ax[2].legend(fontsize=6)
    fig.tight_layout()
    fig.savefig(p / "graficos" / f"{p.name}_m1_monitores.png")
    plt.close(fig)


def graf_residuos(eqs, res, p, r):
    """Resíduos do ponto.

    Args:
        eqs (list[str]): Equações.
        res (np.ndarray): Histórico.
        p (Path): Pasta do ponto.
        r (dict): Resultado do ponto.
    """
    if not len(res):
        return
    sel = (res[:, 0] > r["it_ini"]) & (res[:, 0] <= r["it_fim"])
    if not sel.any():
        return
    fig, ax = plt.subplots(figsize=(6, 3.8))
    for i, e in enumerate(eqs, 1):
        ax.semilogy(res[sel, 0], res[sel, i], label=e, lw=1)
    ax.set(xlabel="iteração", ylabel="resíduo escalado", title=f"Resíduos, {p.name}")
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(p / "graficos" / f"{p.name}_m2_residuos.png")
    plt.close(fig)


def ler_csv(f):
    """Lê um CSV exportado pelo Fluent (cabeçalho com espaços).

    Args:
        f (Path): Arquivo.

    Returns:
        dict[str, np.ndarray]: Coluna -> valores.
    """
    with open(f, newline="") as fh:
        rd = csv.reader(fh)
        cab = [c.strip() for c in next(rd)]
        dados = np.array([[float(x) for x in row] for row in rd if row])
    return {c: dados[:, i] for i, c in enumerate(cab)}


def graf_yplus(p):
    """Histograma de y⁺ por peça e estatísticas.

    Args:
        p (Path): Pasta do ponto.

    Returns:
        dict[str, dict]: Peça -> {média, mediana, p5, p95, máx}.
    """
    est = {}
    fig, ax = plt.subplots(figsize=(6, 3.8))
    for peca in ("bomba", "turbina", "estator"):
        f = p / "dados" / f"{p.name}_parede_{peca}.csv"
        if not f.exists():
            continue
        d = ler_csv(f)
        y = d.get("y-plus")
        if y is None:
            continue
        est[peca] = {"media": float(y.mean()), "mediana": float(np.median(y)), "p5": float(np.percentile(y, 5)),
                     "p95": float(np.percentile(y, 95)), "max": float(y.max())}
        ax.hist(y, bins=np.logspace(-1, 3.5, 60), histtype="step", lw=1.5, color=CORES[peca], label=peca)
    ax.set_xscale("log")
    ax.axvline(30, ls="--", color="gray")
    ax.set(xlabel="y⁺", ylabel="faces", title=f"y⁺ nas pás e cascas, {p.name}")
    ax.legend()
    fig.tight_layout()
    fig.savefig(p / "graficos" / f"{p.name}_m3_yplus_hist.png")
    plt.close(fig)
    return est


def graf_perfis(p):
    """Perfis radiais (médias por faixa de raio) de c_ax, c_r, c_θ e p em cada lado de cada interface.

    Args:
        p (Path): Pasta do ponto.
    """
    arqs = sorted((p / "dados").glob(f"{p.name}_if_*.csv"))
    if not arqs:
        return
    fig, ax = plt.subplots(2, 2, figsize=(11, 7))
    for f in arqs:
        d = ler_csv(f)
        if "y-coordinate" not in d:
            continue
        r = np.hypot(d["y-coordinate"], d["z-coordinate"])
        bins = np.linspace(r.min(), r.max(), 25)
        k = np.digitize(r, bins)
        rc = [r[k == i].mean() for i in range(1, len(bins)) if (k == i).any()]
        rot = f.stem.replace(p.name + "_if_", "")
        for a, campo in zip(ax.flat, ("axial-velocity", "radial-velocity", "tangential-velocity", "pressure")):
            if campo in d:
                a.plot(rc, [d[campo][k == i].mean() for i in range(1, len(bins)) if (k == i).any()], ".-", label=rot, lw=1)
                a.set(xlabel="r [m]", title=campo)
    ax[0, 0].legend(fontsize=6)
    fig.suptitle(f"Perfis radiais nas interfaces (médias por faixa de raio), {p.name}")
    fig.tight_layout()
    fig.savefig(p / "graficos" / f"{p.name}_m4_perfis_interfaces.png")
    plt.close(fig)


def ponto_acoplamento(tab):
    """SR em que o torque do estator troca de sinal (interpolação linear).

    Args:
        tab (list[dict]): Desempenho por ponto, em SR crescente.

    Returns:
        float | None: SR do ponto de acoplamento, ou None se T_e não trocar de sinal.
    """
    for a, b in zip(tab, tab[1:]):
        if a["T_e"] * b["T_e"] < 0:
            return a["SR"] + (b["SR"] - a["SR"]) * a["T_e"] / (a["T_e"] - b["T_e"])
    return None


def marcar_acoplamento(ax, src):
    """Linha vertical no ponto de acoplamento e sombra na faixa acima (estator fixo não é físico ali).

    Args:
        ax (Iterable[matplotlib.axes.Axes]): Eixos.
        src (float | None): SR do ponto de acoplamento.
    """
    if src is None:
        return
    for a in ax:
        a.axvline(src, color="gray", ls=":", lw=1)
        a.axvspan(src, 1.0, color="gray", alpha=0.08)


def curvas(tab):
    """Curvas globais contra SR.

    Args:
        tab (list[dict]): Desempenho por ponto.
    """
    g = DIR / "graficos"
    sr = np.array([t["SR"] for t in tab])
    fig, ax = plt.subplots(1, 3, figsize=(13, 3.8))
    ax[0].errorbar(sr, [t["TR"] for t in tab], yerr=[t.get("dp_TR", 0) for t in tab], fmt="o-", capsize=3)
    ax[0].axhline(1, color="gray", ls="--")
    ax[0].set(xlabel="SR = ω_t/ω_b", ylabel="TR = −T_t/T_b", title="Razão de torque")
    ax[1].plot(sr, [t["eta"] for t in tab], "o-", color="C1")
    ax[1].set(xlabel="SR", ylabel="η = SR·TR", title="Eficiência")
    ax[2].plot(sr, [t["K"] for t in tab], "o-", color="C2")
    ax[2].set(xlabel="SR", ylabel="K = n_b/√|T_b| [rpm/√(N·m)]", title="Fator de capacidade")
    src = ponto_acoplamento(tab)
    marcar_acoplamento(ax, src)
    if src is not None:
        ax[0].text(src, ax[0].get_ylim()[1], f" acoplamento SR≈{src:.2f}" + chr(10) + " (acima: estator fixo, não físico)",
                   va="top", fontsize=7, color="gray")
    fig.tight_layout()
    fig.savefig(g / "g1_TR_eta_K.png")
    plt.close(fig)

    fig, ax = plt.subplots(1, 3, figsize=(13, 3.8))
    for n, k in (("bomba", "T_b"), ("turbina", "T_t"), ("estator", "T_e")):
        ax[0].plot(sr, [t[k] for t in tab], "o-", color=CORES[n], label=k)
    ax[0].axhline(0, color="gray", lw=0.8)
    ax[0].set(xlabel="SR", ylabel="torque [N·m]", title="Torques do fluido nas peças")
    ax[0].legend()
    ax[1].plot(sr, [t["mdot_bt"] for t in tab], "o-", color="C4")
    ax[1].set(xlabel="SR", ylabel="ṁ [kg/s]", title="Vazão de circulação (if bomba→turbina)")
    ax[2].plot(sr, [t["balanco_pct"] for t in tab], "o-", color="k")
    ax[2].set(xlabel="SR", ylabel="(T_b+T_t+T_e)/|T_b| [%]", title="Erro do balanço de torque")
    marcar_acoplamento(ax, ponto_acoplamento(tab))
    fig.tight_layout()
    fig.savefig(g / "g2_torques_vazao_balanco.png")
    plt.close(fig)

    if any(("euler_fluido_bomba" in t) or ("euler_csv_fluido_bomba" in t) for t in tab):
        fig, ax = plt.subplots(1, 3, figsize=(13, 3.8))
        for a, zona in zip(ax, LADOS):
            a.plot(sr, [t.get("euler_" + zona, np.nan) for t in tab], "s--", label="Euler, report do Fluent")
            a.plot(sr, [t.get("euler_csv_" + zona, np.nan) for t in tab], "^:", label="Euler, CSV (valores de célula, aprox.)")
            a.plot(sr, [t.get("parede_" + zona, np.nan) for t in tab], "o-", label="paredes da zona")
            a.set(xlabel="SR", ylabel="torque [N·m]", title=zona)
            a.legend(fontsize=7)
        fig.suptitle("Torque nas paredes × Euler (fluxo líquido de momento angular nas interfaces), por zona")
        fig.tight_layout()
        fig.savefig(g / "g3_torque_vs_euler.png")
        plt.close(fig)


def nota(tab, yplus, figs):
    """Escreve a nota única do Obsidian e copia as imagens.

    Args:
        tab (list[dict]): Desempenho por ponto.
        yplus (dict): SR -> estatísticas de y⁺.
        figs (list[Path]): Imagens a incluir.
    """
    pasta = VAULT / "05_Conversor"
    pasta.mkdir(parents=True, exist_ok=True)
    img = VAULT / "img"
    img.mkdir(exist_ok=True)
    ef = json.loads((DIR / "dados" / "setup.json").read_text()) if (DIR / "dados" / "setup.json").exists() else {}
    L = ["# CT_V1_MRF: conversor de torque, MRF, varredura de SR", "",
         "Gerada por `tools/analise_conversor.py`. Histórico e decisões em [[C_TG]].", "",
         "## Parâmetros de entrada",
         f"- Malha `{ef.get('malha', 'CT_V1_3zonas.msh.h5')}`: 1 900 736 células poly-hexcore sem prismas, 360°, eixo x pela origem.",
         f"- Bomba a {ef.get('rpm_b', 2000):.0f} rpm (ω_b = {ef.get('omega_b', 209.44):.2f} rad/s), sentido do eixo {ef.get('sinal', '?')}·x (teste de sentido em SR = 0). Turbina a SR·ω_b; estator parado.",
         f"- Óleo: ρ = {ef.get('rho', 830)} kg/m³, μ = {ef.get('mu', 0.007)} Pa·s (ATF a ~90 °C; valor típico, não do Kotwicki).",
         "- k-ω SST com correção de curvatura (TG 2.7); PRESTO! (TG 4.3); SIMPLE com URFs padrão: o Coupled pedia ~17,6 GB "
         "nesta malha, acima dos 16 GB de RAM, e o TG 4.4 o recomenda só quando a memória permite; "
         "SOU em momento e turbulência, com 200 iterações iniciais de primeira ordem (TG 2.8); formulação de velocidade relativa (TG 2.6); "
         "inicialização standard absoluta, k = 1 m²/s², ω = 500 1/s (TG 2.8).",
         "- Interfaces `if_*` → 3 mesh interfaces não conformes one-to-one; faces sem par viram paredes `sobra_*`, que giram com a peça do outro lado.",
         "- Paredes `parede_*_lado_*` e `sobra_*`: Moving Wall rotacional absoluta com ω da peça dona.",
         "- Parada: a média dos torques nas últimas 200 iterações muda < 0,5 % em relação às 200 anteriores (T_e relativo a |T_b|). "
         "O critério de máx−mín < 0,5 % foi abandonado porque os torques oscilam ~1–2 % sem deriva (regime estatisticamente estacionário). "
         "Resíduos só monitorados (TG 2.8, 4.4).",
         "", "## Desempenho", "",
         "Valores: médias das últimas 200 iterações de cada ponto; σ = desvio-padrão nessa janela (incerteza iterativa).", "",
         "| SR | T_b [N·m] | σ T_b | T_t [N·m] | σ T_t | T_e [N·m] | TR | η | K [rpm/√(N·m)] | balanço [%] | ṁ [kg/s] | it | média estável |",
         "|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for t in tab:
        L.append(f"| {t['SR']:.2f} | {t['T_b']:.1f} | {t['dp_T_b_pct']:.1f} % | {t['T_t']:.1f} | {t['dp_T_t_pct']:.1f} % | {t['T_e']:.1f} | {t['TR']:.3f} | {t['eta']:.3f} | "
                 f"{t['K']:.1f} | {t['balanco_pct']:.2f} | {t['mdot_bt'] or float('nan'):.1f} | {t['iter']} | {'sim' if t['conv'] else 'não'} |")
    if yplus:
        L += ["", "## y⁺ nas paredes", "", "| SR | peça | média | mediana | p5 | p95 | máx |", "|---|---|---|---|---|---|---|"]
        for sr, e in yplus.items():
            for peca, v in e.items():
                L.append(f"| {sr} | {peca} | {v['media']:.1f} | {v['mediana']:.1f} | {v['p5']:.1f} | {v['p95']:.1f} | {v['max']:.0f} |")
    src = ponto_acoplamento(tab)
    if src is not None:
        L += ["", f"**Ponto de acoplamento (T_e = 0): SR ≈ {src:.3f}.** Acima dele o estator real giraria livre (roda-livre); "
              "com o estator fixo do MRF, os pontos dessa faixa (TR < 1) não representam o conversor e ficam só como referência."]
    L += ["", "## Figuras", ""]
    for f in figs:
        destino = img / f"CT_V1_MRF_{f.name}" if not f.name.startswith("CT_V1_MRF") else img / f.name
        shutil.copy2(f, destino)
        L.append(f"- {f.stem}: ![[{destino.name}]]")
    (pasta / "CT_V1_MRF.md").write_text("\n".join(L) + "\n", encoding="utf-8")


def main(com_nota=True):
    """Roda a análise completa.

    Args:
        com_nota (bool, optional): Se True, escreve a nota do Obsidian. Padrão True.
    """
    nomes, m = monitores()
    logs = sorted((DIR / "logs").glob("varredura*.log")) + sorted((DIR / "logs").glob("sentido*.log"))
    eqs, res = [], np.zeros((0, 0))
    for lg in logs:
        e, r_ = residuos(lg)
        if len(r_):
            eqs = e
            res = r_ if not len(res) else np.vstack([res, r_]) if r_.shape[1] == res.shape[1] else res
    tab, yplus, figs = [], {}, []
    for p, r in pontos():
        if len(m):
            graf_monitores(nomes, m, p, r)
        graf_residuos(eqs, res, p, r)
        yplus[p.name] = graf_yplus(p)
        graf_perfis(p)
        d = desempenho(r)
        for zona, v in euler_csv(p, r.get("sinal", 1)).items():
            d[("euler_csv_" if zona.startswith("fluido") else "csv_") + zona] = v
        tab.append(d)
        figs += sorted((p / "graficos").glob("*.png"))
    if tab:
        curvas(tab)
        with open(DIR / "dados" / "desempenho.csv", "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(tab[0].keys()) + [k for t in tab for k in t if k not in tab[0]],
                               delimiter=";", extrasaction="ignore")
            w.writeheader()
            w.writerows(tab)
    figs = sorted((DIR / "graficos").glob("*.png")) + figs
    for t in tab:
        print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in t.items()})
    if com_nota:
        nota(tab, yplus, figs)
        print("nota:", VAULT / "05_Conversor" / "CT_V1_MRF.md", "imagens:", len(figs))


if __name__ == "__main__":
    main("--sem-nota" not in sys.argv)
