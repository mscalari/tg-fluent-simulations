"""Teste de sinais e monotonicidade das curvas de desempenho do conversor CT_V1 (MRF).

Uso: python tools/tendencia_sinais_conversor.py

Pergunta respondida: os SINAIS das variações entre pontos vizinhos de SR batem com o esperado para um conversor
de torque, e a variação é maior do que a incerteza iterativa dos dois pontos?

Critérios por intervalo [SR_i, SR_i+1]:
  * TR cai sempre (ΔTR < 0);
  * η = SR·TR sobe até o SR de η máximo (da própria série) e cai depois (Δη > 0 antes, Δη < 0 depois);
  * K = n_b/√|T_b| sobe com SR (ΔK > 0);
  * TR = 1 é cruzado uma única vez (ponto de acoplamento).

Classificação de cada intervalo da simulação:
  * "confirmado": sinal certo e |Δ| > U_Δ = √(U_i² + U_j²);
  * "sinal certo, dentro da incerteza": sinal certo, mas |Δ| ≤ U_Δ;
  * "sinal errado": sinal contrário ao esperado.
A incerteza U de cada ponto é U = 2·SEM (fator de abrangência k = 2, ≈ 95 %), com o SEM por lotes (batch means)
da série por iteração do report file (``gci_conversor.serie_janela`` e ``grandezas_ponto``). Só entra a incerteza
iterativa (ruído na janela de médias); erro de discretização e de modelo NÃO entram em U.

As mesmas regras valem para cada curva da literatura (``refs/conversor_literatura.csv``, lida por
``comparacao_literatura_conversor.ler_literatura``). Lá não há SEM; usa-se como U a tolerância de digitalização
(``TOL_LIT``: TR ±0,02, η ±0,01, K ±2 % do valor; U_Δ = √2·u), de modo que patamares dentro do erro de leitura caem
em "dentro da incerteza" e não em "sinal errado".

Saídas:
  * ``results/conversor/CT_V1_MRF/dados/sinais_tendencia.csv`` (separador ``;``);
  * ``results/conversor/CT_V1_MRF/graficos/g9_sinais_tendencia.png``.
"""
import csv
import json
import math
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

sys.path.insert(0, str(Path(__file__).resolve().parent))
import comparacao_literatura_conversor as lit  # noqa: E402
import gci_conversor as gci  # noqa: E402

RAIZ = Path(__file__).resolve().parents[1]
BASE = RAIZ / "results" / "conversor" / "CT_V1_MRF"
DES = BASE / "dados" / "desempenho.csv"
SAIDA_CSV = BASE / "dados" / "sinais_tendencia.csv"
SAIDA_PNG = BASE / "graficos" / "g9_sinais_tendencia.png"
K_ABRANGENCIA = 2.0                       # U = 2·SEM
# Tolerância de digitalização da literatura (leitura de gráfico): u de TR e η em valor absoluto, de K em fração do
# valor (±0,05 em CC ≈ 1 %, ±0,07 ≈ 1,2 %; adota-se 2 % por segurança). U_Δ = √2·u.
TOL_LIT = {"TR": 0.02, "eta": 0.01, "K": 0.02}
GRANDEZAS = ("TR", "eta", "K")
ROT_G = {"TR": "TR", "eta": "η", "K": "K"}
NOME_FONTE = {"zhang2023": "Zhang 2023", "guo2021": "Guo 2021", "ran2021": "Ran 2021",
              "nadeem2020": "Nadeem 2020", "kotwicki_forma_nadeem2020": "Kotwicki (modelo)"}
NOME_ORIGEM = {"exp": "exp.", "cfd_rans": "RANS", "cfd_les": "LES", "modelo": "modelo"}
COR = {"confirmado": "#4caf50", "dentro": "#ffd54f", "errado": "#e57373", "nd": "#e0e0e0"}
CLASSES = {"confirmado": "confirmado", "dentro": "sinal certo, dentro da incerteza", "errado": "sinal errado",
           "nd": "não definido"}
plt.rcParams.update({"font.size": 9})


def br(x, casas=3, sinal=False):
    """Formata um número com vírgula decimal e sinal "−" tipográfico.

    Args:
        x (float): Valor.
        casas (int): Casas decimais.
        sinal (bool): Se True, mostra "+" nos positivos.

    Returns:
        str: Texto em pt-BR (``"n/d"`` se ``x`` não for finito).
    """
    if x is None or not np.isfinite(x):
        return "n/d"
    s = f"{x:{'+' if sinal else ''}.{casas}f}".replace(".", ",")
    return s.replace("-", "−")


def sinal_esperado(grandeza, sr_i, sr_j, sr_etamax):
    """Sinal esperado de Δ em um intervalo.

    Args:
        grandeza (str): ``"TR"``, ``"eta"`` ou ``"K"``.
        sr_i (float): SR do início do intervalo.
        sr_j (float): SR do fim do intervalo.
        sr_etamax (float): SR de η máximo da série.

    Returns:
        int: -1 (cai) ou +1 (sobe).
    """
    if grandeza == "TR":
        return -1
    if grandeza == "K":
        return +1
    return +1 if sr_j <= sr_etamax + 1e-12 else -1


def classificar(delta, esperado, u_delta):
    """Classifica a variação de um intervalo.

    Args:
        delta (float): Δ = valor_j − valor_i.
        esperado (int): Sinal esperado (+1 ou -1).
        u_delta (float): Incerteza U_Δ da variação (NaN se não houver).

    Returns:
        str: ``"confirmado"``, ``"dentro"``, ``"errado"`` ou ``"nd"``; sem incerteza, o melhor que sai é
        ``"certo"`` (sinal certo sem teste de incerteza).
    """
    if not np.isfinite(delta):
        return "nd"
    if np.sign(delta) != esperado:
        return "errado"
    if not np.isfinite(u_delta):
        return "certo"
    return "confirmado" if abs(delta) > u_delta else "dentro"


def cruzamentos(sr, tr):
    """Conta e localiza os cruzamentos de TR = 1.

    Args:
        sr (np.ndarray): SR crescente.
        tr (np.ndarray): TR de cada ponto.

    Returns:
        list[float]: SR de cada cruzamento por interpolação linear (vazia se não houver).
    """
    out = []
    for i in range(len(sr) - 1):
        a, b = tr[i] - 1.0, tr[i + 1] - 1.0
        if np.isfinite(a) and np.isfinite(b) and a * b < 0:
            out.append(float(sr[i] + (1.0 - tr[i]) * (sr[i + 1] - sr[i]) / (tr[i + 1] - tr[i])))
    return out


def ler_simulacao():
    """Lê os 8 pontos da simulação com a incerteza U = 2·SEM de TR, η e K.

    Os valores vêm de ``desempenho.csv``; o SEM, de ``gci_conversor.grandezas_ponto`` com a série por iteração
    (``serie_janela``). Se a série não for achada, cai no SEM conservador por σ da janela (e marca o aviso).

    Returns:
        tuple[dict, list[str]]: ``{"SR", "TR", "eta", "K", "U_TR", "U_eta", "U_K"}`` (arrays) e a lista de avisos.

    Raises:
        FileNotFoundError: Se ``desempenho.csv`` ou algum ``resultado.json`` não existir.
    """
    with open(DES, encoding="utf-8") as f:
        linhas = sorted(csv.DictReader(f, delimiter=";"), key=lambda l: float(l["SR"]))
    d = {k: [] for k in ("SR", "TR", "eta", "K", "U_TR", "U_eta", "U_K")}
    avisos = []
    for l in linhas:
        sr = float(l["SR"])
        arq = gci.pasta_sr(gci.PASTAS["G0"], sr) / "dados" / "resultado.json"
        r = json.loads(arq.read_text(encoding="utf-8"))
        serie, info = gci.serie_janela("G0", r, sr)
        if serie is None:
            avisos.append(f"SR {br(sr, 2)}: série por iteração não achada ({info}); SEM conservador por σ da janela")
        g = gci.grandezas_ponto(r, serie)
        for q in GRANDEZAS:
            v_csv, v_json = float(l[q]), g[q]["v"]
            if abs(v_csv - v_json) > 1e-6 * max(abs(v_csv), 1e-12):
                avisos.append(f"SR {br(sr, 2)}: {q} do desempenho.csv ({v_csv:.6g}) difere do resultado.json ({v_json:.6g})")
            d[q].append(v_csv)
            d["U_" + q].append(K_ABRANGENCIA * g[q]["sem"])
        d["SR"].append(sr)
    return {k: np.array(v) for k, v in d.items()}, avisos


def intervalos_simulacao(d):
    """Aplica o teste de sinais aos intervalos da simulação.

    Args:
        d (dict): Saída de ``ler_simulacao``.

    Returns:
        tuple[list[dict], dict]: Uma linha por intervalo e grandeza, e o resumo (SR de η máximo, acoplamento,
        confirmação do cruzamento).
    """
    sr = d["SR"]
    i_max = int(np.argmax(d["eta"]))
    sr_max = float(sr[i_max])
    linhas = []
    for i in range(len(sr) - 1):
        for q in GRANDEZAS:
            delta = float(d[q][i + 1] - d[q][i])
            u = math.hypot(d["U_" + q][i], d["U_" + q][i + 1])
            esp = sinal_esperado(q, sr[i], sr[i + 1], sr_max)
            cl = classificar(delta, esp, u)
            linhas.append({"fonte": "CT_V1 MRF", "origem": "sim", "grandeza": q, "SR_i": sr[i], "SR_j": sr[i + 1],
                           "valor_i": d[q][i], "valor_j": d[q][i + 1], "delta": delta, "U_delta": u,
                           "delta_sobre_U": abs(delta) / u if u > 0 else float("nan"),
                           "sinal_esperado": "+" if esp > 0 else "−", "classe": cl})
    cr = cruzamentos(sr, d["TR"])
    resumo = {"sr_etamax": sr_max, "eta_max": float(d["eta"][i_max]), "pico_no_limite": i_max in (0, len(sr) - 1),
              "cruzamentos": cr}
    if len(cr) == 1:
        i = next(k for k in range(len(sr) - 1) if (d["TR"][k] - 1) * (d["TR"][k + 1] - 1) < 0)
        resumo["cruz_confirmado"] = bool(abs(d["TR"][i] - 1) > d["U_TR"][i] and abs(d["TR"][i + 1] - 1) > d["U_TR"][i + 1])
        resumo["cruz_margem"] = (float(abs(d["TR"][i] - 1) / d["U_TR"][i]), float(abs(d["TR"][i + 1] - 1) / d["U_TR"][i + 1]))
    return linhas, resumo


def intervalos_literatura():
    """Aplica o teste de sinais a cada curva da literatura (com a tolerância de digitalização como U).

    Tolerância de digitalização ``TOL_LIT`` entra como U_Δ. Para o K usa-se K/K_estol (``KK``) quando existir; o sinal de ΔK não depende da normalização.

    Returns:
        tuple[list[dict], list[dict]]: Linhas por intervalo e grandeza, e o resumo por série (rótulo, SR de η
        máximo, cruzamentos, contagens ``[confirmado, dentro da tolerância, errado, avaliados]``).
    """
    linhas, resumos = [], []
    for s in lit.ler_literatura():
        sr = s["SR"]
        eta = s["eta"]
        ok_eta = np.isfinite(eta)
        sr_max = float(sr[ok_eta][int(np.argmax(eta[ok_eta]))]) if ok_eta.any() else float("nan")
        val = {"TR": s["TR"], "eta": eta, "K": s["KK"]}
        rot = f"{NOME_FONTE.get(s['fonte'], s['fonte'])} {NOME_ORIGEM.get(s['origem'], s['origem'])}"
        if s["rotulo"].endswith(")") and s["origem"] != "exp":
            rot += " " + s["rotulo"][s["rotulo"].index("("):]
        cont = {q: [0, 0, 0, 0] for q in GRANDEZAS}  # [confirmado, dentro da tolerância, errado, avaliados]
        for q in GRANDEZAS:
            for i in range(len(sr) - 1):
                delta = float(val[q][i + 1] - val[q][i])
                if not np.isfinite(delta):
                    continue
                esp = sinal_esperado(q, sr[i], sr[i + 1], sr_max)
                u_rel = TOL_LIT[q] * (abs(val[q][i]) + abs(val[q][i + 1])) / 2 if q == "K" else TOL_LIT[q]
                u_d = math.sqrt(2) * u_rel
                cl = classificar(delta, esp, u_d)
                if cl == "errado" and abs(delta) <= u_d:
                    cl = "dentro"                   # patamar: variação dentro do erro de digitalização
                cont[q][3] += 1
                cont[q][("confirmado", "dentro", "errado").index(cl)] += 1
                linhas.append({"fonte": s["fonte"], "origem": s["origem"], "grandeza": q, "SR_i": float(sr[i]),
                               "SR_j": float(sr[i + 1]), "valor_i": float(val[q][i]), "valor_j": float(val[q][i + 1]),
                               "delta": delta, "U_delta": u_d, "delta_sobre_U": abs(delta) / u_d,
                               "sinal_esperado": "+" if esp > 0 else "−", "classe": cl})
        ok_tr = np.isfinite(s["TR"])
        cr = cruzamentos(sr[ok_tr], s["TR"][ok_tr])
        resumos.append({"rotulo": rot, "fonte": s["fonte"], "origem": s["origem"], "sr_etamax": sr_max,
                        "pico_no_limite": bool(ok_eta.any() and sr_max in (float(sr[ok_eta][0]), float(sr[ok_eta][-1]))),
                        "cruzamentos": cr, "cont": cont})
    ordem = {"zhang2023": 0, "guo2021": 1, "ran2021": 2, "nadeem2020": 3, "kotwicki_forma_nadeem2020": 4}
    resumos.sort(key=lambda r: (ordem.get(r["fonte"], 9), r["origem"]))
    return linhas, resumos


def gravar_csv(linhas):
    """Grava ``sinais_tendencia.csv`` (separador ``;``, ponto decimal).

    Args:
        linhas (list[dict]): Linhas por intervalo e grandeza (simulação e literatura).

    Raises:
        OSError: Se a pasta de saída não for gravável.
    """
    cols = ["fonte", "origem", "grandeza", "SR_i", "SR_j", "valor_i", "valor_j", "delta", "U_delta",
            "delta_sobre_U", "sinal_esperado", "classe"]
    SAIDA_CSV.parent.mkdir(parents=True, exist_ok=True)
    with open(SAIDA_CSV, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter=";")
        w.writerow(cols)
        for l in linhas:
            w.writerow([l[c] if not isinstance(l[c], float) else f"{l[c]:.6g}" for c in cols])


def desenhar(d, lin_sim, res_sim, res_lit):
    """Desenha a figura g9: grade de intervalos da simulação e linha por curva da literatura.

    Args:
        d (dict): Pontos da simulação (``ler_simulacao``).
        lin_sim (list[dict]): Intervalos da simulação.
        res_sim (dict): Resumo da simulação.
        res_lit (list[dict]): Resumo por série da literatura.

    Raises:
        OSError: Se a pasta de saída não for gravável.
    """
    sr = d["SR"]
    n_int = len(sr) - 1
    n_lit = len(res_lit)
    fig = plt.figure(figsize=(12.5, 1.9 + 0.62 * n_int + 0.5 * n_lit + 1.6))
    gs = fig.add_gridspec(2, 1, height_ratios=[n_int + 1.6, n_lit + 1.6], hspace=0.22,
                          left=0.03, right=0.99, top=0.93, bottom=0.10)

    # --- painel A: simulação ---
    ax = fig.add_subplot(gs[0])
    ax.set_xlim(0, 1 + 3 * 3.0 + 0.0)
    ax.set_xlim(-1.6, 10.6)
    ax.set_ylim(n_int + 1.7, -1.0)
    ax.axis("off")
    ax.set_title("A. CT_V1 (MRF): intervalo a intervalo, Δ e Δ/U_Δ  (U = 2·SEM por lotes, só incerteza iterativa)",
                 loc="left", fontsize=10, fontweight="bold")
    x0, larg = 1.6, 3.0
    for c, q in enumerate(GRANDEZAS):
        txt = {"TR": "TR  (esperado ▼)", "eta": f"η  (esperado ▲ até SR = {br(res_sim['sr_etamax'], 2)})",
               "K": "K  (esperado ▲)"}[q]
        ax.text(x0 + c * larg + larg / 2, -0.55, txt, ha="center", va="center", fontweight="bold")
    por = {(l["grandeza"], l["SR_i"]): l for l in lin_sim}
    for i in range(n_int):
        ax.text(x0 - 0.1, i + 0.5, f"{br(sr[i], 2)} → {br(sr[i + 1], 2)}", ha="right", va="center")
        for c, q in enumerate(GRANDEZAS):
            l = por[(q, sr[i])]
            ax.add_patch(Rectangle((x0 + c * larg, i), larg - 0.06, 0.94, fc=COR[l["classe"]], ec="white"))
            seta = "▲" if l["delta"] > 0 else "▼"
            ax.text(x0 + c * larg + larg / 2, i + 0.47, f"{seta} Δ = {br(l['delta'], 3, True)}   Δ/U = {br(l['delta_sobre_U'], 1)}",
                    ha="center", va="center", fontsize=8.5)
    n_tot = {q: {k: sum(1 for l in lin_sim if l["grandeza"] == q and l["classe"] == k) for k in COR} for q in GRANDEZAS}
    ax.text(x0 - 0.1, n_int + 0.45, "Contagem", ha="right", va="center", fontweight="bold")
    for c, q in enumerate(GRANDEZAS):
        t = n_tot[q]
        ax.text(x0 + c * larg + larg / 2, n_int + 0.45,
                f"{t['confirmado']} confirmados · {t['dentro']} na incerteza · {t['errado']} errados",
                ha="center", va="center", fontsize=8.5, fontweight="bold")
    cr = res_sim["cruzamentos"]
    if len(cr) == 1:
        tag = "confirmado" if res_sim.get("cruz_confirmado") else "não confirmado pela incerteza"
        cruz = f"TR = 1 cruzado 1 vez, em SR ≈ {br(cr[0], 2)} ({tag})"
    else:
        cruz = f"TR = 1 cruzado {len(cr)} vezes"
    lim = " (no limite da faixa: pico não resolvido)" if res_sim["pico_no_limite"] else ""
    ax.text(x0, n_int + 1.25, f"{cruz}\nη máx = {br(res_sim['eta_max'], 3)} em SR = {br(res_sim['sr_etamax'], 2)}{lim}",
            ha="left", va="center", fontsize=9)

    # --- painel B: literatura ---
    bx = fig.add_subplot(gs[1])
    bx.set_xlim(-1.6, 10.6)
    bx.set_ylim(n_lit + 0.9, -1.0)
    bx.axis("off")
    bx.set_title("B. Literatura (tolerância de digitalização como U): confirmados[+na tolerância]/avaliados; vermelho = algum errado",
                 loc="left", fontsize=10, fontweight="bold")
    cab = [("TR ▼", 2.0), ("η ▲/▼", 2.0), ("K ▲", 2.0), ("η máx em SR", 1.6), ("cruza TR = 1", 1.6)]
    x = x0
    xs = []
    for nome, w in cab:
        bx.text(x + w / 2, -0.55, nome, ha="center", va="center", fontweight="bold")
        xs.append((x, w))
        x += w
    for k, r in enumerate(res_lit):
        bx.text(x0 - 0.1, k + 0.5, r["rotulo"], ha="right", va="center")
        for c, q in enumerate(GRANDEZAS):
            ok, dn, er, tot = r["cont"][q]
            cl = "nd" if tot == 0 else ("errado" if er else ("confirmado" if ok == tot else "dentro"))
            bx.add_patch(Rectangle((xs[c][0], k), xs[c][1] - 0.05, 0.9, fc=COR[cl], ec="white"))
            bx.text(xs[c][0] + xs[c][1] / 2, k + 0.45, "n/d" if tot == 0 else f"{ok}+{dn}/{tot}" if dn else f"{ok}/{tot}", ha="center", va="center")
        bx.add_patch(Rectangle((xs[3][0], k), xs[3][1] - 0.05, 0.9, fc="#eceff1", ec="white"))
        bx.text(xs[3][0] + xs[3][1] / 2, k + 0.45,
                "n/d" if not np.isfinite(r["sr_etamax"]) else br(r["sr_etamax"], 2) + ("*" if r["pico_no_limite"] else ""),
                ha="center", va="center")
        n_c = len(r["cruzamentos"])
        bx.add_patch(Rectangle((xs[4][0], k), xs[4][1] - 0.05, 0.9, fc="#eceff1", ec="white"))
        bx.text(xs[4][0] + xs[4][1] / 2, k + 0.45,
                "não" if n_c == 0 else (f"{br(r['cruzamentos'][0], 2)}" if n_c == 1 else f"{n_c} vezes"),
                ha="center", va="center")
    # legenda
    handles = [Rectangle((0, 0), 1, 1, fc=COR[k]) for k in ("confirmado", "dentro", "errado", "nd")]
    nomes = ["confirmado (sinal certo e |Δ| > U_Δ; na literatura: |Δ| > √2·tolerância)",
             "sinal certo (ou patamar), dentro da incerteza/tolerância", "sinal errado",
             "sem dados"]
    fig.legend(handles, nomes, loc="lower center", bbox_to_anchor=(0.5, 0.02), ncol=2, frameon=False, fontsize=8.5)
    fig.text(0.99, 0.004, "* pico de η no extremo da faixa digitalizada (não resolvido).  "
             "K da literatura: K/K_estol ou fator de capacidade, só o sinal importa.",
             ha="right", va="bottom", fontsize=7.5, color="#555555")
    fig.suptitle("Teste de sinais e monotonicidade das curvas do conversor CT_V1 (MRF) frente à literatura",
                 fontsize=12, fontweight="bold", y=0.985)
    SAIDA_PNG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(SAIDA_PNG, dpi=150)
    plt.close(fig)


def main():
    """Roda o teste, grava o CSV e a figura e imprime o veredito por grandeza.

    Raises:
        FileNotFoundError: Se faltarem os dados de entrada.
    """
    d, avisos = ler_simulacao()
    lin_sim, res_sim = intervalos_simulacao(d)
    lin_lit, res_lit = intervalos_literatura()
    gravar_csv(lin_sim + lin_lit)
    desenhar(d, lin_sim, res_sim, res_lit)
    for a in avisos:
        print("AVISO:", a)
    print("U = 2·SEM (batch means); pontos SR:", ", ".join(br(s, 2) for s in d["SR"]))
    for q in GRANDEZAS:
        c = {k: sum(1 for l in lin_sim if l["grandeza"] == q and l["classe"] == k) for k in COR}
        print(f"{ROT_G[q]}: {c['confirmado']} confirmados, {c['dentro']} sinal certo dentro da incerteza, "
              f"{c['errado']} sinal errado (de {sum(c.values())})")
        for l in lin_sim:
            if l["grandeza"] == q and l["classe"] != "confirmado":
                print(f"   SR {br(l['SR_i'], 2)}→{br(l['SR_j'], 2)}: Δ = {br(l['delta'], 4, True)}, "
                      f"U_Δ = {br(l['U_delta'], 4)}, {l['classe']}")
    print(f"η máx = {br(res_sim['eta_max'], 3)} em SR = {br(res_sim['sr_etamax'], 2)} "
          f"(no limite da faixa: {res_sim['pico_no_limite']})")
    print("Cruzamentos de TR = 1:", [br(x, 3) for x in res_sim["cruzamentos"]],
          "| confirmado pela incerteza:", res_sim.get("cruz_confirmado"), "| |TR−1|/U nos vizinhos:",
          [br(x, 1) for x in res_sim.get("cruz_margem", ())])
    print("Literatura ([confirmado, dentro, errado]/avaliados; η máx; cruzamentos):")
    for r in res_lit:
        c = r["cont"]
        print(f"  {r['rotulo']:<28} TR {c['TR'][:3]}/{c['TR'][3]}  η {c['eta'][:3]}/{c['eta'][3]}  "
              f"K {c['K'][:3]}/{c['K'][3]}  η máx SR={br(r['sr_etamax'], 2)}{'*' if r['pico_no_limite'] else ''}  "
              f"cruza: {[br(x, 2) for x in r['cruzamentos']]}")
    print("Saídas:", SAIDA_CSV, SAIDA_PNG, sep="\n  ")


if __name__ == "__main__":
    main()
