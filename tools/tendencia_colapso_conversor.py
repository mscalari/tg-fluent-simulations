"""Colapso adimensional das curvas de desempenho do CT_V1 (MRF) contra a literatura (só a FORMA).

Uso (``.venv241`` e ``PYTHONIOENCODING=utf-8``)::

    python tools/tendencia_colapso_conversor.py

Entradas (o Fluent NÃO é executado)::

    results/conversor/CT_V1_MRF/dados/desempenho.csv        8 pontos (SR, TR, η)
    results/conversor/CT_V1_MRF/SRxxx/dados/resultado.json  janela de média de cada ponto
    results/conversor/CT_V1_MRF/dados/*.out                 séries por iteração (SEM por lotes)
    refs/conversor_literatura.csv                           curvas digitalizadas (ver refs/conversor_literatura.md)

Saídas::

    results/conversor/CT_V1_MRF/graficos/g10_colapso_adimensional.png   2 painéis (TR e η colapsados)
    results/conversor/CT_V1_MRF/dados/colapso_adimensional.csv          1 linha por ponto (separador ;)
    results/conversor/CT_V1_MRF/dados/colapso_adimensional_resumo.csv   1 linha por série e painel

Normalização::

    painel TR : x = SR/SR_acopl,  y = (TR − 1)/(TR_estol − 1)    (toda curva vai de (0, 1) a (1, 0))
    painel η  : x = SR/SR*,       y = η/η_máx                    (SR*, η_máx: vértice da parábola nos 3 pontos
                                                                  em torno do máximo discreto)

* SR_acopl (TR = 1) por interpolação linear entre os pontos vizinhos (é o valor usado) e, como verificação, por
  PCHIP monotônico (Fritsch-Carlson, implementado aqui; o ``.venv241`` não tem scipy). A diferença entra como
  incerteza de interpolação.
* Incerteza do CT_V1: SEM por lotes de cada ponto (``gci_conversor.grandezas_ponto`` + ``serie_janela``),
  propagada por Monte Carlo (TR_i ~ N(TR_i, SEM_i), independentes entre pontos) para SR_acopl, x e y. As barras
  são o intervalo de 95 % (percentis 2,5 e 97,5); em x soma-se, em quadratura, a diferença linear × PCHIP.
  η = SR·TR, então η e TR de um ponto são amostrados juntos.
* Literatura: mesma normalização por série. TR_estol é o TR do menor SR da série (SR = 0,01 em Zhang). Entram no
  envelope (faixa mín–máx) do painel TR só as séries que cruzam TR = 1 (Zhang ×3, Nadeem). Guo e Ran (SR ≤ 0,8,
  TR > 1) aparecem tracejadas, com SR_acopl EXTRAPOLADO linearmente pelos dois últimos pontos, e ficam fora do
  envelope. No painel η entram as séries cujo máximo é interior (Zhang ×3 com η = SR·TR no LES, Guo CFD ×2,
  Ran ×3); Guo exp (máximo no último ponto, SR = 0,8) e Nadeem (η cresce até o fim, modo acoplamento) ficam de
  fora. A forma de Kotwicki (ajuste de Nadeem, não medição, sem ponto de estol) fica fora dos dois painéis.
* η_máx do CT_V1: η cresce até o último ponto (SR = 0,9), então o máximo NÃO está na faixa simulada. O painel η
  normaliza o CT_V1 pelo último ponto (SR* ≥ 0,9 e η_máx ≥ η(0,9)); os x e y desenhados são limites superiores.
  O vértice da parábola nos três últimos pontos é só diagnóstico (extrapolação longa) e não é usado.
* Métrica: para cada ponto do CT_V1 não fixado pela normalização e dentro do domínio do envelope, diz se o valor
  central está dentro da faixa e se a barra de 95 % a toca; fora, o desvio vertical até a borda mais próxima.
  Forma de cada curva: índice A = ∫₀¹ (y − corda) dx (corda de (0, 1) a (1, 0) no TR; de (0, 0) a (1, 1) no η);
  A > 0 é côncava (curva acima da corda), A < 0 convexa. Também se conta o sinal das segundas diferenças divididas.

Isto é validação de TENDÊNCIA (geometrias diferentes); a validação forte está nos blocos A–H e nas verificações do
conversor (balanço de torque, Euler nas interfaces, y+, sensibilidade, GCI).
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
from matplotlib.ticker import FuncFormatter

sys.path.insert(0, str(Path(__file__).resolve().parent))
import comparacao_literatura_conversor as lit_mod  # noqa: E402  (leitura da literatura e estilo do g5)
import gci_conversor as gci  # noqa: E402  (SEM por lotes de cada ponto)

RAIZ = Path(__file__).resolve().parents[1]
DIR = RAIZ / "results" / "conversor" / "CT_V1_MRF"
DES = DIR / "dados" / "desempenho.csv"
SAIDA_PNG = DIR / "graficos" / "g10_colapso_adimensional.png"
SAIDA_CSV = DIR / "dados" / "colapso_adimensional.csv"
SAIDA_RES = DIR / "dados" / "colapso_adimensional_resumo.csv"
N_MC = 20000                 # amostras do Monte Carlo
SEMENTE = 20261121
TOL_FORMA = 0.02             # |A| abaixo disto: "quase linear"
NAN = float("nan")


# ---------------------------------------------------------------------------------------------------------------
# Numérica
# ---------------------------------------------------------------------------------------------------------------
def pchip_derivadas(x, y):
    """Derivadas nos nós do interpolador PCHIP (Fritsch-Carlson, como ``scipy.interpolate.PchipInterpolator``).

    Args:
        x (np.ndarray): Abscissas estritamente crescentes (n ≥ 2).
        y (np.ndarray): Ordenadas.

    Returns:
        np.ndarray: Derivada em cada nó; zero onde a inclinação troca de sinal (preserva a monotonicidade).
    """
    h = np.diff(x)
    dl = np.diff(y) / h
    n = len(x)
    d = np.zeros(n)
    if n == 2:
        return np.full(2, dl[0])
    for k in range(1, n - 1):
        if dl[k - 1] * dl[k] > 0:
            w1, w2 = 2 * h[k] + h[k - 1], h[k] + 2 * h[k - 1]
            d[k] = (w1 + w2) / (w1 / dl[k - 1] + w2 / dl[k])

    def borda(h0, h1, m0, m1):
        """Derivada na borda pela fórmula de três pontos não centrada, com as correções de forma do PCHIP.

        Args:
            h0 (float): Passo do intervalo da borda.
            h1 (float): Passo do intervalo vizinho.
            m0 (float): Inclinação do intervalo da borda.
            m1 (float): Inclinação do intervalo vizinho.

        Returns:
            float: Derivada no nó da borda.
        """
        dd = ((2 * h0 + h1) * m0 - h0 * m1) / (h0 + h1)
        if np.sign(dd) != np.sign(m0):
            return 0.0
        if np.sign(m0) != np.sign(m1) and abs(dd) > abs(3 * m0):
            return 3 * m0
        return dd

    d[0] = borda(h[0], h[1], dl[0], dl[1])
    d[-1] = borda(h[-1], h[-2], dl[-1], dl[-2])
    return d


def pchip_avaliar(x, y, d, xq):
    """Avalia o PCHIP (base de Hermite cúbica) num ponto.

    Args:
        x (np.ndarray): Nós crescentes.
        y (np.ndarray): Valores nos nós.
        d (np.ndarray): Derivadas nos nós (``pchip_derivadas``).
        xq (float): Ponto, dentro de [x[0], x[-1]].

    Returns:
        float: Valor interpolado.
    """
    k = int(np.clip(np.searchsorted(x, xq) - 1, 0, len(x) - 2))
    h = x[k + 1] - x[k]
    t = (xq - x[k]) / h
    h00, h10, h01, h11 = 2 * t**3 - 3 * t**2 + 1, t**3 - 2 * t**2 + t, -2 * t**3 + 3 * t**2, t**3 - t**2
    return float(h00 * y[k] + h10 * h * d[k] + h01 * y[k + 1] + h11 * h * d[k + 1])


def cruzamento_linear(sr, tr, alvo=1.0):
    """Primeiro SR em que TR cruza ``alvo``, por interpolação linear entre os pontos vizinhos.

    Args:
        sr (np.ndarray): SR crescentes.
        tr (np.ndarray): TR.
        alvo (float): Valor de TR procurado.

    Returns:
        tuple[float, int]: (SR do cruzamento, índice i do intervalo [i, i+1]); (nan, −1) se não cruzar.
    """
    for i in range(len(sr) - 1):
        a, b = tr[i] - alvo, tr[i + 1] - alvo
        if a * b <= 0 and tr[i] != tr[i + 1]:
            return float(sr[i] + a * (sr[i + 1] - sr[i]) / (tr[i] - tr[i + 1])), i
    return NAN, -1


def cruzamento_pchip(sr, tr, alvo=1.0):
    """SR em que o PCHIP de TR(SR) vale ``alvo``, por bisseção no intervalo do cruzamento linear.

    Args:
        sr (np.ndarray): SR crescentes.
        tr (np.ndarray): TR.
        alvo (float): Valor de TR procurado.

    Returns:
        float: SR do cruzamento, ou nan se TR não cruzar ``alvo``.
    """
    _, i = cruzamento_linear(sr, tr, alvo)
    if i < 0:
        return NAN
    d = pchip_derivadas(sr, tr)
    a, b = float(sr[i]), float(sr[i + 1])
    fa = pchip_avaliar(sr, tr, d, a) - alvo
    for _ in range(80):
        m = 0.5 * (a + b)
        fm = pchip_avaliar(sr, tr, d, m) - alvo
        if fa * fm <= 0:
            b = m
        else:
            a, fa = m, fm
    return 0.5 * (a + b)


def cruzamento_extrapolado(sr, tr, alvo=1.0):
    """SR de TR = ``alvo`` extrapolando a reta dos dois últimos pontos (para séries que não chegam ao acoplamento).

    Args:
        sr (np.ndarray): SR crescentes.
        tr (np.ndarray): TR.
        alvo (float): Valor de TR procurado.

    Returns:
        float: SR extrapolado, ou nan se a reta não descer.
    """
    s = (tr[-1] - tr[-2]) / (sr[-1] - sr[-2])
    return float(sr[-1] + (alvo - tr[-1]) / s) if s < 0 else NAN


def maximo_parabola(x, y):
    """Máximo local por parábola nos três pontos em torno do máximo discreto.

    Args:
        x (np.ndarray): Abscissas crescentes.
        y (np.ndarray): Ordenadas.

    Returns:
        dict: ``x``, ``y`` (vértice; ou o ponto extremo se o máximo discreto estiver na borda), ``interior``
        (bool) e ``x_vertice_diag`` (vértice da parábola nos três últimos pontos, só diagnóstico, quando o
        máximo está na borda; nan se a parábola não tiver máximo).
    """
    j = int(np.nanargmax(y))

    def vertice(i0):
        """Vértice da parábola que passa por (x, y)[i0:i0+3].

        Args:
            i0 (int): Índice do primeiro dos três pontos.

        Returns:
            tuple[float, float]: (x, y) do vértice; (nan, nan) se a parábola for convexa ou degenerada.
        """
        a, b, c = np.polyfit(x[i0:i0 + 3], y[i0:i0 + 3], 2)
        if a >= 0:
            return NAN, NAN
        xv = -b / (2 * a)
        return float(xv), float(np.polyval([a, b, c], xv))

    if 0 < j < len(x) - 1:
        xv, yv = vertice(j - 1)
        return {"x": xv, "y": yv, "interior": True, "x_vertice_diag": NAN}
    xd = vertice(len(x) - 3)[0] if j == len(x) - 1 else NAN
    return {"x": float(x[j]), "y": float(y[j]), "interior": False, "x_vertice_diag": xd}


def indice_forma(x, y, painel):
    """Índice de forma A = ∫₀¹ (y − corda) dx e classificação pelas segundas diferenças divididas.

    Args:
        x (np.ndarray): Abscissas normalizadas crescentes (começam em ~0).
        y (np.ndarray): Ordenadas normalizadas.
        painel (str): ``"TR"`` (corda 1 − x) ou ``"eta"`` (corda x).

    Returns:
        tuple[float, str, str]: (A, forma pelo sinal de A, forma pelo sinal das segundas diferenças em x ≤ 1).
    """
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    corda = (1 - x) if painel == "TR" else x
    m = x <= 1.0 + 1e-9
    xs, ys = list(x[m]), list((y - corda)[m])
    if xs[-1] < 1.0 and (~m).any():                      # fecha a integral em x = 1 por interpolação linear
        k = int(np.argmax(~m))
        xs.append(1.0)
        ys.append(float(np.interp(1.0, [x[k - 1], x[k]], [(y - corda)[k - 1], (y - corda)[k]])))
    A = float(np.trapezoid(ys, xs)) if hasattr(np, "trapezoid") else float(np.trapz(ys, xs))
    f_a = "quase linear" if abs(A) < TOL_FORMA else ("côncava" if A > 0 else "convexa")
    xm, ym = x[m], y[m]
    if len(xm) < 3:
        return A, f_a, "indefinida"
    d1 = np.diff(ym) / np.diff(xm)
    d2 = np.diff(d1) / (0.5 * (xm[2:] - xm[:-2]))
    neg, pos = int((d2 < 0).sum()), int((d2 > 0).sum())
    f_d = "côncava" if pos == 0 else ("convexa" if neg == 0 else f"mista ({neg}−/{pos}+)")
    return A, f_a, f_d


def virgula(v, nd=3):
    """Formata um número com vírgula decimal (textos da figura e do console).

    Args:
        v (float): Valor.
        nd (int): Casas decimais.

    Returns:
        str: Texto, ou "—" se não finito.
    """
    return f"{v:.{nd}f}".replace(".", ",") if np.isfinite(v) else "—"


# ---------------------------------------------------------------------------------------------------------------
# Dados
# ---------------------------------------------------------------------------------------------------------------
def ler_ct():
    """Lê o desempenho do CT_V1 e o SEM por lotes de TR de cada ponto.

    Returns:
        dict: ``SR``, ``TR``, ``eta``, ``sem_TR`` (np.ndarray, ordenados por SR) e ``metodo`` (lista de str).

    Raises:
        FileNotFoundError: Se faltar ``desempenho.csv`` ou o ``resultado.json`` de algum ponto.
        RuntimeError: Se o TR do ``resultado.json`` não conferir com o de ``desempenho.csv``.
    """
    with open(DES, encoding="utf-8") as f:
        linhas = sorted(csv.DictReader(f, delimiter=";"), key=lambda l: float(l["SR"]))
    sr = np.array([float(l["SR"]) for l in linhas])
    tr = np.array([float(l["TR"]) for l in linhas])
    sem, metodo = [], []
    for s, t in zip(sr, tr):
        arq = gci.pasta_sr(DIR, s) / "dados" / "resultado.json"
        r = json.loads(arq.read_text(encoding="utf-8"))
        serie, _ = gci.serie_janela("G0", r, s)
        q = gci.grandezas_ponto(r, serie)
        if abs(q["TR"]["v"] - t) > 1e-6 * abs(t):
            raise RuntimeError(f"SR={s}: TR do resultado.json ({q['TR']['v']}) ≠ desempenho.csv ({t})")
        sem.append(q["TR"]["sem"])
        metodo.append(q["_meta"]["metodo"])
    return {"SR": sr, "TR": tr, "eta": sr * tr, "sem_TR": np.array(sem), "metodo": metodo}


def normalizar_ct(ct):
    """Normaliza o CT_V1 nos dois painéis e propaga o SEM por Monte Carlo.

    Args:
        ct (dict): Saída de ``ler_ct``.

    Returns:
        dict: ``SRc_lin``, ``SRc_pchip``, ``U_SRc_iter`` (meia largura de 95 % do MC), ``d_interp``
        (|linear − PCHIP|), ``U_SRc`` (quadratura das duas), ``eta_max`` (dict de ``maximo_parabola``),
        ``TR`` e ``eta`` (dicts com ``x``, ``y``, ``xlo``, ``xhi``, ``ylo``, ``yhi`` = limites de 95 %).
    """
    sr, tr, sem = ct["SR"], ct["TR"], ct["sem_TR"]
    srl, _ = cruzamento_linear(sr, tr)
    srp = cruzamento_pchip(sr, tr)
    d_int = abs(srl - srp)
    emax = maximo_parabola(sr, ct["eta"])
    rng = np.random.default_rng(SEMENTE)
    amostras = tr + rng.standard_normal((N_MC, len(sr))) * sem
    src = np.array([cruzamento_linear(sr, a)[0] for a in amostras])
    x_tr = sr[None, :] / src[:, None]
    y_tr = (amostras - 1) / (amostras[:, [0]] - 1)
    eta_a = sr * amostras
    nrm = [maximo_parabola(sr, e) for e in eta_a]
    x_e = sr[None, :] / np.array([n["x"] for n in nrm])[:, None]
    y_e = eta_a / np.array([n["y"] for n in nrm])[:, None]

    def faixa(v):
        """Percentis 2,5 e 97,5 por coluna, ignorando amostras inválidas.

        Args:
            v (np.ndarray): Amostras (N_MC × n_pontos).

        Returns:
            tuple[np.ndarray, np.ndarray]: Limites inferior e superior.
        """
        return np.nanpercentile(v, 2.5, axis=0), np.nanpercentile(v, 97.5, axis=0)

    lo_s, hi_s = np.nanpercentile(src, [2.5, 97.5])
    u_iter = 0.5 * (hi_s - lo_s)
    u_tot = math.hypot(u_iter, d_int)
    xt = sr / srl
    xlo, xhi = faixa(x_tr)
    # barra em x: MC (iterativo) e, em quadratura, a diferença linear × PCHIP (interpolação)
    rel_int = d_int / srl
    xlo = xt - np.hypot(xt - xlo, xt * rel_int)
    xhi = xt + np.hypot(xhi - xt, xt * rel_int)
    ylo, yhi = faixa(y_tr)
    out_tr = {"x": xt, "y": (tr - 1) / (tr[0] - 1), "xlo": xlo, "xhi": xhi, "ylo": ylo, "yhi": yhi}
    elo, ehi = faixa(y_e)
    exlo, exhi = faixa(x_e)
    out_eta = {"x": sr / emax["x"], "y": ct["eta"] / emax["y"], "xlo": exlo, "xhi": exhi, "ylo": elo, "yhi": ehi}
    return {"SRc_lin": srl, "SRc_pchip": srp, "U_SRc_iter": u_iter, "d_interp": d_int, "U_SRc": u_tot,
            "SRc_lo": lo_s, "SRc_hi": hi_s, "eta_max": emax, "TR": out_tr, "eta": out_eta}


def normalizar_literatura():
    """Normaliza cada série da literatura nos dois painéis e diz se entra no envelope, tracejada ou fora.

    Returns:
        list[dict]: Uma entrada por série × painel com ``painel``, ``serie``, ``fonte``, ``origem``, ``situacao``
        (``envelope``, ``tracejada_extrapolada`` ou ``fora``), ``motivo``, ``norm_x`` (SR_acopl ou SR*),
        ``norm_x_pchip`` (só TR), ``norm_y`` (TR_estol ou η_máx), ``SR``, ``x``, ``y`` (np.ndarray).
    """
    out = []
    for s in lit_mod.ler_literatura():
        sr, tr, eta = s["SR"], s["TR"], s["eta"]
        base = {"serie": s["rotulo"], "fonte": s["fonte"], "origem": s["origem"], "SR": sr}
        if s["origem"] == "modelo":
            for p in ("TR", "eta"):
                out.append(dict(base, painel=p, situacao="fora", norm_x=NAN, norm_x_pchip=NAN, norm_y=NAN,
                                x=np.array([]), y=np.array([]),
                                motivo="forma de Kotwicki ajustada por Nadeem (não é medição; sem estol)"))
            continue
        # --- TR
        tr0 = tr[0]
        srl, _ = cruzamento_linear(sr, tr)
        if np.isfinite(srl):
            sit, mot, srp = "envelope", f"cruza TR = 1 (TR_estol em SR = {virgula(sr[0], 2)})", cruzamento_pchip(sr, tr)
        else:
            srl = cruzamento_extrapolado(sr, tr)
            sit, srp = "tracejada_extrapolada", NAN
            mot = f"SR máx = {virgula(sr[-1], 2)} com TR = {virgula(tr[-1], 3)} > 1: SR_acopl extrapolado (2 últimos pontos)"
        out.append(dict(base, painel="TR", situacao=sit, motivo=mot, norm_x=srl, norm_x_pchip=srp, norm_y=tr0,
                        x=sr / srl, y=(tr - 1) / (tr0 - 1)))
        # --- η (no LES de Zhang o CSV não traz η: usa η = SR·TR, a definição)
        e = sr * tr if np.isnan(eta).all() else eta
        mx = maximo_parabola(sr, e)
        if mx["interior"]:
            sit, mot = "envelope", "máximo interior" + (" (η = SR·TR)" if np.isnan(eta).all() else "")
        else:
            sit = "fora"
            mot = (f"máximo no último ponto (SR = {virgula(sr[-1], 2)}): η_máx não delimitado"
                   if s["fonte"] != "nadeem2020" else "η cresce até o fim (passa ao modo acoplamento): sem máximo interior")
        out.append(dict(base, painel="eta", situacao=sit, motivo=mot, norm_x=mx["x"], norm_x_pchip=NAN,
                        norm_y=mx["y"], x=sr / mx["x"], y=e / mx["y"]))
    return out


def envelope(series, n=301):
    """Faixa mín–máx das séries de envelope, no domínio comum a todas (0 até o menor x máximo).

    Args:
        series (list[dict]): Séries normalizadas com ``situacao == "envelope"`` de um painel.
        n (int): Pontos da grade.

    Returns:
        tuple[np.ndarray, np.ndarray, np.ndarray]: (grade x, y mínimo, y máximo), interpolação linear em x.
    """
    xmax = min(np.nanmax(s["x"]) for s in series)
    g = np.linspace(0.0, xmax, n)
    ys = np.array([np.interp(g, s["x"], s["y"]) for s in series])
    return g, ys.min(axis=0), ys.max(axis=0)


def avaliar_ct(norm, painel, env):
    """Compara os pontos do CT_V1 de um painel com o envelope.

    Args:
        norm (dict): ``normalizar_ct(...)[painel]``.
        painel (str): ``"TR"`` ou ``"eta"``.
        env (tuple): Saída de ``envelope``.

    Returns:
        list[dict]: Por ponto: ``env_min``, ``env_max``, ``situacao`` (``fixado pela normalização``,
        ``fora do domínio do envelope``, ``dentro``, ``fora (barra toca)``, ``fora``) e ``desvio`` (distância
        vertical do valor central à borda mais próxima; 0 dentro; nan sem comparação).
    """
    g, lo, hi = env
    res = []
    for i, (x, y) in enumerate(zip(norm["x"], norm["y"])):
        fixo = (i == 0) or (painel == "eta" and abs(x - 1) < 1e-12 and abs(y - 1) < 1e-12)
        if fixo:
            res.append({"env_min": NAN, "env_max": NAN, "situacao": "fixado pela normalização", "desvio": NAN})
            continue
        if x > g[-1]:
            res.append({"env_min": NAN, "env_max": NAN, "situacao": "fora do domínio do envelope", "desvio": NAN})
            continue
        a, b = float(np.interp(x, g, lo)), float(np.interp(x, g, hi))
        if a <= y <= b:
            sit, dv = "dentro", 0.0
        else:
            dv = (a - y) if y < a else (y - b)
            toca = norm["yhi"][i] >= a and norm["ylo"][i] <= b
            sit = "fora (barra toca)" if toca else "fora"
        res.append({"env_min": a, "env_max": b, "situacao": sit, "desvio": float(dv) if y < a or y > b else 0.0})
    return res


# ---------------------------------------------------------------------------------------------------------------
# Figura e tabelas
# ---------------------------------------------------------------------------------------------------------------
def estilo_serie(s):
    """Estilo de linha de uma série da literatura (reaproveita ``ESTILO`` e ``CORES`` do g5).

    Args:
        s (dict): Série normalizada.

    Returns:
        dict: kwargs de ``plot``.
    """
    st = dict(lit_mod.ESTILO.get(s["origem"], {}), color=lit_mod.CORES.get(s["fonte"], "k"), alpha=0.85)
    st["ms"] = 3.0
    if s["situacao"] == "tracejada_extrapolada":
        st.update(ls="--", lw=0.8, alpha=0.55, marker=st.get("marker"), mfc="none")
    return st


def grafico(ct, nct, lit, envs, aval):
    """Desenha os dois painéis do colapso e grava o PNG.

    Args:
        ct (dict): Saída de ``ler_ct``.
        nct (dict): Saída de ``normalizar_ct``.
        lit (list[dict]): Saída de ``normalizar_literatura``.
        envs (dict): ``{"TR": envelope, "eta": envelope}``.
        aval (dict): ``{"TR": avaliar_ct, "eta": avaliar_ct}``.
    """
    fmt = FuncFormatter(lambda v, _: f"{v:g}".replace(".", ","))
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(13.5, 9.0))
    est_ct = dict(color="k", marker="D", ms=5, lw=1.8, capsize=2.5, elinewidth=0.9, zorder=6)

    def desenha_lit(ax, painel):
        """Desenha o envelope e as séries da literatura de um painel.

        Args:
            ax (matplotlib.axes.Axes): Eixo.
            painel (str): ``"TR"`` ou ``"eta"``.
        """
        g, lo, hi = envs[painel]
        ax.fill_between(g, lo, hi, color="#e67e22", alpha=0.30, lw=0, label="envelope mín–máx da literatura", zorder=1)
        for s in lit:
            if s["painel"] != painel or s["situacao"] == "fora":
                continue
            rot = s["serie"] + (" [SR_acopl extrapolado]" if s["situacao"] == "tracejada_extrapolada" else "")
            ax.plot(s["x"], s["y"], label=rot, zorder=3, **estilo_serie(s))
            if s["situacao"] == "tracejada_extrapolada":
                ax.plot([s["x"][-1], 1.0], [s["y"][-1], 0.0], ls=":", lw=0.8, color=estilo_serie(s)["color"], alpha=0.6)
                ax.plot(1.0, 0.0, marker="x", ms=5, color=estilo_serie(s)["color"], alpha=0.7)

    def barras(n):
        """Meias larguras assimétricas (95 %) a partir dos limites do Monte Carlo.

        Args:
            n (dict): ``nct["TR"]`` ou ``nct["eta"]``.

        Returns:
            tuple[np.ndarray, np.ndarray]: ``xerr`` e ``yerr`` (2 × n, não negativos).
        """
        xe = np.vstack([np.clip(n["x"] - n["xlo"], 0, None), np.clip(n["xhi"] - n["x"], 0, None)])
        ye = np.vstack([np.clip(n["y"] - n["ylo"], 0, None), np.clip(n["yhi"] - n["y"], 0, None)])
        return xe, ye

    # ---------------- painel TR
    desenha_lit(a1, "TR")
    n = nct["TR"]
    xe, ye = barras(n)
    a1.errorbar(n["x"], n["y"], xerr=xe, yerr=ye, label="CT_V1 MRF (este trabalho), barras de 95 %", **est_ct)
    for xi, yi, s in zip(n["x"], n["y"], ct["SR"]):
        if xi <= 1.33 and yi >= -0.62:
            a1.annotate(f"SR {virgula(s, 2)}", (xi, yi), xytext=(5, 4), textcoords="offset points", fontsize=7)
    fora = [(xi, yi, s) for xi, yi, s in zip(n["x"], n["y"], ct["SR"]) if xi > 1.33 or yi < -0.62]
    if fora:
        txt = "; ".join(f"SR {virgula(s, 2)}: ({virgula(xi, 2)}; {virgula(yi, 2)})" for xi, yi, s in fora)
        txt = txt.replace("; SR", "\nSR")
        a1.annotate("CT_V1 fora do quadro (ver detalhe):\n" + txt, (1.31, 1.09), ha="right", va="top",
                    fontsize=6.5, color="k")
    a1.axhline(0, color="gray", lw=0.7, ls=":")
    a1.axvline(1, color="gray", lw=0.7, ls=":")
    a1.set(xlim=(-0.03, 1.33), ylim=(-0.62, 1.12), xlabel="SR / SR_acopl",
           ylabel="(TR − 1) / (TR_estol − 1)",
           title=(f"(a) Razão de torque colapsada — CT_V1: TR_estol = {virgula(ct['TR'][0], 4)}, "
                  f"SR_acopl = {virgula(nct['SRc_lin'], 3)} ± {virgula(nct['U_SRc'], 3)}"))
    a1.title.set_fontsize(9)
    ins = a1.inset_axes([0.07, 0.07, 0.40, 0.33])
    g, lo, hi = envs["TR"]
    ins.fill_between(g, lo, hi, color="#e67e22", alpha=0.35, lw=0)
    ins.errorbar(n["x"], n["y"], xerr=xe, yerr=ye, color="k", marker="D", ms=3, lw=1.2, capsize=1.5, elinewidth=0.7)
    ins.axhline(0, color="gray", lw=0.6, ls=":")
    ins.axvline(1, color="gray", lw=0.6, ls=":")
    ins.set_title("detalhe: CT_V1 inteiro", fontsize=7)
    ins.tick_params(labelsize=6)
    ins.xaxis.set_major_formatter(fmt)
    ins.yaxis.set_major_formatter(fmt)

    # ---------------- painel η
    desenha_lit(a2, "eta")
    n = nct["eta"]
    xe, ye = barras(n)
    a2.errorbar(n["x"], n["y"], xerr=xe, yerr=ye, mfc="white",
                label="CT_V1 MRF: η_máx não atingido; normalizado por SR = 0,9 (x e y são limites superiores)",
                **est_ct)
    a2.annotate("CT_V1: η ainda cresce em SR = 0,9\n→ SR* ≥ 0,9 e η_máx ≥ " + virgula(ct["eta"][-1], 3)
                + "\n(pontos verdadeiros ficam à esquerda e abaixo)", (1.0, 1.0), xytext=(0.52, 0.30),
                textcoords="data", fontsize=7, arrowprops=dict(arrowstyle="->", lw=0.7))
    a2.plot([0, 1], [0, 1], color="gray", lw=0.7, ls=":", label="corda (forma linear)")
    a2.set(xlim=(-0.02, 1.2), ylim=(-0.02, 1.1), xlabel="SR / SR(η_máx)", ylabel="η / η_máx",
           title="(b) Eficiência colapsada (SR* e η_máx: parábola nos 3 pontos em torno do máximo)")
    a2.title.set_fontsize(9)
    for a in (a1, a2):
        a.grid(alpha=0.3)
        a.xaxis.set_major_formatter(fmt)
        a.yaxis.set_major_formatter(fmt)
    for a in (a1, a2):
        a.legend(fontsize=6.5, loc="upper center", bbox_to_anchor=(0.5, -0.09), ncol=2, frameon=False)

    fora_eta = sorted({s["serie"] for s in lit if s["painel"] == "eta" and s["situacao"] == "fora"
                       and s["origem"] != "modelo"})
    na = {p: sum(r["situacao"] == "dentro" for r in aval[p]) for p in aval}
    nc = {p: sum(r["situacao"] not in ("fixado pela normalização", "fora do domínio do envelope") for r in aval[p])
          for p in aval}
    rod = ("Validação de TENDÊNCIA (geometrias diferentes; só a forma é comparada). A validação forte está nos blocos A–H "
           "e nas verificações do conversor (balanço de torque, Euler, y+, sensibilidade, GCI).\n"
           f"(a) Envelope: Zhang 2023 (exp, SST, LES) e Nadeem 2020. Guo 2021 e Ran 2021 não chegam a TR = 1 (SR ≤ 0,8): "
           "tracejadas com SR_acopl extrapolado (x no fim da reta pontilhada), fora do envelope. "
           f"CT_V1 no envelope: {na['TR']} de {nc['TR']} pontos comparáveis.\n"
           f"(b) Envelope: séries com máximo interior. Fora: {', '.join(fora_eta)} (sem máximo interior). "
           f"CT_V1 no envelope: {na['eta']} de {nc['eta']} (com η_máx provisório). Kotwicki (ajuste, sem estol) fora dos dois.\n"
           f"Barras do CT_V1: 95 % do Monte Carlo do SEM por lotes (gci_conversor); em x soma a diferença linear × PCHIP do "
           f"SR_acopl ({virgula(nct['SRc_lin'], 4)} × {virgula(nct['SRc_pchip'], 4)}). "
           "Estol de Zhang em SR = 0,01.")
    fig.text(0.01, 0.005, rod, fontsize=7, va="bottom", ha="left", wrap=True)
    fig.suptitle("Colapso adimensional das curvas do conversor: CT_V1 (MRF, k-ω SST) × literatura", fontsize=11)
    fig.tight_layout(rect=(0, 0.075, 1, 0.97))
    SAIDA_PNG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(SAIDA_PNG, dpi=150)
    plt.close(fig)


def gravar_tabelas(ct, nct, lit, aval):
    """Grava a tabela por ponto e o resumo por série (separador ;, ponto decimal como os demais CSV do repositório).

    Args:
        ct (dict): Saída de ``ler_ct``.
        nct (dict): Saída de ``normalizar_ct``.
        lit (list[dict]): Saída de ``normalizar_literatura``.
        aval (dict): ``{"TR": avaliar_ct, "eta": avaliar_ct}``.

    Returns:
        list[dict]: Linhas do resumo (para o console).
    """
    col = ["painel", "serie", "fonte", "origem", "situacao_serie", "SR", "x", "y", "x_lo95", "x_hi95",
           "y_lo95", "y_hi95", "env_min", "env_max", "situacao_ponto", "desvio_fora"]
    linhas = []
    for p in ("TR", "eta"):
        n = nct[p]
        for i, sr in enumerate(ct["SR"]):
            a = aval[p][i]
            linhas.append({"painel": p, "serie": "CT_V1 MRF", "fonte": "este_trabalho", "origem": "cfd_rans_mrf",
                           "situacao_serie": "simulacao" + (" (eta_max provisorio: SR=0,9)" if p == "eta" else ""),
                           "SR": sr, "x": n["x"][i], "y": n["y"][i], "x_lo95": n["xlo"][i], "x_hi95": n["xhi"][i],
                           "y_lo95": n["ylo"][i], "y_hi95": n["yhi"][i], "env_min": a["env_min"],
                           "env_max": a["env_max"], "situacao_ponto": a["situacao"], "desvio_fora": a["desvio"]})
        for s in lit:
            if s["painel"] != p or s["situacao"] == "fora":
                continue
            for sr, x, y in zip(s["SR"], s["x"], s["y"]):
                linhas.append({"painel": p, "serie": s["serie"], "fonte": s["fonte"], "origem": s["origem"],
                               "situacao_serie": s["situacao"], "SR": sr, "x": x, "y": y})
    SAIDA_CSV.parent.mkdir(parents=True, exist_ok=True)
    with open(SAIDA_CSV, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=col, delimiter=";", restval="")
        w.writeheader()
        for l in linhas:
            w.writerow({k: (f"{v:.6g}" if isinstance(v, (float, np.floating)) and np.isfinite(v) else
                            ("" if isinstance(v, (float, np.floating)) else v)) for k, v in l.items()})

    res = []
    em = nct["eta_max"]
    for p in ("TR", "eta"):
        n = nct[p]
        A, fa, fd = indice_forma(n["x"], n["y"], p)
        comp = [a for a in aval[p] if a["situacao"] not in ("fixado pela normalização", "fora do domínio do envelope")]
        dev = [a["desvio"] for a in comp if a["desvio"] > 0]
        res.append({"painel": p, "serie": "CT_V1 MRF", "origem": "cfd_rans_mrf", "situacao": "simulacao",
                    "motivo": ("SR_acopl linear; PCHIP como verificação" if p == "TR" else
                               "eta cresce até SR=0,9: normalizado pelo último ponto (limite)"),
                    "norm_x": nct["SRc_lin"] if p == "TR" else em["x"],
                    "norm_x_pchip": nct["SRc_pchip"] if p == "TR" else NAN,
                    "U95_norm_x": nct["U_SRc"] if p == "TR" else NAN,
                    "norm_y": ct["TR"][0] if p == "TR" else em["y"],
                    "indice_forma_A": A, "forma_A": fa, "forma_d2": fd,
                    "pontos_comparaveis": len(comp), "pontos_dentro": sum(a["situacao"] == "dentro" for a in comp),
                    "pontos_barra_toca": sum(a["situacao"] == "fora (barra toca)" for a in comp),
                    "desvio_max_fora": max(dev) if dev else 0.0})
        for s in lit:
            if s["painel"] != p:
                continue
            A, fa, fd = (indice_forma(s["x"], s["y"], p) if len(s["x"]) else (NAN, "", ""))
            res.append({"painel": p, "serie": s["serie"], "origem": s["origem"], "situacao": s["situacao"],
                        "motivo": s["motivo"], "norm_x": s["norm_x"], "norm_x_pchip": s["norm_x_pchip"],
                        "U95_norm_x": NAN, "norm_y": s["norm_y"], "indice_forma_A": A, "forma_A": fa, "forma_d2": fd})
    colr = ["painel", "serie", "origem", "situacao", "motivo", "norm_x", "norm_x_pchip", "U95_norm_x", "norm_y",
            "indice_forma_A", "forma_A", "forma_d2", "pontos_comparaveis", "pontos_dentro", "pontos_barra_toca",
            "desvio_max_fora"]
    with open(SAIDA_RES, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=colr, delimiter=";", restval="")
        w.writeheader()
        for l in res:
            w.writerow({k: (f"{v:.6g}" if isinstance(v, (float, np.floating)) and np.isfinite(v) else
                            ("" if isinstance(v, (float, np.floating)) else v)) for k, v in l.items()})
    return res


def main():
    """Calcula o colapso, grava figura e tabelas e imprime o resumo no console."""
    ct = ler_ct()
    nct = normalizar_ct(ct)
    lit = normalizar_literatura()
    envs = {p: envelope([s for s in lit if s["painel"] == p and s["situacao"] == "envelope"]) for p in ("TR", "eta")}
    aval = {p: avaliar_ct(nct[p], p, envs[p]) for p in ("TR", "eta")}
    grafico(ct, nct, lit, envs, aval)
    res = gravar_tabelas(ct, nct, lit, aval)

    em = nct["eta_max"]
    print(f"SEM por ponto: {', '.join(sorted(set(ct['metodo'])))}")
    print(f"SR_acopl CT_V1: linear = {virgula(nct['SRc_lin'], 4)}  PCHIP = {virgula(nct['SRc_pchip'], 4)}  "
          f"|Δ| = {virgula(nct['d_interp'], 4)}")
    print(f"  U95 iterativo (MC) = ±{virgula(nct['U_SRc_iter'], 4)} [{virgula(nct['SRc_lo'], 4)}; "
          f"{virgula(nct['SRc_hi'], 4)}]   U95 total (com interpolação) = ±{virgula(nct['U_SRc'], 4)}")
    print(f"TR_estol CT_V1 = {virgula(ct['TR'][0], 4)} ± {virgula(2 * ct['sem_TR'][0], 4)} (2·SEM)")
    print(f"η_máx CT_V1: máximo interior = {em['interior']}; usado SR* = {virgula(em['x'], 3)}, "
          f"η = {virgula(em['y'], 4)}; vértice da parábola nos 3 últimos pontos (diagnóstico) SR ≈ "
          f"{virgula(em['x_vertice_diag'], 2)}")
    for p in ("TR", "eta"):
        print(f"\nPainel {p}: CT_V1 x envelope")
        for sr, a, x, y in zip(ct["SR"], aval[p], nct[p]["x"], nct[p]["y"]):
            print(f"  SR {virgula(sr, 2):>5}  x = {virgula(x, 3):>6}  y = {virgula(y, 3):>7}  "
                  f"env [{virgula(a['env_min'], 3)}; {virgula(a['env_max'], 3)}]  {a['situacao']}"
                  + (f"  desvio {virgula(a['desvio'], 3)}" if a["desvio"] and np.isfinite(a["desvio"]) else ""))
    print("\nResumo por série:")
    for r in res:
        print(f"  [{r['painel']:>3}] {r['serie']:<42} {r['situacao']:<22} norm_x = {virgula(r['norm_x'], 4):>6}"
              f"  norm_y = {virgula(r['norm_y'], 4):>6}  A = {virgula(r['indice_forma_A'], 3):>6} {r['forma_A']:<12}"
              f" d2: {r['forma_d2']}" + (f"   ({r['motivo']})" if r["situacao"] != "envelope" else ""))
    print("\ngravado:", SAIDA_PNG)
    print("gravado:", SAIDA_CSV)
    print("gravado:", SAIDA_RES)


if __name__ == "__main__":
    main()
