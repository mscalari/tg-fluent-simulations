"""GCI (Grid Convergence Index) do conversor CT_V1 (MRF) pelo procedimento de Celik et al. (2008).

Celik, Ghia, Roache, Freitas, Coleman, Raad. "Procedure for Estimation and Reporting of Uncertainty Due to
Discretization in CFD Applications". J. Fluids Eng. 130(7), 078001, 2008.

Uso (rodar com ``PYTHONIOENCODING=utf-8`` e o ``.venv241``)::

    python tools/gci_conversor.py                       calcula o GCI com G0, G1 e G2 (diz o que falta, se faltar)
    python tools/gci_conversor.py --n0 1900736 --v0 0.08576847
                                                        força N e V da G0 (padrão: transcript do solver da G0)
    python tools/gci_conversor.py --p-formal 2          ordem formal usada como reserva e no GCI conservador
    python tools/gci_conversor.py --mc 0                desliga o Monte Carlo do SEM (padrão: 300 amostras)
    python tools/gci_conversor.py --teste               valida o código com dados sintéticos e o exemplo do artigo

Entradas (nenhuma é gerada aqui; o Fluent NÃO é executado)::

    G0 (fina)   results/conversor/CT_V1_MRF/SR000|SR040/dados/resultado.json
    G1 (média)  results/conversor/CT_V1_GCI_G1/SR000|SR040/dados/resultado.json
    G2 (grossa) results/conversor/CT_V1_GCI_G2/SR000|SR040/dados/resultado.json
    N e V       preferência: transcript do solver de cada nível (``<pasta do nível>/fluent/*.trn``: linhas
                "<N> cells, 3 cell zones" e "total volume (m3): <V>", na leitura da malha do nível; para a G0 é
                ``CT_V1_MRF/fluent/fluent-20261120-230926-10972.trn``, linhas 56 e 128). Alternativas: ``--n0/--v0``
                (G0), ``CT_V1_GCI/logs/<nivel>_parametros.json`` (chaves ``estatisticas.celulas`` e
                ``estatisticas.volume_m3``, de ``malha_gci_conversor.py``) e, para a G0, as constantes ``G0_PADRAO``
                (as mesmas do transcript).
    Séries      ``<pasta do nível>/dados/*.out`` (Report File do Fluent; ``torque`` t_bomba, t_turbina, t_estator),
                para o ruído iterativo. O arquivo é escolhido pelas médias da janela (conferem com ``media_janela``).

Saídas::

    results/conversor/CT_V1_GCI/dados/gci.csv         (separador ``;``, uma linha por SR x grandeza)
    results/conversor/CT_V1_GCI/graficos/gci_tr.png   (φ x h com φ_ext e barra de GCI na malha fina)

Método (Celik 2008, três malhas, 1 = fina, 3 = grossa)::

    h_i = (V_i / N_i)^(1/3)         r21 = h2/h1        r32 = h3/h2
    ε32 = φ3 - φ2                   ε21 = φ2 - φ1      s = sinal(ε32/ε21)
    p   : raiz de  p·ln r21 = ln|ε32/ε21| + q(p),   q(p) = ln((r21^p - s)/(r32^p - s))   (eqs. 3a-3c)
    φ_ext21 = (r21^p φ1 - φ2)/(r21^p - 1)
    e_a21 = |(φ1-φ2)/φ1|    e_ext21 = |(φ_ext21-φ1)/φ_ext21|    GCI_fina21 = 1,25 e_a21/(r21^p - 1)
    GCI32 análogo, e razão assintótica GCI32/(r21^p GCI21) (≈ |φ1/φ2| para a solução exata; ≈ 1 se φ varia pouco).

Escolhas além do artigo (conferir)::

    * p é a MENOR raiz positiva de (3a-3c) em [0,05; 20] por varredura + bisseção (o ponto fixo do artigo falha
      quando r32 ≫ r21). Não há raiz com s > 0  ->  "divergente". s < 0 -> "oscilatória" (|R| < 1) ou
      "oscilatória divergente" (|R| > 1). R = ε21/ε32 é só informativo: o critério 0<R<1 (monotônica) vale
      apenas com r21 = r32; com r21 ≠ r32 a existência da raiz é que decide.
    * p fora de [0,5; 10], sem raiz ou ε = 0: usa p formal (``--p-formal``, padrão 2) e sinaliza em ``p_origem``.
      Divergentes não têm φ_ext.
    * Coluna conservadora (``GCI21_cons_pct``): Fs = 1,25 só se |p − p_f|/p_f < 10 %; senão Fs = 3 e p limitado
      a min(max(p; 0,5); p_f). O GCI de Celik (Fs = 1,25 com o p calculado) é mantido em ``GCI21_pct``.
    * Ruído iterativo: o erro-padrão da média (SEM) vem da série por iteração da janela (últimas ``janela``
      iterações do ``resultado.json``): TR, K e η são calculados iteração a iteração e o SEM é por lotes (5 lotes
      de 40 iterações); o tempo de autocorrelação integral τ é calculado como conferência (alerta se o SEM de
      τ for muito maior). A diferença é "dentro do ruído iterativo" se |ε| < U_iter = 2·√(SEM_i² + SEM_j²).
      ``Uiter21_sobre_eps21`` (= U_iter/|ε21|; ≥ 1 é ruído) e ``Uiter_sobre_Udisc`` (= 2·SEM1/(GCI21·|φ1|),
      incerteza iterativa da malha fina contra a discretização) quantificam. Sem o .out cobrindo a janela, cai
      no método antigo (σ por iteração como SEM, muito conservador) com aviso.
    * Monte Carlo opcional (``--mc``): φ_i ~ N(φ_i, SEM_i) independentes, 300 amostras com semente fixa; reporta
      os percentis 16 e 84 de p e o 84 do GCI21 e a fração de amostras monotônicas.
    * Sanidade da malha: marca INVÁLIDO se |V_i/V_0 − 1| > 2 % ou se r_ij se afastar mais de 15 % do nominal
      1,35 (``FATORES`` de ``malha_gci_conversor.py``); nesse caso o script sai com código 3.
    * T_e é pequeno e troca de sinal perto do acoplamento: ``e_a21`` relativo a ele explode, então também se
      reporta o GCI como fração de |T_b| (``GCI21_rel_Tb_pct``).
"""
import argparse
import csv
import json
import math
import re
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from analise_conversor import CM_PARA_NM, fator_torque  # noqa: E402,F401  (mesma conversão Cm -> N·m da análise)

RAIZ = Path(__file__).resolve().parents[1]
CONV = RAIZ / "results" / "conversor"
DIR_GCI = CONV / "CT_V1_GCI"
LOGS = DIR_GCI / "logs"
PASTAS = {"G0": CONV / "CT_V1_MRF", "G1": CONV / "CT_V1_GCI_G1", "G2": CONV / "CT_V1_GCI_G2"}
NIVEIS = ("G0", "G1", "G2")
SRS = (0.0, 0.4)
# G0: valores do transcript do solver (CT_V1_MRF/fluent/fluent-20261120-230926-10972.trn, linhas 56 e 128), que
# também aparecem em tools/malha_gci_conversor.py (CELULAS_G0, VOLUME_G0_M3). Só usados se o transcript faltar.
G0_PADRAO = {"celulas": 1_900_736, "volume_m3": 8.576847e-2,
             "fonte": "constantes G0_PADRAO (iguais ao transcript fluent-20261120-230926-10972.trn)"}
R_NOMINAL = 1.35                   # FATORES de malha_gci_conversor.py: G1 = 1,35; G2 = 1,35² (r entre níveis = 1,35)
TOL_R, TOL_V = 0.15, 0.02          # desvios máximos de r (relativo ao nominal) e de V (relativo à G0)
GRANDEZAS = ("TR", "T_b", "T_t", "T_e", "K", "eta")
ROTULOS = {"TR": "TR = −T_t/T_b", "T_b": "T_b [N·m]", "T_t": "T_t [N·m]", "T_e": "T_e [N·m]",
           "K": "K [rpm/√(N·m)]", "eta": "η = SR·TR"}
FS = 1.25                          # fator de segurança de Celik para 3 malhas
FS_CONS = 3.0                      # fator de segurança conservador (p fora de ±10 % do formal)
P_MIN, P_MAX = 0.5, 10.0           # faixa aceita para a ordem aparente
P_BUSCA = (0.05, 20.0)             # intervalo da busca de raiz
N_LOTES = 5                        # lotes do batch means (5 x 40 iterações com janela 200)
COLUNAS = ["SR", "grandeza", "phi1", "phi2", "phi3", "N1", "N2", "N3", "h1_mm", "h2_mm", "h3_mm", "r21", "r32",
           "p", "p_origem", "R", "tipo", "phi_ext21", "e_a21_pct", "e_ext21_pct", "GCI21_pct", "GCI32_pct",
           "razao_assintotica", "GCI21_rel_Tb_pct", "Fs_cons", "p_cons", "GCI21_cons_pct",
           "sem1", "sem2", "sem3", "tau1", "tau2", "tau3", "metodo_sem", "Uiter21_sobre_eps21",
           "Uiter32_sobre_eps32", "Uiter_sobre_Udisc", "dentro_do_ruido", "p_mc_16", "p_mc_84",
           "GCI21_mc_84_pct", "frac_monotonica_mc", "malha_valida", "conclusivo", "observacoes"]

plt.rcParams.update({"figure.dpi": 130, "axes.grid": True, "grid.alpha": 0.3, "font.size": 8})
NAN = float("nan")


# ---------------------------------------------------------------------------------------------------------------
# Núcleo numérico (Celik 2008)
# ---------------------------------------------------------------------------------------------------------------
def residuo_p(p, e21, e32, r21, r32):
    """Resíduo g(p) = p·ln r21 − ln|ε32/ε21| − q(p) das eqs. 3a-3c de Celik (g = 0 define p).

    Args:
        p (float): Ordem aparente tentativa.
        e21 (float): ε21 = φ2 − φ1.
        e32 (float): ε32 = φ3 − φ2.
        r21 (float): h2/h1.
        r32 (float): h3/h2.

    Returns:
        float: g(p); ``nan`` se o argumento do logaritmo de q(p) não for positivo.
    """
    s = 1.0 if e32 / e21 > 0 else -1.0
    num, den = r21 ** p - s, r32 ** p - s
    if num <= 0 or den <= 0:
        return NAN
    return p * math.log(r21) - math.log(abs(e32 / e21)) - math.log(num / den)


def raiz_p(e21, e32, r21, r32):
    """Menor raiz positiva de ``residuo_p`` em ``P_BUSCA`` (varredura em escala log + bisseção).

    Substitui o ponto fixo do artigo, que não converge, por exemplo, com r32 muito maior que r21.

    Args:
        e21 (float): ε21 = φ2 − φ1 (≠ 0).
        e32 (float): ε32 = φ3 − φ2 (≠ 0).
        r21 (float): h2/h1 (> 1).
        r32 (float): h3/h2 (> 1).

    Returns:
        float: A menor raiz, ou ``nan`` se não houver mudança de sinal no intervalo.
    """
    ps = np.geomspace(P_BUSCA[0], P_BUSCA[1], 400)
    g = np.array([residuo_p(p, e21, e32, r21, r32) for p in ps])
    for i in range(len(ps) - 1):
        a, b, ga, gb = ps[i], ps[i + 1], g[i], g[i + 1]
        if not (np.isfinite(ga) and np.isfinite(gb)):
            continue
        if ga == 0:
            return float(a)
        if ga * gb < 0:
            for _ in range(80):
                m = 0.5 * (a + b)
                gm = residuo_p(m, e21, e32, r21, r32)
                if not np.isfinite(gm):
                    break
                if ga * gm <= 0:
                    b = m
                else:
                    a, ga = m, gm
            return float(0.5 * (a + b))
    return NAN


def tipo_convergencia(e21, e32, p_raiz):
    """Classifica a convergência.

    Com s = sinal(ε32/ε21) < 0 a solução alterna: "oscilatória" se |R| < 1 (R = ε21/ε32) e "oscilatória
    divergente" se |R| > 1. Com s > 0 é "monotônica" se existir raiz de (3a-3c) e "divergente" se não existir.
    O critério 0 < R < 1 do artigo só vale com r21 = r32; aqui R é apenas informativo.

    Args:
        e21 (float): ε21 = φ2 − φ1.
        e32 (float): ε32 = φ3 − φ2.
        p_raiz (float): Raiz de ``raiz_p`` (``nan`` se não houver).

    Returns:
        tuple[float, str]: ``(R, tipo)``; tipo ``"indefinido"`` se ε = 0 ou valores não finitos.
    """
    if e21 == 0 or e32 == 0 or not np.isfinite(e21) or not np.isfinite(e32):
        return NAN, "indefinido"
    R = e21 / e32
    if R < 0:
        return R, "oscilatória" if abs(R) < 1 else "oscilatória divergente"
    return R, "monotônica" if np.isfinite(p_raiz) else "divergente"


def ordem_aparente(p_raiz, tipo, p_formal=2.0):
    """Decide o p a usar: a raiz, se válida, ou a ordem formal com o motivo.

    A raiz mais de 10 % acima da ordem formal é limitada a ela (conservador): com ε pequeno ou no ruído
    iterativo, um p aparente grande dá um GCI otimista sem respaldo no esquema de 2ª ordem.

    Args:
        p_raiz (float): Raiz de ``raiz_p`` (``nan`` se não houver).
        tipo (str): Saída de ``tipo_convergencia``.
        p_formal (float): Ordem formal de reserva.

    Returns:
        tuple[float, str]: ``(p, origem)``; origem ``"aparente"`` ou o motivo do p formal.
    """
    if tipo in ("indefinido", "divergente", "oscilatória divergente"):
        return p_formal, tipo
    if not np.isfinite(p_raiz):
        return p_formal, "sem raiz"
    if not P_MIN <= p_raiz <= P_MAX:
        return p_formal, "p fora de [0,5;10] (aparente %.3g)" % p_raiz
    if p_raiz > 1.1 * p_formal:                 # mesma tolerância de 10 % do Fs conservador
        return p_formal, "p aparente %.3g > formal: limitado ao formal" % p_raiz
    return p_raiz, "aparente"


def gci_conservador(phi1, phi2, r21, p_ap, p_formal=2.0):
    """GCI conservador: Fs = 1,25 só se p aparente estiver a menos de 10 % do formal, senão Fs = 3.

    Args:
        phi1 (float): Valor na malha fina.
        phi2 (float): Valor na malha média.
        r21 (float): h2/h1.
        p_ap (float): Ordem aparente bruta (``nan`` se inexistente ou não confiável).
        p_formal (float): Ordem formal.

    Returns:
        tuple[float, float, float]: ``(Fs, p_usado, GCI21)`` (GCI como fração; ``nan`` se φ1 = 0).
    """
    if np.isfinite(p_ap) and abs(p_ap - p_formal) / p_formal < 0.1:
        fs, p = FS, p_ap
    else:
        fs = FS_CONS
        p = min(max(p_ap, P_MIN), p_formal) if np.isfinite(p_ap) else p_formal
    if phi1 == 0 or r21 <= 1:
        return fs, p, NAN
    return fs, p, fs * abs((phi1 - phi2) / phi1) / (r21 ** p - 1)


def celik(phi1, phi2, phi3, r21, r32, p_formal=2.0):
    """Aplica o procedimento de Celik et al. (2008) a uma grandeza em três malhas.

    Args:
        phi1 (float): Valor na malha fina.
        phi2 (float): Valor na malha média.
        phi3 (float): Valor na malha grossa.
        r21 (float): h2/h1 (> 1).
        r32 (float): h3/h2 (> 1).
        p_formal (float): Ordem formal de reserva.

    Returns:
        dict: ``eps21``, ``eps32``, ``R``, ``tipo``, ``p``, ``p_origem``, ``p_raiz`` (raiz bruta), ``phi_ext21``,
        ``e_a21``, ``e_ext21``, ``GCI21``, ``GCI32``, ``razao_assintotica`` (frações; ``nan`` se indefinido, p. ex.
        com ε21 = 0 ou ε32 = 0) e ``Fs_cons``, ``p_cons``, ``GCI21_cons``.
    """
    e21, e32 = phi2 - phi1, phi3 - phi2
    out = {"eps21": e21, "eps32": e32, "R": NAN, "tipo": "indefinido", "p": p_formal, "p_origem": "indefinido",
           "p_raiz": NAN, "phi_ext21": NAN, "e_a21": NAN, "e_ext21": NAN, "GCI21": NAN, "GCI32": NAN,
           "razao_assintotica": NAN, "Fs_cons": NAN, "p_cons": NAN, "GCI21_cons": NAN}
    if e21 == 0 or e32 == 0 or not np.isfinite(e21) or not np.isfinite(e32):
        return out
    p_raiz = raiz_p(e21, e32, r21, r32)
    R, tipo = tipo_convergencia(e21, e32, p_raiz)
    p, origem = ordem_aparente(p_raiz, tipo, p_formal)
    out.update({"R": R, "tipo": tipo, "p": p, "p_origem": origem, "p_raiz": p_raiz})
    p_ap = p_raiz if tipo in ("monotônica", "oscilatória") else NAN
    out["Fs_cons"], out["p_cons"], out["GCI21_cons"] = gci_conservador(phi1, phi2, r21, p_ap, p_formal)
    if phi1 == 0 or phi2 == 0:
        return out
    a21 = r21 ** p
    out["e_a21"] = abs((phi1 - phi2) / phi1)
    e_a32 = abs((phi2 - phi3) / phi2)
    out["GCI21"] = FS * out["e_a21"] / (a21 - 1)
    out["GCI32"] = FS * e_a32 / (r32 ** p - 1)
    if out["GCI21"] > 0:
        out["razao_assintotica"] = out["GCI32"] / (a21 * out["GCI21"])
    if tipo in ("monotônica", "oscilatória"):         # extrapolar dados divergentes não faz sentido
        ext = (a21 * phi1 - phi2) / (a21 - 1)
        out["phi_ext21"] = ext
        out["e_ext21"] = abs((ext - phi1) / ext) if ext != 0 else NAN
    return out


# ---------------------------------------------------------------------------------------------------------------
# Ruído iterativo: SEM por lotes, τ integral e Monte Carlo
# ---------------------------------------------------------------------------------------------------------------
def sem_serie(x):
    """Erro-padrão da média de uma série autocorrelacionada, por lotes, com o τ integral como conferência.

    Batch means: descarta o resto, divide em ``N_LOTES`` lotes de mesmo tamanho (40 iterações com janela 200) e
    usa o desvio-padrão das médias dos lotes dividido por √N_LOTES. O τ integral (janelamento automático de
    Sokal, M ≥ 5τ) dá SEM_τ = σ·√(2τ/n).

    Args:
        x (np.ndarray): Série por iteração na janela.

    Returns:
        tuple[float, float, float]: ``(sem_lotes, tau, sem_tau)`` (τ em iterações, ≥ 0,5).
    """
    x = np.asarray(x, float)
    n = len(x)
    L = n // N_LOTES
    lotes = x[n - N_LOTES * L:].reshape(N_LOTES, L).mean(axis=1)
    sem = float(lotes.std(ddof=1) / math.sqrt(N_LOTES))
    xc = x - x.mean()
    c0 = float(np.dot(xc, xc)) / n
    tau = 0.5
    if c0 > 0:
        for k in range(1, n):
            tau += float(np.dot(xc[:-k], xc[k:])) / n / c0
            if k >= 5 * tau:
                break
        tau = max(tau, 0.5)
    return sem, tau, float(x.std(ddof=1) * math.sqrt(2 * tau / n))


def monte_carlo(vals, sems, r21, r32, p_formal, n_amostras, semente=12345):
    """Propaga o SEM de cada malha (φ_i ~ N(φ_i, SEM_i), independentes) para p, GCI21 e o tipo de convergência.

    Args:
        vals (list[float]): φ1, φ2, φ3.
        sems (list[float]): SEM1, SEM2, SEM3 (mesma unidade de φ).
        r21 (float): h2/h1.
        r32 (float): h3/h2.
        p_formal (float): Ordem formal de reserva.
        n_amostras (int): Número de amostras (0 desliga).
        semente (int): Semente do gerador.

    Returns:
        dict: ``p_mc_16``, ``p_mc_84`` (percentis do p aparente, entre as amostras com p aparente válido),
        ``GCI21_mc_84_pct`` (percentil 84 do GCI21) e ``frac_monotonica_mc``; ``nan`` se desligado ou sem amostras.
    """
    vazio = {"p_mc_16": NAN, "p_mc_84": NAN, "GCI21_mc_84_pct": NAN, "frac_monotonica_mc": NAN}
    if n_amostras <= 0 or not all(np.isfinite(s) for s in sems):
        return vazio
    rng = np.random.default_rng(semente)
    ps, gcis, mono = [], [], 0
    for _ in range(n_amostras):
        a, b, c = (rng.normal(v, s) for v, s in zip(vals, sems))
        r = celik(a, b, c, r21, r32, p_formal)
        mono += r["tipo"] == "monotônica"
        if r["p_origem"] == "aparente":
            ps.append(r["p"])
        if np.isfinite(r["GCI21"]):
            gcis.append(r["GCI21"])
    out = dict(vazio, frac_monotonica_mc=mono / n_amostras)
    if len(ps) >= 20:
        out["p_mc_16"], out["p_mc_84"] = (float(np.percentile(ps, q)) for q in (16, 84))
    if len(gcis) >= 20:
        out["GCI21_mc_84_pct"] = 100 * float(np.percentile(gcis, 84))
    return out


# ---------------------------------------------------------------------------------------------------------------
# Leitura dos dados do repositório
# ---------------------------------------------------------------------------------------------------------------
def pasta_sr(base, sr):
    """Pasta de um ponto de operação dentro da pasta de um nível.

    Args:
        base (Path): Pasta do nível (``PASTAS[nivel]``).
        sr (float): Razão de velocidades (ω_t/ω_b).

    Returns:
        Path: ``base/SRxxx`` (xxx = 100·SR com três dígitos).
    """
    return base / f"SR{int(round(sr * 100)):03d}"


def ler_tabela_out(arq):
    """Lê um Report File do Fluent (mesmo formato de ``conversor_mrf.ler_tabela``, reimplementado aqui).

    Args:
        arq (Path): Arquivo ``.out``.

    Returns:
        tuple[list[str], np.ndarray]: Nomes das colunas do último cabeçalho e matriz só com as linhas numéricas que
        têm o mesmo número de colunas (vazia se não houver).
    """
    nomes, linhas = [], []
    for ln in arq.read_text(errors="replace").splitlines():
        if '"Iteration"' in ln:
            nomes = re.findall(r'"([^"]+)"', ln)
        p = ln.split()
        if p and p[0].isdigit():
            try:
                linhas.append([float(x) for x in p])
            except ValueError:
                pass
    linhas = [l for l in linhas if len(l) == len(nomes)]
    return nomes, (np.array(linhas) if linhas else np.zeros((0, len(nomes))))


_CACHE_OUT = {}


def serie_janela(nivel, r, sr):
    """Acha a série por iteração (em N·m) da janela de um ponto, conferindo com a média do ``resultado.json``.

    Procura em ``<pasta do nível>/dados/*.out`` e ``<pasta do nível>/SRxxx/dados/*.out`` o arquivo cujas
    ``janela`` linhas com it ∈ (it_ini, it_fim] (as últimas) têm médias de t_bomba, t_turbina e t_estator iguais
    (rel. 1e-6) às de ``media_janela``; assim o .out certo é escolhido mesmo com vários na pasta.

    Args:
        nivel (str): ``"G0"``, ``"G1"`` ou ``"G2"``.
        r (dict): ``resultado.json`` do ponto.
        sr (float): Razão de velocidades.

    Returns:
        tuple[dict | None, str]: ``({"t_b", "t_t", "t_e"}: arrays em N·m, ou None)`` e o nome do arquivo ou o motivo
        de não ter achado.
    """
    jan, i0, i1 = r.get("janela"), r.get("it_ini"), r.get("it_fim")
    if not jan or i0 is None or i1 is None or "media_janela" not in r:
        return None, "resultado.json sem janela/it_ini/it_fim"
    c = fator_torque(r)
    chaves = ("t_bomba", "t_turbina", "t_estator")
    cands = sorted(PASTAS[nivel].glob("dados/*.out")) + sorted(pasta_sr(PASTAS[nivel], sr).glob("dados/*.out"))
    for arq in cands:
        if arq not in _CACHE_OUT:
            _CACHE_OUT[arq] = ler_tabela_out(arq)
        nomes, m = _CACHE_OUT[arq]
        if not all(k in nomes for k in chaves) or not len(m):
            continue
        w = m[(m[:, 0] > i0) & (m[:, 0] <= i1)][-jan:]
        if len(w) != jan:
            continue
        cols = [nomes.index(k) for k in chaves]
        if all(abs(w[:, j].mean() - r["media_janela"][k]) <= 1e-6 * max(abs(r["media_janela"][k]), 1e-12)
               for j, k in zip(cols, chaves)):
            return {"t_b": w[:, cols[0]] * c, "t_t": w[:, cols[1]] * c, "t_e": w[:, cols[2]] * c}, arq.name
    return None, f"nenhum .out em {PASTAS[nivel] / 'dados'} cobre a janela ({i0}, {i1}] com as médias do resultado.json"


def grandezas_ponto(r, serie=None):
    """Calcula TR, T_b, T_t, T_e, K e η de um ``resultado.json`` com a incerteza iterativa (SEM) de cada uma.

    Os valores φ vêm de ``media_janela`` (média das últimas ``janela`` iterações), convertidos de Cm para N·m com
    ``fator_torque`` quando ``unidade_torque`` ≠ "N·m". As definições são as de ``analise_conversor.desempenho``:
    TR = −T_t/T_b, η = SR·TR, K = n_b/√|T_b| (n_b em rpm).

    Com ``serie`` (de ``serie_janela``), TR, K e η são calculados iteração a iteração e o SEM sai de ``sem_serie``;
    isso já inclui a correlação entre T_b e T_t (|T_b| e |T_t| são positivamente correlacionados: corr ≈ −0,7 a
    −0,8 com os sinais do Fluent; as flutuações se cancelam em parte na razão).
    Sem série, cai no método antigo, CONSERVADOR e com aviso: SEM := σ por iteração (``desvio_janela``), e para TR
    σ_TR = |TR|·√((σ_b/T_b)² + (σ_t/T_t)²) sem correlação (superestima ~2×), σ_η = SR·σ_TR e σ_K = K·σ_b/(2|T_b|).

    Args:
        r (dict): ``resultado.json`` de um ponto.
        serie (dict | None): Saída de ``serie_janela`` (séries em N·m) ou ``None``.

    Returns:
        dict: ``{grandeza: {"v", "sem", "tau"}}`` para as chaves de ``GRANDEZAS`` (``tau`` é ``nan`` sem série) e
        ``"_meta"`` com ``SR``, ``convergiu``, ``unidade`` e ``metodo`` do SEM.

    Raises:
        KeyError: Se faltar ``media_janela`` ou algum de seus torques.
        ZeroDivisionError: Se T_b ou T_t for zero.
    """
    c = fator_torque(r)
    m, d = r["media_janela"], r.get("desvio_janela", {})
    tb, tt, te = (m[k] * c for k in ("t_bomba", "t_turbina", "t_estator"))
    sr = r.get("sr")
    if sr is None:
        sr = r["omega_t"] / r["omega_b"]
    rpm = r.get("rpm_b") or r["omega_b"] * 30.0 / math.pi
    tr = -tt / tb
    k = rpm / math.sqrt(abs(tb))
    val = {"TR": tr, "T_b": tb, "T_t": tt, "T_e": te, "K": k, "eta": sr * tr}
    out = {g: {"v": val[g], "sem": NAN, "tau": NAN} for g in GRANDEZAS}
    if serie is not None:
        x = {"T_b": serie["t_b"], "T_t": serie["t_t"], "T_e": serie["t_e"]}
        x["TR"] = -serie["t_t"] / serie["t_b"]
        x["K"] = rpm / np.sqrt(np.abs(serie["t_b"]))
        x["eta"] = sr * x["TR"]
        for g in GRANDEZAS:
            out[g]["sem"], out[g]["tau"], _ = sem_serie(x[g])
        metodo = f"lotes {N_LOTES}x{len(x['T_b']) // N_LOTES} it"
    else:
        sb, st, se = (abs(d.get(q, NAN) * c) for q in ("t_bomba", "t_turbina", "t_estator"))
        s_tr = abs(tr) * math.hypot(sb / tb, st / tt)
        for g, s in (("TR", s_tr), ("T_b", sb), ("T_t", st), ("T_e", se), ("K", k * sb / (2 * abs(tb))),
                     ("eta", sr * s_tr)):
            out[g]["sem"] = s
        metodo = "σ por iteração (SEM indisponível)"
    out["_meta"] = {"SR": sr, "convergiu": r.get("convergiu"), "unidade": "N·m" if c == 1.0 else "Cm×0,6125",
                    "metodo": metodo}
    return out


def ler_transcript(nivel):
    """Lê N e V da malha de um nível no transcript do solver (``<pasta do nível>/fluent/*.trn``).

    Considera só transcripts cuja linha "Reading from ..." leia a malha do nível (``CT_V1_3zonas.msh.h5`` para a
    G0; ``CT_V1_3zonas_<nivel>.msh.h5`` para as demais) e pega as linhas "<N> cells, <k> cell zones" e
    "total volume (m3): <V>". Usa o primeiro transcript (em ordem de nome) e avisa se outro divergir em N.

    Args:
        nivel (str): ``"G0"``, ``"G1"`` ou ``"G2"``.

    Returns:
        dict | None: ``{"N", "V", "fonte", "avisos"}`` ou ``None`` se não houver transcript utilizável.
    """
    esperado = "CT_V1_3zonas.msh.h5" if nivel == "G0" else f"CT_V1_3zonas_{nivel}.msh.h5"
    achados = []
    for arq in sorted((PASTAS[nivel] / "fluent").glob("*.trn")):
        txt = arq.read_text(encoding="utf-8", errors="replace")
        if not re.search(r'Reading from [^\n]*[\\/]' + re.escape(esperado) + '"', txt):
            continue
        mn = re.search(r"(\d+)\s+cells,\s+\d+\s+cell zones", txt)
        mv = re.search(r"total volume \(m3\):\s*([\d.eE+-]+)", txt)
        if mn and mv:
            achados.append((arq.name, int(mn.group(1)), float(mv.group(1))))
    if not achados:
        return None
    avisos = [f"{a[0]}: N = {a[1]} difere do primeiro transcript" for a in achados[1:] if a[1] != achados[0][1]]
    return {"N": achados[0][1], "V": achados[0][2], "fonte": f"transcript do solver {achados[0][0]}",
            "avisos": avisos}


def ler_malha(nivel, n0=None, v0=None):
    """Obtém a contagem de células e o volume de um nível.

    Prioridade: (G0) ``n0``/``v0`` da linha de comando; transcript do solver (``ler_transcript``);
    ``LOGS/<nivel>_parametros.json`` (``estatisticas.celulas`` e ``estatisticas.volume_m3``); (G0) ``G0_PADRAO``.
    Se só o volume do G1/G2 faltar, usa o da G0 e avisa.

    Args:
        nivel (str): ``"G0"``, ``"G1"`` ou ``"G2"``.
        n0 (int | None): Células da G0 (linha de comando).
        v0 (float | None): Volume da G0 em m³ (linha de comando).

    Returns:
        dict | None: ``{"N", "V", "h", "fonte", "avisos"}`` (``h`` em m), ou ``None`` se não houver dados.

    Raises:
        ValueError: Se o JSON de parâmetros existir mas for inválido (``json.JSONDecodeError`` é subclasse).
    """
    def monta(n, v, fonte, avisos):
        """Monta o dicionário de retorno com h = (V/N)^(1/3).

        Args:
            n (int): Células.
            v (float): Volume em m³.
            fonte (str): Origem dos dados.
            avisos (list[str]): Avisos.

        Returns:
            dict: ``{"N", "V", "h", "fonte", "avisos"}``.
        """
        return {"N": int(n), "V": float(v), "h": (v / n) ** (1 / 3), "fonte": fonte, "avisos": avisos}

    if nivel == "G0" and (n0 or v0):
        n, v = n0 or G0_PADRAO["celulas"], v0 or G0_PADRAO["volume_m3"]
        return monta(n, v, "linha de comando" + ("" if (n0 and v0) else " (parcial; o resto do padrão)"), [])
    t = ler_transcript(nivel)
    if t:
        return monta(t["N"], t["V"], t["fonte"], t["avisos"])
    arq = LOGS / f"{nivel}_parametros.json"
    if arq.exists():
        j = json.loads(arq.read_text(encoding="utf-8"))
        est = j.get("estatisticas") or {}
        n, v = est.get("celulas", j.get("celulas")), est.get("volume_m3", j.get("volume_m3"))
        if n:
            avisos = []
            if not v:
                v = v0 or G0_PADRAO["volume_m3"]
                avisos.append("volume não lido em %s: usei o da G0" % arq.name)
            return monta(n, v, arq.name, avisos)
    if nivel == "G0":
        return monta(G0_PADRAO["celulas"], G0_PADRAO["volume_m3"], G0_PADRAO["fonte"], [])
    return None


def carregar_nivel(nivel, n0=None, v0=None):
    """Carrega malha e resultados (um por SR) de um nível, listando o que falta.

    Args:
        nivel (str): ``"G0"``, ``"G1"`` ou ``"G2"``.
        n0 (int | None): Células da G0 (só vale para a G0).
        v0 (float | None): Volume da G0 em m³ (só vale para a G0).

    Returns:
        tuple[dict, list[str]]: ``({"malha": dict | None, "pontos": {SR: grandezas_ponto(...)}, "avisos_serie":
        {SR: str}}, faltas)`` com ``faltas`` descrevendo cada arquivo ou dado ausente.

    Raises:
        ValueError: Se um JSON estiver corrompido (propagado de ``ler_malha`` e da leitura do ``resultado.json``).
        KeyError: Se o ``resultado.json`` não tiver ``media_janela`` (propagado de ``grandezas_ponto``).
    """
    faltas, pontos, avisos = [], {}, {}
    malha = ler_malha(nivel, n0, v0)
    if malha is None:
        faltas.append(f"{nivel}: sem contagem de células/volume (transcript .trn em {PASTAS[nivel] / 'fluent'} ou "
                      f"{LOGS / (nivel + '_parametros.json')}; gerar com tools/malha_gci_conversor.py {nivel})")
    for sr in SRS:
        arq = pasta_sr(PASTAS[nivel], sr) / "dados" / "resultado.json"
        if not arq.exists():
            faltas.append(f"{nivel}: falta {arq}")
            continue
        r = json.loads(arq.read_text(encoding="utf-8"))
        serie, info = serie_janela(nivel, r, sr)
        pontos[sr] = grandezas_ponto(r, serie)
        avisos[sr] = "" if serie is not None else f"{nivel} SR {sr:g}: sem série ({info}); SEM := σ por iteração"
    return {"malha": malha, "pontos": pontos, "avisos_serie": avisos}, faltas


def validar_malhas(malhas):
    """Confere a sanidade de V e de r entre os níveis (aborta o uso dos resultados se falhar).

    Args:
        malhas (list[dict]): ``ler_malha`` de G0, G1 e G2.

    Returns:
        list[str]: Problemas encontrados (vazia se tudo certo): |V_i/V_0 − 1| > 2 % e r_ij a mais de 15 % do
        nominal 1,35.
    """
    pb = []
    for n, m in zip(NIVEIS[1:], malhas[1:]):
        dv = m["V"] / malhas[0]["V"] - 1
        if abs(dv) > TOL_V:
            pb.append(f"V({n}) difere {100 * dv:+.1f} % da G0 (limite 2 %)")
    for nome, a, b in (("r21", malhas[1], malhas[0]), ("r32", malhas[2], malhas[1])):
        r = a["h"] / b["h"]
        if abs(r / R_NOMINAL - 1) > TOL_R:
            pb.append(f"{nome} = {r:.3f} está a {100 * (r / R_NOMINAL - 1):+.0f} % do nominal {R_NOMINAL} (limite 15 %)")
    return pb


# ---------------------------------------------------------------------------------------------------------------
# Cálculo para o conjunto de grandezas
# ---------------------------------------------------------------------------------------------------------------
def pct(x):
    """Converte fração em porcentagem preservando ``nan``.

    Args:
        x (float): Fração.

    Returns:
        float: 100·x, ou ``nan`` se x não for finito.
    """
    return 100 * x if np.isfinite(x) else NAN


def calcular_grandeza(g, sr, pts, malhas, tb1, p_formal, n_mc, problemas, avisos_serie):
    """Calcula o GCI de uma grandeza num SR e monta a linha do CSV.

    Args:
        g (str): Nome da grandeza (chave de ``GRANDEZAS``).
        sr (float): Razão de velocidades.
        pts (list[dict]): ``grandezas_ponto`` de G0, G1 e G2 neste SR.
        malhas (list[dict]): ``ler_malha`` de G0, G1 e G2.
        tb1 (float): T_b da malha fina (N·m), para ``GCI21_rel_Tb_pct``.
        p_formal (float): Ordem formal de reserva.
        n_mc (int): Amostras do Monte Carlo (0 desliga).
        problemas (list[str]): Saída de ``validar_malhas``.
        avisos_serie (list[str]): Avisos de séries ausentes (um por nível/SR; vazios são ignorados).

    Returns:
        dict: Linha com as chaves de ``COLUNAS`` (e ``"_h"`` com h1, h2, h3 e ``"_c"`` com o resultado cru de
        ``celik``, usados no gráfico).
    """
    vals = [p[g]["v"] for p in pts]
    sems = [p[g]["sem"] for p in pts]
    taus = [p[g]["tau"] for p in pts]
    h = [m["h"] for m in malhas]
    r21, r32 = h[1] / h[0], h[2] / h[1]
    obs = []
    ordenadas = r21 > 1 and r32 > 1
    if problemas:
        obs.append("MALHA INVÁLIDA: " + " | ".join(problemas))
    if not ordenadas:
        obs.append("h não cresce da fina para a grossa (r ≤ 1): sem GCI")
    elif min(r21, r32) < 1.3:
        obs.append("r < 1,3 (recomendação de Celik)")
    for m in malhas:
        obs += m["avisos"]
    obs += [a for a in avisos_serie if a]
    for n, p in zip(NIVEIS, pts):
        if p["_meta"]["convergiu"] is not True:
            obs.append(f"{n}: ponto {'não convergido' if p['_meta']['convergiu'] is False else 'sem registro de convergência'}")
    phi1, phi2, phi3 = vals
    if g == "eta" and sr == 0:                # η = SR·TR é identicamente zero em SR = 0
        c = celik(1.0, 1.0, 1.0, 1.5, 1.5, p_formal)
        obs.append("η = SR·TR = 0 em SR = 0 (sem sentido de GCI)")
    elif ordenadas:
        c = celik(phi1, phi2, phi3, r21, r32, p_formal)
        if c["tipo"] == "indefinido":
            obs.append("ε21 = 0 ou ε32 = 0: GCI indefinido")
    else:
        c = celik(1.0, 1.0, 1.0, 1.5, 1.5, p_formal)
    u21 = 2 * math.hypot(sems[0], sems[1])
    u32 = 2 * math.hypot(sems[1], sems[2])
    q21 = u21 / abs(c["eps21"]) if c["eps21"] else NAN
    q32 = u32 / abs(c["eps32"]) if c["eps32"] else NAN
    ruido = bool((np.isfinite(q21) and q21 >= 1) or (np.isfinite(q32) and q32 >= 1))
    if ruido:
        obs.append("diferença dentro do ruído iterativo (|ε| < U_iter = 2·√(SEM_i²+SEM_j²)): GCI não conclusivo")
    if any(np.isfinite(t) and t > N_LOTES * 8 for t in taus):
        obs.append("τ integral > 40 it: lotes de 40 iterações menores que τ; SEM por lotes pode estar subestimado")
    if c["p_origem"] not in ("aparente", "indefinido") and c["tipo"] != "indefinido":
        obs.append("p formal = %g usado (%s)" % (p_formal, c["p_origem"]))
    if c["tipo"] == "oscilatória":
        obs.append("convergência oscilatória: p e GCI com cautela (3 malhas não bastam para confirmar)")
    if c["tipo"] in ("divergente", "oscilatória divergente"):
        obs.append(f"{c['tipo']}: GCI só de referência (p formal); sem φ_ext")
    ra = c["razao_assintotica"]
    if c["p_origem"] == "aparente":   # p ajustado aos 3 pontos: a razão sai igual a |φ1/φ2| por construção
        obs.append("razão assintótica = |φ1/φ2| por construção (p aparente): não testa a faixa assintótica")
        if abs(c["p"] - p_formal) / p_formal >= 0.1 and np.isfinite(c["GCI21_cons"]):
            obs.append("p aparente %.2f longe do formal %g: use o GCI conservador (Fs = %g) = %.3g %%"
                       % (c["p"], p_formal, c["Fs_cons"], 100 * c["GCI21_cons"]))
    elif np.isfinite(ra) and phi2 != 0 and not 0.8 <= ra / abs(phi1 / phi2) <= 1.2:
        obs.append("razão assintótica %.2f longe de |φ1/φ2| = %.2f: fora da faixa assintótica" % (ra, abs(phi1 / phi2)))
    rel_tb = NAN
    if g in ("T_b", "T_t", "T_e") and np.isfinite(c["GCI21"]):
        rel_tb = c["GCI21"] * abs(phi1) / abs(tb1)
        if g == "T_e":
            obs.append("T_e pequeno: use GCI21_rel_Tb_pct (erro relativo a |T_b|), não GCI21_pct")
    udisc = c["GCI21"] * abs(phi1) if np.isfinite(c["GCI21"]) else NAN
    mc = (monte_carlo(vals, sems, r21, r32, p_formal, n_mc) if ordenadas and c["tipo"] != "indefinido"
          else monte_carlo(vals, sems, r21, r32, p_formal, 0))
    valida = not problemas and ordenadas
    conclusivo = bool(valida and not ruido and c["tipo"] == "monotônica" and c["p_origem"] == "aparente"
                      and np.isfinite(c["GCI21"]))
    lin = {"SR": sr, "grandeza": g, "phi1": phi1, "phi2": phi2, "phi3": phi3,
           "N1": malhas[0]["N"], "N2": malhas[1]["N"], "N3": malhas[2]["N"],
           "h1_mm": 1e3 * h[0], "h2_mm": 1e3 * h[1], "h3_mm": 1e3 * h[2], "r21": r21, "r32": r32,
           "p": c["p"], "p_origem": c["p_origem"], "R": c["R"], "tipo": c["tipo"], "phi_ext21": c["phi_ext21"],
           "e_a21_pct": pct(c["e_a21"]), "e_ext21_pct": pct(c["e_ext21"]), "GCI21_pct": pct(c["GCI21"]),
           "GCI32_pct": pct(c["GCI32"]), "razao_assintotica": ra, "GCI21_rel_Tb_pct": pct(rel_tb),
           "Fs_cons": c["Fs_cons"], "p_cons": c["p_cons"], "GCI21_cons_pct": pct(c["GCI21_cons"]),
           "sem1": sems[0], "sem2": sems[1], "sem3": sems[2], "tau1": taus[0], "tau2": taus[1], "tau3": taus[2],
           "metodo_sem": pts[0]["_meta"]["metodo"] + (" (G0)" if len({p["_meta"]["metodo"] for p in pts}) > 1 else ""),
           "Uiter21_sobre_eps21": q21, "Uiter32_sobre_eps32": q32,
           "Uiter_sobre_Udisc": 2 * sems[0] / udisc if udisc and np.isfinite(udisc) else NAN,
           "dentro_do_ruido": ruido, **mc, "malha_valida": valida, "conclusivo": conclusivo,
           "observacoes": "; ".join(obs), "_h": h, "_c": c}
    return lin


def calcular_tudo(niveis, p_formal=2.0, n_mc=300):
    """Calcula o GCI de todas as grandezas em todos os SR disponíveis nos três níveis.

    Args:
        niveis (dict): ``{"G0": ..., "G1": ..., "G2": ...}`` com a saída de ``carregar_nivel`` (primeiro item).
        p_formal (float): Ordem formal de reserva.
        n_mc (int): Amostras do Monte Carlo (0 desliga).

    Returns:
        list[dict]: Linhas (ver ``calcular_grandeza``), uma por SR x grandeza, só para SR presentes nos 3 níveis.
    """
    malhas = [niveis[n]["malha"] for n in NIVEIS]
    problemas = validar_malhas(malhas)
    linhas = []
    for sr in SRS:
        if not all(sr in niveis[n]["pontos"] for n in NIVEIS):
            continue
        pts = [niveis[n]["pontos"][sr] for n in NIVEIS]
        av = [niveis[n]["avisos_serie"].get(sr, "") for n in NIVEIS]
        for g in GRANDEZAS:
            linhas.append(calcular_grandeza(g, sr, pts, malhas, pts[0]["T_b"]["v"], p_formal, n_mc, problemas, av))
    return linhas


# ---------------------------------------------------------------------------------------------------------------
# Saídas
# ---------------------------------------------------------------------------------------------------------------
def gravar_csv(linhas, arq):
    """Grava o CSV (separador ``;``, ponto decimal).

    Args:
        linhas (list[dict]): Saída de ``calcular_tudo``.
        arq (Path): Arquivo de saída.
    """
    arq.parent.mkdir(parents=True, exist_ok=True)
    with arq.open("w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f, delimiter=";")
        w.writerow(COLUNAS)
        for lin in linhas:
            w.writerow([f"{lin[c]:.8g}" if isinstance(lin[c], float) else lin[c] for c in COLUNAS])


def fmt(x, formato, largura=9):
    """Formata um número para a tabela do console (``-`` se não for finito).

    Args:
        x (float): Valor.
        formato (str): Especificação de formato (``".4g"``).
        largura (int): Largura da coluna.

    Returns:
        str: Texto alinhado à direita.
    """
    return (format(x, formato) if np.isfinite(x) else "-").rjust(largura)


def imprimir_resumo(linhas):
    """Imprime a tabela-resumo no console.

    Args:
        linhas (list[dict]): Saída de ``calcular_tudo``.
    """
    print(f"{'SR':>4} {'grand.':<6}{'phi1':>10}{'phi2':>10}{'phi3':>10} {'tipo':<12}{'p':>6}{'phi_ext':>11}"
          f"{'GCI21[%]':>9}{'cons[%]':>9}{'raz.ass':>8}{'Uit/e21':>8}{'Uit/e32':>8}  concl.")
    for lin in linhas:
        print(f"{lin['SR']:>4.1f} {lin['grandeza']:<6}{lin['phi1']:>10.4g}{lin['phi2']:>10.4g}{lin['phi3']:>10.4g} "
              f"{lin['tipo'][:11]:<12}{lin['p']:>6.2f}{fmt(lin['phi_ext21'], '.4g', 11)}{fmt(lin['GCI21_pct'], '.2f')}"
              f"{fmt(lin['GCI21_cons_pct'], '.2f')}{fmt(lin['razao_assintotica'], '.2f', 8)}"
              f"{fmt(lin['Uiter21_sobre_eps21'], '.2f', 8)}{fmt(lin['Uiter32_sobre_eps32'], '.2f', 8)}"
              f"  {'sim' if lin['conclusivo'] else 'NÃO'}")
        if lin["grandeza"] in ("T_b", "T_t", "T_e") and np.isfinite(lin["GCI21_rel_Tb_pct"]):
            print(f"{'':>11}GCI21 relativo a |T_b| = {lin['GCI21_rel_Tb_pct']:.2f} %")
        if lin["observacoes"]:
            print(f"{'':>11}obs.: {lin['observacoes']}")


def grafico(linhas, arq):
    """Gera φ × h por SR (linhas) e grandeza (colunas), com φ_ext e barra de GCI na malha fina.

    Args:
        linhas (list[dict]): Saída de ``calcular_tudo``.
        arq (Path): PNG de saída.
    """
    srs = sorted({lin["SR"] for lin in linhas})
    fig, axs = plt.subplots(len(srs), len(GRANDEZAS), figsize=(3.1 * len(GRANDEZAS), 3.1 * len(srs)),
                            squeeze=False)
    for i, sr in enumerate(srs):
        for j, g in enumerate(GRANDEZAS):
            ax = axs[i][j]
            lin = next(x for x in linhas if x["SR"] == sr and x["grandeza"] == g)
            h = np.array(lin["_h"]) * 1e3
            phi = np.array([lin["phi1"], lin["phi2"], lin["phi3"]])
            ax.set(title=f"{ROTULOS[g]}, SR = {sr:g}", xlabel="h = (V/N)^(1/3) [mm]")
            ax.plot(h, phi, "o", color="C0", label="malhas G0, G1, G2")
            for hh, ph, tx in zip(h, phi, NIVEIS):
                ax.annotate(tx, (hh, ph), textcoords="offset points", xytext=(3, 4), fontsize=7)
            c = lin["_c"]
            if np.isfinite(c["phi_ext21"]) and np.isfinite(lin["p"]):
                hh = np.linspace(0, h[2] * 1.05, 60)
                C = (phi[0] - c["phi_ext21"]) / h[0] ** lin["p"]
                ax.plot(hh, c["phi_ext21"] + C * hh ** lin["p"], "-", color="C0", alpha=0.5, lw=1,
                        label=f"φ_ext + C·h^{lin['p']:.2f}")
                ax.plot([0], [c["phi_ext21"]], "s", color="C3", label="φ_ext")
                if np.isfinite(lin["GCI21_pct"]):
                    ax.errorbar([h[0]], [phi[0]], yerr=[lin["GCI21_pct"] / 100 * abs(phi[0])], fmt="none",
                                ecolor="C1", capsize=4, lw=1.6, label=f"GCI21 = {lin['GCI21_pct']:.2f} %")
            elif g == "eta" and sr == 0:
                ax.text(0.5, 0.5, "η = 0 em SR = 0", transform=ax.transAxes, ha="center")
            ax.set_xlim(left=0)
            if lin["dentro_do_ruido"]:
                ax.text(0.03, 0.04, "dentro do ruído iterativo", transform=ax.transAxes, color="C3", fontsize=7)
            if not lin["malha_valida"]:
                ax.text(0.03, 0.90, "MALHA INVÁLIDA", transform=ax.transAxes, color="C3", fontsize=8, weight="bold")
            if i == 0 and j == 0:
                ax.legend(fontsize=6)
            ax.set_ylabel(ROTULOS[g].split(" [")[0])
    fig.suptitle("Estudo de independência de malha (GCI, Celik et al. 2008) do conversor CT_V1, MRF")
    fig.tight_layout()
    arq.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(arq)
    plt.close(fig)


# ---------------------------------------------------------------------------------------------------------------
# Autoteste
# ---------------------------------------------------------------------------------------------------------------
def confere(nome, cond, detalhe=""):
    """Imprime o resultado de uma conferência do autoteste.

    Args:
        nome (str): Descrição.
        cond (bool): Resultado.
        detalhe (str): Texto extra.

    Returns:
        bool: ``bool(cond)``.
    """
    print(f"[{'OK' if cond else 'FALHOU'}] {nome} {detalhe}")
    return bool(cond)


def autoteste():
    """Valida o código com soluções conhecidas.

    Casos: (1) φ = φ_ext + C·h^p com r uniforme e não uniforme, inclusive r21 = 1,35 e r32 = 1,25 com p = 1 (que o
    critério antigo 0<R<1 chamava de divergente) e r32 ≫ r21 (onde o ponto fixo falha); (2) exemplo numérico da
    Tabela 1 de Celik et al. (2008) (entradas e p, φ_ext21 de memória: CONFERIR no artigo; e_a21, e_ext21 e GCI21 são
    verificados contra as fórmulas, não contra o valor publicado, que pode diferir do 2,17 % calculado aqui);
    (3) classificação oscilatória, oscilatória divergente, divergente e ε = 0; (4) pipeline completo.

    Returns:
        bool: True se todas as conferências passarem.
    """
    oks = []
    casos = ((1.8, 100.0, 5.0, (1.0, 1.5, 2.25)), (2.5, 30.0, 0.7, (1.0, 1.3, 2.08)),
             (0.9, 12.0, -2.0, (1.0, 1.4, 2.0)), (1.0, 50.0, 3.0, (1.0, 1.35, 1.6875)),
             (2.0, 7.0, 0.4, (1.0, 1.1, 2.2)))
    for p_ver, ext, C, hs in casos:
        phi = [ext + C * h ** p_ver for h in hs]
        r21, r32 = hs[1] / hs[0], hs[2] / hs[1]
        c = celik(*phi, r21, r32, p_formal=max(2.0, p_ver))   # p_formal ≥ p_ver: testa só a raiz
        oks.append(confere(f"sintético p={p_ver}, r21={r21:.2f}, r32={r32:.2f}: monotônica, recupera p e φ_ext",
                           abs(c["p"] - p_ver) < 1e-6 and abs(c["phi_ext21"] - ext) < 1e-6
                           and c["p_origem"] == "aparente" and c["tipo"] == "monotônica",
                           f"(p={c['p']:.8f}, φ_ext={c['phi_ext21']:.8f}, tipo={c['tipo']}, R={c['R']:.3f})"))
        e_abs = abs(phi[1] - phi[0]) * FS / (r21 ** p_ver - 1)
        oks.append(confere("   GCI21 = 1,25·|ε21|/(r21^p−1)/|φ1| e razão assintótica = |φ1/φ2| (solução exata)",
                           abs(c["GCI21"] - e_abs / abs(phi[0])) < 1e-12
                           and abs(c["razao_assintotica"] / abs(phi[0] / phi[1]) - 1) < 1e-6,
                           f"(GCI21={100 * c['GCI21']:.4f} %, razão={c['razao_assintotica']:.4f})"))
    phi = [30.0 + 0.7 * h ** 2.5 for h in (1.0, 1.3, 2.08)]
    c = celik(*phi, 1.3, 2.08 / 1.3)
    oks.append(confere("p aparente 2,5 > 1,1·p_formal (2): limitado ao formal",
                       c["p"] == 2.0 and c["p_origem"].startswith("p aparente"), f"(p={c['p']}, {c['p_origem']})"))
    r21, r32 = math.sqrt(18000 / 8000), math.sqrt(8000 / 4500)       # Celik 2008, Tab. 1 (2D)
    f1, f2, f3 = 6.063, 5.972, 5.863
    c = celik(f1, f2, f3, r21, r32)
    a = r21 ** c["p"]
    oks.append(confere("Celik 2008 Tab. 1: p ≈ 1,53 e φ_ext21 ≈ 6,168 (de memória; conferir no artigo)",
                       abs(c["p"] - 1.53) < 0.02 and abs(c["phi_ext21"] - 6.168) < 0.01,
                       f"(r21={r21:.3f}, r32={r32:.3f}, p={c['p']:.3f}, φ_ext21={c['phi_ext21']:.4f})"))
    oks.append(confere("Celik Tab. 1: e_a21, e_ext21 e GCI21 batem com as fórmulas",
                       abs(c["e_a21"] - abs((f1 - f2) / f1)) < 1e-12
                       and abs(c["e_ext21"] - abs((c["phi_ext21"] - f1) / c["phi_ext21"])) < 1e-12
                       and abs(c["GCI21"] - 1.25 * abs((f1 - f2) / f1) / (a - 1)) < 1e-12,
                       f"(e_a21={100 * c['e_a21']:.2f} %, e_ext21={100 * c['e_ext21']:.2f} %, "
                       f"GCI21={100 * c['GCI21']:.2f} %; valor publicado do GCI21: conferir no artigo)"))
    oks.append(confere("oscilatória (|R|<1)", celik(10.0, 10.1, 9.95, 1.5, 1.5)["tipo"] == "oscilatória"))
    co = celik(10.0, 10.3, 10.1, 1.5, 1.5)
    oks.append(confere("oscilatória divergente (|R|>1): p formal, sem φ_ext",
                       co["tipo"] == "oscilatória divergente" and not np.isfinite(co["phi_ext21"]) and co["p"] == 2.0))
    cd = celik(1.0, 1.1, 1.15, 1.5, 1.5)
    oks.append(confere("divergente (r constante, R>1): sem φ_ext e com p formal",
                       cd["tipo"] == "divergente" and not np.isfinite(cd["phi_ext21"]) and cd["p"] == 2.0))
    ce = celik(1.0, 1.0, 1.2, 1.5, 1.5)
    oks.append(confere("ε21 = 0: indefinido e GCI NaN (não 0 %)",
                       ce["tipo"] == "indefinido" and not np.isfinite(ce["GCI21"])))
    oks.append(confere("p aparente fora de [0,5;10] cai em p formal",
                       celik(0.0, 1.0e-9, 1.0e-9 + 1.0e-3, 1.5, 1.5)["p_origem"] != "aparente"))
    fs, pu, g = gci_conservador(1.0, 0.99, 1.5, 1.95)
    fs2, pu2, g2 = gci_conservador(1.0, 0.99, 1.5, 0.7)
    oks.append(confere("GCI conservador: Fs=1,25 perto do formal; Fs=3 e p limitado fora dele",
                       fs == 1.25 and pu == 1.95 and fs2 == 3.0 and pu2 == 0.7 and g2 > g))
    x = np.random.default_rng(1).normal(size=200)
    sem, tau, sem_tau = sem_serie(x)
    oks.append(confere("SEM de ruído branco: ≈ σ/√n e τ ≈ 0,5", abs(sem / (1 / math.sqrt(200)) - 1) < 0.6
                       and tau < 1.5, f"(sem={sem:.4f}, τ={tau:.2f}, sem_τ={sem_tau:.4f})"))
    oks.append(autoteste_pipeline())
    return all(oks)


def escrever_ponto(niv, sr, hh, rng, com_out):
    """Grava um ``resultado.json`` (e opcionalmente o ``.out``) sintético de um ponto para o autoteste.

    Args:
        niv (str): Nível (``"G0"``, ``"G1"`` ou ``"G2"``).
        sr (float): Razão de velocidades.
        hh (float): "h" relativo do nível (1, 1,5, 2,25); T_b e T_t variam com hh² e T_e linearmente.
        rng (np.random.Generator): Gerador de números aleatórios.
        com_out (bool): Se True, grava o Report File com ruído AR(1) (φ = 0,9; σ = 30 N·m).
    """
    n, jan = 300, 200
    base = {"t_bomba": -40000 * (1 + 0.01 * hh ** 2), "t_turbina": (42000 + 5 * sr * 1000) * (1 - 0.012 * hh ** 2),
            "t_estator": -2500 + 40 * hh}
    it = np.arange(1, n + 1)
    ser = {}
    for k, b in base.items():
        e = np.zeros(n)
        for i in range(1, n):
            e[i] = 0.9 * e[i - 1] + rng.normal(scale=30 * math.sqrt(1 - 0.81))
        ser[k] = b + e
    d = pasta_sr(PASTAS[niv], sr) / "dados"
    d.mkdir(parents=True, exist_ok=True)
    w = slice(n - jan, n)
    (d / "resultado.json").write_text(json.dumps({
        "sr": sr, "omega_b": 209.4395, "omega_t": 209.4395 * sr, "rpm_b": 2000.0, "convergiu": True,
        "unidade_torque": "N·m", "janela": jan, "it_ini": 0, "it_fim": n,
        "media_janela": {k: float(v[w].mean()) for k, v in ser.items()},
        "desvio_janela": {k: float(v[w].std(ddof=1)) for k, v in ser.items()}}), encoding="utf-8")
    if com_out:
        (PASTAS[niv] / "dados").mkdir(parents=True, exist_ok=True)
        out = PASTAS[niv] / "dados" / f"rep_SR{int(round(sr * 100)):03d}.out"
        linhas = ['"mon"', '"Iteration" "t_bomba etc.."', '("Iteration" "t_bomba" "t_turbina" "t_estator")']
        linhas += [f"{i} {float(ser['t_bomba'][i - 1])!r} {float(ser['t_turbina'][i - 1])!r} {float(ser['t_estator'][i - 1])!r}" for i in it]
        out.write_text("\n".join(linhas), encoding="utf-8")


def autoteste_pipeline():
    """Roda o pipeline completo (arquivos, SEM por lotes, fallback, CSV e gráfico) com dados sintéticos.

    Usa uma pasta temporária (nada é gravado no repositório). G0 e G1 têm .out; a G2 não (testa o fallback).

    Returns:
        bool: True se todas as conferências desta etapa passarem.
    """
    import tempfile
    global CONV, DIR_GCI, LOGS, PASTAS
    salvo = (CONV, DIR_GCI, LOGS, dict(PASTAS))
    oks = []
    try:
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            CONV, DIR_GCI = tmp, tmp / "CT_V1_GCI"
            LOGS = DIR_GCI / "logs"
            PASTAS = {"G0": tmp / "G0", "G1": tmp / "G1", "G2": tmp / "G2"}
            LOGS.mkdir(parents=True)
            _CACHE_OUT.clear()
            rng = np.random.default_rng(7)
            Ns = {"G1": 1_900_736 / 1.35 ** 3, "G2": 1_900_736 / 1.35 ** 6}
            for niv, hh in (("G0", 1.0), ("G1", 1.35), ("G2", 1.8225)):
                if niv != "G0":
                    (LOGS / f"{niv}_parametros.json").write_text(json.dumps(
                        {"estatisticas": {"celulas": int(Ns[niv]), "volume_m3": 8.576847e-2}}), encoding="utf-8")
                for sr in SRS:
                    escrever_ponto(niv, sr, hh, rng, com_out=niv != "G2")
            niv, faltas = {}, []
            for n in NIVEIS:
                niv[n], f = carregar_nivel(n)
                faltas += f
            oks.append(confere("pipeline: lê G0 (constantes), G1 e G2 sem faltas", not faltas, str(faltas)))
            oks.append(confere("pipeline: série da janela achada na G0/G1 e ausente na G2",
                               niv["G0"]["avisos_serie"][0.4] == "" and "sem série" in niv["G2"]["avisos_serie"][0.4]))
            linhas = calcular_tudo(niv, n_mc=60)
            oks.append(confere("pipeline: 2 SR x 6 grandezas", len(linhas) == 12))
            tb = next(x for x in linhas if x["SR"] == 0.4 and x["grandeza"] == "T_b")
            oks.append(confere("pipeline: T_b monotônica, p ≈ 2, fora do ruído, malha válida",
                               tb["tipo"] == "monotônica" and abs(tb["p"] - 2) < 0.15 and not tb["dentro_do_ruido"]
                               and tb["malha_valida"],
                               f"(p={tb['p']:.3f}, GCI21={tb['GCI21_pct']:.3f} %, cons={tb['GCI21_cons_pct']:.3f} %, "
                               f"Uiter/ε21={tb['Uiter21_sobre_eps21']:.3f}, τ1={tb['tau1']:.1f}, "
                               f"p_mc=[{tb['p_mc_16']:.2f}; {tb['p_mc_84']:.2f}], {tb['metodo_sem']})"))
            te = next(x for x in linhas if x["SR"] == 0.0 and x["grandeza"] == "T_e")
            oks.append(confere("pipeline: T_e com ε ~ SEM marcado como ruído",
                               te["dentro_do_ruido"] and not te["conclusivo"],
                               f"(Uiter/ε21={te['Uiter21_sobre_eps21']:.2f}, Uiter/ε32={te['Uiter32_sobre_eps32']:.2f})"))
            oks.append(confere("pipeline: η em SR=0 sinalizado", "η = SR·TR = 0" in
                               next(x for x in linhas if x["SR"] == 0.0 and x["grandeza"] == "eta")["observacoes"]))
            gravar_csv(linhas, DIR_GCI / "dados" / "gci.csv")
            grafico(linhas, DIR_GCI / "graficos" / "gci_tr.png")
            oks.append(confere("pipeline: CSV e PNG gravados", (DIR_GCI / "dados" / "gci.csv").stat().st_size > 500
                               and (DIR_GCI / "graficos" / "gci_tr.png").stat().st_size > 5000))
            imprimir_resumo(linhas)
            m2 = [dict(m) for m in (niv[n]["malha"] for n in NIVEIS)]
            m2[1]["V"] *= 1.05
            oks.append(confere("sanidade: V 5 % diferente da G0 é marcado inválido", bool(validar_malhas(m2))))
            m3 = [dict(m) for m in (niv[n]["malha"] for n in NIVEIS)]
            m3[2]["h"] *= 1.3
            oks.append(confere("sanidade: r fora de ±15 % do nominal é marcado inválido", bool(validar_malhas(m3))))
            oks.append(confere("sanidade: malhas do teste são válidas", not validar_malhas([niv[n]["malha"] for n in NIVEIS])))
    finally:
        CONV, DIR_GCI, LOGS, PASTAS = salvo
        _CACHE_OUT.clear()
    return all(oks)


# ---------------------------------------------------------------------------------------------------------------
# Principal
# ---------------------------------------------------------------------------------------------------------------
def main(argv=None):
    """Interpreta a linha de comando, calcula o GCI e grava CSV e gráfico.

    Args:
        argv (list[str] | None): Argumentos; ``None`` usa ``sys.argv``.

    Returns:
        int: 0 se calculou; 1 se o autoteste falhou; 2 se faltam dados (G1/G2 ou algum ponto); 3 se a malha foi
        marcada inválida (V ou r fora dos limites; os arquivos são gravados mesmo assim, com a marca).
    """
    ap = argparse.ArgumentParser(description="GCI (Celik et al. 2008) do conversor CT_V1, MRF, com G0/G1/G2.")
    ap.add_argument("--n0", type=int, help="células da G0 (padrão: transcript do solver)", default=None)
    ap.add_argument("--v0", type=float, help="volume da G0 em m³ (padrão: transcript do solver)", default=None)
    ap.add_argument("--p-formal", type=float, default=2.0, help="ordem formal de reserva e do GCI conservador")
    ap.add_argument("--mc", type=int, default=300, help="amostras do Monte Carlo do SEM (0 desliga)")
    ap.add_argument("--teste", action="store_true", help="roda o autoteste com dados sintéticos e sai")
    a = ap.parse_args(argv)
    if a.teste:
        ok = autoteste()
        print("AUTOTESTE:", "tudo certo" if ok else "FALHOU")
        return 0 if ok else 1
    niveis, faltas = {}, []
    for n in NIVEIS:
        niveis[n], f = carregar_nivel(n, a.n0, a.v0)
        faltas += f
        m = niveis[n]["malha"]
        if m:
            print(f"{n}: N = {m['N']:,}, V = {m['V']:.6e} m³, h = {1e3 * m['h']:.3f} mm ({m['fonte']})")
    srs_ok = [sr for sr in SRS if all(sr in niveis[n]["pontos"] for n in NIVEIS)]
    if faltas or not srs_ok or any(niveis[n]["malha"] is None for n in NIVEIS):
        print("\nNão dá para calcular o GCI ainda. Falta:")
        for f in faltas:
            print("  -", f)
        print("\nRode as malhas (tools/malha_gci_conversor.py G1 G2) e os pontos SR = 0 e 0,4 em G1/G2 antes.")
        if not srs_ok or any(niveis[n]["malha"] is None for n in NIVEIS):
            return 2
        print(f"(Seguindo só com os SR completos: {srs_ok})")
    linhas = calcular_tudo(niveis, a.p_formal, a.mc)
    imprimir_resumo(linhas)
    gravar_csv(linhas, DIR_GCI / "dados" / "gci.csv")
    grafico(linhas, DIR_GCI / "graficos" / "gci_tr.png")
    print("\nGravado:", DIR_GCI / "dados" / "gci.csv", "e", DIR_GCI / "graficos" / "gci_tr.png")
    if linhas and not linhas[0]["malha_valida"]:
        print("\n*** MALHA INVÁLIDA: " + linhas[0]["observacoes"].split("MALHA INVÁLIDA: ")[-1].split("; ")[0] + " ***")
        return 3
    return 0 if len(srs_ok) == len(SRS) and not faltas else 2


if __name__ == "__main__":
    sys.exit(main())
