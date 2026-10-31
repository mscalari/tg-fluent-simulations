"""Análise das execuções: métricas e gráficos de cada rodada.

Uso: python tools/analise_caso.py R02 [R03 ...]   (analisa as execuções e refaz registro e notas)
     python tools/analise_caso.py --todos          (analisa todas as que têm resultado)
     python tools/analise_caso.py --selftest       (confere as métricas no CSV de referência do guia)
"""
import json, math, re, sys
from pathlib import Path
import numpy as np
import pandas as pd
from casos import ROOT, CASES, TIPOS, caso

GHIA_CENTRO = (0.6172, 0.7344)  # Ghia et al. (1982), Tabela III: centro do vórtice principal, Re = 100
GHIA = pd.read_csv(ROOT / "refs" / "ghia1982.csv", comment="#")
DRIVER = pd.read_csv(ROOT / "refs" / "bfs_driver1985.csv", comment="#")
XR_H_EXP, XR_H_INC = 6.26, 0.10  # Driver e Seegmiller (1985), via NPARC


# ---------------------------------------------------------------- leitura

def ler_campo(path):
    """Lê o ``_campo.csv`` exportado pelo Fluent (Cell Center, vírgula).

    Args:
        path (str | Path): Caminho do CSV.

    Na cavidade rotativa (``cavr``), as colunas ``x-velocity``, ``y-velocity`` e ``velocity-magnitude``
    passam a ser as do referencial rotativo; as absolutas ficam com o sufixo ``-abs``.

    Returns:
        pd.DataFrame: Uma linha por célula, com os nomes de coluna sem espaços.
    """
    df = pd.read_csv(path, skipinitialspace=True)
    df.columns = df.columns.str.strip()
    if "relative-x-velocity" in df:
        for k, r in (("x-velocity", "relative-x-velocity"), ("y-velocity", "relative-y-velocity"),
                     ("velocity-magnitude", "rel-velocity-magnitude")):
            df[k + "-abs"] = df[k]
            df[k] = df[r]
    return df


def ler_rep(path):
    """Lê as linhas numéricas de um Report File do Fluent.

    Args:
        path (str | Path): Caminho do ``_rep.out``.

    Returns:
        np.ndarray: Matriz (n, 2) com iteração e valor.
    """
    rows = [ln.split() for ln in Path(path).read_text().splitlines()]
    return np.array([(int(p[0]), float(p[1])) for p in rows if len(p) >= 2 and p[0].isdigit()])


def ler_residuos(csv_path, log_path):
    """Histórico de resíduos a partir das linhas que o Fluent imprime no log, gravado em ``csv_path``.

    Args:
        csv_path (Path): ``_residuos.csv`` a (re)gravar.
        log_path (Path): Log da solução (stdout do runner); os ``continuacao*.log`` da mesma
            pasta entram em seguida, em ordem de nome.

    Returns:
        pd.DataFrame: Coluna ``iter`` e uma coluna por equação, com os nomes do cabeçalho do Fluent.
    """
    pat = re.compile(r"^\s+(\d+)((?:\s+[\d.]+e[-+]\d+)+)")
    cab = re.compile(r"^\s*iter\s+(continuity.*?)\s+time/iter")
    logs = [log_path] + sorted(log_path.parent.glob("continuacao*.log"))  # continuações (run_caso --continuar)
    cols, rows = ["continuity", "x-velocity", "y-velocity"], []
    for lp in logs:
        for ln in lp.read_text(errors="replace").splitlines():
            if m := cab.match(ln):
                cols = m.group(1).split()
            elif (m := pat.match(ln)) and len(v := m.group(2).split()) == len(cols):
                rows.append([m.group(1)] + v)
    res = pd.DataFrame(rows, columns=["iter"] + cols).astype(float)
    res["iter"] = res["iter"].astype(int)
    res = res.drop_duplicates("iter", keep="last")
    res.to_csv(csv_path, index=False)
    return res


def arquivos(c):
    """Caminhos dos arquivos de uma execução.

    Args:
        c (dict): Parâmetros completos do caso (``casos.caso``).

    Returns:
        dict[str, Path]: ``campo``, ``rep``, ``res``, ``setup``, ``met``, ``log`` e ``graf`` (pasta).
    """
    d, n = c["dir"], c["nome"]
    return dict(campo=d / "dados" / f"{n}_campo.csv", rep=d / "dados" / f"{n}_rep.out",
                res=d / "dados" / f"{n}_residuos.csv", setup=d / "dados" / f"{n}_setup.json",
                met=d / "dados" / f"{n}_metricas.json", log=d / "logs" / "run.log", graf=d / "graficos")


# ---------------------------------------------------------------- métricas

def u_analitica(d):
    """Perfil analítico de Couette-Poiseuille do caso-base (U=2, P=10).

    Args:
        d (float | np.ndarray): Distância à parede parada, d = −y [m].

    Returns:
        float | np.ndarray: u(d) = 2d + 20d(1−d) [m/s].
    """
    return 2 * d + 20 * d * (1 - d)


def perto(serie, valor, tol=1e-3):
    """Máscara das células cuja coordenada está a menos de ``tol`` de ``valor``.

    Os centros exportados têm ruído de ~1e-6 m (nós da malha fora da grade ideal), então
    comparação exata ou ``np.isclose`` padrão perde células. 1 mm é menor que meia célula
    em todas as malhas (a menor célula, N128, tem 7,8 mm).

    Args:
        serie (pd.Series): Coordenada das células.
        valor (float): Coordenada alvo [m].
        tol (float, optional): Tolerância [m]. Padrão 1e-3.

    Returns:
        pd.Series: Máscara booleana.
    """
    return (serie - valor).abs() < tol


def coluna(df, eixo, valor):
    """Células da coluna (ou linha) de centros mais próxima de uma coordenada.

    Args:
        df (pd.DataFrame): Campo lido por :func:`ler_campo`.
        eixo (str): ``"x-coordinate"`` ou ``"y-coordinate"``.
        valor (float): Coordenada desejada [m].

    Returns:
        tuple[float, pd.DataFrame]: Coordenada efetiva e as células, ordenadas pelo outro eixo.
    """
    vals = df[eixo].unique()
    v0 = vals[np.argmin(abs(vals - valor))]
    outro = "y-coordinate" if eixo == "x-coordinate" else "x-coordinate"
    return v0, df[perto(df[eixo], v0)].sort_values(outro)


def metricas_cp(df):
    """Perfil em x ≈ 2,5 contra a analítica, nos centros de célula (d = −y).

    Args:
        df (pd.DataFrame): Campo lido por :func:`ler_campo`.

    Returns:
        dict: ``x0``, ``ubar`` (média das células), ``eL2`` (‖u−uₐ‖₂/‖uₐ‖₂, em %), ``einf`` (max|u−uₐ|).
    """
    x0, col = coluna(df, "x-coordinate", 2.5)
    d, u = -col["y-coordinate"].to_numpy(), col["x-velocity"].to_numpy()
    ua = u_analitica(d)
    return dict(x0=float(x0), ubar=float(u.mean()), eL2=float(100 * np.sqrt(((u - ua) ** 2).sum() / (ua ** 2).sum())),
                einf=float(abs(u - ua).max()))


def max_du(df, cid_ref, detalhe=False):
    """Maior |Δu| célula a célula contra outra execução na mesma malha (casamento por coordenada).

    Args:
        df (pd.DataFrame): Campo desta execução.
        cid_ref (str): ID da execução de referência.
        detalhe (bool, optional): Se True, devolve também o RMS de |Δu| e o máximo fora da
            fileira de células junto à tampa da cavidade (onde ficam as singularidades dos cantos).

    Returns:
        float | dict | None: max|u − u_ref| [m/s] (ou dict ``max``, ``rms``, ``max_sem_tampa``
        com ``detalhe``), ou None se a referência ainda não existir.
    """
    p = arquivos(caso(cid_ref))["campo"]
    if not p.exists():
        return None
    ref = ler_campo(p)
    k = ["x-coordinate", "y-coordinate"]
    m = df[k + ["x-velocity"]].round({k[0]: 9, k[1]: 9}).merge(ref[k + ["x-velocity"]].round({k[0]: 9, k[1]: 9}),
                                                                on=k, suffixes=("", "_ref"))
    assert len(m) == len(df), "malhas diferentes"
    du = abs(m["x-velocity"] - m["x-velocity_ref"])
    if not detalhe:
        return float(du.max())
    tampa = perto(m[k[1]], m[k[1]].max())  # tolerância de 1 mm: as coordenadas têm ruído de ~1e-6 m
    return dict(max=float(du.max()), rms=float(np.sqrt((du ** 2).mean())), max_sem_tampa=float(du[~tampa].max()))


def metricas_ldc(df):
    """Centro do vórtice principal e extremos nas linhas centrais da cavidade.

    Args:
        df (pd.DataFrame): Campo lido por :func:`ler_campo` (com ``stream-function``).

    Returns:
        dict: ``xc``, ``yc`` (centro do vórtice), ``dist_ghia``, ``u_min_vert`` (x ≈ 0,5),
        ``v_max_hor``, ``v_min_hor`` (y ≈ 0,5) e ``u_max`` (max|u| no campo).
    """
    # O ψ do Fluent tem constante arbitrária: mede-se contra o valor médio na fileira de células junto às paredes.
    # Centro = célula de max|ψ − ψ_parede| (resolução de uma célula); interpolar se o bloco C precisar de mais.
    xs, ys = df["x-coordinate"], df["y-coordinate"]
    borda = perto(xs, xs.min()) | perto(xs, xs.max()) | perto(ys, ys.min()) | perto(ys, ys.max())
    i = (df["stream-function"] - df.loc[borda, "stream-function"].mean()).abs().idxmax()
    xc, yc = float(df.at[i, "x-coordinate"]), float(df.at[i, "y-coordinate"])
    _, vert = coluna(df, "x-coordinate", 0.5)
    _, hor = coluna(df, "y-coordinate", 0.5)
    g = comparar_ghia(df)
    return dict(xc=xc, yc=yc, dist_ghia=float(math.dist((xc, yc), GHIA_CENTRO)),
                u_min_vert=float(vert["x-velocity"].min()), v_max_hor=float(hor["y-velocity"].max()),
                v_min_hor=float(hor["y-velocity"].min()), u_max=float(df["x-velocity"].abs().max()),
                ghia_du_max=float(g["I"]["erro"].abs().max()), ghia_du_rms=float(np.sqrt((g["I"]["erro"] ** 2).mean())),
                ghia_dv_max=float(g["II"]["erro"].abs().max()), ghia_dv_rms=float(np.sqrt((g["II"]["erro"] ** 2).mean())))


def linha_central(df, eixo, comp):
    """Perfil de ``comp`` sobre a linha central (``eixo`` = 0,5), com os valores de parede nas pontas.

    Com N par não há centro de célula em 0,5: faz-se a média das duas colunas vizinhas (interpolação linear).

    Args:
        df (pd.DataFrame): Campo lido por :func:`ler_campo`.
        eixo (str): Coordenada fixa em 0,5 (``"x-coordinate"`` ou ``"y-coordinate"``).
        comp (str): ``"x-velocity"`` ou ``"y-velocity"``.

    Returns:
        tuple[np.ndarray, np.ndarray]: Posição ao longo da linha (0 a 1) e a componente; na parede
        vale 0, exceto u = 1 na tampa (ponta y = 1 da vertical x = 0,5).

    Raises:
        AssertionError: Se as duas colunas vizinhas não tiverem o mesmo número de células.
    """
    # Supõe malha uniforme com N par (N040, N128); com N ímpar usar a coluna em 0,5 direto.
    h = 1 / math.sqrt(len(df))
    _, a = coluna(df, eixo, 0.5 - h / 2)
    _, b = coluna(df, eixo, 0.5 + h / 2)
    assert len(a) == len(b), "colunas vizinhas de 0,5 diferentes"
    outro = "y-coordinate" if eixo == "x-coordinate" else "x-coordinate"
    s = a[outro].to_numpy()
    f = (a[comp].to_numpy() + b[comp].to_numpy()) / 2
    topo = 1.0 if (eixo, comp) == ("x-coordinate", "x-velocity") else 0.0
    return np.r_[0.0, s, 1.0], np.r_[0.0, f, topo]


def comparar_ghia(df):
    """Interpola o perfil do Fluent nos pontos de Ghia et al. (1982) e calcula o erro ponto a ponto.

    Args:
        df (pd.DataFrame): Campo da cavidade lido por :func:`ler_campo`.

    Returns:
        dict[str, pd.DataFrame]: ``"I"`` (u na vertical x = 0,5) e ``"II"`` (v na horizontal y = 0,5),
        com ``pos``, ``vel`` (Ghia), ``fluent`` e ``erro`` = fluent − Ghia, só nos pontos internos
        (as pontas são condição de contorno).
    """
    out = {}
    for tab, eixo, comp in (("I", "x-coordinate", "x-velocity"), ("II", "y-coordinate", "y-velocity")):
        s, f = linha_central(df, eixo, comp)
        g = GHIA[(GHIA["tabela"] == tab) & (GHIA["pos"] > 0) & (GHIA["pos"] < 1)].copy()
        g["fluent"] = np.interp(g["pos"], s, f)
        g["erro"] = g["fluent"] - g["vel"]
        out[tab] = g
    return out


def metricas_cavr(df, cid):
    """Métricas do bloco G (Adendo, 3.3) contra a referência R26, no referencial rotativo.

    e_L2 e e∞ usam a referência interpolada linearmente (triangulação) nos centros desta malha.
    Overshoot e subvalor: max u − max u_ref e min u − min u_ref na linha central horizontal (y = 0,5).
    CB = max|p(i+1) − 2p(i) + p(i−1)| / (max p − min p) na fileira de células mais próxima de y = 0,5.

    Args:
        df (pd.DataFrame): Campo lido por :func:`ler_campo`.
        cid (str): ID da execução.

    Returns:
        dict: ``u_max``, ``cb``, ``u_max_lin``, ``u_min_lin`` e, se R26 existir e não for a própria,
        ``eL2_u``, ``eL2_v`` [%], ``einf_u``, ``einf_v`` [m/s], ``overshoot``, ``subvalor`` [m/s];
        com sufixo ``_R27``, os mesmos erros contra o R27 (nas rodadas N040); no R27 também ``ref_dif``
        (max|u − u_R26| / max|V|_R26, em %), ``ref_dif_sem_tampa`` e ``ref_dif_rms``.
    """
    import matplotlib.tri as mtri
    _, lin = coluna(df, "y-coordinate", 0.5)
    p = lin["pressure"].to_numpy()
    _, ul = linha_central(df, "y-coordinate", "x-velocity")
    m = dict(u_max=float(np.hypot(df["x-velocity"], df["y-velocity"]).max()),
             cb=float(np.abs(p[2:] - 2 * p[1:-1] + p[:-2]).max() / np.ptp(p)),
             u_max_lin=float(ul.max()), u_min_lin=float(ul.min()))
    # R26 (PRESTO!) é a referência do Adendo; contra R27 (Second Order) a comparação não favorece o PRESTO!.
    for rid, suf in (("R26", ""), ("R27", "_R27")):
        pr = arquivos(caso(rid))["campo"]
        if cid in ("R26", rid) or not pr.exists():
            continue
        ref = ler_campo(pr)
        tri = mtri.Triangulation(ref["x-coordinate"], ref["y-coordinate"])
        for k, n in (("x-velocity", "u"), ("y-velocity", "v")):
            r = mtri.LinearTriInterpolator(tri, ref[k])(df["x-coordinate"], df["y-coordinate"])
            assert not np.ma.is_masked(r), "ponto fora da malha de referência"
            e = df[k].to_numpy() - np.asarray(r)
            m[f"eL2_{n}{suf}"] = float(100 * np.sqrt((e ** 2).sum() / (np.asarray(r) ** 2).sum()))
            m[f"einf_{n}{suf}"] = float(np.abs(e).max())
    pr = arquivos(caso("R26"))["campo"]
    if cid == "R26" or not pr.exists():
        return m
    ref = ler_campo(pr)
    _, ur = linha_central(ref, "y-coordinate", "x-velocity")
    m.update(overshoot=float(ul.max() - ur.max()), subvalor=float(ul.min() - ur.min()))
    if cid == "R27":
        vref = np.hypot(ref["x-velocity"], ref["y-velocity"]).max()
        d = max_du(df, "R26", detalhe=True)
        m.update(ref_dif=float(100 * d["max"] / vref), ref_dif_sem_tampa=float(100 * d["max_sem_tampa"] / vref),
                 ref_dif_rms=float(100 * d["rms"] / vref))
    return m


def metricas_bfs(df):
    """Comprimento de recolamento x_r pela troca de sinal de u na 1ª linha de células acima de y = 0.

    Args:
        df (pd.DataFrame): Campo lido por :func:`ler_campo`.

    Returns:
        dict: ``xr`` (medido a partir do degrau, x = 1) e ``xr_h`` (x_r / h, h = 0,5 m); se a bolha
        atravessa a saída, ``xr`` é None e ``xr_maior_que`` traz a cota inferior.
    """
    y0 = df["y-coordinate"][df["y-coordinate"] > 0].min()
    lin = df[perto(df["y-coordinate"], y0) & (df["x-coordinate"] > 1)].sort_values("x-coordinate")
    x, u = lin["x-coordinate"].to_numpy(), lin["x-velocity"].to_numpy()
    j = np.where((u[:-1] < 0) & (u[1:] >= 0))[0]
    if not len(j):
        # Sem troca de sinal com u < 0 na última célula: a bolha atravessa a saída e x_r só tem cota inferior.
        return dict(xr=None, xr_h=None, xr_maior_que=float(x[-1] - 1) if u[-1] < 0 else None)
    j = j[0]
    xs = x[j] - u[j] * (x[j + 1] - x[j]) / (u[j + 1] - u[j])
    return dict(xr=float(xs - 1), xr_h=float((xs - 1) / 0.5))


def ler_parede(c, zona):
    """Lê o CSV de uma parede exportada pelo runner.

    Args:
        c (dict): Parâmetros completos do caso.
        zona (str): Nome da zona de parede.

    Returns:
        pd.DataFrame: Uma linha por face, com os nomes de coluna sem espaços.
    """
    return ler_campo(c["dir"] / "dados" / f"{c['nome']}_parede_{zona}.csv")


def metricas_duto(c, df):
    """Métricas do duto: |V_s|max/U_b no plano xz, Re_b, Cf e y+.

    Args:
        c (dict): Parâmetros completos do caso.
        df (pd.DataFrame): Campo lido por :func:`ler_campo`.

    Returns:
        dict: ``ub`` (média de v por volume), ``vs_max``, ``vs_ub`` [%], ``re_b``, ``cf`` (2τ/ρU_b², com τ = 1 Pa
        pelo balanço de forças G·A = τ·P), ``yplus_max``, ``yplus_med``.
    """
    vol = df["cell-volume"].to_numpy()
    ub = float((df["y-velocity"] * vol).sum() / vol.sum())
    vs = float(np.hypot(df["x-velocity"], df["z-velocity"]).max())
    yp = ler_parede(c, "wall")["y-plus"]
    return dict(ub=ub, vs_max=vs, vs_ub=100 * vs / ub, re_b=c["rho"] * ub / c["mu"], cf=2 / (c["rho"] * ub ** 2),
                yplus_max=float(yp.max()), yplus_med=float(yp.mean()))


def paredes_bft(c):
    """Cf e Cp ao longo das paredes do degrau, com sinal e referência iguais aos de Driver e Seegmiller.

    O sinal de τ_x é fixado pela região recolada junto à saída, onde o escoamento segue em +x.
    O Cp é referido à pressão em x/H = −4 na parede inferior de entrada, como nos dados originais.

    Args:
        c (dict): Parâmetros completos do caso.

    Returns:
        tuple[pd.DataFrame, pd.DataFrame]: Paredes do lado do degrau e oposta, com ``x``, ``cf``, ``cp``, ``yplus``.
    """
    q = 0.5 * c["rho"] * c["u_ref"] ** 2
    ent = ler_parede(c, "wall_in_bottom").sort_values("x-coordinate")
    bai = ler_parede(c, "wall_down_bottom").sort_values("x-coordinate")
    top = ler_parede(c, "wall_top").sort_values("x-coordinate")
    sinal = np.sign(bai["x-wall-shear"].iloc[-1])
    pref = float(np.interp(-4, ent["x-coordinate"], ent["pressure"]))
    lado = pd.concat([ent, bai])
    fr = [pd.DataFrame(dict(x=d["x-coordinate"].to_numpy(), cf=(sinal * d["x-wall-shear"] / q).to_numpy(),
                            cp=((d["pressure"] - pref) / q).to_numpy(), yplus=d["y-plus"].to_numpy()))
          for d in (lado, top)]
    return fr[0], fr[1]


def metricas_bft(c, df):
    """Métricas do degrau turbulento: x_r por τ_w, μt na camada cisalhante, y+ e erros contra Driver.

    Args:
        c (dict): Parâmetros completos do caso.
        df (pd.DataFrame): Campo lido por :func:`ler_campo`.

    Returns:
        dict: ``xr``, ``xr_h``, ``erro_xr`` [%] contra 6,26; ``mut_max``, ``mut_mu_max``; ``yplus_max``, ``yplus_med``
        no fundo a jusante; ``cf_rms``, ``cp_step_rms``, ``cp_top_rms``.
    """
    lado, top = paredes_bft(c)
    b = lado[lado["x"] > 0]
    x, cf = b["x"].to_numpy(), b["cf"].to_numpy()
    j = np.where((cf[:-1] < 0) & (cf[1:] >= 0))[0]  # último negativo→positivo: ignora o vórtice de canto
    m = dict(xr=None, xr_h=None, erro_xr=None)
    if len(j):
        j = j[-1]
        xr = float(x[j] - cf[j] * (x[j + 1] - x[j]) / (cf[j + 1] - cf[j]))
        m.update(xr=xr, xr_h=xr, erro_xr=100 * (xr - XR_H_EXP) / XR_H_EXP)
    cam = df[(df["x-coordinate"] > 0) & (df["x-coordinate"] < 12) & (df["y-coordinate"] > 0.2) & (df["y-coordinate"] < 2.5)]
    m.update(mut_max=float(cam["viscosity-turb"].max()), mut_mu_max=float(cam["viscosity-turb"].max() / c["mu"]),
             yplus_max=float(b["yplus"].max()), yplus_med=float(b["yplus"].mean()))
    for nome, serie, d, col in (("cf_rms", "cf", lado, "cf"), ("cp_step_rms", "cp_step", lado, "cp"),
                                ("cp_top_rms", "cp_top", top, "cp")):
        e = DRIVER[(DRIVER["serie"] == serie) & (DRIVER["pos"] > -8) & (DRIVER["pos"] < 29)]
        d = d.sort_values("x")
        m[nome] = float(np.sqrt(((np.interp(e["pos"], d["x"], d[col]) - e["val"]) ** 2).mean()))
    return m


def calcular(cid):
    """Calcula as métricas de uma execução e grava ``_metricas.json``.

    Args:
        cid (str): ID da execução.

    Returns:
        dict: Métricas, iterações, convergência, tempo e comparações com outras execuções.
    """
    c = caso(cid)
    a = arquivos(c)
    df = ler_campo(a["campo"])
    rep = ler_rep(a["rep"])
    res = ler_residuos(a["res"], a["log"])
    st = json.loads(a["setup"].read_text())
    ult = res.iloc[-1]
    eqs = [k for k in res.columns if k != "iter"]
    t = c["tipo"]
    # No duto a continuidade fica fora do critério: o fator de escala é ~1e-16 (ver run_caso.setup).
    crit = [k for k in eqs if not (t == "duto" and k == "continuity")]
    m = dict(n_celulas=len(df), iteracoes=int(res["iter"].max()), tempo_s=st.get("tempo_parede_s"),
             convergiu=bool(all(ult[k] < 1e-6 for k in crit)), residuos={k: float(ult[k]) for k in eqs},
             var_monitor_50=float(np.ptp(rep[-50:, 1])), monitor_final=float(rep[-1, 1]))
    if t == "cp_per":
        m.update(metricas_cp(df))
        f = st["reports_finais"]
        m.update(fx_top=f["fx_top"], fx_bottom=f["fx_bottom"], balanco=f["fx_top"] + f["fx_bottom"])
        if cid == "R05":
            m["max_du_R02"] = max_du(df, "R02")
    elif t == "cp_io":
        m.update(metricas_cp(df))
        f = st["reports_finais"]
        m.update(p_in=f["p_in"], p_out=f["p_out"], dp=f["p_in"] - f["p_out"], fx_top=f["fx_top"],
                 fx_bottom=f["fx_bottom"])
    elif t == "ldc":
        m.update(metricas_ldc(df))
        if cid in ("R12", "R13", "R14", "R18", "R19", "R20"):
            m["max_du_R08"] = max_du(df, "R08")
        if cid in ("R18", "R20"):
            m["max_du_R19"] = max_du(df, "R19")
            m["du_R19"] = max_du(df, "R19", detalhe=True)
    elif t == "cavr":
        m.update(metricas_cavr(df, cid))
    elif t == "bfs":
        m.update(metricas_bfs(df))
    elif t == "duto":
        m.update(metricas_duto(c, df))
    elif t == "bft":
        m.update(metricas_bft(c, df))
    a["met"].write_text(json.dumps(m, indent=1))
    return m


def ordem_observada():
    """Ordem observada p = ln(e_k / e_{k+1}) / ln 2 na série M1–M4 (R01–R04), quando disponível.

    Returns:
        dict[str, float]: Chave ``"R01→R02"`` etc. com o p de e_L2.
    """
    serie = ["R01", "R02", "R03", "R04"]
    e = {}
    for k in serie:
        p = arquivos(caso(k))["met"]
        if p.exists():
            e[k] = json.loads(p.read_text())["eL2"]
    return {f"{a}→{b}": math.log(e[a] / e[b]) / math.log(2) for a, b in zip(serie, serie[1:]) if a in e and b in e}


# ---------------------------------------------------------------- gráficos

def graficos(cid, m):
    """Gera os gráficos de uma execução em ``graficos/``.

    ``campos.png``: |V| em cores com vetores (normalizados, para ver a direção mesmo onde a
    velocidade é baixa) e, na cavidade, linhas de ψ; e a pressão em cores.
    ``convergencia.png``: resíduos e monitor. ``perfis.png``: perfis do tipo de caso.

    Args:
        cid (str): ID da execução.
        m (dict): Métricas de :func:`calcular`.

    Returns:
        list[Path]: PNGs gravados.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.tri as mtri

    c = caso(cid)
    a = arquivos(c)
    a["graf"].mkdir(parents=True, exist_ok=True)
    df = ler_campo(a["campo"])
    if c["tipo"] == "duto":  # a seção fica no plano xz: desenha z no lugar de y e o escoamento secundário como vetores
        df = df.assign(**{"y-coordinate": df["z-coordinate"], "v-axial": df["y-velocity"], "y-velocity": df["z-velocity"]})
    x, y = df["x-coordinate"].to_numpy(), df["y-coordinate"].to_numpy()
    u, v = df["x-velocity"].to_numpy(), df["y-velocity"].to_numpy()
    tri = mtri.Triangulation(x, y)
    if c["tipo"] in ("bfs", "bft"):  # remove triângulos que atravessam o degrau (a região sólida)
        cx, cy = x[tri.triangles].mean(1), y[tri.triangles].mean(1)
        tri.set_mask((cx < 1) & (cy < 0.5) if c["tipo"] == "bfs" else (cx < 0) & (cy < 1))
    lx, ly = np.ptp(x), np.ptp(y)
    larg = 11 if lx / ly < 3 else 14
    alt = larg * ly / lx * 2.4 + 1.2 if lx / ly < 3 else larg * ly / lx * 2 + 2
    out = []

    fig, ax = plt.subplots(2, 1, figsize=(larg, min(alt, 13)))
    cf = ax[0].tricontourf(tri, df["velocity-magnitude"], levels=40, cmap="viridis")
    fig.colorbar(cf, ax=ax[0], label="|V| [m/s]")
    passo = max(1, len(df) // 900)
    s = slice(None, None, passo)
    mag = np.hypot(u, v) + 1e-30
    ax[0].quiver(x[s], y[s], (u / mag)[s], (v / mag)[s], color="w", alpha=0.75, scale=60 if lx / ly < 3 else 110,
                 width=0.0018)
    if "stream-function" in df and c["tipo"] in ("ldc", "bfs"):
        ax[0].tricontour(tri, df["stream-function"], levels=25, colors="k", linewidths=0.5, alpha=0.6)
    ax[0].set(title="Módulo da velocidade, vetores (direção) e linhas de corrente" if c["tipo"] in ("ldc", "bfs")
              else "Módulo da velocidade e vetores (direção)", xlabel="x [m]", ylabel="y [m]", aspect="equal")
    cp = ax[1].tricontourf(tri, df["pressure"], levels=40, cmap="coolwarm")
    fig.colorbar(cp, ax=ax[1], label="p [Pa]")
    tit = ("Pressão periódica (sem a parcela linear −G·x; o que sobra é ruído numérico)" if c["tipo"] == "cp_per"
           else "Pressão estática (manométrica)")
    ax[1].set(title=tit, xlabel="x [m]", ylabel="y [m]", aspect="equal")
    fig.suptitle(c["nome"])
    fig.tight_layout()
    out.append(a["graf"] / f"{c['nome']}_campos.png")
    fig.savefig(out[-1], dpi=120)
    plt.close(fig)

    res = ler_residuos(a["res"], a["log"])
    rep = ler_rep(a["rep"])
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
    for k in res.columns.drop("iter"):
        ax[0].semilogy(res["iter"], res[k].clip(lower=1e-20), label=k)
    ax[0].axhline(1e-6, color="k", ls="--", lw=0.8, label="critério 1e-6")
    ax[0].set(xlabel="iteração", ylabel="resíduo absoluto", title="Resíduos")
    ax[0].legend()
    ax[1].plot(rep[:, 0], rep[:, 1])
    mon = {"duto": "v médio na linha z = 0 [m/s]", "bft": "τ_x médio em wall_down_bottom [Pa]"}.get(c["tipo"])
    ax[1].set(xlabel="iteração", ylabel=mon or "u médio na linha [m/s]",
              title="Monitor" if mon else "Monitor (média de u na linha de medição)")
    if c["tipo"] == "cp_per":
        ax[1].axhline(4.3333, color="k", ls="--", lw=0.8, label="ū analítico 4,3333")
        ax[1].legend()
    fig.suptitle(c["nome"])
    fig.tight_layout()
    out.append(a["graf"] / f"{c['nome']}_convergencia.png")
    fig.savefig(out[-1], dpi=120)
    plt.close(fig)

    fig, ax = plt.subplots(1, 2, figsize=(12, 4.8))
    if c["tipo"].startswith("cp"):
        x0, col = coluna(df, "x-coordinate", 2.5)
        d = -col["y-coordinate"]
        dd = np.linspace(0, 1, 200)
        ax[0].plot(u_analitica(dd), dd, "k-", lw=1, label="analítica")
        ax[0].plot(col["x-velocity"], d, "o", mfc="none", label="Fluent")
        ax[0].set(xlabel="u [m/s]", ylabel="d = −y [m]", title=f"Perfil em x = {x0:.3f} m")
        ax[1].plot(col["x-velocity"] - u_analitica(d), d, "o-")
        ax[1].axvline(0, color="k", lw=0.5)
        ax[1].set(xlabel="u − u_analítica [m/s]", ylabel="d [m]", title="Erro pontual")
        if c["tipo"] == "cp_io":
            for xx in (0.05, 0.5, 4.95):
                xq, cq = coluna(df, "x-coordinate", xx)
                ax[0].plot(cq["x-velocity"], -cq["y-coordinate"], ".-", lw=0.7, label=f"x = {xq:.2f}")
        ax[0].legend()
    elif c["tipo"] == "ldc":
        s, f = linha_central(df, "x-coordinate", "x-velocity")
        ax[0].plot(f, s, "-", label="Fluent (x = 0,5, interpolado)")
        ax[1].plot(*linha_central(df, "y-coordinate", "y-velocity"), "-", label="Fluent (y = 0,5, interpolado)")
        gi, gii = GHIA[GHIA["tabela"] == "I"], GHIA[GHIA["tabela"] == "II"]
        ax[0].plot(gi["vel"], gi["pos"], "ks", mfc="none", label="Ghia et al. (1982)")
        ax[1].plot(gii["pos"], gii["vel"], "ks", mfc="none", label="Ghia et al. (1982)")
        ax[0].set(xlabel="u [m/s]", ylabel="y [m]",
                  title=f"u na vertical x = 0,5 (max|Δu| Ghia = {m['ghia_du_max']:.4f})")
        ax[1].set(xlabel="x [m]", ylabel="v [m/s]",
                  title=f"v na horizontal y = 0,5 (max|Δv| Ghia = {m['ghia_dv_max']:.4f})")
        for a_ in ax:
            a_.grid(alpha=0.3)
            a_.legend()
    elif c["tipo"] == "cavr":
        rc = arquivos(caso("R26"))["campo"]
        for d_, lab, est in ((df, "Fluent", "-"),) + (((ler_campo(rc), "R26 (referência)", "k--"),)
                                                      if cid != "R26" and rc.exists() else ()):
            s, f = linha_central(d_, "y-coordinate", "x-velocity")
            ax[0].plot(s, f, est, label=lab)
            s, f = linha_central(d_, "x-coordinate", "y-velocity")
            ax[1].plot(f, s, est, label=lab)
        ax[0].set(xlabel="x [m]", ylabel="u relativo [m/s]", title="u na horizontal y = 0,5 (overshoot/subvalor)")
        ax[1].set(xlabel="v relativo [m/s]", ylabel="y [m]", title="v na vertical x = 0,5")
        for a_ in ax:
            a_.grid(alpha=0.3)
            a_.legend()
    elif c["tipo"] == "duto":
        lin = df[perto(df["y-coordinate"], df["y-coordinate"].abs().min(), 1e-6 + 1e-3)].sort_values("x-coordinate")
        ax[0].plot(lin["x-coordinate"], lin["v-axial"], ".-")
        ax[0].set(xlabel="x [m]", ylabel="v [m/s]", title=f"Velocidade axial em z ≈ 0 (U_b = {m['ub']:.3f} m/s)")
        dg = df.assign(s=(df["x-coordinate"] - df["y-coordinate"]).abs()).nsmallest(40, "s").sort_values("x-coordinate")
        ax[1].plot(dg["x-coordinate"], 100 * np.hypot(dg["x-velocity"], dg["y-velocity"]) / m["ub"], ".-")
        ax[1].set(xlabel="x na diagonal [m]", ylabel="|V_s|/U_b [%]",
                  title=f"Escoamento secundário na diagonal (max no campo = {m['vs_ub']:.3g} %)")
        for a_ in ax:
            a_.grid(alpha=0.3)
    elif c["tipo"] == "bft":
        lado, top = paredes_bft(c)
        lado, top = lado.sort_values("x"), top.sort_values("x")
        ex = DRIVER.set_index("serie")
        ax[0].plot(lado["x"], lado["cf"], label="Fluent")
        ax[0].errorbar(ex.loc["cf", "pos"], ex.loc["cf", "val"], yerr=abs(ex.loc["cf", "val"] * ex.loc["cf", "dcf_rel"]),
                       fmt="ks", mfc="none", ms=4, label="Driver e Seegmiller (1985)")
        ax[0].axhline(0, color="k", lw=0.5)
        ax[0].set(xlim=(-8, 30), ylim=(-0.004, 0.006), xlabel="x/H", ylabel="Cf",
                  title=f"Cf na parede do degrau (x_r/H = {m['xr_h'] if m['xr_h'] is None else round(m['xr_h'], 3)}; exp. 6,26)")
        ax[1].plot(lado["x"], lado["cp"], label="Fluent, lado do degrau")
        ax[1].plot(top["x"], top["cp"], label="Fluent, parede oposta")
        ax[1].plot(ex.loc["cp_step", "pos"], ex.loc["cp_step", "val"], "ks", mfc="none", ms=4, label="Exp., lado do degrau")
        ax[1].plot(ex.loc["cp_top", "pos"], ex.loc["cp_top", "val"], "k^", mfc="none", ms=4, label="Exp., parede oposta")
        ax[1].set(xlim=(-8, 30), xlabel="x/H", ylabel="Cp (ref. x/H = −4)", title="Coeficiente de pressão")
        for a_ in ax:
            a_.grid(alpha=0.3)
            a_.legend()
    else:
        y0 = df["y-coordinate"][df["y-coordinate"] > 0].min()
        lin = df[perto(df["y-coordinate"], y0) & (df["x-coordinate"] > 1)].sort_values("x-coordinate")
        ax[0].plot(lin["x-coordinate"], lin["x-velocity"])
        ax[0].axhline(0, color="k", lw=0.5)
        if m.get("xr") is not None:
            ax[0].axvline(1 + m["xr"], color="r", ls="--", label=f"recolamento x = {1 + m['xr']:.3f} m")
            ax[0].legend()
        ax[0].set(xlabel="x [m]", ylabel="u [m/s]", title=f"u na 1ª célula acima do fundo (y = {y0:.4f} m)")
        for xx in (0.0125, 1.5, 2.45, df["x-coordinate"].max()):
            xq, cq = coluna(df, "x-coordinate", xx)
            ax[1].plot(cq["x-velocity"], cq["y-coordinate"], label=f"x = {xq:.3f}")
        ax[1].set(xlabel="u [m/s]", ylabel="y [m]", title="Perfis de u")
        ax[1].legend()
    fig.suptitle(c["nome"])
    fig.tight_layout()
    out.append(a["graf"] / f"{c['nome']}_perfis.png")
    fig.savefig(out[-1], dpi=120)
    plt.close(fig)
    return out


def animacao(cid):
    """Monta o GIF da evolução do campo ao longo das iterações a partir de ``dados/anim/*.csv``.

    Cada quadro: |V| em cores (escala fixa, a do campo final) e vetores normalizados (direção).

    Args:
        cid (str): ID da execução.

    Returns:
        Path | None: GIF gravado em ``graficos/``, ou None se a rodada não tem quadros.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.tri as mtri
    from matplotlib.animation import FuncAnimation, PillowWriter

    c = caso(cid)
    a = arquivos(c)
    quadros = sorted((c["dir"] / "dados" / "anim").glob("*_it*.csv"))
    if not quadros:
        return None
    dfs = [ler_campo(q) for q in quadros]
    x, y = dfs[0]["x-coordinate"].to_numpy(), dfs[0]["y-coordinate"].to_numpy()
    tri = mtri.Triangulation(x, y)
    vmax = dfs[-1]["velocity-magnitude"].max()
    s = slice(None, None, max(1, len(x) // 900))
    fig, ax = plt.subplots(figsize=(6.4, 5.6))

    def quadro(i):
        """Desenha o quadro ``i``.

        Args:
            i (int): Índice do quadro.
        """
        ax.clear()
        d = dfs[i]
        ax.tricontourf(tri, d["velocity-magnitude"].clip(upper=vmax), levels=np.linspace(0, vmax, 31), cmap="viridis")
        u, v = d["x-velocity"].to_numpy(), d["y-velocity"].to_numpy()
        mag = np.hypot(u, v) + 1e-30
        ax.quiver(x[s], y[s], (u / mag)[s], (v / mag)[s], color="w", alpha=0.75, scale=60, width=0.0025)
        it = int(re.search(r"_it(\d+)", quadros[i].stem).group(1))
        ax.set(title=f"{c['nome']}\niteração {it}", xlabel="x [m]", ylabel="y [m]", aspect="equal")

    anim = FuncAnimation(fig, quadro, frames=len(dfs))
    out = a["graf"] / f"{c['nome']}_animacao.gif"
    anim.save(out, writer=PillowWriter(fps=6))
    plt.close(fig)
    return out


# ---------------------------------------------------------------- resumo

def resumo_metricas(cid, m):
    """Texto ``nome=valor`` para a coluna ``Metrica_integral`` do registro.

    Args:
        cid (str): ID da execução.
        m (dict): Métricas de :func:`calcular`.

    Returns:
        str: Métricas principais do tipo de caso.
    """
    t = CASES[cid]["tipo"]
    if t == "cp_per":
        s = f"ubar={m['ubar']:.5f} m/s; Fx_sup(movel)={m['fx_top']:.4f} N/m; Fx_inf(parada)={m['fx_bottom']:.4f} N/m"
        if m.get("max_du_R02") is not None:
            s += f"; max|u_R05-u_R02|={m['max_du_R02']:.2e} m/s"
        return s
    if t == "cp_io":
        return f"dp={m['dp']:.4f} Pa; ubar(x=2.5)={m['ubar']:.5f} m/s"
    if t == "ldc":
        s = (f"centro=({m['xc']:.4f}, {m['yc']:.4f}) m; u_min(x=0.5)={m['u_min_vert']:.4f}; v_max/min(y=0.5)={m['v_max_hor']:.4f}/{m['v_min_hor']:.4f}"
             f"; Ghia max|du|={m['ghia_du_max']:.4f} max|dv|={m['ghia_dv_max']:.4f}")
        for k in ("max_du_R08", "max_du_R19"):
            if m.get(k) is not None:
                s += f"; {k}={m[k]:.2e} m/s ({100 * m[k] / m['u_max']:.3f}% de u_max)"
        return s
    if t == "cavr":
        s = f"cb={m['cb']:.3e}"
        for k in ("overshoot", "subvalor", "eL2_u", "eL2_v", "ref_dif"):
            if k in m:
                s += f"; {k}={m[k]:.4g}"
        return s
    if t == "duto":
        return (f"ub={m['ub']:.4f} m/s; vs_ub={m['vs_ub']:.4g} %; Re_b={m['re_b']:.0f}; Cf={m['cf']:.3e}; "
                f"y+max={m['yplus_max']:.2f}")
    if t == "bft" and m.get("xr") is not None:
        return (f"xr_h={m['xr_h']:.3f}; erro={m['erro_xr']:+.1f}%; mut/mu_max={m['mut_mu_max']:.0f}; y+max={m['yplus_max']:.2f}; "
                f"Cf_rms={m['cf_rms']:.2e}; Cp_rms(degrau/oposta)={m['cp_step_rms']:.3f}/{m['cp_top_rms']:.3f}")
    if m.get("xr") is not None:
        return f"xr={m['xr']:.4f} m; xr_h={m['xr_h']:.3f}"
    return f"xr>{m['xr_maior_que']:.3f} m (bolha atravessa a saida)" if m.get("xr_maior_que") else "xr=nao encontrado"


def analisar(cid):
    """Analisa uma execução: métricas, gráficos e animação.

    Args:
        cid (str): ID da execução.

    Returns:
        dict: Métricas de :func:`calcular`.
    """
    m = calcular(cid)
    graficos(cid, m)
    animacao(cid)
    print(cid, resumo_metricas(cid, m), "| convergiu:", m["convergiu"], "| it:", m["iteracoes"])
    return m


def selftest():
    """Confere as métricas no CSV de referência do guia (R02 pela GUI).

    Raises:
        AssertionError: Se ū, e_L2 ou e∞ divergirem do guia.
    """
    assert GHIA.groupby("tabela").size().to_dict() == {"I": 17, "II": 17}, "ghia1982.csv incompleto"
    ref = ROOT / "refs" / "r02cpp_m2_sim_sou_campo.csv"
    if not ref.exists():
        print("selftest parcial: CSV de referência do R02 ausente")
        return
    m = metricas_cp(ler_campo(ref))
    assert abs(m["ubar"] - 4.3499) < 5e-5 and abs(m["eL2"] - 0.2693) < 5e-5 and abs(m["einf"] - 0.01249) < 5e-6, m
    print("selftest ok:", m)


if __name__ == "__main__":
    if sys.argv[1] == "--selftest":
        selftest()
    else:
        ids = [k for k in CASES if arquivos(caso(k))["setup"].exists()] if sys.argv[1] == "--todos" else sys.argv[1:]
        # A série de malhas precisa de todas as métricas antes das notas (ordem observada).
        for k in ids:
            calcular(k)
        for k in ids:
            analisar(k)
