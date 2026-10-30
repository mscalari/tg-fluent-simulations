"""Análise das execuções: métricas, gráficos, registro_pyfluent.csv e nota de cada rodada.

Uso: python tools/analise_caso.py R02 [R03 ...]   (analisa as execuções e refaz registro e notas)
     python tools/analise_caso.py --todos          (analisa todas as que têm resultado)
     python tools/analise_caso.py --selftest       (confere as métricas no CSV de referência do guia)
"""
import csv, json, math, re, shutil, sys
from pathlib import Path
import numpy as np
import pandas as pd
from casos import ROOT, CASES, TIPOS, caso

ADMIN = ROOT / "00_admin"
REG_BASE = ADMIN / "registro_execucoes_v2.csv"
REG_PYF = ADMIN / "registro_pyfluent.csv"
VAULT = ROOT / "notas"
FEITO, PENDENTE = "✅", "⏳"
# Previsões do FVM 1D de 2ª ordem (guia, seção 2): ū, e_L2 [%], e∞ [m/s]
PREVISTO = {"M1": (4.4000, 1.082, 0.0500), "M2": (4.3500, 0.271, 0.0125),
            "M3": (4.3375, 0.068, 0.0031), "M4": (4.3344, 0.017, 0.0008)}
GHIA_CENTRO = (0.6172, 0.7344)  # Ghia et al. (1982), Tabela III: centro do vórtice principal, Re = 100
GHIA = pd.read_csv(ROOT / "refs" / "ghia1982.csv", comment="#")


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
        pd.DataFrame: Colunas ``iter``, ``continuity``, ``x-velocity`` e ``y-velocity``.
    """
    pat =re.compile(r"^\s+(\d+)\s+([\d.]+e[-+]\d+)\s+([\d.]+e[-+]\d+)\s+([\d.]+e[-+]\d+)\s")
    logs = [log_path] + sorted(log_path.parent.glob("continuacao*.log"))  # continuações (run_caso --continuar)
    rows = [m.groups() for lp in logs for ln in lp.read_text(errors="replace").splitlines() if (m := pat.match(ln))]
    res = pd.DataFrame(rows, columns=["iter", "continuity", "x-velocity", "y-velocity"]).astype(float)
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
    m = dict(n_celulas=len(df), iteracoes=int(res["iter"].max()), tempo_s=st.get("tempo_parede_s"),
             convergiu=bool(all(ult[k] < 1e-6 for k in ("continuity", "x-velocity", "y-velocity"))),
             residuos={k: float(ult[k]) for k in ("continuity", "x-velocity", "y-velocity")},
             var_monitor_50=float(np.ptp(rep[-50:, 1])), monitor_final=float(rep[-1, 1]))
    t = c["tipo"]
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
    x, y = df["x-coordinate"].to_numpy(), df["y-coordinate"].to_numpy()
    u, v = df["x-velocity"].to_numpy(), df["y-velocity"].to_numpy()
    tri = mtri.Triangulation(x, y)
    if c["tipo"] == "bfs":  # remove triângulos que atravessam o degrau (região sólida x < 1, y < 0,5)
        cx, cy = x[tri.triangles].mean(1), y[tri.triangles].mean(1)
        tri.set_mask((cx < 1) & (cy < 0.5))
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
    for k in ("continuity", "x-velocity", "y-velocity"):
        ax[0].semilogy(res["iter"], res[k], label=k)
    ax[0].axhline(1e-6, color="k", ls="--", lw=0.8, label="critério 1e-6")
    ax[0].set(xlabel="iteração", ylabel="resíduo absoluto", title="Resíduos")
    ax[0].legend()
    ax[1].plot(rep[:, 0], rep[:, 1])
    ax[1].set(xlabel="iteração", ylabel="u médio na linha [m/s]", title="Monitor (média de u na linha de medição)")
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


# ---------------------------------------------------------------- registro e notas

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
    if m.get("xr") is not None:
        return f"xr={m['xr']:.4f} m; xr_h={m['xr_h']:.3f}"
    return f"xr>{m['xr_maior_que']:.3f} m (bolha atravessa a saida)" if m.get("xr_maior_que") else "xr=nao encontrado"


def registrar():
    """Reescreve ``registro_pyfluent.csv`` com as 43 execuções e o status na 1ª coluna.

    Parte das linhas planejadas em ``registro_execucoes_v2.csv`` (o original não é alterado)
    e preenche os resultados das execuções que têm ``_metricas.json``. Status: ✅ feito, ⏳ pendente.
    """
    linhas = list(csv.DictReader(open(REG_BASE, encoding="utf-8"), delimiter=";"))
    campos = ["Status"] + list(linhas[0].keys())
    for r in linhas:
        cid = r["ID"]
        met = arquivos(caso(cid))["met"] if cid in CASES else None
        if not (met and met.exists()):
            r["Status"] = PENDENTE
            continue
        m = json.loads(met.read_text())
        c = caso(cid)
        st = json.loads(arquivos(c)["setup"].read_text())
        urf = st.get("urf", {})
        r.update(Status=FEITO, Convergiu_S_N="S" if m["convergiu"] else "N", Iteracoes=str(m["iteracoes"]),
                 Tempo_parede_s=str(m["tempo_s"]), URF_p=str(urf.get("pressure", "padrao Coupled")),
                 URF_mom=str(urf.get("mom", "padrao Coupled")),
                 Arquivos_gerados=f"results/{c['dir'].parent.name}/{c['nome']}/",
                 eL2_rel=f"{m['eL2']:.4f}%" if "eL2" in m else "", eInf=f"{m['einf']:.5f} m/s" if "einf" in m else "",
                 Metrica_integral=resumo_metricas(cid, m),
                 Obs=f"PyFluent {st.get('pyfluent')} + Fluent 2024 R1, serial 2ddp sem GUI. Malha {c['malha']}. "
                     "Residuos finais: " + ", ".join(f"{k}={v:.3e}" for k, v in m["residuos"].items()))
    with open(REG_PYF, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, campos, delimiter=";")
        w.writeheader()
        w.writerows(linhas)


def tabela_resultados(cid, m):
    """Linhas Markdown de resultado × esperado × critério para a nota da rodada.

    Args:
        cid (str): ID da execução.
        m (dict): Métricas de :func:`calcular`.

    Returns:
        list[str]: Linhas da tabela (sem cabeçalho).
    """
    c = caso(cid)
    t = c["tipo"]
    L = []
    if t == "cp_per":
        pu, pe, pi = PREVISTO[c["malha"]]
        L += [f"| ū (média das células em x = {m['x0']:.3f}) | {m['ubar']:.5f} m/s | FVM 1D: {pu} (analítico 4,3333) | ±0,5 % do analítico |",
              f"| e_L2 = ‖u−uₐ‖₂/‖uₐ‖₂ | {m['eL2']:.4f} % | FVM 1D: {pe} % | ≤ 0,5 % na M2 |",
              f"| e∞ = max|u−uₐ| | {m['einf']:.5f} m/s | FVM 1D: {pi} | — |",
              f"| Fx na parede móvel (`wall_top`) | {m['fx_top']:.4f} N/m | ±9 N/m | — |",
              f"| Fx na parede parada (`wall_bottom`) | {m['fx_bottom']:.4f} N/m | ±11 N/m | — |",
              f"| Soma das forças nas paredes | {m['balanco']:.4f} N/m | ±20 N/m (= G·L·H) | — |"]
        if "max_du_R02" in m:
            L.append(f"| max|u_R05 − u_R02| (célula a célula) | {m['max_du_R02']:.2e} m/s | ≈ 0 | < 1e-4 m/s |")
        if cid in ("R01", "R02", "R03", "R04"):
            for k, p in ordem_observada().items():
                if cid in k:
                    L.append(f"| Ordem observada p ({k}) | {p:.3f} | 2 | entre 1,8 e 2,2 |")
    elif t == "cp_io":
        L += [f"| Δp = p_entrada − p_saída (média por área) | {m['dp']:.4f} Pa | 20 Pa | ±1 % |",
              f"| ū em x = {m['x0']:.3f} | {m['ubar']:.5f} m/s | 4,3333 | — |",
              f"| e_L2 em x = {m['x0']:.3f} | {m['eL2']:.4f} % | — | — |",
              f"| e∞ em x = {m['x0']:.3f} | {m['einf']:.5f} m/s | — | — |"]
    elif t == "ldc":
        L += [f"| Centro do vórtice (célula de max|ψ − ψ_parede|) | ({m['xc']:.4f}; {m['yc']:.4f}) m | ≈ (0,617; 0,734) | distância {m['dist_ghia']:.4f} m (≤ 1 célula = {1 / math.sqrt(m['n_celulas']):.4f} m) |",
              f"| u mínimo na coluna mais próxima de x = 0,5 | {m['u_min_vert']:.4f} m/s | −0,2109 (Ghia, ponto y = 0,4531) | — |",
              f"| v máx / mín na linha mais próxima de y = 0,5 | {m['v_max_hor']:.4f} / {m['v_min_hor']:.4f} m/s | 0,1753 / −0,2453 (Ghia, nos pontos da tabela) | — |",
              f"| max|Δu| nos 15 pontos internos da Tabela I de Ghia (RMS) | {m['ghia_du_max']:.4f} m/s ({m['ghia_du_rms']:.4f}) | 0 | ≤ 0,02 m/s |",
              f"| max|Δv| nos 15 pontos internos da Tabela II de Ghia (RMS) | {m['ghia_dv_max']:.4f} m/s ({m['ghia_dv_rms']:.4f}) | 0 | ≤ 0,02 m/s |"]
        if m.get("max_du_R08") is not None:
            crit = "< 1e-4 m/s" if cid in ("R12", "R13", "R14") else "< 1 % de u_max (G1)"
            L.append(f"| max|u − u_R08| | {m['max_du_R08']:.2e} m/s ({100 * m['max_du_R08'] / m['u_max']:.3f} % de u_max) | ≈ 0 | {crit} |")
        if m.get("max_du_R19") is not None:
            L.append(f"| max|u − u_R19| (G1) | {m['max_du_R19']:.2e} m/s ({100 * m['max_du_R19'] / m['u_max']:.3f} % de u_max) | ≈ 0 | < 1 % de u_max |")
            d = m["du_R19"]
            L.append(f"| max|u − u_R19| sem a fileira junto à tampa (cantos singulares) | {d['max_sem_tampa']:.2e} m/s ({100 * d['max_sem_tampa'] / m['u_max']:.3f} %) | — | diagnóstico |")
            L.append(f"| RMS de |u − u_R19| | {d['rms']:.2e} m/s ({100 * d['rms'] / m['u_max']:.3f} %) | — | diagnóstico |")
    elif t == "cavr":
        L += [f"| max|V| relativo no campo | {m['u_max']:.4f} m/s | — | — |",
              f"| u máx / mín na horizontal y = 0,5 | {m['u_max_lin']:.4f} / {m['u_min_lin']:.4f} m/s | — | — |",
              f"| Índice de checkerboard CB (fileira y ≈ 0,5) | {m['cb']:.3e} | — | comparar entre esquemas |"]
        if "eL2_u" in m:
            L += [f"| e_L2(u) / e_L2(v) contra R26 | {m['eL2_u']:.3f} % / {m['eL2_v']:.3f} % | — | G3: R23 < R21, R22; G4: R24, R25 ≥ 5 % abaixo de R23 |",
                  f"| e∞(u) / e∞(v) contra R26 | {m['einf_u']:.4f} / {m['einf_v']:.4f} m/s | — | — |",
                  f"| Overshoot / subvalor (max u − max u_R26; min u − min u_R26) | {m['overshoot']:+.4f} / {m['subvalor']:+.4f} m/s | — | G2: R22 acima dos demais |"]
        if "eL2_u_R27" in m:
            L.append(f"| e_L2(u) / e_L2(v) contra R27 (Second Order N128, não favorece o PRESTO!) | {m['eL2_u_R27']:.3f} % / {m['eL2_v_R27']:.3f} % | — | diagnóstico |")
        if "ref_dif" in m:
            L += [f"| ref_dif = max|u − u_R26| / max|V|_R26 | {m['ref_dif']:.3f} % | — | < 1 % (senão a referência não vale) |",
                  f"| ref_dif sem a fileira junto à tampa (cantos singulares) | {m['ref_dif_sem_tampa']:.3f} % | — | diagnóstico |",
                  f"| RMS de |u − u_R26| / max|V|_R26 | {m['ref_dif_rms']:.3f} % | — | diagnóstico |"]
    else:
        xr = (f"{m['xr']:.4f} m (x_r/h = {m['xr_h']:.3f})" if m.get("xr") is not None else
              f"> {m['xr_maior_que']:.3f} m: a bolha atravessa a saída (domínio curto)" if m.get("xr_maior_que")
              else "não encontrado")
        L.append(f"| Comprimento de recolamento x_r | {xr} | o Roteiro não informa | — |")
    # Barras dentro das células (|Δu|, max|ψ|) quebram a tabela Markdown: escapa célula a célula.
    return ["| " + " | ".join(x.replace("|", "\\|") for x in l[2:-2].split(" | ")) + " |" for l in L]


def nota(cid, m, pngs):
    """Escreve a nota Markdown da rodada em ``notas/`` e copia os gráficos.

    Args:
        cid (str): ID da execução.
        m (dict): Métricas de :func:`calcular`.
        pngs (list[Path]): Gráficos gerados por :func:`graficos`.

    Returns:
        Path: Caminho da nota.
    """
    c = caso(cid)
    tp = TIPOS[c["tipo"]]
    st = json.loads(arquivos(c)["setup"].read_text())
    img = VAULT / "img"
    img.mkdir(parents=True, exist_ok=True)
    for p in pngs:
        shutil.copy(p, img / p.name)
    nomes_esq = {"first-order-upwind": "First Order Upwind (FOU)", "second-order-upwind": "Second Order Upwind (SOU)",
                 "quick": "QUICK", "third-order-muscl": "MUSCL (3ª ordem)", "second-order": "Second Order",
                 "standard": "Standard", "presto!": "PRESTO!"}
    urf = st.get("urf")
    urf_txt = f"pressão {urf['pressure']}, momento {urf['mom']}" if urf else "não se aplica (Coupled usa os padrões próprios)"
    mov = {k: v for k, v in st["paredes"].items() if v.get("motion_bc") == "Moving Wall"}
    cont = "".join(
        f"\n- **Desvio do protocolo (continuação):** +{x['iter_extra']} iterações a partir do .cas/.dat"
        + (f", com o corte por resíduo desligado (1e-12) até o monitor variar menos de {x['tol_monitor']:g} m/s em "
           "50 iterações (2ª metade do critério do guia: monitor estável)" if x.get("tol_monitor") else
           f", com o mesmo critério de resíduo ({x['residuo_abs']:g})")
        for x in st.get("continuacoes", []))
    txt = f"""# {cid} · {c['nome']}

**Status:** {FEITO} feito · **Bloco:** {c['dir'].parent.name} · **Convergiu:** {'sim' if m['convergiu'] else 'NÃO'} em {m['iteracoes']} iterações · **Tempo de parede:** {m['tempo_s']} s

## O que esta rodada testa
{c['muda']}

## Parâmetros de entrada
**Geometria.** {tp['geometria']}

**Malha.** `mesh/{c['malha']}.msh`, {m['n_celulas']} células quadriláteras.

**Fluido.** {tp['fluido']} Material `fluido-tg` (cópia do ar com ρ e μ trocados), atribuído à zona `{st['zona_fluido']}`.

**Condições de contorno.** {tp['contorno']}

**Modelo.** Pressure-based, regime permanente, 2D planar, laminar, sem equação de energia.

**Solver.**
- Acoplamento pressão-velocidade: **{st['methods']['p_v']}**
- Gradiente: Least Squares Cell Based
- Interpolação de pressão: **{nomes_esq.get(st['methods']['disc']['pressure'], st['methods']['disc']['pressure'])}**
- Momento: **{nomes_esq.get(st['methods']['disc']['mom'], st['methods']['disc']['mom'])}**
- Sub-relaxação (URF): {urf_txt}
- Critério de parada: resíduos absolutos < 1e-6 em continuity, x-velocity e y-velocity; teto de {st['iter_max']} iterações
- Inicialização: standard, velocidade e pressão nulas (igual em todas as rodadas){cont}

**Paredes móveis (lidas de volta do Fluent).** {', '.join(f"`{k}`: {v['vmag']['value']} m/s em {v['wall_translation']}, relativa={v['relative']}" for k, v in mov.items()) or 'nenhuma'}

**Monitor de convergência.** {', '.join(st['monitores'])} (média de u por área na linha de medição).

## O que se espera
{tp['esperado']}

## Resultados
| Grandeza | Fluent | Esperado | Critério |
|---|---|---|---|
""" + "\n".join(tabela_resultados(cid, m)) + f"""

Resíduos finais: {', '.join(f"{k} = {v:.3e}" for k, v in m['residuos'].items())}. Variação do monitor nas últimas 50 iterações: {m['var_monitor_50']:.2e} m/s.

## Gráficos
""" + "\n".join(f"![[{p.name}]]" for p in pngs) + f"""

## Arquivos
`results/{c['dir'].parent.name}/{c['nome']}/` com `fluent/` (.cas.h5, .dat.h5), `dados/` (campo, rep.out, resíduos, setup.json, métricas.json), `graficos/` e `logs/`.

Gerado por `tools/analise_caso.py`.
"""
    p = VAULT / f"{c['nome']}.md"
    p.write_text(txt, encoding="utf-8")
    return p


def analisar(cid):
    """Analisa uma execução: métricas, gráficos e nota.

    Args:
        cid (str): ID da execução.

    Returns:
        dict: Métricas de :func:`calcular`.
    """
    m = calcular(cid)
    pngs = graficos(cid, m)
    gif = animacao(cid)
    if gif:
        pngs.append(gif)
    nota(cid, m, pngs)
    print(cid, resumo_metricas(cid, m), "| convergiu:", m["convergiu"], "| it:", m["iteracoes"])
    return m


def selftest():
    """Confere as métricas no CSV de referência do guia (R02 pela GUI).

    Raises:
        AssertionError: Se ū, e_L2 ou e∞ divergirem do guia.
    """
    m = metricas_cp(ler_campo(ROOT / "refs" / "r02cpp_m2_sim_sou_campo.csv"))
    assert abs(m["ubar"] - 4.3499) < 5e-5 and abs(m["eL2"] - 0.2693) < 5e-5 and abs(m["einf"] - 0.01249) < 5e-6, m
    assert GHIA.groupby("tabela").size().to_dict() == {"I": 17, "II": 17}, "ghia1982.csv incompleto"
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
        registrar()
        print("registro:", REG_PYF)
