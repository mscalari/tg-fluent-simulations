"""Balanço de momento angular (Euler) por zona do conversor CT_V1 MRF, com termos convectivo, difusivo e de pressão.

Uso: PYTHONIOENCODING=utf-8 python tools/euler_interfaces_conversor.py            (todos os SRxxx com dados)
     PYTHONIOENCODING=utf-8 python tools/euler_interfaces_conversor.py 0 20 60     (só alguns SR, em %)
     ... acrescente --diag-vorticidade para o diagnóstico (inconclusivo) de gradientes absolutos/relativos

Física (regime permanente, eixo x, origem 0, n = normal para FORA da zona, dA = vetor área da face):
    T_parede_sobre_fluido + T_interfaces_sobre_fluido = ∮ ρ (r × u)_x (u·n) dA      (fluxo que SAI)
    T_interfaces_sobre_fluido = ∫ r × (−p n + τ·n) dA
logo, o torque do fluido nas paredes da zona (convenção de ``analise_conversor.euler_csv``) é
    T_euler = −sinal · [ H_saída − (D + P) ],   com
    H = Σ ρ (u·dA) (y w − z v)                       (convectivo, u ABSOLUTA, fluxo ρ u·dA: variante principal)
    H_rel = Σ ρ (u − ω×r)·dA (y w − z v)             (convectivo com fluxo RELATIVO: variante extra)
    D = Σ (y t_z − z t_y),  t_i dA = μ_eff Σ_j (∂u_i/∂x_j + ∂u_j/∂x_i) dA_j      (difusivo: viscoso + turbulento)
    P = −Σ p (y dA_z − z dA_y)                       (pressão; nula em superfície de revolução)
Em zona rotativa (MRF) a superfície da interface está presa ao referencial da zona, então o fluxo que a atravessa é
ρ (u − ω×r)·dA (ω = sinal·ω_zona no eixo x; estator: 0). Isso só coincide com ρ u·dA em superfície de revolução,
e a if_bt não é (Σ ρ ω (y dA_z − z dA_y) ≈ 7 % de ṁ no lado da bomba). Por isso as colunas ``*_rel`` trazem o
convectivo recalculado assim; as colunas sem sufixo continuam sendo o convectivo de ``euler_csv``.
H, D e P são torques/fluxos sobre/da zona vistos de cada lado da interface; como em ``euler_csv``, o valor na
interface é a média dos dois lados (F = (X_A − X_B)/2 para a zona A e −F para a zona B, X orientado para fora da
própria zona). A diferença entre os lados (X_A + X_B, que seria 0 numa solução conservativa) é reportada por interface.
Pressão: além da média dos dois lados (``T_pres``), grava-se ``T_pres_prop``, com P só do lado próprio de cada
zona: somado às paredes da mesma zona ele é invariante ao nível de pressão (∮ p0 n dA = 0 na superfície fechada da
zona), e ``dTdp0_*`` dá a sensibilidade a um deslocamento p0 do nível de pressão [N·m/bar]: −sinal·ΣI·1e5, com
I = Σ (y dA_z − z dA_y) (``dTdp0_parede`` = lado próprio = efeito sobre o torque de parede; ``dTdp0_media`` = efeito
sobre o termo da média; a diferença é o que a média deixa de cancelar).
O resultado é comparado com o torque de parede por zona de ``resultado.json`` (``desempenho()['parede_<zona>']``,
que inclui as paredes de sobra do lado da zona).

Aproximações e incertezas (os termos aqui são ESTIMATIVAS, não correções que fechem o balanço):
    * O export usa location="cell-center": cada linha é o valor da CÉLULA adjacente à face (ρ, u, μ_eff, ∇u, p),
      colocado na posição e com a área da face. Nesse pós-processamento por centro de célula, as vazões de um lado
      diferem até ±18 % da do solver e a reconstrução célula→face com o gradiente muda a razão Euler/parede em
      ±3–4 %. O resíduo de 1–3 % da razão está dentro dessa incerteza e é comparável ao desequilíbrio global de torque
      (até 1,2 %). Difusivo (~1 %) e pressão (~1,5 %, das faces inclinadas junto às arestas das pás na if_bt) são da
      mesma ordem e não se conservam entre os lados; não se deve ler a variação da razão como melhora do fechamento.
    * τ = μ_eff (∇u + ∇uᵀ), incompressível (∇·u = 0 anula o termo −(2/3) μ ∇·u δ). O termo −(2/3) ρ k δ_ij da
      hipótese de Boussinesq NÃO foi incluído porque k não foi exportado; ele é isotrópico e age como uma pressão
      (−(2/3)ρk n) nas faces inclinadas, logo da mesma ordem de P e não desprezível a priori.
    * μ_eff = μ_lam + μ_t vem do CSV (``viscosity-eff``), por célula.
    * Pressão: manométrica como no CSV antigo (a pressão operacional é a mesma das paredes do ``paredes_zona``).
    * Os CSVs novo (``_euler_``) e antigo têm o mesmo conjunto de pontos mas ordens diferentes e ``cellnumber`` é só
      índice de linha: casa-se por (x, y, z) exatas (junção 1:1 verificada; CSV truncado/vazio ignora o ponto).
    * A vazão média dos dois lados é comparada com os reports ``m_*`` de ``resultado.json`` (``media_janela``); aviso
      se a diferença passar de 5 %.
    * Rotação de corpo rígido não gera deformação: τ é o mesmo com u absoluta ou relativa. O diagnóstico de
      vorticidade dos gradientes (``--diag-vorticidade``) é INCONCLUSIVO e fica desligado por padrão.

Saídas: CT_V1_MRF/dados/euler_interfaces.csv (por SR e zona), CT_V1_MRF/dados/euler_interfaces_if.csv (por SR e
interface: conservação entre os lados), CT_V1_MRF/graficos/g7_euler_interfaces.png e tabela no console.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import analise_conversor as ac  # importa matplotlib (Agg); sem efeitos além de rcParams
import matplotlib.pyplot as plt

DIR = ac.DIR
PARES = {"bomba_turbina": ("fluido_bomba", "fluido_turbina"),
         "turbina_estator": ("fluido_turbina", "fluido_estator"),
         "estator_bomba": ("fluido_estator", "fluido_bomba")}   # mesma tabela de ``analise_conversor.euler_csv``
ZONAS = list(ac.LADOS)
COORD = ["x-coordinate", "y-coordinate", "z-coordinate"]


def ler_lado(p, par, peca):
    """Lê o CSV de Euler e o CSV antigo (pressão) de um lado e junta por coordenada.

    Args:
        p (Path): Pasta do ponto (``SRxxx``).
        par (str): Nome do par de interface, p. ex. ``bomba_turbina``.
        peca (str): Peça do lado (``bomba``, ``turbina`` ou ``estator``).

    Returns:
        pd.DataFrame: Dados do lado com a coluna ``pressure`` acrescentada.

    Raises:
        FileNotFoundError: Se faltar um dos dois CSVs.
        pd.errors.EmptyDataError: Se um dos CSVs estiver vazio.
        ValueError: Se os CSVs tiverem números de linhas diferentes (truncado), a junção por coordenada
            não for 1:1 ou houver NaN nas colunas usadas.
    """
    nome = f"if_{par}_lado_{peca}"
    fe, fo = p / "dados" / f"{p.name}_euler_{nome}.csv", p / "dados" / f"{p.name}_{nome}.csv"
    for f in (fe, fo):
        if not f.exists():
            raise FileNotFoundError(f)
    e = pd.read_csv(fe, skipinitialspace=True)
    o = pd.read_csv(fo, skipinitialspace=True)
    e.columns, o.columns = e.columns.str.strip(), o.columns.str.strip()
    if len(e) != len(o) or len(e) == 0:
        raise ValueError(f"{nome}: CSVs com tamanhos diferentes ou truncados ({len(e)} x {len(o)} linhas)")
    m = e.merge(o[COORD + ["pressure"]], on=COORD, how="inner", validate="one_to_one")
    if len(m) != len(e):
        raise ValueError(f"{nome}: junção por coordenada não é 1:1 ({len(m)} de {len(e)} linhas casaram)")
    if m.isna().any().any():
        raise ValueError(f"{nome}: há NaN nas colunas usadas (CSV truncado?)")
    return m


def tracao_axial(d):
    """Torque axial da tração viscosa+turbulenta τ·n dA nas faces de um lado.

    Args:
        d (pd.DataFrame): Dados do lado (gradientes ``d<a>-velocity-d<b>`` = ∂u_a/∂x_b, ``viscosity-eff``, área).

    Returns:
        float: D = Σ (y t_z − z t_y), com t_i dA = μ_eff Σ_j (∂u_i/∂x_j + ∂u_j/∂x_i) dA_j; área para fora da zona.
    """
    ax = ("x", "y", "z")
    A = {a: d[f"{a}-face-area"].to_numpy() for a in ax}
    t = {}
    for i in ax:
        t[i] = sum((d[f"d{i}-velocity-d{j}"].to_numpy() + d[f"d{j}-velocity-d{i}"].to_numpy()) * A[j] for j in ax)
        t[i] = d["viscosity-eff"].to_numpy() * t[i]
    y, z = d["y-coordinate"].to_numpy(), d["z-coordinate"].to_numpy()
    return float(np.sum(y * t["z"] - z * t["y"]))


def vorticidade(d, sinal, w_zona):
    """Diagnóstico: os gradientes exportados são da velocidade absoluta ou da relativa?

    Compara a vorticidade axial média (ponderada pela área) dos gradientes, ω_x = ∂w/∂y − ∂v/∂z, com a do perfil
    radial de v_θ absoluta, (1/r) d(r v_θ)/dr (médias por faixa de raio; vale para escoamento quase axissimétrico),
    e com 2ω da zona. Se ω_grad − ω_perfil ≈ 0 os gradientes são da absoluta; se ≈ −2ω, da relativa.

    Args:
        d (pd.DataFrame): Dados do lado.
        sinal (int): Sentido do eixo de ω_b (+1/−1).
        w_zona (float): Velocidade angular da zona (módulo) [rad/s].

    Returns:
        dict[str, float]: ``vort_grad`` (média ponderada), ``vort_perfil`` (mesma ponderação), ``dif_sobre_2w``
        ((vort_grad − vort_perfil)/(2 sinal ω) ou NaN se ω = 0) e ``vort_sobre_2w`` (vort_grad/(2 sinal ω)).
    """
    y, z = d["y-coordinate"].to_numpy(), d["z-coordinate"].to_numpy()
    v, w = d["y-velocity"].to_numpy(), d["z-velocity"].to_numpy()
    A = np.sqrt(d["x-face-area"] ** 2 + d["y-face-area"] ** 2 + d["z-face-area"] ** 2).to_numpy()
    r = np.hypot(y, z)
    vt = (y * w - z * v) / np.maximum(r, 1e-12)
    om = (d["dz-velocity-dy"] - d["dy-velocity-dz"]).to_numpy()
    nb = 16
    bins = np.linspace(r.min(), r.max() * (1 + 1e-9), nb + 1)
    k = np.digitize(r, bins) - 1
    ok = [i for i in range(nb) if (k == i).sum() >= 5]
    if len(ok) < 4:
        return {"vort_grad": np.nan, "vort_perfil": np.nan, "dif_sobre_2w": np.nan, "vort_sobre_2w": np.nan}
    rc = np.array([np.average(r[k == i], weights=A[k == i]) for i in ok])
    vtc = np.array([np.average(vt[k == i], weights=A[k == i]) for i in ok])
    omc = np.array([np.average(om[k == i], weights=A[k == i]) for i in ok])
    wa = np.array([A[k == i].sum() for i in ok])
    perfil = np.gradient(rc * vtc, rc) / rc
    vg, vp = float(np.average(omc, weights=wa)), float(np.average(perfil, weights=wa))
    n2 = 2 * sinal * w_zona
    return {"vort_grad": vg, "vort_perfil": vp,
            "dif_sobre_2w": (vg - vp) / n2 if w_zona > 0 else np.nan,
            "vort_sobre_2w": vg / n2 if w_zona > 0 else np.nan}


def termos_lado(d, w):
    """Fluxos de um lado de interface: massa, momento angular convectivo, tração difusiva e pressão.

    Args:
        d (pd.DataFrame): Dados do lado (valores da célula adjacente; área orientada para fora da zona).
        w (float): ω da zona do lado no eixo x, já com o sinal (sinal·ω_zona) [rad/s]; 0 no estator.

    Returns:
        dict[str, float]: ``mdot`` (Σ ρ u·A), ``H`` (Σ ρ u·A (y w − z v)), ``mdot_rel`` e ``H_rel`` (mesmos, com o
        fluxo relativo ρ (u − ω×r)·A, ω×r = ω (0, −z, y)), ``D`` (ver ``tracao_axial``), ``P`` (−Σ p (y A_z − z A_y)),
        ``I`` (Σ (y A_z − z A_y), geometria, para a sensibilidade ao nível de pressão) e ``dA`` (área total).
    """
    A = d[["x-face-area", "y-face-area", "z-face-area"]].to_numpy()
    u = d[["x-velocity", "y-velocity", "z-velocity"]].to_numpy()
    y, z = d["y-coordinate"].to_numpy(), d["z-coordinate"].to_numpy()
    rho = d["density"].to_numpy()
    vA = np.einsum("ij,ij->i", u, A)
    gA = y * A[:, 2] - z * A[:, 1]                      # (ω×r)·A / ω
    vA_rel = vA - w * gA
    lx = y * u[:, 2] - z * u[:, 1]
    return {"mdot": float(np.sum(rho * vA)), "H": float(np.sum(rho * vA * lx)),
            "mdot_rel": float(np.sum(rho * vA_rel)), "H_rel": float(np.sum(rho * vA_rel * lx)),
            "D": tracao_axial(d), "P": float(-np.sum(d["pressure"].to_numpy() * gA)), "I": float(np.sum(gA)),
            "dA": float(np.sum(np.linalg.norm(A, axis=1)))}


def processar_ponto(p, r, diag_vort=False):
    """Calcula os balanços de Euler de um ponto.

    Args:
        p (Path): Pasta do ponto.
        r (dict): ``resultado.json`` do ponto.
        diag_vort (bool): Se True, calcula também o diagnóstico (inconclusivo) de vorticidade.

    Returns:
        tuple[list[dict], list[dict], list[dict]]: Linhas por zona, por interface e de diagnóstico de vorticidade
        (vazia se ``diag_vort`` for False); listas vazias, com aviso, se faltarem CSVs ou estiverem truncados.

    Raises:
        KeyError: Se ``resultado.json`` não tiver ``omega_b``/``omega_t``.
    """
    try:
        lados = {(par, peca): ler_lado(p, par, peca) for par, zz in PARES.items() for peca in (z.split("_")[1] for z in zz)}
    except (FileNotFoundError, pd.errors.EmptyDataError, pd.errors.ParserError, ValueError) as ex:
        print(f"AVISO {p.name}: ponto ignorado ({type(ex).__name__}: {ex})")
        return [], [], []
    sinal = r.get("sinal", 1)
    desemp = ac.desempenho(r)
    omega = {"bomba": r["omega_b"], "turbina": r["omega_t"], "estator": 0.0}
    rep = r.get("media_janela") or r.get("reports_inst") or {}
    chaves = ("H", "H_rel", "D", "P", "I")
    sai = {z: {k: 0.0 for k in chaves} for z in ZONAS}
    prop = {z: {k: 0.0 for k in chaves} for z in ZONAS}     # lado próprio de cada zona (sem média)
    ifs, vort = [], []
    for par, (za, zb) in PARES.items():
        pa, pb = za.split("_")[1], zb.split("_")[1]
        ta = termos_lado(lados[(par, pa)], sinal * omega[pa])
        tb = termos_lado(lados[(par, pb)], sinal * omega[pb])
        linha = {"SR": r["sr"], "interface": par, "mdot_a": ta["mdot"], "mdot_b": tb["mdot"],
                 "mdot_rel_a": ta["mdot_rel"], "mdot_rel_b": tb["mdot_rel"]}
        mf = rep.get(f"m_{par}_{pa}")        # Fluent: positivo ENTRANDO na zona; para fora = −m
        if mf:
            mm = (ta["mdot"] - tb["mdot"]) / 2
            dev = (mm + mf) / abs(mf)
            linha["mdot_fluent"] = -mf
            linha["mdot_dev_pct"] = 100 * dev
            if abs(dev) > 0.05:
                print(f"AVISO {p.name} {par}: ṁ médio dos lados ({mm:.1f} kg/s) difere "
                      f"{100 * dev:+.1f} % do ṁ do Fluent ({-mf:.1f} kg/s)")
        for k in chaves:
            F = (ta[k] - tb[k]) / 2
            sai[za][k] += F
            sai[zb][k] -= F
            prop[za][k] += ta[k]
            prop[zb][k] += tb[k]
            linha.update({f"{k}_{pa}": ta[k], f"{k}_{pb}": tb[k], f"{k}_soma": ta[k] + tb[k], f"{k}_media": F})
        ifs.append(linha)
        if diag_vort:
            for peca in (pa, pb):
                vort.append({"SR": r["sr"], "interface": par, "lado": peca,
                             **vorticidade(lados[(par, peca)], sinal, omega[peca])})
    zonas = []
    for z in ZONAS:
        # D e P entram com sinal oposto a H: T_euler = −sinal·(H_saída − D − P)
        conv, conv_rel = -sinal * sai[z]["H"], -sinal * sai[z]["H_rel"]
        dif, pres, pres_p = sinal * sai[z]["D"], sinal * sai[z]["P"], sinal * prop[z]["P"]
        par_ = desemp.get("parede_" + z)
        lin = {"SR": r["sr"], "zona": z, "T_conv": conv, "T_dif": dif, "T_pres": pres,
               "T_conv_dif": conv + dif, "T_soma": conv + dif + pres, "T_parede": par_,
               "T_conv_rel": conv_rel, "T_soma_rel": conv_rel + dif + pres,
               "T_pres_prop": pres_p, "T_soma_prop": conv + dif + pres_p, "T_soma_rel_prop": conv_rel + dif + pres_p,
               "dTdp0_parede_Nm_bar": -sinal * prop[z]["I"] * 1e5, "dTdp0_media_Nm_bar": -sinal * sai[z]["I"] * 1e5}
        lin["dTdp0_dif_Nm_bar"] = lin["dTdp0_media_Nm_bar"] - lin["dTdp0_parede_Nm_bar"]
        if par_:
            for nome, t in (("conv", conv), ("conv_dif", conv + dif), ("soma", conv + dif + pres),
                            ("conv_rel", conv_rel), ("soma_rel", conv_rel + dif + pres),
                            ("soma_prop", conv + dif + pres_p), ("soma_rel_prop", conv_rel + dif + pres_p)):
                lin["razao_" + nome] = t / par_
            lin["residuo_pct"] = 100 * (par_ - (conv + dif + pres)) / abs(par_)
        zonas.append(lin)
    return zonas, ifs, vort


def verificar_contra_euler_csv(p, sinal, zonas):
    """Confere que o termo convectivo reproduz ``analise_conversor.euler_csv`` (mesmas convenções).

    Args:
        p (Path): Pasta do ponto.
        sinal (int): Sentido do eixo de ω_b.
        zonas (list[dict]): Linhas por zona de ``processar_ponto``.

    Returns:
        float: Maior diferença relativa entre os dois cálculos do convectivo.
    """
    ref = ac.euler_csv(p, sinal)
    return max(abs(l["T_conv"] - ref[l["zona"]]) / max(abs(ref[l["zona"]]), 1e-9) for l in zonas)


def grafico(tab, saida):
    """Razão Euler/parede contra SR por zona, com e sem o difusivo e com a variante de fluxo relativo.

    Args:
        tab (pd.DataFrame): Tabela por zona (colunas ``SR``, ``zona``, ``razao_*``).
        saida (Path): Arquivo PNG.
    """
    fig, ax = plt.subplots(1, 3, figsize=(13, 4.2), sharex=True)
    for a, z in zip(ax, ZONAS):
        t = tab[tab["zona"] == z].sort_values("SR")
        a.plot(t["SR"], t["razao_conv"], "o-", color="#888888", label="convectivo (ρ u·A)")
        a.plot(t["SR"], t["razao_conv_dif"], "s-", color="#d62728", label="convectivo + difusivo")
        a.plot(t["SR"], t["razao_soma"], "^--", color="#1f77b4", label="conv + dif + p (média dos lados)")
        a.plot(t["SR"], t["razao_conv_rel"], "o:", color="#2ca02c", label="conv. fluxo relativo")
        a.plot(t["SR"], t["razao_soma_rel_prop"], "v:", color="#9467bd", label="rel + dif + p (lado próprio)")
        a.axhline(1, color="k", lw=0.8, ls=":")
        a.set(xlabel="SR = ω_t/ω_b", title=z.replace("fluido_", "Zona "), ylabel="T_Euler / T_parede")
    ax[0].legend(fontsize=6)
    fig.suptitle("Razão entre o torque por momento angular nas interfaces e o torque nas paredes "
                 "(variações de 1–3 % estão dentro da incerteza do pós-processamento por centro de célula)", fontsize=9)
    fig.tight_layout()
    fig.savefig(saida)
    plt.close(fig)


def main(srs=None, diag_vort=False):
    """Processa os pontos existentes, grava os CSVs, o gráfico e imprime a tabela.

    Args:
        srs (list[float] | None): SR (fração) a processar; None = todos com ``resultado.json`` e CSVs.
        diag_vort (bool): Liga o diagnóstico (inconclusivo) de vorticidade.
    """
    zonas, ifs, vorts = [], [], []
    for p, r in ac.pontos():
        if srs is not None and not any(abs(r["sr"] - s) < 1e-6 for s in srs):
            continue
        z, i, v = processar_ponto(p, r, diag_vort)
        if not z:
            print(f"{p.name}: sem CSVs de interface utilizáveis, ignorado")
            continue
        dif = verificar_contra_euler_csv(p, r.get("sinal", 1), z)
        print(f"{p.name}: convectivo confere com euler_csv (dif. relativa máx. {dif:.1e})")
        zonas += z
        ifs += i
        vorts += v
    if not zonas:
        print("nenhum ponto processado")
        return
    tz, ti = pd.DataFrame(zonas), pd.DataFrame(ifs)
    (DIR / "dados").mkdir(exist_ok=True)
    tz.to_csv(DIR / "dados" / "euler_interfaces.csv", sep=";", index=False)
    ti.to_csv(DIR / "dados" / "euler_interfaces_if.csv", sep=";", index=False)
    grafico(tz, DIR / "graficos" / "g7_euler_interfaces.png")
    pd.set_option("display.width", 250, "display.max_columns", 40, "display.float_format", "{:.4g}".format)
    print("\nPor zona (torques em N·m; razão = T_Euler/T_parede; *_rel = fluxo relativo; *_prop = pressão do lado próprio):")
    print(tz[["SR", "zona", "T_conv", "T_dif", "T_pres", "T_pres_prop", "T_parede", "razao_conv", "razao_conv_dif",
              "razao_soma", "razao_soma_prop", "razao_conv_rel", "razao_soma_rel", "razao_soma_rel_prop"]].to_string(index=False))
    print("\nSensibilidade ao nível de pressão p0 [N·m/bar] (parede = lado próprio; média = termo usado em T_pres):")
    print(tz[["SR", "zona", "dTdp0_parede_Nm_bar", "dTdp0_media_Nm_bar", "dTdp0_dif_Nm_bar"]].to_string(index=False))
    print("\nPor interface (X_a + X_b = 0 se os dois lados forem conservativos; media = fluxo usado):")
    cols = ["SR", "interface"] + [c for c in ti.columns if c.endswith(("_soma", "_media")) and not c.startswith("I_")]
    print(ti[cols].to_string(index=False))
    print("\nNota: difusivo (~1 %) e pressão (~1,5 %, faces inclinadas junto às arestas das pás na if_bt) são estimativas da")
    print("mesma ordem da incerteza do pós-processamento por centro de célula (vazão de um lado até ±18 % da do solver;")
    print("reconstrução célula→face muda a razão ±3–4 %) e do desequilíbrio global de torque (até 1,2 %); não se conservam")
    print("entre os lados. O resíduo de 1–3 % não deve ser lido como melhora ou piora do fechamento.")
    if diag_vort:
        print("\nDiagnóstico de vorticidade axial dos gradientes: INCONCLUSIVO (perfil radial grosseiro, Δr pequeno na te).")
        print(pd.DataFrame(vorts).to_string(index=False))


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    main([float(a) / 100 for a in args] or None, "--diag-vorticidade" in sys.argv)
