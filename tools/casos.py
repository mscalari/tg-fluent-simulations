"""Tabela das execuções que já podem rodar e textos que explicam os parâmetros de entrada.

Fonte: guia validado (seções 2–5), registro_execucoes_v2.csv e Roteiro PDF (fichas A–G).
Quando o guia e o Roteiro divergem, vale o guia.
"""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"

BLOCOS = {
    "A": "A_CP_periodico", "B": "B_CP_entrada_saida", "C": "C_LDC_esquemas", "D": "D_LDC_algoritmos",
    "E": "E_degrau_laminar", "F": "F_CP_entrada_uniforme", "G": "G_interp_pressao", "H": "H_turbulencia",
}

# Física e contorno de cada tipo de caso. Os valores numéricos são os únicos usados pelo runner.
TIPOS = {
    "cp_per": dict(
        rho=1.0, mu=0.1, u_parede=2.0, dpdx=4.0, iter_max=3000,
        geometria="Canal plano 5 × 1 m, desenhado em y ∈ [−1, 0]. `wall_top` (móvel) fica em y = −1 e "
                  "`wall_bottom` (parada) em y = 0; por isso a análise usa d = −y (distância à parede parada).",
        fluido="ρ = 1 kg/m³, μ = 0,1 Pa·s. Re = ρ·U·H/μ = 1·2·1/0,1 = 20. P = G·H²/(μ·U) = 10.",
        contorno="`periodic_a` (x = 0) e `periodic_b` (x = 5) viram um par periódico (make-periodic, deltas 5 e 0). "
                 "Gradiente de pressão imposto de 4 Pa/m com direção de escoamento (−1, 0) — com (1, 0) o ū sai "
                 "negativo (guia, seção 3). `wall_top` é parede móvel absoluta a 2 m/s na direção (1, 0); "
                 "`wall_bottom` é parede parada. Sem escorregamento nas duas.",
        esperado="Escoamento de Couette-Poiseuille desenvolvido: u(d) = 2d + 20d(1−d), ū = 4,3333 m/s. "
                 "Balanço de forças: G·L·H = 20 N/m = |Fx parede parada| (11 N/m) + |Fx parede móvel| (9 N/m)."),
    "cp_io": dict(
        rho=1.0, mu=0.1, u_parede=2.0, iter_max=3000,
        geometria="Mesmo canal do CP periódico (5 × 1 m, y ∈ [−1, 0], parede móvel em y = −1), malha M2_io "
                  "(50 × 20) com as bordas renomeadas: `inlet` em x = 0 e `outlet` em x = 5.",
        fluido="ρ = 1 kg/m³, μ = 0,1 Pa·s, Re = 20.",
        contorno="`inlet`: Velocity Inlet, módulo normal à fronteira. `outlet`: Pressure Outlet com 0 Pa manométrico. "
                 "`wall_top` móvel absoluta a 2 m/s em (1, 0); `wall_bottom` parada. Não há gradiente imposto: "
                 "a queda de pressão aparece sozinha.",
        esperado="Δp = p_entrada − p_saída ≈ G·L = 4 Pa/m · 5 m = 20 Pa; em x = 2,5 o perfil deve ser o de Couette-Poiseuille."),
    "ldc": dict(
        rho=1.0, mu=0.01, u_parede=1.0, iter_max=3000,
        geometria="Cavidade quadrada 1 × 1 m, 0 ≤ x, y ≤ 1. `lid` (tampa) em y = 1; `wall_left`, `wall_right` e "
                  "`wall_bottom` paradas.",
        fluido="ρ = 1 kg/m³, μ = 0,01 Pa·s. Re = ρ·U·L/μ = 1·1·1/0,01 = 100 (caso clássico de Ghia et al., 1982). "
               "Valores REC/FIG do Roteiro; confirmar no código MATLAB.",
        contorno="`lid`: parede móvel absoluta a 1 m/s na direção (1, 0). Demais paredes paradas, sem escorregamento. "
                 "Não há entrada nem saída: a pressão fica definida a menos de uma constante.",
        esperado="Um vórtice principal com centro perto de (0,617; 0,734) m (Roteiro; conferir na tabela de Ghia) e "
                 "vórtices secundários pequenos nos cantos inferiores."),
    "cavr": dict(
        rho=1.0, mu=0.01, u_parede=1.0, omega=1.0, eixo=[0.5, 0.5], iter_max=3000, anim_cada=20,
        geometria="Mesma cavidade 1 × 1 m do lid-driven (`lid` em y = 1), resolvida em 2D planar: o Fluent 2024 R1 "
                  "aceita Frame Motion em 2D (testado em `tools/teste_frame_motion.py`), então não foi preciso "
                  "a versão 3D de uma célula com faces Symmetry prevista no Adendo.",
        fluido="ρ = 1 kg/m³, μ = 0,01 Pa·s (Re = 100, igual ao lid-driven).",
        contorno="Zona fluida com Frame Motion: Ω = 1 rad/s em torno de z, eixo "
                 "no centro da cavidade (0,5; 0,5) — o Adendo não define a origem; escolha registrada como hipótese. "
                 "Ro = U/(ΩL) ≈ 1: Coriolis (≈ 2ΩU = 2 m/s²) comparável à inércia. `lid`: Moving Wall 1 m/s em (1, 0), "
                 "relativa à zona adjacente. Demais paredes: Stationary Wall, que no Fluent fica parada em relação à "
                 "zona (gira com o referencial). Velocidades analisadas no referencial rotativo (relative-x/y-velocity).",
        esperado="Vórtice da cavidade deformado pela rotação. Hipóteses G2–G5 do Adendo: Standard com sobrevalores, "
                 "Second Order melhor que Linear e Standard, BFW e PRESTO! ≥ 5 % melhores que Second Order, "
                 "diferenças caindo com o refino. Referência: R26 (N128, PRESTO!), checada pelo R27."),
    "bfs": dict(
        rho=1.0, mu=0.01, iter_max=3000,
        geometria="Degrau (backward-facing step): canal de entrada 0 ≤ x ≤ 1, 0,5 ≤ y ≤ 1,5; degrau de altura "
                  "h = 0,5 m em x = 1; trecho a jusante 1 ≤ x ≤ x_saída, 0 ≤ y ≤ 1,5. Geometria FIG (Roteiro, ficha E): "
                  "confirmar no código MATLAB.",
        fluido="ρ = 1 kg/m³, μ = 0,01 Pa·s (REC, a substituir pelo valor do código). Re_H = ρ·ū·H_in/μ = 1·1·1/0,01 = 100 "
               "(H_in = 1 m, altura do canal de entrada); Re_h = 50 com a altura do degrau.",
        contorno="`inlet` (x = 0): Velocity Inlet com o perfil parabólico de `refs/bfs_inlet.prof` (ū = 1 m/s, "
                 "u_max = 1,5 m/s). `outlet`: Pressure Outlet 0 Pa. `wall_top`, `wall_in_bottom`, `wall_step` e "
                 "`wall_down_bottom`: paredes paradas.",
        esperado="Bolha de recirculação atrás do degrau; x_r = comprimento de recolamento medido a partir do degrau "
                 "(x = 1) pela troca de sinal de u na primeira linha de células acima de `wall_down_bottom`. "
                 "O Roteiro não traz o x_r esperado."),
}

# Padrões comuns: SIMPLE, LSQ cell based, pressão de 2ª ordem, momento SOU, URF 0,3/0,7, resíduos 1e-6.
BASE = dict(alg="SIMPLE", mom="second-order-upwind", pres="second-order", urf_p=0.3, urf_mom=0.7)

CASES = {
    "R01": dict(tipo="cp_per", malha="M1", muda="Malha M1 (25 × 10, 250 células), a mais grossa da série de refino."),
    "R02": dict(tipo="cp_per", malha="M2", muda="Nada: é o caso-base validado (M2, 50 × 20)."),
    "R03": dict(tipo="cp_per", malha="M3", muda="Malha M3 (100 × 40, 4000 células)."),
    "R04": dict(tipo="cp_per", malha="M4", muda="Malha M4 (200 × 80, 16 000 células), a mais fina."),
    "R05": dict(tipo="cp_per", malha="M2", mom="first-order-upwind",
                muda="Só o esquema de momento: First Order Upwind em vez de SOU. Em escoamento desenvolvido não há "
                     "convecção, então o perfil deve ser igual ao do R02 (|Δu| < 1e-4 m/s)."),
    "R06": dict(tipo="cp_io", malha="M2_io", perfil="cp_inlet_ymenos.prof",
                muda="Periodicidade trocada por entrada e saída: perfil analítico de Couette-Poiseuille na entrada "
                     "e 0 Pa na saída. Mede-se a queda de pressão. O `refs/cp_inlet.prof` original está em y ∈ [0, 1] "
                     "(u = 0 em y = 0, u = 2 em y = 1), ao contrário do que diz o guia; na malha (y ∈ [−1, 0]) o Fluent "
                     "usava o ponto mais próximo (u = 0) em toda a entrada. Usa-se a cópia espelhada "
                     "`refs/cp_inlet_ymenos.prof` (y → −y), conferida ponto a ponto contra u(d) = 2d + 20d(1−d)."),
    "R07": dict(tipo="ldc", malha="N040", mom="first-order-upwind", muda="Esquema de momento First Order Upwind."),
    "R08": dict(tipo="ldc", malha="N040", muda="Nada: base da cavidade (SIMPLE + SOU, N040 40 × 40)."),
    "R09": dict(tipo="ldc", malha="N040", mom="quick", muda="Esquema de momento QUICK."),
    "R10": dict(tipo="ldc", malha="N040", mom="third-order-muscl", muda="Esquema de momento MUSCL (3ª ordem)."),
    "R11": dict(tipo="ldc", malha="N128", iter_max=8000,
                muda="Malha fina N128 (128 × 128, 16 384 células), referência para comparar com Ghia (1982). "
                     "Teto de 8000 iterações."),
    "R12": dict(tipo="ldc", malha="N040", alg="SIMPLEC", urf_p=1.0,
                muda="Só o acoplamento pressão-velocidade: SIMPLEC, com URF de pressão 1,0."),
    "R13": dict(tipo="ldc", malha="N040", alg="PISO", muda="Só o acoplamento pressão-velocidade: PISO."),
    "R14": dict(tipo="ldc", malha="N040", alg="Coupled",
                muda="Só o acoplamento pressão-velocidade: Coupled (os URFs do SIMPLE não se aplicam; ficam os "
                     "padrões do Coupled)."),
    "R15": dict(tipo="bfs", malha="DEG15", perfil="bfs_inlet.prof",
                muda="Base do degrau: saída em x = 2,5 m (malha DEG15, Δ = 0,025 m, 5200 células)."),
    "R16": dict(tipo="bfs", malha="DEG16", perfil="bfs_inlet.prof",
                muda="Só o comprimento da saída: x = 11 m (malha DEG16, 25 600 células), para ver a influência da saída."),
    "R17": dict(tipo="cp_io", malha="M2_io", uin=4.3333,
                muda="Entrada uniforme ū = 4,3333 m/s em vez do perfil (região de desenvolvimento). Opcional; o Roteiro "
                     "não define a saída: usada a mesma do R06 (0 Pa) — hipótese."),
    "R18": dict(tipo="ldc", malha="N040", pres="presto!", bloco="G",
                muda="Interpolação de pressão PRESTO! (controle Ω = 0 do bloco G, hipótese G1). "
                     "G1 avaliado pelo mesmo critério da referência do bloco G, sem a fileira de células junto à tampa; R18 (1,10 %) e R20 (2,45 %) ficam acima de 1 %, então G1 não se confirma (o RMS fica abaixo de 1 %: 0,55 % e 0,21 %)."),
    "R19": dict(tipo="ldc", malha="N040", bloco="G",
                muda="Interpolação de pressão Second Order (controle G1). Rodado em 2D (N040): com Ω = 0 a versão 3D "
                     "de uma célula é equivalente (pré-requisitos, decisão 3). Configuração igual à do R08."),
    "R20": dict(tipo="ldc", malha="N040", pres="standard", bloco="G",
                muda="Interpolação de pressão Standard (controle G1). Rodado em 2D (N040), Ω = 0. "
                     "G1 avaliado pelo mesmo critério da referência do bloco G, sem a fileira de células junto à tampa; R18 (1,10 %) e R20 (2,45 %) ficam acima de 1 %, então G1 não se confirma (o RMS fica abaixo de 1 %: 0,55 % e 0,21 %)."),
    "R21": dict(tipo="cavr", malha="N040", pres="linear", bloco="G", muda="Interpolação de pressão Linear (G3)."),
    "R22": dict(tipo="cavr", malha="N040", pres="standard", bloco="G",
                muda="Interpolação de pressão Standard (G2: espera-se sobrevalor)."),
    "R23": dict(tipo="cavr", malha="N040", bloco="G", muda="Interpolação de pressão Second Order (base do bloco G)."),
    "R24": dict(tipo="cavr", malha="N040", pres="body-force-weighted", bloco="G",
                muda="Interpolação de pressão Body Force Weighted (G4)."),
    "R25": dict(tipo="cavr", malha="N040", pres="presto!", bloco="G", muda="Interpolação de pressão PRESTO! (G4)."),
    "R26": dict(tipo="cavr", malha="N128", anim_cada=100, pres="presto!", bloco="G", iter_max=8000,
                muda="Malha fina N128 com PRESTO!: referência do bloco G. Teto de 8000 iterações."),
    "R27": dict(tipo="cavr", malha="N128", anim_cada=100, bloco="G", iter_max=8000,
                muda="Malha fina N128 com Second Order: checagem da referência (|R26 − R27| < 1 % de u_max). "
                     "A referência é aceita pelo critério sem a fileira de células junto à "
                     "tampa (0,92 %), porque o máximo de 14,1 % vem dos cantos singulares da tampa, como no G1."),
}

BLOCO_ID = {"R01": "A", "R02": "A", "R03": "A", "R04": "A", "R05": "A", "R06": "B", "R07": "C", "R08": "C",
            "R09": "C", "R10": "C", "R11": "C", "R12": "D", "R13": "D", "R14": "D", "R15": "E", "R16": "E",
            "R17": "F"}


def caso(cid):
    """Junta a linha de ``CASES`` com os padrões, a física do tipo e o nome-base do registro.

    Args:
        cid (str): ID da execução, por exemplo ``"R02"``.

    Returns:
        dict: Parâmetros completos (``id``, ``nome``, ``bloco``, ``dir``, física do tipo e métodos).
    """
    import csv
    with open(ROOT / "00_admin" / "registro_execucoes_v2.csv", encoding="utf-8") as f:
        reg = {r["ID"]: r for r in csv.DictReader(f, delimiter=";")}
    c = {**TIPOS[CASES[cid]["tipo"]], **BASE, **CASES[cid]}
    c["id"] = cid
    c["nome"] = reg[cid]["Arquivo_base"]
    c["bloco"] = c.get("bloco") or BLOCO_ID[cid]
    c["dir"] = RESULTS / BLOCOS[c["bloco"]] / c["nome"]
    return c


if __name__ == "__main__":
    for k in CASES:
        c = caso(k)
        assert c["dir"].parent.name.startswith(c["bloco"]), c
        assert (ROOT / "mesh" / f"{c['malha']}.msh").exists(), c["malha"]
    print("casos ok:", len(CASES))
