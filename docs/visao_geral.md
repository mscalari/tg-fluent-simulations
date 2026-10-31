# Visão geral

## O que é o projeto

O projeto de graduação visa realizar um estudo acerca de métodos aplicados no software Ansys Fluent, fazendo também comparações com códigos MATLAB de desenvolvimento próprio do autor. Este repositório cuida da parte do Fluent: cada **execução** (uma simulação,
identificada por `R01` … `R43`) é montada, resolvida e exportada por um script Python que controla
o Fluent pelo **PyFluent** (pacote `ansys-fluent-core`). 

Para quem nunca usou o Fluent: o Fluent é um programa de CFD (dinâmica dos fluidos computacional).
Ele lê uma **malha** (o domínio dividido em células), recebe a física (fluido, paredes, entradas e
saídas) e os **métodos numéricos**, e resolve as equações iterativamente até os **resíduos** (o
quanto as equações ainda não estão satisfeitas) ficarem abaixo de um critério. Normalmente isso é
feito clicando na interface; aqui, o script faz os mesmos passos por código.

A geometria e as malhas **não** são feitas por código: foram geradas manualmente uma vez no Workbench e estão em
`mesh/`.

## Os casos

| Caso | O que é | Com o que se compara |
|---|---|---|
| **Couette-Poiseuille periódico** (CP) | Canal plano 5 × 1 m entre uma parede parada e outra móvel (2 m/s), com gradiente de pressão imposto de 4 Pa/m. As pontas do canal são ligadas (periodicidade), então o escoamento é totalmente desenvolvido. | Solução analítica u(d) = 2d + 20d(1−d), ū = 4,3333 m/s |
| **CP entrada/saída** | O mesmo canal, mas com entrada (perfil de velocidade) e saída (pressão 0). | Queda de pressão esperada Δp = 20 Pa e o mesmo perfil analítico |
| **Cavidade (lid-driven)**, Re = 100 | Caixa quadrada 1 × 1 m fechada, só a tampa se move (1 m/s). Forma um vórtice. | Ghia et al. (1982) e comparações entre execuções |
| **Degrau laminar** (backward-facing step) | Canal que se alarga de repente (degrau de 0,5 m). Atrás do degrau forma-se uma bolha de recirculação. | Comprimento de recolamento x_r |
| **Duto quadrado 3D** | Seção de 1 × 1 m no plano xz, uma célula em y, periódico com gradiente de 4 Pa/m (Re_τ = 150). Modelos de turbulência e, com SST, esquemas de pressão. | Pirozzoli et al. (2018), só ordem de grandeza |
| **Degrau turbulento** | Driver e Seegmiller (1985), Re_H = 36 000, expansão 9/8, entrada lei 1/7. | x_r/H = 6,26 ± 0,10, Cf e Cp de `refs/bfs_driver1985.csv` |

## Blocos e execuções

| Bloco | Pasta em `results/` | Execuções | O que muda |
|---|---|---|---|
| A | `A_CP_periodico` | R01–R05 | Malha (M1–M4) e esquema de momento (R05: FOU) |
| B | `B_CP_entrada_saida` | R06 | Entrada com perfil e saída de pressão |
| C | `C_LDC_esquemas` | R07–R11 | Esquema de momento (FOU, SOU, QUICK, MUSCL) e malha fina N128 |
| D | `D_LDC_algoritmos` | R12–R14 | Acoplamento pressão-velocidade (SIMPLEC, PISO, Coupled) |
| E | `E_degrau_laminar` | R15–R16 | Comprimento da saída do degrau |
| F | `F_CP_entrada_uniforme` | R17 | Entrada uniforme | 
| G | `G_interp_pressao` | R18–R30 | Interpolação da pressão; R18–R20 sem rotação | 
| H | `H_turbulencia` | R31–R43 | Modelos de turbulência: duto 3D (R31–R34) e degrau turbulento (R35–R43; R41 na BFYW, demais na BFY1) |

A lista exata do que roda é o dicionário `CASES` em `tools/casos.py` (R01–R43). 

Convenção de nomes: `<ID>_<CASO>_<MALHA>_<ALG>_<ESQ>[_<VAR>]`, por exemplo `R02_CPP_M2_SIM_SOU`.

## Organização das pastas

```
mesh/            malhas de entrada (.msh): M1–M4, M2_io (CP), N040, N128 (cavidade), DEG15, DEG16 (degrau), DU040 (duto), BFY1, BFYW (degrau turbulento)
refs/            perfis de entrada (cp_inlet.prof original em y ∈ [0, 1]; cp_inlet_ymenos.prof espelhado, usado no R06; bfs_inlet.prof; bfy_inlet.prof, gerado por tools/gerar_bfy_inlet.py), bfs_driver1985.csv, duto_secundario.csv, o CSV de referência do R02 e ghia1982.csv (Ghia et al., Tabelas I e II, Re = 100)
tools/           código (ver docs/codigo.md)
results/<Bloco>/<nome da execução>/
  fluent/        .cas.h5 e .dat.h5 (caso e solução do Fluent)
  dados/         _campo.csv, _parede_<zona>.csv (R28–R43), _rep.out, _residuos.csv, _setup.json, _metricas.json
  graficos/      _campos.png, _convergencia.png, _perfis.png
  logs/          run.log (saída do script) e transcrições .trn do Fluent
docs/            esta documentação
```

O conteúdo de cada arquivo de saída está em [como_rodar.md](como_rodar.md#o-que-cada-arquivo-de-saída-contém).