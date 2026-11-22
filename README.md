# tg-fluent-simulations

Automação, com **PyFluent**, das simulações de validação do Trabalho de Graduação (TG) no
**Ansys Fluent 2024 R1** (Windows). Os scripts montam o caso no Fluent, resolvem, exportam os
resultados e calculam as métricas em Python, sem abrir a interface gráfica.

Casos cobertos:

- **Couette-Poiseuille (CP)** periódico e com entrada/saída, comparado com a solução analítica;
- **cavidade com tampa móvel** (lid-driven cavity), Re = 100: esquemas de momento, algoritmos de
  acoplamento pressão-velocidade e interpolação da pressão;
- **degrau** (backward-facing step) laminar;
- **duto quadrado 3D** periódico, com modelos RANS (escoamento secundário);
- **degrau turbulento** de Driver e Seegmiller (1985), com modelos RANS;
- **conversor de torque CT_V1** (3D, MRF), com varredura da razão de velocidades, sensibilidades e GCI (seção abaixo).

O plano completo tem 43 execuções (R01–R43) em 8 blocos (A–H), todas definidas em `tools/casos.py`
(a cavidade rotativa R21–R27 roda em 2D, com Frame Motion).

## Instalação

Requer Ansys Fluent 2024 R1 com licença válida e Python 3.13. A versão 0.37.2 do PyFluent é a última
compatível com a 2024 R1.

```bash
python -m venv .venv241
.venv241/Scripts/python.exe -m pip install -r requirements.txt
```

## Começo rápido

No Git Bash, na raiz do repositório:

```bash
export PYTHONIOENCODING=utf-8
.venv241/Scripts/python.exe tools/casos.py                     # confere a tabela de casos e as malhas
.venv241/Scripts/python.exe tools/analise_caso.py --selftest   # confere as métricas (não abre o Fluent)
.venv241/Scripts/python.exe tools/gerar_bfy_inlet.py           # gera o perfil de entrada do degrau turbulento
bash tools/rodar_lote.sh R02                                   # roda o R02 no Fluent
MAX=4 bash tools/rodar_fila.sh R35 R36 R37 R38                 # fila com até 4 Fluents simultâneos
.venv241/Scripts/python.exe tools/analise_caso.py R02          # métricas e gráficos
```

Cada Fluent compromete cerca de 4 GB de memória virtual: dimensione `MAX` pela memória da máquina.
Casos turbulentos gravam `.cas/.dat` a cada 2000 iterações. Na fila, `Rxx:c` retoma até os resíduos
convergirem e `Rxx:t` retoma até o monitor estabilizar.

## Resultados dos modelos RANS

**Duto quadrado (R28–R34), Re_τ = 150.** Referência DNS: |V_s|/U_b de 1–2 %.

| Rodada | Modelo | U_b [m/s] | Re_b | \|V_s\|/U_b |
|---|---|---|---|---|
| R28–R30, R32 | SST (três esquemas de pressão) | 14,112 | 4234 | zero numérico |
| R31 | k-ε Standard com EWT | 13,618 | 4085 | zero numérico |
| R33 | RSM linear com EWT | 12,681 | 3804 | 2,28 % |
| R34 | SST com Corner Flow Correction | 13,925 | 4178 | 1,54 % |

**Degrau turbulento (R35–R43), Re_H = 36 000.** Experimento: x_r/H = 6,26 ± 0,10.

| Rodada | Modelo | x_r/H | Erro | Atingiu 1e-6 |
|---|---|---|---|---|
| R35 | k-ε Standard | 5,858 | −6,4 % | sim |
| R36 | k-ε Realizable | 6,812 | +8,8 % | sim |
| R37 | k-ω BSL | 6,201 | −0,9 % | não |
| R38 | k-ω SST | 6,871 | +9,8 % | não |
| R39 | SST com Curvature Correction | 7,038 | +12,4 % | não |
| R40 | RSM linear (τ_w oscila 4 %) | 5,348 | −14,6 % | não |
| R41 | SST, malha y+ ≈ 20 | 6,777 | +8,3 % | sim |
| R42 | k-ε RNG | 6,631 | +5,9 % | sim |
| R43 | GEKO, C_SEP = 2,5 | 7,392 | +18,1 % | não |

O perfil de entrada é uma lei 1/7, não o perfil medido; por isso a comparação entre modelos vale mais
que o valor absoluto de x_r. Nas rodadas que não atingiram 1e-6, a continuidade estagnou; o monitor de
τ_w oscila menos de 0,2 %, exceto no R40. O valor de x_r do R40 é um instante dessa oscilação.

## Conversor de torque CT_V1 (MRF)

Simulação 3D do conversor (bomba, turbina e estator; malha G0 de 1 900 736 células, `mesh/CT_V1_3zonas.msh.h5`,
fora do git) com referenciais múltiplos (MRF). A geometria CT_V1 é dada: o TG valida o Fluent e o setup, não o
conversor. A bomba gira a 2000 rpm em torno de +x (sentido escolhido pelo teste de SR = 0: só nele TR > 1) e a
turbina a SR·ω_b; o estator fica parado. Óleo com ρ = 830 kg/m³ e μ = 0,007 Pa·s.

**Setup:** k-ω SST com correção de curvatura, PRESTO!, SIMPLE (o Coupled não cabe nos 16 GB de RAM desta malha),
segunda ordem, formulação de velocidade relativa, inicialização standard no referencial absoluto. Cada ponto de
operação para quando a média dos torques em janelas de 200 iterações varia menos que 0,5 % (TOL) entre janelas
vizinhas; o valor oficial é essa média, e o desvio dentro da janela é a incerteza iterativa. Os torques saem em N·m:
no report de momento da 24.1, a opção `Drag Force` já é o momento dimensional (`Drag Coefficient` = Cm; M = 0,6125·Cm
com os valores de referência padrão).

```bash
export PYTHONIOENCODING=utf-8
.venv241/Scripts/python.exe tools/conversor_mrf.py setup                     # monta o caso e grava o .cas de setup
.venv241/Scripts/python.exe tools/conversor_mrf.py sentido 500               # SR = 0 nos dois sentidos de rotação
.venv241/Scripts/python.exe tools/conversor_mrf.py varredura 1               # SR = 0 … 0,9, cada ponto parte do anterior
.venv241/Scripts/python.exe tools/conversor_mrf.py varredura 1 0.9 0.95 1.0 --ini <.dat.h5 do ponto anterior>
.venv241/Scripts/python.exe tools/conversor_mrf.py sensibilidade             # sst_sem_cc ke_realizable fou meio_passo (coupled à parte)
.venv241/Scripts/python.exe tools/conversor_mrf.py refazer_figuras 0 0.4     # relê .cas/.dat e refaz as imagens do Fluent
# malhas do GCI: --malha mesh/CT_V1_3zonas_G1.msh.h5 --saida CT_V1_GCI_G1 (vale para qualquer comando)
```

Saídas em `results/conversor/CT_V1_MRF/SRxxx/{fluent,dados,graficos}`: cerca de 60 imagens do Fluent por ponto
(plano meridional, cortes axiais, raios constantes, pás, interfaces, linhas de corrente, critério Q), monitores,
resíduos, y⁺ e perfis nas interfaces (matplotlib). Resultados globais em `CT_V1_MRF/{dados,graficos}`. Na
varredura as imagens do Fluent saem em branco (a janela gráfica fica presa aos plots de monitor), por isso
`refazer_figuras` as refaz depois, numa sessão nova. Os `.cas/.dat` e os exports por face ficam fora do git.

### Análises sem Fluent

Rodam sobre o que já está em `results/conversor/` (todas com `.venv241` e `PYTHONIOENCODING=utf-8`). Detalhes de
entradas e saídas em [docs/codigo.md](docs/codigo.md).

| Script | O que faz | Saídas principais |
|---|---|---|
| `analise_conversor.py` | Curvas TR, η e K, torques, vazão e balanço; nota do Obsidian | `g1`–`g3`, `dados/desempenho.csv` |
| `euler_interfaces_conversor.py` | Balanço de momento angular por zona contra o torque das paredes | `g7`, `dados/euler_interfaces*.csv` |
| `analise_yplus_conversor.py` | y⁺ por faixas e peça, ponderado por área | `g4`, `dados/yplus_faixas.csv` |
| `analise_sensibilidade_conversor.py` | Δ% de TR, T_b, T_t e K de cada variante contra o estol-base | `g6`, `dados/sensibilidade.csv` |
| `comparacao_literatura_conversor.py` | Curvas adimensionais contra a literatura | `g5` |
| `tendencia_derivada_conversor.py` | dTR/dSR, dη/dSR e dK/dSR com incerteza | `g8`, `dados/derivada_tendencia.csv` |
| `tendencia_sinais_conversor.py` | Sinais e monotonicidade por intervalo de SR | `g9`, `dados/sinais_tendencia.csv` |
| `tendencia_colapso_conversor.py` | Colapso adimensional da forma das curvas | `g10`, `dados/colapso_adimensional*.csv` |
| `animacao_meio_passo.py` | GIF da turbina girada de meio passo (10°) | `graficos/animacao_meio_passo.gif` |
| `gci_conversor.py` | GCI (Celik et al., 2008) com G0, G1 e G2 | `CT_V1_GCI/dados/gci.csv`, `graficos/gci_tr.png` |
| `malha_gci_conversor.py` | Gera as malhas G1 e G2 (usa o Fluent Meshing) | `mesh/CT_V1_3zonas_G{1,2}.msh.h5` |
| `imagens_malha_conversor.py` | Imagens das malhas G0, G1 e G2 e montagens lado a lado | `CT_V1_GCI/graficos/malha_*.png` |
| `comparacao_blocos.py` | Figuras-resumo de R01–R43 contra as referências | `results/comparacao_blocos/` |

### Independência de malha (GCI)

Malhas G1 e G2 por semelhança geométrica com a G0, com fatores nominais 1,35 e 1,8225. A calibração pela
superfície da G1 deu r efetivo = 1,308 (842 090 faces contra 1 440 968 da G0). Rodadas em SR = 0 e SR = 0,4 em
`CT_V1_GCI_G1` e `CT_V1_GCI_G2`; `gci_conversor.py` usa p por raiz das equações de Celik, Fs conservador e o erro-padrão
da média por batch means (p aparente acima de 1,1·p formal é limitado ao formal). A G1 e a G2 rodaram com URF
menores (`--urf`). Resultado (22/11/2026): em SR 0 o TR não depende da malha dentro do ruído iterativo (1,065 / 1,062 /
1,069; GCI21 0,49 %); em SR 0,4 os torques absolutos caem ~5 % por nível, fora da faixa assintótica, e o TR da G0
(1,027) difere do da G1/G2 (1,015–1,016), diferença confundida com o protocolo (URF e caminho de SR). Os testes
`CT_V1_GCI_G0_urf_reinicio` e `CT_V1_GCI_G0_urf` separam essas causas. Os `.cas/.dat` antigos da G0 estão em
`D:\TG_backup_fluent\CT_V1_MRF\` (o C: encheu).

### Monitor local

A pasta `monitor/` (ignorada pelo git) tem um painel somente leitura em `http://127.0.0.1:8765/`, com processos do
Fluent, RAM, logs, resultados e filas: `powershell -ExecutionPolicy Bypass -File monitor\abrir_monitor.ps1`.

## Documentação

| Arquivo | Conteúdo |
|---|---|
| [docs/visao_geral.md](docs/visao_geral.md) | O que o projeto valida, blocos e execuções, organização das pastas |
| [docs/como_rodar.md](docs/como_rodar.md) | Ambiente, comandos, arquivos de saída, como adicionar uma execução, armadilhas |
| [docs/codigo.md](docs/codigo.md) | Cada módulo de `tools/` (casos, runner, análises e o conversor), em linguagem simples |
| [docs/handoff_proxima_sessao.md](docs/handoff_proxima_sessao.md) | Estado atual, o que está rodando, próximos passos e pendências |
