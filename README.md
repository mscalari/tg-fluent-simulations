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
- **degrau turbulento** de Driver e Seegmiller (1985), com modelos RANS.

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

Simulação 3D do conversor (bomba, turbina e estator; 1,9 M células; malha `mesh/CT_V1_3zonas.msh.h5`,
fora do git) com referenciais múltiplos (MRF), k-ω SST com correção de curvatura, PRESTO! e SIMPLE.
A bomba gira a 2000 rpm em torno de +x e a turbina a SR·ω_b; o estator fica parado.

```bash
export PYTHONIOENCODING=utf-8
.venv241/Scripts/python.exe tools/conversor_mrf.py sentido 500     # SR = 0 nos dois sentidos de rotação
.venv241/Scripts/python.exe tools/conversor_mrf.py varredura 1     # SR = 0 … 0,9, cada ponto parte do anterior
.venv241/Scripts/python.exe tools/conversor_mrf.py refazer_figuras # relê o .cas/.dat de cada ponto e refaz as imagens
.venv241/Scripts/python.exe tools/analise_conversor.py             # gráficos, curvas TR/η/K, Euler e nota do Obsidian
```

Saídas em `results/conversor/CT_V1_MRF/SRxxx/{fluent,dados,graficos}`: cerca de 60 imagens do Fluent por
ponto (plano meridional, cortes axiais, raios constantes, pás, interfaces, linhas de corrente, critério Q)
mais monitores, resíduos, y⁺ e perfis nas interfaces (matplotlib); curvas globais em `CT_V1_MRF/graficos`.
O Coupled não cabe nos 16 GB de RAM desta malha (swap); por isso SIMPLE. Os torques oscilam ~1,5 % sem
deriva: o ponto para quando a média de 200 iterações estabiliza, e os valores são essas médias.
Na varredura as imagens do Fluent saem em branco (só a legenda): a janela gráfica fica presa aos plots
de monitor. Por isso as figuras são refeitas depois, relendo o caso de cada ponto numa sessão nova.

## Documentação

| Arquivo | Conteúdo |
|---|---|
| [docs/visao_geral.md](docs/visao_geral.md) | O que o projeto valida, blocos e execuções, organização das pastas |
| [docs/como_rodar.md](docs/como_rodar.md) | Ambiente, comandos, arquivos de saída, como adicionar uma execução, armadilhas |
| [docs/codigo.md](docs/codigo.md) | Cada módulo e cada função de `tools/`, em linguagem simples |
