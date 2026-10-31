# O código (`tools/`)

Visão do conjunto:

```
casos.py ──► run_caso.py ──► results/<bloco>/<nome>/{fluent,dados,logs}
   │           ▲                         │
   │      rodar_lote.sh, rodar_fila.sh       ▼
   └─────────────────────────► analise_caso.py ──► _metricas.json, graficos/
gerar_bfy_inlet.py: refs/bfy_inlet.prof
```

| Arquivo | Papel | Abre o Fluent? |
|---|---|---|
| `casos.py` | Tabela das execuções e dos parâmetros físicos | não |
| `run_caso.py` | Monta, resolve e exporta uma execução no Fluent | **sim** |
| `rodar_lote.sh` | Roda várias execuções, em sequência ou com `PARALELO=N`, com log | **sim** (via `run_caso.py`) |
| `rodar_fila.sh` | Roda uma fila com teto de Fluents simultâneos e retoma casos de `.cas/.dat` | **sim** (via `run_caso.py`) |
| `analise_caso.py` | Métricas e gráficos | não |
| `gerar_bfy_inlet.py` | Gera `refs/bfy_inlet.prof` (lei 1/7 nas duas paredes, δ = 1,4H, x = −10) | não |

---

## `casos.py`

Define **o que** cada execução é. Não faz nada no Fluent; é importado pelos outros scripts.

### Constantes

| Nome | O que guarda |
|---|---|
| `ROOT` | Raiz do repositório (pasta acima de `tools/`). |
| `RESULTS` | `ROOT/results`. |
| `BLOCOS` | Letra do bloco → nome da pasta em `results/` (ex.: `"A"` → `"A_CP_periodico"`). |
| `TIPOS` | Um dicionário por **tipo de caso** (`cp_per`, `cp_io`, `ldc`, `cavr`, `bfs`, `duto`, `bft`). Cada um tem os números usados pelo runner (`rho`, `mu`, `u_parede`, `dpdx` no CP periódico, `omega` e `eixo` na cavidade rotativa, `anim_cada` quando há animação, `iter_max`; no duto e no degrau turbulento, `dpdy`, `dim`, `I`, `mut_mu`, `u_ref`, `ini` e `paredes`, conforme o tipo) e quatro textos explicativos: `geometria`, `fluido`, `contorno`, `esperado`. |
| `BASE` | Métodos-padrão de todas as execuções: `alg="SIMPLE"`, `mom="second-order-upwind"`, `pres="second-order"`, `urf_p=0.3`, `urf_mom=0.7`. |
| `KES`, `REAL`, `RNG`, `SST`, `BSL`, `GEKO`, `RSM` | Atributos de `setup.models.viscous` de cada modelo de turbulência; o RSM é `rsm-linear-pressure-strain`. |
| `CASES` | Uma linha por execução (R01–R43). Cada linha diz o `tipo`, a `malha`, o texto `muda` (o que a execução testa) e só o que difere do padrão (`mom`, `pres`, `alg`, `urf_p`, `iter_max`, `perfil`, `uin`, `bloco`, `anim_cada`, `turb`, `canto`, `curvatura`, `csep`). |
| `BLOCO_ID` | ID → letra do bloco, para R01–R17. |
| `NOMES` | ID → nome-base da rodada: pasta e prefixo dos arquivos em `results/`. |

### Funções

| Função | O que faz | Entrada | Saída |
|---|---|---|---|
| `caso(cid)` | Junta tudo o que se sabe de uma execução: a física do tipo (`TIPOS`), os padrões (`BASE`) e a linha de `CASES` (nessa ordem; a linha de `CASES` vence). Pega o nome-base em `NOMES`, decide o bloco (campo `bloco` da linha, senão `BLOCO_ID`, senão `"G"`) e monta a pasta de resultados. | `cid`: ID, ex. `"R02"` | `dict` com todos os parâmetros, mais `id`, `nome` (ex. `R02_CPP_M2_SIM_SOU`), `bloco` e `dir` (`results/<bloco>/<nome>`) |

**Executado direto** (`python tools/casos.py`): para cada execução de `CASES`, confere que a pasta
cai no bloco certo e que `mesh/<malha>.msh` existe; imprime `casos ok: <n>`.

---

## `run_caso.py`

O **runner**: abre o Fluent, monta o caso a partir de `caso(cid)`, resolve, grava e exporta.

Uso:

```
python tools/run_caso.py R02             # resolve e exporta
python tools/run_caso.py --exportar R02  # relê .cas/.dat e refaz só a exportação
python tools/run_caso.py --continuar R04 20000 1e-6  # relê .cas/.dat e itera mais (ver continuar)
```

Constantes: `ITER_SINAL = 50`, iterações feitas no CP periódico e no duto antes de conferir o sinal
de ū; `BLOCO_SALVA = 2000`, iterações entre gravações de `.cas/.dat` nos casos turbulentos.

| Função | O que faz | Entrada | Saída |
|---|---|---|---|
| `ler_rep(path)` | Lê o Report File do Fluent, pulando o cabeçalho: só as linhas que começam com um número. | caminho do `_rep.out` | lista de pares `(iteração, valor)` |
| `parede_movel(st, zona, u, relativa=False)` | Transforma uma parede em parede móvel, translação na direção (1, 0); velocidade absoluta, ou relativa à zona adjacente com `relativa=True` (cavidade rotativa). | `st` (= `session.settings`), nome da zona, velocidade [m/s], `relativa` | nada (altera a sessão) |
| `linha_monitor(st, nome, p0, p1, campo)` | Cria uma linha de medição entre dois pontos e uma report definition `rd_<nome>` com a média por área de `campo` (padrão `x-velocity`) nessa linha. | `st`, nome, ponto inicial e final `[x, y, z]`, campo | nada (altera a sessão) |
| `modelo_turb(s, c)` | Aplica os atributos de `c["turb"]`; liga Enhanced Wall Treatment em k-ε e RSM; liga `canto` (Corner Flow Correction) e `curvatura` (Curvature Correction); define `csep` com `rpsetvar 'geko-csep` e confere o valor lido. | sessão, parâmetros | estado de `viscous`, mais `geko_csep` quando houver |
| `setup(s, c)` | Monta o caso inteiro (ver abaixo) e lê de volta da sessão o que ficou configurado. | sessão do solver, parâmetros do caso | `dict` com o setup efetivo (vira o `_setup.json`) |
| `exportar_campo(st, c, arquivo)` | Exporta todas as células (centro de célula, vírgula): pressão, pressão total, \|V\|, u, v, função de corrente em 2D ou w em 3D; na cavidade rotativa, `relative-x-velocity`, `relative-y-velocity`, `rel-velocity-magnitude`; nos casos turbulentos, k, μt, μt/μ, volume da célula e, no RSM, as tensões de Reynolds. | `st`, parâmetros, caminho relativo a `fluent/` | nada |
| `iterar_com_quadros(s, c, rep, it)` | Itera em blocos de `anim_cada` e exporta um quadro do campo (`dados/anim/<nome>_itNNNNN.csv`) ao fim de cada bloco; para quando o Fluent converge dentro do bloco ou no teto `it`. | sessão, parâmetros, Report File, teto | nada |
| `iterar_salvando(s, c, rep, it, ef)` | Itera em blocos de `BLOCO_SALVA`, gravando `.cas/.dat` e `_setup.json` ao fim de cada bloco e somando o tempo de parede; para quando o Fluent converge dentro do bloco ou no teto `it`. Permite retomar com `--continuar` depois de um reinício. | sessão, parâmetros, Report File, teto, setup efetivo | nada |
| `exportar(s, c)` | Exporta o campo para `_campo.csv` com `exportar_campo` e cada zona de `c["paredes"]` para `_parede_<zona>.csv`; calcula os valores finais dos reports `fx_top`, `fx_bottom`, `p_in`, `p_out` que existirem. | sessão com solução, parâmetros do caso | `dict` com os valores finais dos reports (vazio se não houver) |
| `encerrar(s)` | Fecha o Fluent de verdade. Manda `/exit yes` numa thread (a chamada nunca retorna), chama `s.exit()`, espera até 30 s e, se os processos continuarem vivos, mata a árvore do Fluent pelo PID com `taskkill`. | sessão aberta com `cleanup_on_exit=False` | `True` se não sobrou processo da sessão |
| `lancar(c)` | Cria `results/<bloco>/<nome>/{fluent,dados,graficos,logs}` e abre o Fluent 2024 R1: 2D ou 3D (`c["dim"]`), dupla precisão, `c["procs"]` processos (padrão 1), sem interface, diretório de trabalho `fluent/`. | parâmetros do caso | sessão do solver conectada |
| `arrumar(c)` | Depois de fechar: move os `.trn` e `*error.log` de `fluent/` para `logs/` e apaga os `cleanup-fluent*` que o Fluent deixa. | parâmetros do caso | nada |
| `run(cid)` | Execução completa: `lancar` → `setup` → (no CP periódico e no duto, 50 iterações e checagem ū > 0) → iterar até `iter_max` (o Fluent para sozinho ao atingir 1e-6; com `turb`, via `iterar_salvando`; com `anim_cada`, via `iterar_com_quadros`) → gravar `.cas.h5`/`.dat.h5` → `exportar` → gravar `_setup.json` com tempo de parede, versão do PyFluent e reports finais. Sempre termina com `encerrar` e `arrumar`, mesmo com erro. | ID | nada; lança `RuntimeError` se ū ≤ 0 após 50 iterações |
| `reexportar(cid)` | Abre o Fluent, lê `.cas.h5`/`.dat.h5` gravados e chama só `exportar`. Útil quando a exportação mudou. Não regrava o `_setup.json`. | ID | nada |
| `continuar(cid, n_iter, tol_mon)` | Lê `.cas.h5`/`.dat.h5` e segue iterando. Sem `tol_mon`: até `n_iter` a mais (para se os resíduos convergirem). Com `tol_mon`: aplica a 2ª metade do critério do guia ("monitor estável por 50 iterações"): baixa o critério de resíduo para 1e-12 e itera em blocos de 250 até a variação do monitor nas últimas 50 iterações ficar abaixo de `tol_mon` [m/s] ou somar `n_iter`. Registra o desvio em `_setup.json` (`continuacoes`, `iter_extra`, tempo somado), regrava `.cas`/`.dat` e reexporta. | ID, máximo de iterações extras, tolerância do monitor | nada |

### O que `setup()` faz, em ordem

1. Lê `mesh/<malha>.msh`.
2. CP periódico e duto: `make-periodic` de `periodic_a` com `periodic_b` (pela TUI) e confere que as
   duas zonas viraram periódicas.
3. Viscoso **laminar** ou `modelo_turb` quando há `turb`; energia **desligada**.
4. Material `fluido-tg` = cópia do ar com ρ e μ do caso, atribuído à zona de fluido.
5. Contorno por tipo:

   | Tipo | Contorno | Linha de medição | Reports extras |
   |---|---|---|---|
   | `cp_per` | `wall_top` móvel; gradiente periódico `dpdx` na direção (−1, 0) via TUI, conferido por `rpgetvar` | `line_x2p5`, (2,5; −1) → (2,5; 0) | `fx_top`, `fx_bottom` |
   | `cp_io` | `wall_top` móvel; `inlet` com o perfil `cp_inlet_ymenos.prof` (R06) ou com velocidade uniforme `uin` (R17); `outlet` a 0 Pa | `line_x2p5` | `p_in`, `p_out` (pressão média), `fx_top`, `fx_bottom` |
   | `ldc` | `lid` móvel | `line_x0p5`, (0,5; 0) → (0,5; 1) | — |
   | `cavr` | Frame Motion na zona fluida (`reference_frame.frame_motion`, `mrf_omega` = Ω, `reference_frame_axis_origin` = eixo); `lid` móvel relativa; demais paredes Stationary (paradas em relação à zona) | `line_x0p5` com `relative-x-velocity` | — |
   | `bfs` | `inlet` com o perfil `bfs_inlet.prof`; `outlet` a 0 Pa | `line_x2`, (2; 0) → (2; 1,5) | — |
   | `duto` | gradiente periódico `dpdy` via TUI, `4 0 -1 0` (direção oposta ao escoamento, como no CP), conferido por `rpgetvar` | `line_ub` (v) e `line_vs` (u na bissetriz do canto) | — |
   | `bft` | `inlet` com `bfy_inlet.prof`; I e μt/μ em `inlet` e `outlet`; `outlet` a 0 Pa | nenhuma | `rd_tau`, τ_x médio em `wall_down_bottom` |

6. Report File (`../dados/<nome>_rep.out`) e Report Plot com todos os reports `rd_*`.
7. Métodos: algoritmo (`alg`), gradiente Least Squares Cell Based, pressão (`pres`), momento (`mom`),
   sem High Order Term Relaxation; as equações de turbulência usam o mesmo esquema do momento (SOU),
   porque o padrão do Fluent muda com o modelo. URFs de pressão e momento (exceto no Coupled); as de
   turbulência ficam no padrão.
8. Resíduos **absolutos** 1e-6 em todas as equações. No duto, `continuity` fica fora do critério: o
   campo inicial já a satisfaz e o fator de escala fica em ~1e-16.
9. Inicialização standard com pressão e velocidades nulas, sobrepostas pelos valores de `ini` do tipo.

---

## `rodar_lote.sh`

Script de Git Bash. Para cada ID passado (ou, sem argumentos, a lista padrão R01–R20 numa ordem
fixa), descobre a pasta da execução com `casos.caso`, cria `logs/` e roda
`.venv241/Scripts/python.exe -u tools/run_caso.py <ID>`, mandando toda a saída para
`logs/run.log`. Imprime no terminal o início e o `exit=` de cada execução e `LOTE FIM` no final.
Uma execução que falha não interrompe as seguintes. Já define `PYTHONIOENCODING=utf-8`.
Com `PARALELO=N`, roda até N IDs ao mesmo tempo (`xargs -P`); o padrão é 1.

## `rodar_fila.sh`

Fila com teto de Fluents simultâneos: `MAX=4 bash tools/rodar_fila.sh R38:c R40 R33:t`. Conta as
sessões vivas de qualquer lote e só lança a próxima quando há vaga. Cada Fluent compromete cerca de
4 GB de memória virtual; dimensione `MAX` pela memória da máquina.

| Forma | O que faz |
|---|---|
| `Rxx` | Roda do zero, com a saída em `logs/run.log`. |
| `Rxx:c` | Continua de `.cas/.dat` até os resíduos convergirem, com até 16 000 iterações a mais. |
| `Rxx:t` | Continua até o monitor variar menos que 1e-4 em 50 iterações, com até 5000 iterações a mais. |

As continuações gravam em `logs/continuacao1.log`.

---

## `analise_caso.py`

Lê os arquivos de `dados/`, calcula as métricas e faz os gráficos. Não abre o Fluent.

Uso:

```
python tools/analise_caso.py R02 [R03 ...]   # analisa estas
python tools/analise_caso.py --todos         # todas as de CASES que têm _setup.json
python tools/analise_caso.py --selftest      # confere as métricas no CSV de referência do R02
```

Com IDs ou `--todos`, primeiro calcula as métricas de todas (para a ordem observada da série de
malhas estar disponível), depois gera os gráficos de cada uma.

### Constantes

| Nome | O que guarda |
|---|---|
| `GHIA_CENTRO` | Centro do vórtice principal da cavidade, Re = 100 (Ghia, Tabela III): (0,6172; 0,7344) |
| `GHIA` | `refs/ghia1982.csv` lido com pandas: colunas `tabela` (I = u na vertical x = 0,5; II = v na horizontal y = 0,5), `pos`, `vel` |
| `DRIVER` | `refs/bfs_driver1985.csv` lido com pandas: colunas `serie` (`xr`, `cp_step`, `cp_top`, `cf`), `pos`, `val`, `dcf_rel` |
| `XR_H_EXP`, `XR_H_INC` | x_r/H experimental e sua incerteza: 6,26 e 0,10 |

### Leitura

| Função | O que faz | Entrada | Saída |
|---|---|---|---|
| `ler_campo(path)` | Lê o `_campo.csv` e tira os espaços dos nomes das colunas. | caminho do CSV | `DataFrame`, uma linha por célula |
| `ler_rep(path)` | Lê o `_rep.out`, só as linhas numéricas. | caminho | matriz `numpy` (n × 2): iteração, valor |
| `ler_residuos(csv_path, log_path)` | Extrai os resíduos das linhas de iteração que o Fluent imprime no `run.log` e nos `continuacao*.log` (o log é completo; o histórico do serviço de monitores vem com lacunas depois de uma continuação) e grava o `_residuos.csv`. Os nomes das equações vêm do cabeçalho `iter ... time/iter` do log, então qualquer modelo funciona. | caminho do CSV e do log | `DataFrame` com `iter` e uma coluna por equação |
| `arquivos(c)` | Monta os caminhos dos arquivos de uma execução. | parâmetros do caso | `dict` com `campo`, `rep`, `res`, `setup`, `met`, `log`, `graf` |

### Métricas

| Função | O que faz | Entrada | Saída |
|---|---|---|---|
| `u_analitica(d)` | Perfil analítico do CP: u(d) = 2d + 20d(1−d). | d = −y [m] (número ou vetor) | u [m/s] |
| `perto(serie, valor, tol)` | Máscara das células cuja coordenada está a menos de `tol` (1 mm) de `valor`. Necessária porque os centros exportados têm ruído de ~1e-6 m; comparação exata perde células. | coordenadas, valor, tolerância | máscara booleana |
| `coluna(df, eixo, valor)` | Pega a coluna (ou linha) de centros de célula mais próxima de uma coordenada, ordenada pelo outro eixo. | campo, `"x-coordinate"` ou `"y-coordinate"`, valor [m] | (coordenada efetiva, células) |
| `metricas_cp(df)` | Compara o perfil na coluna mais próxima de x = 2,5 com a analítica, nos centros de célula, usando d = −y. | campo | `x0`, `ubar` (média simples das células), `eL2` = ‖u−uₐ‖₂/‖uₐ‖₂ em %, `einf` = max\|u−uₐ\| |
| `max_du(df, cid_ref)` | Maior diferença de u célula a célula contra outra execução na mesma malha (casa as células pelas coordenadas). | campo, ID da referência | max\|Δu\| [m/s] (com `detalhe=True`: também o RMS e o máximo sem a fileira junto à tampa), ou `None` se a referência ainda não rodou; falha se as malhas forem diferentes |
| `metricas_ldc(df)` | Centro do vórtice = célula com maior \|ψ − ψ_parede\| (ψ_parede = média de ψ na fileira junto às paredes, porque o ψ do Fluent tem constante arbitrária; precisão de uma célula); distância ao centro de Ghia; extremos de u na vertical x ≈ 0,5 e de v na horizontal y ≈ 0,5; max\|u\|; erro contra Ghia via `comparar_ghia`. | campo (com `stream-function`) | `xc`, `yc`, `dist_ghia`, `u_min_vert`, `v_max_hor`, `v_min_hor`, `u_max`, `ghia_du_max`, `ghia_du_rms`, `ghia_dv_max`, `ghia_dv_rms` |
| `linha_central(df, eixo, comp)` | Perfil na linha central (eixo = 0,5): média das duas colunas vizinhas de 0,5 (N par, sem centro de célula em 0,5) mais os valores de parede nas pontas (u = 1 só na ponta y = 1 da vertical x = 0,5, que é a tampa; 0 nas demais). | campo, eixo, componente | (posição, valor) |
| `comparar_ghia(df)` | Interpola o perfil de `linha_central` nos 15 pontos internos de cada tabela de Ghia e calcula erro = Fluent − Ghia. | campo da cavidade | dict `I`/`II` de DataFrames (`pos`, `vel`, `fluent`, `erro`) |
| `metricas_cavr(df, cid)` | Bloco G (Adendo 3.3), no referencial rotativo: e_L2 e e∞ de u e v contra o R26 interpolado (triangulação linear do matplotlib) nos centros desta malha; overshoot e subvalor de u na horizontal y = 0,5 contra o R26; índice de checkerboard CB da pressão na fileira mais próxima de y = 0,5; `ref_dif` (% de max\|V\|) no R27, também sem a fileira junto à tampa e em RMS. Os erros são calculados também contra o R27 (Second Order, sufixo `_R27`), porque a referência R26 é PRESTO! e favorece o R25. | campo, ID | `u_max`, `cb`, `u_max_lin`, `u_min_lin`, `eL2_u`, `eL2_v`, `einf_u`, `einf_v`, `overshoot`, `subvalor`, os mesmos com `_R27`, `ref_dif`, `ref_dif_sem_tampa`, `ref_dif_rms` |
| `metricas_bfs(df)` | Comprimento de recolamento: na primeira linha de células acima de y = 0 e depois do degrau (x > 1), acha onde u passa de negativo para positivo e interpola linearmente. | campo | `xr` (a partir de x = 1) e `xr_h` = x_r / 0,5; se a bolha atravessa a saída (u < 0 na última célula), `xr = None` e `xr_maior_que` traz a cota inferior |
| `ler_parede(c, zona)` | Lê o `_parede_<zona>.csv`. | parâmetros, zona | `DataFrame`, uma linha por face |
| `metricas_duto(c, df)` | U_b pela média de v ponderada por volume, max\|V_s\| em (u, w), Re_b, Cf = 2τ/ρU_b² com τ = 1 Pa pelo balanço G·A = τ·P, y+ de `wall`. | parâmetros, campo | `ub`, `vs_max`, `vs_ub` [%], `re_b`, `cf`, `yplus_max`, `yplus_med` |
| `paredes_bft(c)` | Junta `wall_in_bottom` e `wall_down_bottom` (lado do degrau) e lê `wall_top`; fixa o sinal de τ_x pela região recolada junto à saída; Cf = τ_x/q e Cp referido à pressão em x/H = −4, como nos dados originais. | parâmetros | dois `DataFrame` com `x`, `cf`, `cp`, `yplus` |
| `metricas_bft(c, df)` | x_r pela última troca negativo→positivo de Cf em x > 0 (ignora o vórtice de canto), erro contra 6,26, μt máximo em 0 < x < 12 e 0,2 < y < 2,5, y+ no fundo a jusante e RMS de Cf e Cp contra Driver em −8 < x/H < 29. | parâmetros, campo | `xr`, `xr_h`, `erro_xr`, `mut_max`, `mut_mu_max`, `yplus_max`, `yplus_med`, `cf_rms`, `cp_step_rms`, `cp_top_rms` |
| `calcular(cid)` | Lê campo, monitor, resíduos e setup; calcula as métricas gerais (células, iterações, convergência = todos os resíduos finais < 1e-6, sem a continuidade no duto, resíduos finais, variação do monitor nas últimas 50 iterações) e as do tipo. Comparações automáticas: R05 contra R02; R12–R14 e R18–R20 contra R08; R18 e R20 contra R19. Grava `_metricas.json`. | ID | `dict` de métricas |
| `ordem_observada()` | Ordem observada p = ln(e_k/e_{k+1}) / ln 2 do e_L2 entre execuções vizinhas da série R01–R04, só para pares que já têm `_metricas.json`. | — | `dict`, ex. `{"R01→R02": 2.0, ...}` |

### Gráficos

| Função | O que faz | Entrada | Saída |
|---|---|---|---|
| `graficos(cid, m)` | Gera três PNGs em `graficos/` a partir do campo exportado (triangulação dos centros de célula; no degrau, mascara a região sólida: x < 1, y < 0,5 no laminar e x < 0, y < 1 no turbulento; no duto, desenha o plano xz com vetores do escoamento secundário). `_campos.png`: \|V\| com vetores de direção (e linhas de corrente na cavidade e no degrau) e a pressão. `_convergencia.png`: resíduos e monitor (no CP periódico, com a linha ū = 4,3333). `_perfis.png`: perfis do tipo de caso (ver `docs/como_rodar.md`). | ID, métricas | lista dos PNGs gravados |

### Resumo

| Função | O que faz | Entrada | Saída |
|---|---|---|---|
| `resumo_metricas(cid, m)` | Texto curto `nome=valor` com as métricas principais, para o terminal. | ID, métricas | `str` |
| `animacao(cid)` | Monta `graficos/<nome>_animacao.gif` com os quadros de `dados/anim/`: \|V\| em cores (escala do campo final) e vetores de direção, um quadro por bloco de iterações. | ID | caminho do GIF ou `None` |
| `analisar(cid)` | `calcular` + `graficos` + `animacao`, e imprime o resumo. | ID | métricas |
| `selftest()` | Confere que `ghia1982.csv` tem 17 + 17 pontos. Se existir o CSV local `refs/r02cpp_m2_sim_sou_campo.csv` (fora do repositório), calcula `metricas_cp` nele e confere ū = 4,3499, e_L2 = 0,2693 %, e∞ = 0,01249 m/s; sem ele, avisa que o teste foi parcial. | — | nada; `AssertionError` se divergir |
