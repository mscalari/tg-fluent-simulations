# O código (`tools/`)

Visão do conjunto:

```
casos.py ──► run_caso.py ──► results/<bloco>/<nome>/{fluent,dados,logs}
   │           ▲                         │
   │      rodar_lote.sh, rodar_fila.sh       ▼
   └─────────────────────────► analise_caso.py ──► _metricas.json, graficos/
gerar_bfy_inlet.py: refs/bfy_inlet.prof
```

O conversor CT_V1 (MRF) e as suas análises estão na segunda metade deste arquivo, a partir de "Conversor de torque CT_V1 (MRF)".

| Arquivo | Papel | Abre o Fluent? |
|---|---|---|
| `casos.py` | Tabela das execuções e dos parâmetros físicos | não |
| `run_caso.py` | Monta, resolve e exporta uma execução no Fluent | **sim** |
| `rodar_lote.sh` | Roda várias execuções, em sequência ou com `PARALELO=N`, com log | **sim** (via `run_caso.py`) |
| `rodar_fila.sh` | Roda uma fila com teto de Fluents simultâneos e retoma casos de `.cas/.dat` | **sim** (via `run_caso.py`) |
| `analise_caso.py` | Métricas e gráficos | não |
| `gerar_bfy_inlet.py` | Gera `refs/bfy_inlet.prof` (lei 1/7 nas duas paredes, δ = 1,4H, x = −10) | não |
| `status.py` | Painel ao vivo do lote, a partir do `run.log` (`python tools/status.py R31 R32 ...`); local, fora do git | não |
| `check_msh.py`, `msh_bordas_qualidade.py` | Conferência de malhas `.msh`; locais, fora do git | não |
| `notas_obsidian.py`, `hipoteses_gh.py` | Registro `registro_pyfluent.csv` e notas do Obsidian por rodada (`Rxx_<nome>.md`) | não (locais, fora do repositório) |

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

---

# Conversor de torque CT_V1 (MRF)

Os scripts do conversor ficam separados do runner dos casos R01–R43: `conversor_mrf.py` e `malha_gci_conversor.py`
são os que abrem o Fluent; as análises só leem arquivos já gravados. Tudo roda com `.venv241` e
`PYTHONIOENCODING=utf-8`. Cada ponto de operação fica em `results/conversor/<saída>/SRxxx/{fluent,dados,graficos}`
(`SR040` = SR 0,4; `SR085` = 0,85); a base é `CT_V1_MRF`, e as malhas do GCI usam `CT_V1_GCI_G1` e `CT_V1_GCI_G2`.

```
conversor_mrf.py ──► CT_V1_MRF/SRxxx/{fluent,dados,graficos}, dados/*.out, sens/<variante>/
        │                         │
        │                         ├─► analise_conversor.py ──► desempenho.csv, g1–g3, nota do Obsidian
        │                         ├─► euler_interfaces_conversor.py (g7)    analise_yplus_conversor.py (g4)
        │                         ├─► analise_sensibilidade_conversor.py (g6)
        │                         ├─► comparacao_literatura_conversor.py (g5) ─► tendencia_{derivada,sinais,colapso} (g8–g10)
        │                         └─► gci_conversor.py ◄── CT_V1_GCI_G1/G2 ◄── conversor_mrf.py --malha ... --saida ...
malha_gci_conversor.py ──► mesh/CT_V1_3zonas_G1/G2.msh.h5   (Fluent Meshing)
```

| Arquivo | Papel | Abre o Fluent? |
|---|---|---|
| `conversor_mrf.py` | Setup MRF, varredura de SR, sensibilidades, figuras | **sim** |
| `malha_gci_conversor.py` | Gera as malhas G1 e G2 do GCI | **sim** (Meshing) |
| `analise_conversor.py` | Curvas globais, torque × Euler, nota do Obsidian | não |
| `euler_interfaces_conversor.py` | Balanço de momento angular por zona | não |
| `analise_yplus_conversor.py` | y⁺ por faixas e peça | não |
| `analise_sensibilidade_conversor.py` | Efeito de cada variante do setup no estol | não |
| `comparacao_literatura_conversor.py` | TR, η e K contra a literatura | não |
| `tendencia_derivada_conversor.py`, `tendencia_sinais_conversor.py`, `tendencia_colapso_conversor.py` | Três testes de tendência | não |
| `animacao_meio_passo.py` | GIF da turbina girada de meio passo | não |
| `gci_conversor.py` | GCI de Celik et al. (2008) | não |
| `imagens_malha_conversor.py` | 10 vistas por malha (G0, G1, G2) e montagens lado a lado | **sim** (2 processos, ≤ 1,7 GB) |
| `comparacao_blocos.py` | Figuras-resumo de R01–R43 contra as referências | não |

## `conversor_mrf.py`

Opções globais, válidas em qualquer comando: `--malha <arquivo .msh.h5>` (padrão `mesh/CT_V1_3zonas.msh.h5`),
`--saida <pasta em results/conversor>` (padrão `CT_V1_MRF`) e `--urf chave=valor,...` (sub-relaxação extra, global
`URF`). Servem para rodar as malhas do GCI sem tocar na base. Funções acrescentadas para as malhas do GCI:
`normalizar_zonas` (tira o sufixo `_NN` das zonas fluid e interior), `renomear_paredes_stl` (paredes `stlb_*` →
`parede_<dono>_lado_<lado>`), `sincronizar_paredes_abs` (refaz `PAREDES_ABS` pelas paredes da malha; `PAREDES_ABS_G0`
guarda as da G0) e `divergiu` (guarda de divergência chamada a cada bloco em `iterar_ponto`).

| Comando | O que faz |
|---|---|
| `setup` | Monta o caso, testa imagens e grava `CT_V1_MRF_setup.cas.h5`. |
| `sentido [N]` | SR = 0, N iterações em cada sentido do eixo (+x e −x). Escolhido: +x (TR > 1 e estator reagindo). |
| `varredura S [SRs] [--ini X.dat.h5]` | S = +1 ou −1; SRs padrão: 0; 0,2; 0,4; 0,6; 0,7; 0,8; 0,85; 0,9. Cada ponto parte do anterior; `--ini` dá a partida do 1º. Grava `resultado.json` por ponto. |
| `sensibilidade [variantes]` | Variantes no SR = 0, em `<saída>/sens/<variante>`: `sst_sem_cc`, `ke_realizable` (funções de parede *scalable*), `fou` (1ª ordem em momento e turbulência), `meio_passo` (gira a zona da turbina 10° com `mesh.modify_zones.rotate_zone`, antes de criar as interfaces, e parte do repouso) e `coupled` (só em malha grossa, por RAM). Cada variante usa o seu próprio report file. |
| `refazer_figuras [SRs]` | Relê o `.cas/.dat` de cada ponto numa sessão nova e refaz as imagens, a verificação Cm → N·m, o torque das paredes `sobra_*` e os exports Euler (os dois lados de cada interface). |
| `figuras_teste <caso>` | Gera as imagens de um caso gravado (diagnóstico). |

**Setup (`montar`).** Lê a malha (em mm, convertida para m), troca para `interface` as zonas das interfaces
(`if_bt` bomba–turbina, `if_te` turbina–estator, `if_eb` estator–bomba; as faces sem par viram paredes
`sobra_<if>_lado_<peça>`), liga o Frame Motion nas zonas de fluido, põe as paredes em movimento absoluto com a peça
vizinha, cria os reports de torque por peça (`report_definitions.moment`) e de vazão, e define o modelo, os
métodos e a inicialização: velocidade relativa (`velocity_formulation = "relative"`), k-ω SST com
`curvature_correction`, PRESTO! na pressão, segunda ordem nas demais equações, SIMPLE, gradiente Least Squares Cell
Based, resíduos só monitorados (a parada é pelos torques) e inicialização standard em referencial absoluto.

**Critério de parada** (`iterar_ponto`, `convergiu`, `media_janela`): blocos de `BLOCO` = 100 iterações; para quando
a média dos torques das últimas `JANELA` = 200 iterações varia menos que `TOL` = 0,5 % em relação às 200 anteriores
(T_e relativo a |T_b|), com mínimo de 1000 (1º ponto) ou 600 iterações e teto de 3000 ou 2000. Os valores oficiais são
essas médias, e o desvio na janela é a incerteza iterativa. Levanta `RuntimeError` se um bloco não acrescentar linhas ao
report file.

**Torque em N·m.** `momento_em_nm` escolhe, em `report_output_type.allowed_values()`, a opção sem "coef"
(`Drag Force`); `corrigir_momentos` converte casos gravados antes da correção (M = 0,6125·Cm). Resultados novos
gravam `unidade_torque = "N·m"`.

**Funções principais.**

| Função | O que faz |
|---|---|
| `abrir()` | Abre o Fluent 3D, dupla precisão, 6 processos, sem interface. |
| `criar_interfaces`, `adjacencia`, `nome_lado` | Cria as três mesh interfaces e renomeia as paredes pela adjacência. |
| `aplicar_rotacao`, `frame`, `girar_parede` | Frame Motion nas zonas e paredes rotativas para o SR dado. |
| `setup_do_caso(st)` | Recupera `sobras` e `grupos_torque` de um caso lido (quando falta o `setup.json`). |
| `montar(s, sinal)` | Setup completo (acima). |
| `reports`, `torques_por_zona`, `finais` | Reports de torque e vazão e leitura dos valores. |
| `camera`, `VISTAS`, `contorno`, `vetores`, `salvar`, `figuras` | Imagens do Fluent (ver armadilhas em `como_rodar.md`). |
| `exportar`, `verificar_ponto` | Exports ASCII (campo, paredes com área, lados das interfaces) e verificação no Fluent. |
| `cmd_*` | Um comando de linha de comando por função (`cmd_varredura`, `cmd_sensibilidade`, `cmd_refazer_figuras`, ...). |
| `opcoes_globais(argv)` | Trata `--malha` e `--saida`. |

Em malhas grossas, as paredes de `PAREDES_ABS` que não existirem (zonas pequenas somem) são removidas com aviso.

## Análises (sem Fluent)

Convenções: CSVs com `;` como separador; figuras `g1`…`g10` em `results/conversor/CT_V1_MRF/graficos/`, tabelas em
`.../dados/`. Torques gravados em Cm são convertidos para N·m por `analise_conversor.fator_torque`.
TR = −T_t/T_b, η = SR·TR, K = n_b/√|T_b| (n_b em rpm).

| Script | Como rodar | Entradas | Saídas |
|---|---|---|---|
| `analise_conversor.py` | `python tools/analise_conversor.py [--sem-nota]` | `SRxxx/dados/resultado.json`, `dados/*.out`, exports por face | `g1_TR_eta_K`, `g2_torques_vazao_balanco`, `g3_torque_vs_euler`; gráficos por ponto em `SRxxx/graficos`; `dados/desempenho.csv`; nota `05_Conversor/CT_V1_MRF.md` no Obsidian (imagens copiadas para `C_TG/img/`). Marca o ponto de acoplamento (TR = 1, SR ≈ 0,51). |
| `euler_interfaces_conversor.py` | `python tools/euler_interfaces_conversor.py [SR em % ...] [--diag-vorticidade]` | Exports `_euler_<lado>.csv` e `_area.csv` por ponto | `g7_euler_interfaces`, `dados/euler_interfaces.csv` e `euler_interfaces_if.csv`. Balanço de momento angular por zona (convectivo com fluxo absoluto e relativo, difusivo, pressão) contra o torque das paredes; razão Euler/parede. |
| `analise_yplus_conversor.py` | `python tools/analise_yplus_conversor.py` | `SRxxx/dados/SRxxx_parede_<peça>[_area].csv` | `g4_yplus_faixas`, `dados/yplus_faixas.csv`. Frações por faixa (< 5, 5–30, 30–300, 300–1000, > 1000), P50 e média; ponderado por área quando há `_area.csv`, senão por face (marcado com `*`). |
| `analise_sensibilidade_conversor.py` | `python tools/analise_sensibilidade_conversor.py [--dir-sens P] [--saida P]` | Base `SR000` e `sens/<variante>/dados/resultado.json` | `g6_sensibilidade`, `dados/sensibilidade.csv`: Δ% de TR, T_b, T_t e K e \|Δ\|/U, com U = √(SEM_var² + SEM_base²) (SEM por batch means). Variante ausente é pulada. |
| `comparacao_literatura_conversor.py` | `python tools/comparacao_literatura_conversor.py` | `dados/desempenho.csv`, `refs/conversor_literatura.csv` (fontes em `refs/conversor_literatura.md`) | `g5_comparacao_literatura` (TR, TR/TR_estol, η e K/K_estol × SR). Compara tendência, não valor: as geometrias são outras. |
| `tendencia_derivada_conversor.py` | `python tools/tendencia_derivada_conversor.py` | `desempenho.csv`, `.out`, literatura | `g8_derivada_tendencia`, `dados/derivada_tendencia.csv`: diferença finita no ponto médio de cada intervalo de SR, com U = 2·SEM. |
| `tendencia_sinais_conversor.py` | `python tools/tendencia_sinais_conversor.py` | idem | `g9_sinais_tendencia`, `dados/sinais_tendencia.csv`: cada intervalo é "confirmado", "sinal certo, dentro da incerteza" ou "sinal errado". |
| `tendencia_colapso_conversor.py` | `python tools/tendencia_colapso_conversor.py` | idem | `g10_colapso_adimensional`, `dados/colapso_adimensional.csv` e `_resumo.csv`: (TR−1)/(TR_estol−1) × SR/SR_acopl e η/η_máx × SR/SR*, com PCHIP implementado no script (o `.venv241` não tem scipy). |
| `animacao_meio_passo.py` | `python tools/animacao_meio_passo.py [--limiar 0.85] [--quadros 48] [--fps 12] [--teste]` | Faces de parede do `SR000` | `animacao_meio_passo.gif` e os PNGs de 0° e 10°. Conta as pás (bomba 20, turbina 18, estator 18) pela harmônica angular. |
| `gci_conversor.py` | `python tools/gci_conversor.py [--n0 N --v0 V] [--p-formal 2] [--mc 0] [--teste]` | `resultado.json` e `.out` de SR000 e SR040 em `CT_V1_MRF`, `CT_V1_GCI_G1`, `CT_V1_GCI_G2`; N e V do transcript do solver | `CT_V1_GCI/dados/gci.csv`, `graficos/gci_tr.png`. Ver abaixo. |
| `malha_gci_conversor.py` | `python tools/malha_gci_conversor.py G1 G2 [--so-superficie] [--retomar-superficie] [--partir-da-superficie] [--fator f]` | `.wft` da G0 e STLs | `mesh/CT_V1_3zonas_G{1,2}.msh.h5`, `CT_V1_GCI/logs/<nível>_parametros.json`. |
| `imagens_malha_conversor.py` | `python tools/imagens_malha_conversor.py G2 [G1 G0] [--processos N]` ou `--so-montagem` | `mesh/CT_V1_3zonas[_G1\|_G2].msh.h5` | `CT_V1_GCI/graficos/malha_<nível>_NN_<vista>.png`, `malha_comparacao_<vista>.png`, `CT_V1_GCI/dados/malha_<nível>_info.json`. Recusa começar com menos de 3000 MB livres de RAM. |
| `comparacao_blocos.py` | `python tools/comparacao_blocos.py [A B ... H]` | `_campo.csv`, `_parede_*.csv`, `_metricas.json`, `refs/` | `results/comparacao_blocos/bloco_{A..G}.png`, `bloco_H_duto.png`, `bloco_H_degrau.png`, `resumo.csv`. Rodadas não convergidas aparecem com `*`, hachura e linha tracejada. Sem referência, só mostra as rodadas. |

### GCI (`gci_conversor.py`, `malha_gci_conversor.py`)

- **Malhas.** G0 = 1 900 736 células e V = 8,576847e-2 m³ (h médio = (V/N)^(1/3) = 3,56 mm). G1 e G2 saem do mesmo
  Fault-tolerant Meshing da G0, escalando os tamanhos por f (G1: 1,35; G2: 1,35² = 1,8225), a taxa de crescimento por
  1 + 0,2·f (limite 1,34), o ângulo de curvatura por f e as células por folga por 1/f (limite 1,775). A calibração pela
  contagem de faces da superfície deu r efetivo = 1,308 na G1 (842 090 faces contra 1 440 968 da G0); decidiu-se manter os
  fatores nominais. O GCI usa o r medido (N e V reais do solver).
- **Rodadas.** SR = 0 e SR = 0,4 em cada malha (`conversor_mrf.py varredura 1 0 0.4 --malha ... --saida ...`).
- **Método (Celik et al., 2008).** p é a menor raiz das equações 3a–3c em [0,05; 20] (varredura + bisseção; o ponto
  fixo do artigo falha com r32 ≫ r21); sem raiz ou fora de [0,5; 10], usa p formal e sinaliza; raiz acima de 1,1·p formal
  é limitada ao formal (evita GCI otimista com ε no ruído). Com p aparente, a razão assintótica é |φ1/φ2| por construção e
  não é usada como teste; se p estiver a ≥ 10 % do formal, a observação manda usar o GCI conservador. A coluna `GCI21_cons_pct`
  usa Fs = 3 e p limitado quando \|p − p_f\|/p_f ≥ 10 %. O ruído iterativo entra pelo SEM por batch means (5 lotes de 40
  iterações), com τ integral como conferência, e há Monte Carlo opcional (300 amostras). Marca a malha INVÁLIDA se
  \|V/V0 − 1\| > 2 % ou se r se afastar mais de 15 % de 1,35 (saída com código 3). `GCI21_rel_Tb_pct` dá o GCI de T_e
  como fração de |T_b|, pois T_e é pequeno e troca de sinal perto do acoplamento.
- **Pendente.** O exemplo da Tabela 1 de Celik dá GCI21 = 2,17 % pelas fórmulas e o revisor lembra 2,09 % publicado:
  o Murilo confere no artigo.
- **Resultado (2026-11-22, revisado por um revisor independente).** SR 0: TR 1,065 / 1,062 / 1,069 (G0/G1/G2), GCI21
  0,49 % com p formal, diferença no ruído iterativo. SR 0,4: torques caem ~5 % por nível (p 0,65, fora da faixa
  assintótica; GCI de T_b 32 %, conservador 78 %); TR da G0 1,027 contra 1,015–1,016 em G1/G2, confundido com o protocolo
  (só a G0 usou URF padrão e passou por SR 0,2). Testes para separar: `CT_V1_GCI_G0_urf_reinicio` (G0 SR 0,4 reiniciado
  com URF reduzida) e `CT_V1_GCI_G0_urf` (G0 SR 0 → 0,4 direto com URF reduzida).

### `monitor/` (fora do git)

Servidor local somente leitura (`servidor.py`, biblioteca padrão, 127.0.0.1:8765), página `index.html` e
`abrir_monitor.ps1`. Mostra processos do Fluent (e órfãos), scripts Python e filas, RAM, logs com o último it+ e alerta
de log parado, resultados da varredura e da sensibilidade. Não inicia nem mata nada. Ver `monitor/README.md`.
