# O código (`tools/`)

Visão do conjunto:

```
casos.py ──► run_caso.py ──► results/<bloco>/<nome>/{fluent,dados,logs}
   │           ▲                         │
   │      rodar_lote.sh (lote)           ▼
   └─────────────────────────► analise_caso.py ──► _metricas.json, graficos/,
                                                    registro_pyfluent.csv, notas
check_msh.py, msh_bordas_qualidade.py: conferência das malhas (independentes do Fluent)
```

| Arquivo | Papel | Abre o Fluent? |
|---|---|---|
| `casos.py` | Tabela das execuções e dos parâmetros físicos | não |
| `run_caso.py` | Monta, resolve e exporta uma execução no Fluent | **sim** |
| `rodar_lote.sh` | Roda várias execuções em sequência, com log | **sim** (via `run_caso.py`) |
| `analise_caso.py` | Métricas, gráficos, registro CSV e notas | não |
| `check_msh.py` | Lê um `.msh` binário 2D e confere contagens, extents e zonas | não |
| `msh_bordas_qualidade.py` | Espaçamento nas fronteiras e qualidade das células de um `.msh` | não |

---

## `casos.py`

Define **o que** cada execução é. Não faz nada no Fluent; é importado pelos outros scripts.

### Constantes

| Nome | O que guarda |
|---|---|
| `ROOT` | Raiz do repositório (pasta acima de `tools/`). |
| `RESULTS` | `ROOT/results`. |
| `BLOCOS` | Letra do bloco → nome da pasta em `results/` (ex.: `"A"` → `"A_CP_periodico"`). |
| `TIPOS` | Um dicionário por **tipo de caso** (`cp_per`, `cp_io`, `ldc`, `cavr`, `bfs`). Cada um tem os números usados pelo runner (`rho`, `mu`, `u_parede`, `dpdx` no CP periódico, `omega` e `eixo` na cavidade rotativa, `anim_cada` quando há animação, `iter_max`) e quatro textos explicativos que vão para a nota: `geometria`, `fluido`, `contorno`, `esperado`. |
| `BASE` | Métodos-padrão de todas as execuções: `alg="SIMPLE"`, `mom="second-order-upwind"`, `pres="second-order"`, `urf_p=0.3`, `urf_mom=0.7`. |
| `CASES` | Uma linha por execução rodável (hoje R01–R27). Cada linha diz o `tipo`, a `malha`, o texto `muda` (o que a execução testa) e só o que difere do padrão (`mom`, `pres`, `alg`, `urf_p`, `iter_max`, `perfil`, `uin`, `bloco`, `anim_cada`). |
| `BLOCO_ID` | ID → letra do bloco, para R01–R17. |

### Funções

| Função | O que faz | Entrada | Saída |
|---|---|---|---|
| `caso(cid)` | Junta tudo o que se sabe de uma execução: a física do tipo (`TIPOS`), os padrões (`BASE`) e a linha de `CASES` (nessa ordem; a linha de `CASES` vence). Lê o nome-base (`Arquivo_base`) em `00_admin/registro_execucoes_v2.csv`, decide o bloco (campo `bloco` da linha, senão `BLOCO_ID`, senão `"G"`) e monta a pasta de resultados. | `cid`: ID, ex. `"R02"` | `dict` com todos os parâmetros, mais `id`, `nome` (ex. `R02_CPP_M2_SIM_SOU`), `bloco` e `dir` (`results/<bloco>/<nome>`) |

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

Constante `ITER_SINAL = 50`: iterações feitas no CP periódico antes de conferir o sinal de ū.

| Função | O que faz | Entrada | Saída |
|---|---|---|---|
| `ler_rep(path)` | Lê o Report File do Fluent, pulando o cabeçalho: só as linhas que começam com um número. | caminho do `_rep.out` | lista de pares `(iteração, valor)` |
| `parede_movel(st, zona, u, relativa=False)` | Transforma uma parede em parede móvel, translação na direção (1, 0); velocidade absoluta, ou relativa à zona adjacente com `relativa=True` (cavidade rotativa). | `st` (= `session.settings`), nome da zona, velocidade [m/s], `relativa` | nada (altera a sessão) |
| `linha_monitor(st, nome, p0, p1, campo)` | Cria uma linha de medição entre dois pontos e uma report definition `rd_<nome>` com a média por área de `campo` (padrão `x-velocity`) nessa linha. | `st`, nome, ponto inicial e final `[x, y, z]`, campo | nada (altera a sessão) |
| `setup(s, c)` | Monta o caso inteiro (ver abaixo) e lê de volta da sessão o que ficou configurado. | sessão do solver, parâmetros do caso | `dict` com o setup efetivo (vira o `_setup.json`) |
| `exportar_campo(st, c, arquivo)` | Exporta todas as células (centro de célula, vírgula): pressão, pressão total, \|V\|, u, v, função de corrente e, na cavidade rotativa, `relative-x-velocity`, `relative-y-velocity`, `rel-velocity-magnitude`. | `st`, parâmetros, caminho relativo a `fluent/` | nada |
| `iterar_com_quadros(s, c, rep, it)` | Itera em blocos de `anim_cada` e exporta um quadro do campo (`dados/anim/<nome>_itNNNNN.csv`) ao fim de cada bloco; para quando o Fluent converge dentro do bloco ou no teto `it`. | sessão, parâmetros, Report File, teto | nada |
| `exportar(s, c)` | Exporta o campo para `_campo.csv` com `exportar_campo`; calcula os valores finais dos reports `fx_top`, `fx_bottom`, `p_in`, `p_out` que existirem. | sessão com solução, parâmetros do caso | `dict` com os valores finais dos reports (vazio se não houver) |
| `encerrar(s)` | Fecha o Fluent de verdade. Manda `/exit yes` numa thread (a chamada nunca retorna), chama `s.exit()`, espera até 30 s e, se os processos continuarem vivos, mata a árvore do Fluent pelo PID com `taskkill`. | sessão aberta com `cleanup_on_exit=False` | `True` se não sobrou processo da sessão |
| `lancar(c)` | Cria `results/<bloco>/<nome>/{fluent,dados,graficos,logs}` e abre o Fluent 2024 R1: 2D, dupla precisão, 1 processo, sem interface, diretório de trabalho `fluent/`. | parâmetros do caso | sessão do solver conectada |
| `arrumar(c)` | Depois de fechar: move os `.trn` e `*error.log` de `fluent/` para `logs/` e apaga os `cleanup-fluent*` que o Fluent deixa. | parâmetros do caso | nada |
| `run(cid)` | Execução completa: `lancar` → `setup` → (no CP periódico, 50 iterações e checagem ū > 0) → iterar até `iter_max` (o Fluent para sozinho ao atingir 1e-6; com `anim_cada`, via `iterar_com_quadros`) → gravar `.cas.h5`/`.dat.h5` → `exportar` → gravar `_setup.json` com tempo de parede, versão do PyFluent e reports finais. Sempre termina com `encerrar` e `arrumar`, mesmo com erro. | ID | nada; lança `RuntimeError` se ū ≤ 0 após 50 iterações |
| `reexportar(cid)` | Abre o Fluent, lê `.cas.h5`/`.dat.h5` gravados e chama só `exportar`. Útil quando a exportação mudou. Não regrava o `_setup.json`. | ID | nada |
| `continuar(cid, n_iter, tol_mon)` | Lê `.cas.h5`/`.dat.h5` e segue iterando. Sem `tol_mon`: até `n_iter` a mais (para se os resíduos convergirem). Com `tol_mon`: aplica a 2ª metade do critério do guia ("monitor estável por 50 iterações"): baixa o critério de resíduo para 1e-12 e itera em blocos de 250 até a variação do monitor nas últimas 50 iterações ficar abaixo de `tol_mon` [m/s] ou somar `n_iter`. Registra o desvio em `_setup.json` (`continuacoes`, `iter_extra`, tempo somado), regrava `.cas`/`.dat` e reexporta. | ID, máximo de iterações extras, tolerância do monitor | nada |

### O que `setup()` faz, em ordem

1. Lê `mesh/<malha>.msh`.
2. CP periódico: `make-periodic` de `periodic_a` com `periodic_b` (pela TUI) e confere que as duas
   zonas viraram periódicas.
3. Viscoso **laminar**, energia **desligada**.
4. Material `fluido-tg` = cópia do ar com ρ e μ do caso, atribuído à zona de fluido.
5. Contorno por tipo:

   | Tipo | Contorno | Linha de medição | Reports extras |
   |---|---|---|---|
   | `cp_per` | `wall_top` móvel; gradiente periódico `dpdx` na direção (−1, 0) via TUI, conferido por `rpgetvar` | `line_x2p5`, (2,5; −1) → (2,5; 0) | `fx_top`, `fx_bottom` |
   | `cp_io` | `wall_top` móvel; `inlet` com o perfil `cp_inlet_ymenos.prof` (R06) ou com velocidade uniforme `uin` (R17); `outlet` a 0 Pa | `line_x2p5` | `p_in`, `p_out` (pressão média), `fx_top`, `fx_bottom` |
   | `ldc` | `lid` móvel | `line_x0p5`, (0,5; 0) → (0,5; 1) | — |
   | `cavr` | Frame Motion na zona fluida (`reference_frame.frame_motion`, `mrf_omega` = Ω, `reference_frame_axis_origin` = eixo); `lid` móvel relativa; demais paredes Stationary (paradas em relação à zona) | `line_x0p5` com `relative-x-velocity` | — |
   | `bfs` | `inlet` com o perfil `bfs_inlet.prof`; `outlet` a 0 Pa | `line_x2`, (2; 0) → (2; 1,5) | — |

6. Report File (`../dados/<nome>_rep.out`) e Report Plot com todos os reports `rd_*`.
7. Métodos: algoritmo (`alg`), gradiente Least Squares Cell Based, pressão (`pres`), momento (`mom`),
   sem High Order Term Relaxation. URFs de pressão e momento (exceto no Coupled).
8. Resíduos **absolutos** 1e-6 em continuity, x-velocity e y-velocity.
9. Inicialização standard com pressão e velocidades nulas.

---

## `rodar_lote.sh`

Script de Git Bash. Para cada ID passado (ou, sem argumentos, a lista padrão R01–R20 numa ordem
fixa), descobre a pasta da execução com `casos.caso`, cria `logs/` e roda
`.venv241/Scripts/python.exe -u tools/run_caso.py <ID>`, mandando toda a saída para
`logs/run.log`. Imprime no terminal o início e o `exit=` de cada execução e `LOTE FIM` no final.
Uma execução que falha não interrompe as seguintes. Já define `PYTHONIOENCODING=utf-8`.

## `status.py`

Painel de acompanhamento ao vivo: `python tools/status.py R26 R27 ...` redesenha a cada 5 s uma
linha por execução, lida do `logs/run.log`: ⏳ na fila (sem log), ▶️ rodando, ✅ feito (o log tem o
tempo de parede) ou ❌ erro, com a última iteração e os resíduos de continuity, u e v. Ctrl+C sai.

| Função | O que faz | Entrada | Saída |
|---|---|---|---|
| `linha(cid)` | Resume o estado de uma execução a partir do `run.log`. O traceback "Channel closed" que o `/exit` do Fluent 24.1 sempre deixa não conta como erro quando o log tem o fim normal. | ID | texto da linha |
| `main()` | Limpa a tela e redesenha o painel a cada 5 s até Ctrl+C. | IDs na linha de comando | nada |

---

## `analise_caso.py`

Lê os arquivos de `dados/`, calcula as métricas, faz os gráficos, escreve a nota e
reescreve o registro. Não abre o Fluent.

Uso:

```
python tools/analise_caso.py R02 [R03 ...]   # analisa estas
python tools/analise_caso.py --todos         # todas as de CASES que têm _setup.json
python tools/analise_caso.py --selftest      # confere as métricas no CSV de referência do R02
```

Com IDs ou `--todos`, primeiro calcula as métricas de todas (para a ordem observada da série de
malhas estar disponível), depois gera gráficos e notas de cada uma, e por fim reescreve o registro.

### Constantes

| Nome | O que guarda |
|---|---|
| `REG_BASE`, `REG_PYF` | `00_admin/registro_execucoes_v2.csv` (plano, só lido) e `00_admin/registro_pyfluent.csv` (reescrito) |
| `VAULT` | Pasta das notas: `notas/` |
| `FEITO`, `PENDENTE` | `✅` e `⏳`, o status da 1ª coluna do registro |
| `PREVISTO` | Previsão do FVM 1D de 2ª ordem para M1–M4: (ū, e_L2 [%], e∞ [m/s]) |
| `GHIA_CENTRO` | Centro do vórtice principal da cavidade, Re = 100 (Ghia, Tabela III): (0,6172; 0,7344) |
| `GHIA` | `refs/ghia1982.csv` lido com pandas: colunas `tabela` (I = u na vertical x = 0,5; II = v na horizontal y = 0,5), `pos`, `vel` |

### Leitura

| Função | O que faz | Entrada | Saída |
|---|---|---|---|
| `ler_campo(path)` | Lê o `_campo.csv` e tira os espaços dos nomes das colunas. | caminho do CSV | `DataFrame`, uma linha por célula |
| `ler_rep(path)` | Lê o `_rep.out`, só as linhas numéricas. | caminho | matriz `numpy` (n × 2): iteração, valor |
| `ler_residuos(csv_path, log_path)` | Extrai os resíduos das linhas de iteração que o Fluent imprime no `run.log` e nos `continuacao*.log` (o log é completo; o histórico do serviço de monitores vem com lacunas depois de uma continuação) e grava o `_residuos.csv`. | caminho do CSV e do log | `DataFrame` com `iter, continuity, x-velocity, y-velocity` |
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
| `calcular(cid)` | Lê campo, monitor, resíduos e setup; calcula as métricas gerais (células, iterações, convergência = todos os resíduos finais < 1e-6, resíduos finais, variação do monitor nas últimas 50 iterações) e as do tipo. Comparações automáticas: R05 contra R02; R12–R14 e R18–R20 contra R08; R18 e R20 contra R19. Grava `_metricas.json`. | ID | `dict` de métricas |
| `ordem_observada()` | Ordem observada p = ln(e_k/e_{k+1}) / ln 2 do e_L2 entre execuções vizinhas da série R01–R04, só para pares que já têm `_metricas.json`. | — | `dict`, ex. `{"R01→R02": 2.0, ...}` |

### Gráficos

| Função | O que faz | Entrada | Saída |
|---|---|---|---|
| `graficos(cid, m)` | Gera três PNGs em `graficos/` a partir do campo exportado (triangulação dos centros de célula; no degrau, mascara a região sólida x < 1, y < 0,5). `_campos.png`: \|V\| com vetores de direção (e linhas de corrente na cavidade e no degrau) e a pressão. `_convergencia.png`: resíduos e monitor (no CP periódico, com a linha ū = 4,3333). `_perfis.png`: perfis do tipo de caso (ver `docs/como_rodar.md`). | ID, métricas | lista dos PNGs gravados |

### Registro e notas

| Função | O que faz | Entrada | Saída |
|---|---|---|---|
| `resumo_metricas(cid, m)` | Texto curto `nome=valor` com as métricas principais, para a coluna `Metrica_integral` do registro e para o terminal. | ID, métricas | `str` |
| `registrar()` | Reescreve `registro_pyfluent.csv` a partir das 43 linhas de `registro_execucoes_v2.csv` (que não é alterado), com a coluna `Status` na frente. Execuções com `_metricas.json` recebem ✅ e os resultados (convergiu, iterações, tempo, URFs, pasta, e_L2, e∞, métricas, observações); as outras recebem ⏳ e ficam como no plano. | — | nada (grava o CSV) |
| `tabela_resultados(cid, m)` | Linhas da tabela Markdown "Grandeza \| Fluent \| Esperado \| Critério" da nota, conforme o tipo de caso (inclui a ordem observada na série R01–R04). | ID, métricas | lista de linhas |
| `nota(cid, m, pngs)` | Escreve a nota `<nome>.md` no vault (`VAULT`): o que a execução testa, parâmetros de entrada (textos de `TIPOS` + setup lido do Fluent), o que se espera, tabela de resultados, resíduos finais, gráficos (copiados para `VAULT/img/`) e onde estão os arquivos. | ID, métricas, PNGs | caminho da nota |
| `animacao(cid)` | Monta `graficos/<nome>_animacao.gif` com os quadros de `dados/anim/`: \|V\| em cores (escala do campo final) e vetores de direção, um quadro por bloco de iterações. | ID | caminho do GIF ou `None` |
| `analisar(cid)` | `calcular` + `graficos` + `animacao` + `nota`, e imprime o resumo. | ID | métricas |
| `selftest()` | Calcula `metricas_cp` no CSV de referência `refs/r02cpp_m2_sim_sou_campo.csv` (R02 feito pela interface) e confere com o guia: ū = 4,3499, e_L2 = 0,2693 %, e∞ = 0,01249 m/s; confere também que `ghia1982.csv` tem 17 + 17 pontos. | — | nada; `AssertionError` se divergir |

---

## `check_msh.py`

Verificador das malhas, independente do Fluent: lê o arquivo `.msh` **binário 2D** exportado pelo
Meshing e imprime um relatório.

Uso: `python -I tools/check_msh.py mesh/M2.msh [mesh/N040.msh ...]`. O tipo de conferência é
escolhido pelo nome do arquivo: M1–M4 e N040/N128 usam `main` (com as divisões esperadas);
DEG15, DEG16, BFY1 e BFYW usam `main_generic` (com nós e células esperados).

| Função | O que faz | Entrada | Saída |
|---|---|---|---|
| `parse(path)` | Decodifica o `.msh`: cabeçalhos globais, nomes das zonas, nós (coordenadas binárias), faces lineares e seções de células. | caminho | `dict` com `nodes`, `faces` (por zona), `cells`, `zones`, `hdr`; `SystemExit` se houver face não linear |
| `main(path, nx, ny)` | Confere uma malha retangular estruturada: nós, faces e células contra nx × ny, extents, número de divisões e espaçamento em x e y, cada zona de contorno (faces, extents, faces com 0 ou 2 células) e, no CP, se `periodic_a` e `periodic_b` têm os mesmos y. | caminho, divisões em x e y | nada (imprime) |
| `main_generic(path, exp_nodes, exp_cells, dxy)` | O mesmo para domínios não retangulares (degrau), sem conferência de divisões: nós e células esperados, extents, espaçamentos contra o alvo `dxy`, zonas. | caminho, nós e células esperados, espaçamento-alvo (0,025) | nada (imprime) |

## `msh_bordas_qualidade.py`

Complemento do `check_msh.py` (importa `parse` dele). Para uma malha 2D, imprime:

- para cada fronteira (menos as interiores): se é vertical ou horizontal, número de divisões,
  espaçamento mínimo, máximo, primeiro e último, e o maior salto entre espaçamentos vizinhos;
- para as células: quantas têm 4 nós, área total, área mínima e máxima, maior razão de aspecto,
  menor ângulo interno e quantas células têm ângulo < 30° e < 15°.

Uso: `python -I msh_bordas_qualidade.py <malha>.msh`, **rodando de dentro de `tools/`** (ele procura
o `check_msh.py` na pasta atual). Funções internas: `zone_nodes(z)` (coordenadas dos nós de uma zona
de faces) e `spacing(vals)` (diferenças entre valores consecutivos após ordenar).
