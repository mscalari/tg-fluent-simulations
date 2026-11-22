# Como rodar

## 1. Ambiente

O conversor CT_V1 (seção 7) usa o mesmo ambiente. Para os casos R01–R43, siga as seções 1 a 6.

| Item | Valor |
|---|---|
| Fluent | Ansys 2024 R1 (versão 24.1), instalado em `D:\ansys\ANSYS Inc\v241` |
| Variável de ambiente | `AWP_ROOT241 = D:\ansys\ANSYS Inc\v241` (o PyFluent a usa para achar o Fluent) |
| Python | 3.13, no ambiente virtual `.venv241` na raiz do repositório |
| PyFluent | `ansys-fluent-core==0.37.2`, a última versão que suporta o Fluent 2024 R1. A 0.43 **não** suporta. |
| Outros pacotes | `pandas`, `matplotlib` (e `numpy`, que vem junto) |
| Codificação | `PYTHONIOENCODING=utf-8`, senão o console do Windows (cp1252) falha ao imprimir "−" |

### Recriar o `.venv241`

No PowerShell, na raiz do repositório:

```powershell
py -V:3.13 -m venv .venv241
.venv241\Scripts\python.exe -m pip install ansys-fluent-core==0.37.2 pandas matplotlib
```

### Conferir sem abrir o Fluent

```bash
export PYTHONIOENCODING=utf-8
.venv241/Scripts/python.exe tools/casos.py                     # "casos ok: 43"
.venv241/Scripts/python.exe tools/analise_caso.py --selftest   # "selftest ok: ..."
```

- `casos.py` confere que cada execução de `CASES` monta os parâmetros, cai na pasta do bloco certo e
  tem a malha em `mesh/`.
- `--selftest` calcula ū, e_L2 e e∞ no CSV de referência do R02
  (`refs/r02cpp_m2_sim_sou_campo.csv`) e confere com a referência (ū = 4,3499 m/s, e_L2 = 0,2693 %,
  e∞ = 0,01249 m/s).

## 2. Fluxo de trabalho

Comandos no Git Bash, na raiz do repositório.

**1. Rodar no Fluent.** Um lote, em sequência:

```bash
bash tools/rodar_lote.sh R02 R05 R03     # só estas
bash tools/rodar_lote.sh                 # lista padrão (R01–R20)
PARALELO=4 bash tools/rodar_lote.sh R31 R32 R33 R34 R28 R29 R30   # até 4 casos ao mesmo tempo
```

Cada Fluent é serial; `PARALELO` só decide quantos rodam juntos. `encerrar()` mata apenas os PIDs da
própria sessão, então um caso não derruba o outro. Dois lotes não podem escrever no mesmo arquivo de log.

O script já usa `.venv241` e `PYTHONIOENCODING=utf-8`. A saída de cada execução vai para
`results/<bloco>/<nome>/logs/run.log`; no terminal aparece só `=== R02 hh:mm:ss ===` e
`R02 exit=0 hh:mm:ss` (exit diferente de 0 = falhou; veja o `run.log`). Uma falha não para o lote.

Para limitar os Fluents simultâneos ou retomar casos de `.cas/.dat`, use a fila (`Rxx:c` continua até
os resíduos convergirem; `Rxx:t`, até o monitor estabilizar):

```bash
MAX=4 bash tools/rodar_fila.sh R38:c R40 R33:t
```

Cada Fluent compromete cerca de 4 GB de memória virtual: oito simultâneos esgotaram os 32 GB da máquina.

Uma execução só, com a saída no terminal:

```bash
.venv241/Scripts/python.exe tools/run_caso.py R02
```

Refazer só a exportação (relê o `.cas.h5`/`.dat.h5` já gravado, sem resolver de novo):

```bash
.venv241/Scripts/python.exe tools/run_caso.py --exportar R02
# relê .cas/.dat e itera mais; com o 3º número, até o monitor variar menos que ele [m/s] em 50 iterações
.venv241/Scripts/python.exe tools/run_caso.py --continuar R04 20000 1e-6
```

Os casos turbulentos gravam `.cas/.dat` e `_setup.json` a cada 2000 iterações (`BLOCO_SALVA`). Se o PC
reiniciar no meio, `--continuar <ID> <iterações>` retoma do último bloco, sem recomeçar.

Antes do degrau turbulento, `refs/bfy_inlet.prof` precisa existir; para recriá-lo:
`.venv241/Scripts/python.exe tools/gerar_bfy_inlet.py`.

**2. Analisar.**

```bash
.venv241/Scripts/python.exe tools/analise_caso.py R02 R05    # estas
.venv241/Scripts/python.exe tools/analise_caso.py --todos    # todas que têm _setup.json
```

Para cada execução: calcula as métricas (`_metricas.json`) e gera os três gráficos.

Antes de rodar o Fluent, confira que não ficaram processos `fluent.exe` / `cx2410.exe` de execuções
anteriores (Gerenciador de Tarefas).

## 3. O que acontece numa execução

`run_caso.py <ID>`:

1. Lê os parâmetros da execução em `tools/casos.py`.
2. Cria `results/<bloco>/<nome>/{fluent,dados,graficos,logs}` e abre o Fluent 2D (3D no duto), dupla
   precisão, serial, sem interface, com diretório de trabalho em `fluent/`.
3. Lê a malha; no CP periódico e no duto liga `periodic_a`/`periodic_b` (make-periodic).
4. Define laminar ou o modelo de turbulência do caso, energia desligada, material `fluido-tg` (cópia
   do ar com ρ e μ do caso).
5. Aplica o contorno do tipo de caso (parede móvel, gradiente periódico, entrada/saída, perfis).
6. Cria a linha de medição e o monitor de u médio nela, os reports de força/pressão, os métodos
   (algoritmo, gradiente Least Squares Cell Based, esquemas de pressão e momento), URFs, resíduos
   absolutos 1e-6 e inicialização com zeros (nos casos turbulentos, com os valores de `ini` do tipo).
7. No CP periódico e no duto roda 50 iterações e confere que ū é positivo (se não for, para com erro).
8. Itera até o teto (`iter_max`); o Fluent para sozinho quando os resíduos ficam abaixo de 1e-6.
   Nos casos turbulentos, grava `.cas/.dat` a cada 2000 iterações.
9. Grava caso e dados, exporta o campo, as paredes e os reports, grava o `_setup.json`.
10. Fecha o Fluent (inclusive os processos que a 24.1 deixa abertos) e move as transcrições para `logs/`.

## 4. O que cada arquivo de saída contém

Em `results/<bloco>/<nome>/` (`<nome>` = por exemplo `R02_CPP_M2_SIM_SOU`):

| Arquivo | Quem gera | Conteúdo |
|---|---|---|
| `fluent/<nome>.cas.h5` | `run_caso.py` | Caso do Fluent: malha + setup completo. Pode ser aberto na interface do Fluent. |
| `fluent/<nome>.dat.h5` | `run_caso.py` | Solução (campos) da última iteração. |
| `dados/<nome>_campo.csv` | `run_caso.py` | Uma linha por célula, valores no **centro da célula**. Colunas: `cellnumber, x-coordinate, y-coordinate, pressure, total-pressure, velocity-magnitude, x-velocity, y-velocity, stream-function`. Unidades SI. As linhas **não** estão em ordem geométrica. Em 3D, `z-coordinate` e `z-velocity` no lugar de `stream-function`. Nos casos turbulentos, também `turb-kinetic-energy`, `viscosity-turb`, `viscosity-ratio`, `cell-volume` e, no RSM, as tensões de Reynolds (`vw` e `uw` só em 3D). |
| `dados/<nome>_parede_<zona>.csv` | `run_caso.py` | Só em `duto` e `bft`, um arquivo por parede: `pressure`, `x-wall-shear`, `y-wall-shear`, `wall-shear`, `y-plus` (e `z-wall-shear` em 3D). |
| `dados/<nome>_rep.out` | Fluent (Report File) | Histórico do monitor: iteração e u médio por área na linha de medição (`rd_line_x2p5` no CP, `rd_line_x0p5` na cavidade, `rd_line_x2` no degrau). No duto, `rd_line_ub` (v médio em z = 0) e `rd_line_vs` (u na bissetriz do canto); no degrau turbulento, `rd_tau` (τ_x médio em `wall_down_bottom`). Três linhas de cabeçalho. |
| `dados/<nome>_residuos.csv` | `analise_caso.py`, a partir dos logs | `iter` e uma coluna por equação, com os nomes do cabeçalho do Fluent: resíduos por iteração. |
| `dados/<nome>_setup.json` | `run_caso.py` | O setup **lido de volta do Fluent** (não o pedido): malha, zona de fluido, solver, modelo viscoso, energia, ρ, μ, paredes, métodos, URFs, critérios, monitores, inicialização, teto de iterações, tempo de parede, versão do PyFluent e valores finais dos reports (`fx_top`, `fx_bottom`, `p_in`, `p_out`, conforme o caso). |
| `dados/<nome>_metricas.json` | `analise_caso.py` | Métricas: nº de células, iterações, convergiu (todos os resíduos < 1e-6; no duto, sem a continuidade), resíduos finais, variação do monitor nas últimas 50 iterações, e as do tipo de caso (ver abaixo). |
| `graficos/<nome>_campos.png` | `analise_caso.py` | Em cima: módulo da velocidade em cores com vetores de direção (e linhas de corrente na cavidade e no degrau). Embaixo: pressão. |
| `graficos/<nome>_convergencia.png` | `analise_caso.py` | Resíduos (escala log) com a linha do critério 1e-6, e o histórico do monitor. |
| `graficos/<nome>_perfis.png` | `analise_caso.py` | CP: perfil em x ≈ 2,5 contra a analítica e o erro ponto a ponto. Cavidade: u na vertical x ≈ 0,5 e v na horizontal y ≈ 0,5. Degrau: u na primeira linha de células acima do fundo (com o ponto de recolamento) e perfis de u em três posições. Duto: v axial em z ≈ 0 e \|V_s\|/U_b na diagonal. Degrau turbulento: Cf e Cp contra Driver e Seegmiller. |
| `graficos/<nome>_animacao.gif` | `analise_caso.py` | Só nas rodadas com `anim_cada` (bloco G rotativo): evolução de \|V\| e dos vetores ao longo das iterações, um quadro a cada `anim_cada` iterações (quadros em `dados/anim/`). |
| `logs/run.log` | `rodar_lote.sh` | Tudo o que o `run_caso.py` imprimiu, incluindo `SETUP EFETIVO`, as iterações do Fluent e erros. |
| `logs/*.trn` | Fluent | Transcrição da sessão do Fluent (o que o console do Fluent mostrou). |

Métricas por tipo de caso (`_metricas.json`):

| Tipo | Métricas |
|---|---|
| CP periódico | `x0` (x da coluna usada), `ubar`, `eL2` [%], `einf` [m/s], `fx_top`, `fx_bottom`, `balanco` (soma das forças nas paredes); no R05, `max_du_R02` |
| CP entrada/saída | as mesmas de perfil, mais `p_in`, `p_out`, `dp`, `fx_top`, `fx_bottom` |
| Cavidade | `xc`, `yc` (centro do vórtice), `dist_ghia`, `u_min_vert`, `v_max_hor`, `v_min_hor`, `u_max`, `ghia_du_max`/`ghia_dv_max` (erro máximo contra Ghia nas linhas centrais); e `max_du_R08` / `max_du_R19` nas comparações entre execuções |
| Cavidade rotativa | no referencial rotativo: `u_max`, `cb` (checkerboard), `u_max_lin`/`u_min_lin` (horizontal y = 0,5) e, contra o R26, `eL2_u`, `eL2_v` [%], `einf_u`, `einf_v`, `overshoot`, `subvalor`; no R27, `ref_dif` [% de max\|V\|] |
| Degrau | `xr` (recolamento medido a partir do degrau, x = 1) e `xr_h` (x_r / h); se a bolha atravessa a saída, `xr_maior_que` (cota inferior) |
| Duto | `ub`, `vs_max`, `vs_ub` [%], `re_b`, `cf` (2τ/ρU_b², τ = 1 Pa pelo balanço de forças), `yplus_max`, `yplus_med` |
| Degrau turbulento | `xr`, `xr_h`, `erro_xr` [%] contra 6,26, `mut_max`, `mut_mu_max` (janela 0 < x < 12, 0,2 < y < 2,5), `yplus_max`, `yplus_med`, `cf_rms`, `cp_step_rms`, `cp_top_rms` |

## 5. Como adicionar uma execução

1. Confira que a malha está em `mesh/<malha>.msh`.
2. Acrescente o nome-base em `NOMES` e uma linha em `CASES`, em `tools/casos.py`, só com o que muda em
   relação aos padrões:

   ```python
   "R44": dict(tipo="ldc", malha="N040", pres="presto!", bloco="G",
               muda="Texto explicando o que esta rodada testa."),
   ```

   Chaves possíveis: `tipo` (obrigatória), `malha` (obrigatória), `muda` (obrigatória, o que a
   rodada testa), `alg`, `mom`, `pres`, `urf_p`, `urf_mom`, `iter_max`, `perfil`, `uin`, `bloco`, `anim_cada` (exporta um quadro a cada N iterações e gera o GIF da evolução), `turb` (dicionário de atributos de `setup.models.viscous`, como `SST` ou `RSM` em `casos.py`), `canto`, `curvatura`, `csep`, `procs`. Os valores
   aceitos em `alg`/`mom`/`pres` são os nomes do Fluent, por exemplo `"SIMPLEC"`, `"quick"`,
   `"presto!"`. Se o ID não estiver em `BLOCO_ID`, informe `bloco`.
3. Se for uma física nova (outro tipo de domínio ou contorno), crie também uma entrada em `TIPOS`
   (ρ, μ, teto de iterações e os textos `geometria`, `fluido`, `contorno`, `esperado`) e trate o novo
   `tipo` em `setup()` de `run_caso.py` e em `calcular()`, `graficos()` e `resumo_metricas()`
   de `analise_caso.py`.
4. Rode `python tools/casos.py` para conferir.
5. Se a execução deve entrar no lote padrão, acrescente o ID na lista `IDS` de `tools/rodar_lote.sh`.

## 6. Armadilhas conhecidas

- **Padrões da sessão nova.** Uma sessão nova do Fluent 2024 R1 aberta pelo PyFluent começa com
  acoplamento Coupled, modelo k-omega, inicialização híbrida e resíduos 1e-3, diferente do que se
  espera pela interface. Por isso o runner define **tudo** explicitamente (laminar, algoritmo,
  inicialização standard, critério 1e-6) e grava o que ficou no `_setup.json`.
- **Gradiente de pressão periódico só pela TUI.** A API de settings da 24.1 não expõe o gradiente
  periódico; o runner usa o comando de texto
  `/define/periodic-conditions/pressure-gradient-specification 4 -1 0` e confere o valor lido de
  volta. A direção é **(−1, 0)**: com (1, 0) o ū sai negativo. No duto vale o mesmo: `4 0 -1 0`
  para escoamento em +y.
- **Continuidade fora do critério no duto.** O campo inicial já satisfaz a continuidade, então o fator
  de escala do resíduo fica em ~1e-16 e o valor escalado não mede nada (no R34 ficaria preso em ~0,1).
  O runner desliga `check_convergence` só nessa equação, só no duto.
- **Tratamento de parede.** k-ε e RSM vêm com função de parede, inválida com y+ ≈ 1: o runner liga
  Enhanced Wall Treatment. Os modelos de ω integram até a parede. No RSM isso é desvio do Adendo.
- **Esquema da turbulência.** O padrão do Fluent muda com o modelo (First Order em k-ε e RSM, SOU em
  SST): o runner força SOU em todas as equações de turbulência.
- **C_SEP do GEKO.** Não aparece em `geko_options`; é a variável `geko-csep`, definida por `rpsetvar`.
- **O Fluent não fecha sozinho.** Na 24.1 com a 0.37.2, `session.exit()` não encerra os processos.
  A função `encerrar()` manda `/exit yes`, chama `exit()`, espera e, se ainda houver processos,
  mata a árvore do Fluent pelo PID. Mesmo assim, confira se não sobraram `fluent.exe`/`cx2410.exe`.
- **Geometria do CP invertida.** O canal está em y ∈ [−1, 0] e a parede **móvel** (`wall_top`) fica
  em **y = −1**. A comparação com a analítica usa d = −y (distância à parede parada). Com y + 1 o
  erro aparece como ~25 %.
- **Exportação de todas as células.** O campo é exportado com `surface_name_list=[]` (lista vazia =
  todas as células do domínio), em Cell Center, separado por vírgula.
- **Ordem das linhas do CSV.** As linhas do `_campo.csv` não seguem a geometria; o código sempre
  filtra e ordena pelas coordenadas. Leia as colunas pelo nome, não pela posição.
- **make-periodic a cada malha.** A periodicidade é refeita pelo runner sempre que lê a malha do CP
  (a malha sai do Workbench com `periodic_a`/`periodic_b` como paredes).
- **Console em cp1252.** Sem `PYTHONIOENCODING=utf-8`, imprimir "−" ou "ū" quebra o script.

## 7. Conversor CT_V1 (MRF)

Comandos no Git Bash, na raiz, com `export PYTHONIOENCODING=utf-8` e `.venv241/Scripts/python.exe`. O conversor é 3D,
com ~1,9 M de células: **uma sessão do Fluent por vez** (~9,4 GB de RAM na G0, numa máquina de 16 GB). Não rode
outro Fluent enquanto uma fila estiver ativa (ver o painel em `monitor/`).

### 7.1 Rodar no Fluent

```bash
python tools/conversor_mrf.py setup                          # monta o caso e testa as imagens
python tools/conversor_mrf.py sentido 500                    # SR = 0 em +x e −x (escolhido: +x)
python tools/conversor_mrf.py varredura 1                    # SR 0 … 0,9, cada ponto parte do anterior
python tools/conversor_mrf.py varredura 1 0.9 0.95 1.0 --ini results/conversor/CT_V1_MRF/SR090/fluent/<...>.dat.h5
python tools/conversor_mrf.py sensibilidade                  # sst_sem_cc ke_realizable fou meio_passo
python tools/conversor_mrf.py sensibilidade coupled --malha mesh/CT_V1_3zonas_G2.msh.h5 --saida CT_V1_GCI_G2
python tools/conversor_mrf.py refazer_figuras 0 0.2          # uma sessão por chamada; use poucos SR por vez
python tools/malha_gci_conversor.py G1 --so-superficie       # calibração da superfície (~30 min)
python tools/malha_gci_conversor.py G1 --retomar-superficie  # volume a partir da superfície gravada
python tools/malha_gci_conversor.py G2                       # do STL
python tools/conversor_mrf.py varredura 1 0 0.4 --malha mesh/CT_V1_3zonas_G1.msh.h5 --saida CT_V1_GCI_G1     --urf mom=0.5,k=0.6,omega=0.6,turb-viscosity=0.8          # URF menores: a G1 divergiu com as padrão
python tools/imagens_malha_conversor.py G2 G1 G0             # 10 vistas por malha + montagens G2 | G1 | G0
python tools/imagens_malha_conversor.py --so-montagem        # só refaz as montagens a partir dos PNGs
```

`--urf chave=valor,...` aplica sub-relaxações extras em `solution.controls.under_relaxation` logo depois de escolher o
SIMPLE (vale para qualquer comando). `iterar_ponto` para o ponto com `RuntimeError` se um torque não for finito, passar
de `T_DIVERGE` (10⁶ N·m) ou se, depois de 2 blocos, |soma dos torques| > 0,5·|T_b| (guarda de divergência; antes a G1
divergiu e gastou horas iterando lixo).

`--malha` e `--saida` valem para qualquer comando e mantêm os resultados da malha base intactos. Cada ponto grava
`SRxxx/dados/resultado.json` com as médias da janela, o desvio e `conv`; os `.cas/.dat` ficam em `SRxxx/fluent/`.

Duração típica na G0: ~4,5 s/iteração com SIMPLE; o teto é de 3000 iterações no 1º ponto e 2000 nos demais. O Fluent
roda sem GUI (`UIMode.NO_GUI`), inclusive para as imagens.

**Lançar filas longas.** Processos iniciados pelo Bash morrem quando o Claude Code reinicia (a fila 3 caiu assim). Para
uma fila que deve sobreviver, lance desacoplada pelo WMI e registre o estado num arquivo de status:

```powershell
Invoke-CimMethod -ClassName Win32_Process -MethodName Create -Arguments @{
  CommandLine = '"C:\Program Files\Git\bin\bash.exe" -l "C:\caminho\da\fila.sh"';
  CurrentDirectory = 'C:\Users\murilo\Repositorios\03 - TG\tg-fluent-simulations' }
```

A fila começa com `cd` para a raiz do repositório. Ao gravar o código de saída no status, guarde `rc=$?` numa
variável **antes** de montar a linha: em `echo "[$(date +%T)] passo rc=$?"` o `$(date)` zera o `$?` e o status
mostra sempre `rc=0` (aconteceu na fila 6; confira o resultado pelo log do passo).

**Disco.** Cada ponto da G0 grava ~700 MB de `.cas/.dat` (G1 ~340 MB, G2 ~190 MB). Em 22/11/2026 o C: (223 GB)
encheu e o Fluent caiu no `write_case_data` (HDF5-DIAG, "Connection reset"), perdendo o resultado da G1. Os casos que
as filas não usam mais foram movidos (não apagados) para `D:\TG_backup_fluent\CT_V1_MRF\` com a mesma estrutura
(`robocopy ... *.cas.h5 *.dat.h5 /MOV`); `--ini` aceita o caminho no D:. As filas conferem o espaço livre antes de cada
passo e pulam o passo com menos de 3 GB.

Entre passos da fila, confira e limpe `fluent.exe`/`cx2410.exe`/`fl_mpi2410.exe` órfãos pela árvore de PIDs da
sessão, nunca por nome.

### 7.2 Analisar

```bash
python tools/analise_conversor.py                  # desempenho.csv, g1–g3, nota do Obsidian (--sem-nota: só gráficos)
python tools/euler_interfaces_conversor.py         # g7 (aceita SR em %: 0 20 60)
python tools/analise_yplus_conversor.py            # g4
python tools/analise_sensibilidade_conversor.py    # g6
python tools/comparacao_literatura_conversor.py    # g5
python tools/tendencia_derivada_conversor.py       # g8
python tools/tendencia_sinais_conversor.py         # g9
python tools/tendencia_colapso_conversor.py        # g10
python tools/animacao_meio_passo.py                # GIF do meio passo
python tools/gci_conversor.py                      # GCI (diz o que falta se G1/G2 não existirem)
python tools/comparacao_blocos.py                  # R01–R43 contra as referências (results/comparacao_blocos/)
```

Entradas e saídas de cada script: `docs/codigo.md`. As análises dependem de `desempenho.csv` (gerado por
`analise_conversor.py`, que deve rodar primeiro) e dos report files `dados/*.out`.

### 7.3 Estrutura das saídas

```
results/conversor/CT_V1_MRF/
  SR000 … SR090/{fluent,dados,graficos}   um ponto de operação por pasta (SR040 = SR 0,4)
  dados/        report files (*.out), desempenho.csv, sensibilidade.csv, yplus_faixas.csv, euler_interfaces*.csv,
                derivada_tendencia.csv, sinais_tendencia.csv, colapso_adimensional*.csv, setup.json
  graficos/     g1–g10, animacao_meio_passo.gif
  sens/<variante>/   sst_sem_cc, ke_realizable, fou, meio_passo
results/conversor/CT_V1_GCI/            logs e parâmetros das malhas G1/G2, dados/gci.csv, graficos/gci_tr.png
results/conversor/CT_V1_GCI_G1|G2/      rodadas de SR 0 e 0,4 nas malhas G1 e G2 (mesma estrutura da base)
results/comparacao_blocos/              figuras-resumo de R01–R43 e resumo.csv
```

`.cas/.dat`, `.h5` e os exports por face (`SR*/dados/*.csv`) ficam fora do git (GB). O `.gitignore` também exclui
`monitor/`, as malhas `mesh/*.msh.h5` e os `*_workflow_files` (exceto o `.wft`).

### 7.4 Lições do Fluent 24.1 sem GUI

- **Câmera.** Só funciona `views.restore_view("front")` seguido de `camera.orbit(right=..., up=...)`. Ajustar
  `camera.position/target` fora do eixo z ou chamar `views.auto_scale()` deixa a imagem em branco, só com a legenda.
  `orbit(right=−90)` leva a câmera a +x (o rótulo `de_mais_x` das imagens segue isso; havia rótulos trocados, já
  corrigidos). `camera(st, direita, cima)` e `VISTAS` em `conversor_mrf.py` aplicam isso.
- **Imagens da varredura em branco.** Durante as iterações os plots de monitor prendem a janela gráfica. Por isso as
  figuras vêm de `refazer_figuras`, que relê o `.cas/.dat` numa sessão nova.
- **`st.file.write_mesh` não existe** nos settings da 24.1. Para girar uma zona (variante `meio_passo`), use
  `mesh.modify_zones.rotate_zone(zone_names, rotation_angle, origin, axis)` logo depois de `read_mesh` e **antes**
  de criar as interfaces (elas alteram a malha), sem gravar malha intermediária.
- **`.dat` não volta num setup recriado do zero.** Criar interfaces altera a malha, e o Fluent recusa ler o `.dat`
  ("The mesh has been modified"). Para retomar entre sessões, use `read_case_data` do caso gravado junto.
- **Nome das cell zones depois do merge no Meshing.** `/mesh/manage/merge "fluido_turbina*"` mantém o nome de uma
  das zonas juntadas: na G0 saiu `fluido_turbina`, nas malhas do GCI saiu `fluido_turbina_18` (e
  `interior--fluido_turbina_18`). `conversor_mrf.normalizar_zonas` renomeia no solver logo depois de `read_mesh`:
  `cell_zone_conditions.fluid.rename(novo, velho)` para as fluid e `boundary_conditions.set_zone_name(zonename=,
  newname=)` para as interior (este último recusa as fluid: "api-rename-obj: Object is invalid").
- **Paredes extras das malhas do GCI.** A G2 (e a G1) sai do Meshing com paredes `stlb_<...>.<n>...fluido_<zona>`
  (pás de uma peça tocando o fluido de outra). `renomear_paredes_stl` as renomeia para `parede_<dono>[_k]_lado_<lado>`
  (dono pelo sufixo `.n`: 1 bomba, 2 estator, 3 turbina) e `sincronizar_paredes_abs` refaz `PAREDES_ABS` a partir das
  paredes da malha/caso lido, para que elas girem com a peça dona e entrem no torque dela.
- **Zonas interior das interfaces** (`if_bt`, `if_te`, `if_eb`) são recusadas como superfície de report e de
  exportação ("Values contain disallowed entries"). Use os dois lados de cada interface (`_euler_<lado>.csv`).
- **O TUI `/define/mesh-interfaces/create` entra em laço infinito.** Use `setup.mesh_interfaces.interface.create(nome)`
  e `interface[nome].zone_names = [lado1, lado2]`.
- **Report file não é sobrescrito** se já tiver linhas além da iteração atual: o Fluent cria `<nome>_<it+1>.out` e o
  script continua lendo o antigo (a sensibilidade travou assim). Use um report file por variante
  (`CT_V1_MRF_sens_<variante>_<data>.out`); `iterar_ponto` levanta erro se um bloco não acrescentar linhas.
- **O `.out` do report file fica aberto pelo Fluent** (`Path.replace` falha com WinError 32): grave o caso antes e
  copie o `.out`.
- **Trocar para k-ε põe k e ε em primeira ordem sozinho.** Depois de mudar o modelo, force `second-order-upwind` em
  tudo, exceto na pressão (a variante `fou` é a exceção).
- **Momento em N·m.** Na 24.1 o report de momento oferece `['Drag Coefficient', 'Drag Force']`; `Drag Force` é o
  momento em N·m e vale 0,6125 × Cm (valores de referência padrão: ρ = 1,225, v = A = L = 1). Escolha a opção sem
  "coef" (`momento_em_nm`).
- **`expr:lx` e velocidade relativa.** O campo `expr:lx` (named expression) foi recusado como campo de report;
  `rel-tangential-velocity` é relativa à zona de referência (ω_b) em todas as zonas.
- **Coupled não cabe na G0** (~17,6 GB de memória privada em 16 GB; swap a ~15 s/iteração). Por isso SIMPLE; o
  Coupled só roda na G2 (variante `coupled`).
- **Travamentos de imagem.** O Fluent já deu erro interno do Scheme (`wta(1st) to string->symbol`) e travou numa
  figura do SR 0,9; e caiu ("disconnected") na exportação Euler depois de 4 pontos na mesma sessão. Grave
  `resultado.json` antes das figuras e use uma sessão por ponto em `refazer_figuras`.
- **Checkpoint nativo.** Criar `check-fluent.txt` no diretório de trabalho faz o Fluent gravar `.cas/.dat` na
  iteração seguinte sem parar.
- **`session.exit()` não encerra o Fluent 24.1**: use `run_caso.encerrar()` e confira que não ficaram
  `fluent.exe`/`cx2410.exe`/`fl_mpi2410.exe`. Rode sempre com `PYTHONIOENCODING=utf-8`.
- **Torque nas paredes `sobra_*`.** As faces sem par das interfaces viram paredes que giram com a peça do outro lado
  e entram no torque dela; no SR000 o maior efeito é 1,3 % de |T_b|.
