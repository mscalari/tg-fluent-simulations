# Como rodar

## 1. Ambiente

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
