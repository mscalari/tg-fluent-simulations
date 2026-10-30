# tg-fluent-simulations

Automação, com **PyFluent**, das simulações de validação do Trabalho de Graduação (TG) no
**Ansys Fluent 2024 R1** (Windows). Os scripts montam o caso no Fluent, resolvem, exportam os
resultados e calculam as métricas em Python, sem abrir a interface gráfica.

Casos cobertos:

- **Couette-Poiseuille (CP)** periódico e com entrada/saída, comparado com a solução analítica;
- **cavidade com tampa móvel** (lid-driven cavity), Re = 100: esquemas de momento, algoritmos de
  acoplamento pressão-velocidade e interpolação da pressão;
- **degrau** (backward-facing step) laminar.

O plano completo tem 43 execuções (R01–R43) em 8 blocos (A–H). Hoje rodam R01–R27 (a cavidade
rotativa R21–R27 roda em 2D, com Frame Motion); o duto (R28–R34) e o degrau turbulento (R35–R43)
aguardam as malhas DU040, BFY1 e BFYW.

## Começo rápido

No Git Bash, na raiz do repositório:

```bash
export PYTHONIOENCODING=utf-8
.venv241/Scripts/python.exe tools/casos.py                     # confere a tabela de casos e as malhas
.venv241/Scripts/python.exe tools/analise_caso.py --selftest   # confere as métricas (não abre o Fluent)
bash tools/rodar_lote.sh R02                                   # roda o R02 no Fluent
.venv241/Scripts/python.exe tools/status.py R02               # acompanha ao vivo (outro terminal)
.venv241/Scripts/python.exe tools/analise_caso.py R02          # métricas, gráficos, registro e nota
```

## Documentação

| Arquivo | Conteúdo |
|---|---|
| [docs/visao_geral.md](docs/visao_geral.md) | O que o projeto valida, blocos e execuções, organização das pastas |
| [docs/como_rodar.md](docs/como_rodar.md) | Ambiente, comandos, arquivos de saída, como adicionar uma execução, armadilhas |
| [docs/codigo.md](docs/codigo.md) | Cada módulo e cada função de `tools/`, em linguagem simples |

Documentos de planejamento (em `00_admin/`): `handoff_pyfluent.md` (o que os scripts devem fazer),
`guia_validado_fluent2024R1.md` (o procedimento que funcionou pela interface gráfica, que prevalece
sobre o Roteiro PDF) e `pre_requisitos_geometria_malhas.md` (como as malhas foram feitas).
