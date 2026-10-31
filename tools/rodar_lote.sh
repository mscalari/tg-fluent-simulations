#!/usr/bin/env bash
# Roda em sequência as execuções passadas (ou todas as rodáveis) e grava o log de cada uma em
# results/<bloco>/<nome>/logs/run.log. Uma falha não interrompe o lote.
# Uso (Git Bash, na raiz do repo): [PARALELO=N] bash tools/rodar_lote.sh [R02 R05 ...]
cd "$(dirname "$0")/.." || exit 1
PY=.venv241/Scripts/python.exe
export PYTHONIOENCODING=utf-8
IDS=${*:-R02 R05 R03 R04 R01 R06 R17 R08 R07 R09 R10 R12 R13 R14 R18 R19 R20 R15 R16 R11}
rodar() {
  dir=$($PY -I -c "import sys; sys.path.insert(0, 'tools'); from casos import caso; print(caso('$1')['dir'])")
  mkdir -p "$dir/logs"
  echo "=== $1 $(date +%H:%M:%S) ==="
  $PY -u tools/run_caso.py "$1" > "$dir/logs/run.log" 2>&1
  echo "$1 exit=$? $(date +%H:%M:%S)"
}
export -f rodar
export PY
# PARALELO=N roda N casos ao mesmo tempo (padrão 1); encerrar() só mata os PIDs da própria sessão.
printf '%s\n' $IDS | xargs -P "${PARALELO:-1}" -I{} bash -c 'rodar {}'
echo "LOTE FIM"
