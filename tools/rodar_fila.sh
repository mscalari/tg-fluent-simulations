#!/usr/bin/env bash
# Fila com teto de Fluents simultâneos: cada Fluent compromete ~4 GB de memória virtual, e 8 ao mesmo tempo
# esgotaram os ~32 GB da máquina (31/10/2026). Conta as sessões vivas de qualquer lote, não só as desta fila.
# Uso (Git Bash, na raiz do repo): MAX=4 bash tools/rodar_fila.sh R38:c R40 R33:t ...
#   Rxx    roda do zero
#   Rxx:c  continua de .cas/.dat até os resíduos convergirem (teto de 16 000 iterações a mais)
#   Rxx:t  continua até o monitor variar menos que 1e-4 em 50 iterações (teto de 5000)
cd "$(dirname "$0")/.." || exit 1
PY=.venv241/Scripts/python.exe
export PYTHONIOENCODING=utf-8
MAX=${MAX:-4}

vivos() {
  tasklist | grep -c "fl_mpi2410.exe"
}

rodar() {
  id=${1%%:*}; modo=${1#*:}
  dir=$($PY -I -c "import sys; sys.path.insert(0, 'tools'); from casos import caso; print(caso('$id')['dir'])")
  mkdir -p "$dir/logs"
  echo "=== $1 $(date +%H:%M:%S) ==="
  case "$1" in
    *:c) $PY -u tools/run_caso.py --continuar "$id" 16000 > "$dir/logs/continuacao1.log" 2>&1 ;;
    *:t) $PY -u tools/run_caso.py --continuar "$id" 5000 1e-4 > "$dir/logs/continuacao1.log" 2>&1 ;;
    *)   $PY -u tools/run_caso.py "$id" > "$dir/logs/run.log" 2>&1 ;;
  esac
  echo "$1 exit=$? $(date +%H:%M:%S)"
}

for job in "$@"; do
  until [ "$(vivos)" -lt "$MAX" ]; do sleep 30; done
  rodar "$job" &
  sleep 90  # tempo para o Fluent novo aparecer na contagem
done
wait
echo "LOTE FIM"
