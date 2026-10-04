#!/bin/bash
# Heavy local job under shared resource management (restart 2026-09-22).
# Usage: tasks/heavy_run.sh <ram_gb> <vram_gb> <command...>
# Acquire /tmp/tenancy/bigmem.lock (all sessions, one heavy job at a time), then /tmp/gpu.lock if vram>0,
# run resource_gate (exit 1 = wait 60 s and try again), thread cap 2, run the command in its own systemd scope with MemoryMax.
set -u
RAM=$1; VRAM=$2; shift 2
mkdir -p /tmp/tenancy
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 NUMEXPR_NUM_THREADS=2 CCX_NPROC_STIFFNESS=2 CCX_NPROC_EQUATION_SOLVER=2
MEMMAX=$(python3 -c "print(int(float('$RAM')*1.5+1))")G
run() {
  until python3 source_repository/scripts/resource_gate.py check --ram-gb "$RAM" --vram-gb "$VRAM" --max-load-per-core 1.0 >/dev/null; do
    echo "[heavy_run] resource_gate: wait" >&2; sleep 60; done
  exec systemd-run --user --scope --quiet -p MemoryMax=$MEMMAX nice -n 10 "$@"
}
export -f run 2>/dev/null
if [ "$(python3 -c "print(float('$VRAM')>0)")" = True ]; then
  exec flock /tmp/tenancy/bigmem.lock flock /tmp/gpu.lock bash -c "$(declare -f run); RAM=$RAM VRAM=$VRAM MEMMAX=$MEMMAX run \"\$@\"" _ "$@"
else
  exec flock /tmp/tenancy/bigmem.lock bash -c "$(declare -f run); RAM=$RAM VRAM=$VRAM MEMMAX=$MEMMAX run \"\$@\"" _ "$@"
fi
