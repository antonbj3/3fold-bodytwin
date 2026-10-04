#!/usr/bin/env bash
set -euo pipefail
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 NUMEXPR_NUM_THREADS=2 PYTHONDONTWRITEBYTECODE=1
export PYTHONPATH="$(dirname "$0")${PYTHONPATH:+:$PYTHONPATH}"
exec python -m surgical_chain "$@"
