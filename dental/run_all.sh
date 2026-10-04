#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "${BASH_SOURCE[0]}")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1
python_bin="${DENTAL_PYTHON:-python3}"
"$python_bin" tools/scrub.py
"$python_bin" -m unittest discover -s tests -v
"$python_bin" -m dental_release.demo "$@"
# Positional arguments retain the original nine-profile interface.
cap_args=()
while (($#)); do
  if [[ "$1" == "--output" ]]; then
    cap_args+=(--output "$2")
    shift 2
  elif [[ "$1" == "--list" ]]; then
    cap_args+=(--list)
    shift
  else
    shift
  fi
done
exec "$python_bin" -m dental_release.capabilities "${cap_args[@]}"
