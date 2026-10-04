#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
python3 code/freeze_r1.py
python3 code/prepare.py
python3 code/run_r1.py
python3 code/compare_r1.py
python3 code/freeze_r2.py
python3 code/run_r2.py
python3 code/freeze_r3.py
python3 code/run_r3.py
python3 code/freeze_r4.py
python3 code/run_r4.py
python3 code/freeze_r5.py
python3 code/run_r5.py
python3 code/freeze_r6.py
python3 code/run_r6.py
python3 code/freeze_r7.py
python3 code/run_r7.py
python3 code/run_r8.py
python3 code/run_r9.py
python3 code/gap_port.py
python3 code/verify.py
PYTHONNOUSERSITE=1 /usr/bin/python3 code/plot.py
