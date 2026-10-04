#!/usr/bin/env bash
set -euo pipefail
python3 code/extract.py
python3 code/r1.py
python3 code/r2.py
python3 code/r3.py
python3 code/r4.py
/usr/bin/python3 -s code/plot.py
python3 code/assemble.py
python3 code/verify.py
