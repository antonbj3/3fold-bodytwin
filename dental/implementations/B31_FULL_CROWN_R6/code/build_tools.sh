#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
D6=@DENTAL_WORK_ROOT@/PROOF_LANE_FULL_CROWN_R6
for part in intersections repair repair_local repair_protected; do
  opt=-O0
  if [[ "$part" == intersections ]]; then opt=-O1; fi
  /usr/bin/time -v g++ "$opt" -std=c++17 -I"$D6/deps/usr/include" -I"$D6/deps/usr/include/x86_64-linux-gnu" -I/usr/include/eigen3 "code/$part.cpp" -o "$D6/$part" -L"$D6/deps/usr/lib/x86_64-linux-gnu" -lmpfr -lgmp 2> "raw/rebuild_${part}.log"
done
