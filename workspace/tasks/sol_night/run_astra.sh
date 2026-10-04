#!/bin/bash
# proof_lane (lane-model_lane, xhigh) on a named question. Own runner so the Sun rotation is not disturbed.
# The model is lane_runner, not coordinator — no coordinator quota is consumed.
set -u
T=$1; D=$2; B=$3
W=
export PATH=local_config_path/bin:local_config_path/bin:/usr/local/bin:/usr/bin:/bin:$PATH
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4 NUMEXPR_NUM_THREADS=4 TOKENIZERS_PARALLELISM=false
mkdir -p "$W/tasks/build_night/logs" "$D"
P="$(cat "$W/tasks/build_night/COMMON.md" 2>/dev/null)

$(cat "$B")

Arbetskatalog: $W. Your catalog for this round: $D (tagg $T). Slutmeddelandet sparas som PROOF_LANE_FINAL_$T.md."
echo "{\"tag\":\"$T\",\"dir\":\"$D\",\"started\":\"$(date -Is)\"}" > "$W/tasks/build_night/logs/$T.start.json"
nice -n 10 lane_runner exec --skip-git-repo-check -m lane-model_lane -c model_reasoning_effort="xhigh" -c tools.web_search=true \
  --dangerously-bypass-approvals-and-sandbox -C "$W" -o "$D/PROOF_LANE_FINAL_$T.md" "$P" \
  > "$W/tasks/build_night/logs/$T.log" 2>&1
echo "{\"tag\":\"$T\",\"finished\":\"$(date -Is)\",\"rc\":$?}" > "$W/tasks/build_night/logs/$T.end.json"
