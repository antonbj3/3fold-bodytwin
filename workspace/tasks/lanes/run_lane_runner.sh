#!/bin/bash
# Drive a BodyTwin-lane with lane runner (lane model, medium effort, full permission according to Anton 24/9). Kostnadsbesparande:
# one exec per round, at most MAX_ROUNDS rounds, continues only if RESULTS.md is missing.
# Usage: tasks/lanes/run_lane_runner.sh <LANE> [results-directory]   (default results/<LANE>)
set -u
LANE=${1:?}; ROOT=.
OUTDIR=${2:-results/$LANE}; OUT=$ROOT/$OUTDIR/RESULTS.md; LOG=$ROOT/tasks/lanes/codex_$LANE.log
MODEL=${MODEL:-lane-model}; EFFORT=${EFFORT:-medium}; MAX_ROUNDS=${MAX_ROUNDS:-3}
mkdir -p "$ROOT/$OUTDIR"; exec 9>"$ROOT/$OUTDIR/.lane_runner.lock"; flock -n 9 || { echo "$LANE is already running"; exit 0; }
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 NUMEXPR_NUM_THREADS=2
done_ok(){ [ -s "$OUT" ] && head -3 "$OUT" | grep -q "$LANE"; }
for i in $(seq 1 $MAX_ROUNDS); do
  done_ok && { echo "[$(date +%T)] $LANE klar" >> "$LOG"; exit 0; }
  echo "[$(date +%T)] $LANE round $i ($MODEL/$EFFORT)" >> "$LOG"
  timeout 7200 nice -n 10 lane_runner exec -m "$MODEL" -c model_reasoning_effort="$EFFORT" --dangerously-bypass-approvals-and-sandbox --skip-git-repo-check \
    -C "$ROOT" "You have a mandate from the coordinator (Anton) to perform the task without asking; choose reasonable defaults and account for them. Read and perform tasks/lanes/$LANE.md in its entirety. First check what already exists in $OUTDIR/ and build on it; do not redo finished steps. Full permission (Anton 24/9): networks, clouds (Modal/OVH via tasks/cloud_run.sh), and writing where the assignment requires it in the work surface (t.ex. results/<id>/-paket, tasks/lanes/bt_queue.txt med >>, notes/RESULTS_INDEX.md med >>). Egen utdata i $OUTDIR/, intermediate results at /media/anton/sdc1-tmp/bodytwin/$LANE/. Read-only: ~/projects/bodytwin and other sessions workspaces. No mail, pushes or publishing; do not touch credentials. Finish with $OUTDIR/RESULTS.md starting with the line '# $LANE'." \
    </dev/null >> "$LOG" 2>&1
done
done_ok && echo "[$(date +%T)] $LANE done" >> "$LOG" || echo "[$(date +%T)] $LANE NOT done after $MAX_ROUNDS" >> "$LOG"
