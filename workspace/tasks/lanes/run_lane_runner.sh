#!/bin/bash
# Run a BodyTwin lane with lane_runner (build-lane model, medium effort, full permissions according to Anton 24/9). Cost-saving:
# one exec per round, at most MAX_ROUNDS rounds, continues only if RESULTS.md is missing.
# Usage: tasks/lanes/run_lane_runner.sh <LANE> [results-directory]   (default results/<LANE>)
set -u
LANE=${1:?}; ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
OUTDIR=${2:-results/$LANE}; OUT=$ROOT/$OUTDIR/RESULTS.md; LOG=$ROOT/tasks/lanes/lane_runner_$LANE.log
MODEL=${MODEL:-build-lane model}; EFFORT=${EFFORT:-medium}; MAX_ROUNDS=${MAX_ROUNDS:-3}
mkdir -p "$ROOT/$OUTDIR"; exec 9>"$ROOT/$OUTDIR/.lane_runner.lock"; flock -n 9 || { echo "$LANE is already running"; exit 0; }
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 NUMEXPR_NUM_THREADS=2
done_ok(){ [ -s "$OUT" ] && head -3 "$OUT" | grep -q "$LANE"; }
for i in $(seq 1 $MAX_ROUNDS); do
  done_ok && { echo "[$(date +%T)] $LANE done" >> "$LOG"; exit 0; }
  echo "[$(date +%T)] $LANE round $i ($MODEL/$EFFORT)" >> "$LOG"
  timeout 7200 nice -n 10 lane_runner exec -m "$MODEL" -c model_reasoning_effort="$EFFORT" --dangerously-bypass-approvals-and-sandbox --skip-git-repo-check \
    -C "$ROOT" "You have a mandate from the coordinator (Anton) to carry out the task without asking; choose reasonable defaults yourself and report them. Read and carry out tasks/lanes/$LANE.md in full. First check what already exists in $OUTDIR/ and build on it; do not redo completed steps. Full permissions (Anton 24/9): network, cloud (Modal/OVH via tasks/cloud_run.sh), and writing where the task requires it in the workspace (e.g. results/<id>/ packages, tasks/lanes/bt_queue.txt with >>, notes/RESULTS_INDEX.md with >>). Own output in $OUTDIR/, intermediate results on external_media$LANE/. Read-only: ~/projects/bodytwin and other sessions' workspaces. No emails, pushes or publishing; do not touch credentials. Finish with $OUTDIR/RESULTS.md starting with the line '# $LANE'." \
    </dev/null >> "$LOG" 2>&1
done
done_ok && echo "[$(date +%T)] $LANE done" >> "$LOG" || echo "[$(date +%T)] $LANE NOT done after $MAX_ROUNDS" >> "$LOG"
