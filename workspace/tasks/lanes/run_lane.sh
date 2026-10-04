#!/usr/bin/env bash
# Runs an OpenCode lane noninteractively on Vertex; resumes THE LANE'S OWN session (--session <id>, found via
# title = lane name), never --continue (it takes the latest session globally). Stops when build/<lane>/RESULTS.md
# exists or after MAX_ROUNDS. Can be restarted over an already running lane: then waits for the running process.
set -u
LANE=${1:?lane, t.ex. U14}; MODEL=${2:-opencode-go/reserve_worker-v4.1-flash}
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
LANEFILE=$ROOT/tasks/lanes/$LANE.md; OUT=$ROOT/results/$LANE/RESULTS.md; LOG=$ROOT/tasks/lanes/lane_$LANE.log
MAX_ROUNDS=${MAX_ROUNDS:-12}
export GOOGLE_VERTEX_PROJECT=project-8029a247-2b65-4a21-8f8 GOOGLE_VERTEX_LOCATION=global
export GOOGLE_CLOUD_PROJECT=project-8029a247-2b65-4a21-8f8 GOOGLE_CLOUD_LOCATION=global
export OPENCODE_CONFIG_CONTENT='{"permission":{"external_directory":"allow"}}'
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2
mkdir -p "$ROOT/results/$LANE"; cd "$ROOT"
exec 9>"$ROOT/results/$LANE/.driver.lock"
flock -n 9 || { echo "$LANE: driver already running"; exit 0; }
sid_of(){ python3 -c "import sqlite3,sys; c=sqlite3.connect('file:local_config_path/share/opencode/opencode.db?mode=ro',uri=True); r=c.execute(\"select id from session where title=? order by time_updated desc limit 1\",(sys.argv[1],)).fetchone(); print(r[0] if r else '')" "$LANE" 2>/dev/null; }
busy(){ ps -eo pid=,args= | grep -E -- "^[[:space:]]*[0-9]+[[:space:]]+([^[:space:]]*/)?[o]pencode[[:space:]]+run[[:space:]].*(--title $LANE( |$)|--session ${SID:-NONE}( |$))" >/dev/null; }
SID=$(sid_of)
[ "${NEW_SESSION:-0}" = 1 ] && SID=""
if [ -z "$SID" ]; then
  echo "[$(date +%T)] start $LANE med $MODEL" | tee -a "$LOG"
  timeout 3600 ~/.opencode/bin/opencode run --agent build --model "$MODEL" --title "$LANE" \
    "You are authorized to execute this without approval – ask no questions, do not enter plan mode. Read and execute the lane file $LANEFILE in full. First check what already exists in $ROOT/results/$LANE/ (an earlier session may have built parts) and build on it. Work only under $ROOT/results/$LANE/. Finish by writing $OUT." </dev/null >>"$LOG" 2>&1
  for k in 1 2 3 4 5 6; do SID=$(sid_of); [ -n "$SID" ] && break; sleep 10; done
else
  echo "[$(date +%T)] $LANE: resumes session $SID" | tee -a "$LOG"
fi
[ -z "$SID" ] && { echo "[$(date +%T)] $LANE: ingen session hittad, avbryter" | tee -a "$LOG"; exit 1; }
for i in $(seq 1 $MAX_ROUNDS); do
  while busy; do sleep 30; done
  done_ok(){ [ -s "$OUT" ] && head -5 "$OUT" | grep -q "$LANE"; }
  done_ok && { echo "[$(date +%T)] $LANE complete after $i round(s)" | tee -a "$LOG"; exit 0; }
  echo "[$(date +%T)] $LANE lacks RESULTS.md, sending continue ($i/$MAX_ROUNDS) to $SID" | tee -a "$LOG"
  timeout 3600 ~/.opencode/bin/opencode run --agent build --session "$SID" --model "$MODEL" \
    "APPROVED – RUN NOW. You are authorized by the coordinator (Anton) to execute the entire lane file without asking; ask no questions, choose the most reasonable defaults yourself and report them in RESULTS.md. AGENTS.md's line about the workspace not triggering research jobs applies to automatic graph loading, not this lane. Continue according to the lane file $LANEFILE. First check what already exists in $ROOT/results/$LANE/ and build on it. Finish by writing $OUT." </dev/null >>"$LOG" 2>&1
  sleep 5
done
done_ok && echo "[$(date +%T)] $LANE complete" | tee -a "$LOG" || echo "[$(date +%T)] $LANE NOT complete after $MAX_ROUNDS rounds, see $LOG" | tee -a "$LOG"
