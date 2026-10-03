#!/bin/bash
# Keeps three build lane lanes (lane_runner exec) going all night. Run by bt-solnight.timer every 5th min and by coordinator cron.
# Each lane = fixed systemd-unit bt-solnight-<LANE>; one round per unit run; new round starts when the previous one ends.
# The castles are in slots.txt (one LANE per row, max 3). The coordinator changes lanes by changing the row.
# Quick error (<180 s) is counted; after 3 in a row, the lane waits for 30 min and writes ALERT.
set -u
ROOT=.
D=$ROOT/tasks/build_night; ST=$D/state; mkdir -p "$ST" "$D/steer"
LOG=$D/ensure.log
MODEL=${MODEL:-lane-model}; EFFORT=${EFFORT:-high}; ROUND_TIMEOUT=${ROUND_TIMEOUT:-14400}
exec 8>"$ST/ensure.lock"; flock -n 8 || exit 0
[ -f "$D/STOP" ] && { echo "[$(date +%F' '%T)] STOP file exists, starting nothing" >> "$LOG"; exit 0; }

n=0
while read -r LANE; do
  LANE=${LANE%%#*}; LANE=$(echo "$LANE" | tr -d '[:space:]'); [ -z "$LANE" ] && continue
  # 2/10 21:25 (anton-5f): roof load out of max_lanes, with 3 as a fallback if the file is missing or is scratched.
  # Anton's directive: send a sol wave. The fallback is deliberately 3 said that a broken or removed
  # file can never stop the driver - it has died silently once before and it cost an entire night.
  MAXL=$(cat "$D/max_lanes" 2>/dev/null); case "$MAXL" in (*[!0-9]*|"") MAXL=3;; esac
  # Count RUNNING units, not read rows: otherwise the ceiling cuts off the bottom slots
  # permanent (2/10: eye and laser were last and never started). Order now only rules
  # who starts first when there is room.
  RUN=$(systemctl --user list-units "bt-solnight-*" --no-legend 2>/dev/null | wc -l)
  case "$RUN" in (*[!0-9]*|"") RUN=0;; esac
  [ "$RUN" -ge "$MAXL" ] && break
  # "A>B" = switch slot from A to B when A's current round ends (A is never interrupted).
  if [[ "$LANE" == *">"* ]]; then
    OLD=${LANE%%>*}; NEW=${LANE##*>}
    systemctl --user is-active --quiet "bt-solnight-$OLD" && continue
    sed -i "s|^$OLD>$NEW\s*$|$NEW|" "$D/slots.txt"; echo "[$(date +%F' '%T)] slotbyte $OLD -> $NEW" >> "$LOG"; LANE=$NEW
  fi
  U=bt-solnight-$LANE
  systemctl --user is-active --quiet "$U" && continue
  { # everything per lane in its own group so one error does not stop the others
    now=$(date +%s)
    back=$(cat "$ST/$LANE.backoff_until" 2>/dev/null || echo 0)
    [ "$now" -lt "$back" ] && continue
    # post the previous round
    if [ -f "$ST/$LANE.started" ]; then
      s=$(cat "$ST/$LANE.started"); dur=$((now - s)); rm -f "$ST/$LANE.started"
      if [ $dur -lt 180 ]; then
        f=$(( $(cat "$ST/$LANE.fastfail" 2>/dev/null || echo 0) + 1 )); echo $f > "$ST/$LANE.fastfail"
        echo "[$(date +%F' '%T)] $LANE snabbfel $f (${dur}s)" >> "$LOG"
        if [ $f -ge 3 ]; then echo $((now + 1800)) > "$ST/$LANE.backoff_until"; echo 0 > "$ST/$LANE.fastfail"
          echo "[$(date +%F' '%T)] ALERT $LANE 3 fast failure, waiting 30 min" >> "$LOG" | tee -a "$D/ALERTS.log" >/dev/null; continue; fi
      else
        # 3/10 17:50: measured 30 of 437 rounds that died on "Selected model is at capacity" and 37 on
        # speed limit, AFTER real work (one of them 143 784 tokens) and without leaving
        # report file. The driver posted them as complete, so about 7 % of the Sol capacity was silently dropped.
        # A round that ends in a provider error and has no report is now counted as aborted and gets
        # rerun, not as completed. All in subshell with || true so the detector can never stop
        # the driver - it has died silently once and it cost a night.
        R_PREV=$(cat "$ST/$LANE.round" 2>/dev/null || echo 0)
        LOGF="$ROOT/results/$LANE/night_rounds/codex_r$R_PREV.log"
        ABORTED=0
        if [ -f "$LOGF" ] && [ ! -f "$ROOT/results/$LANE/night_rounds/r$R_PREV.json" ]; then
          if tail -40 "$LOGF" 2>/dev/null | grep -qiE 'at capacity|rate limit|usage limit'; then ABORTED=1; fi
        fi
        if [ "$ABORTED" = "1" ]; then
          echo 0 > "$ST/$LANE.fastfail"
          echo "[$(date +%F' '%T)] $LANE round $R_PREV ABORTED by provider error after ${dur}s, rerunning" >> "$LOG"
          echo $((R_PREV - 1)) > "$ST/$LANE.round"
        else
          echo 0 > "$ST/$LANE.fastfail"; echo "[$(date +%F' '%T)] $LANE round ended after ${dur}s" >> "$LOG"
        fi
      fi
      # 3/10 11:55 (anton-5f): ONLY remove .go if it's older than the round we're just posting.
      # Otherwise, a review written by the coordinator while the round was ending will be deleted —
      # it happened for DOMAIN_DATA_TO_CELLS at 11:43 and the lane was idle waiting for a new .go.
      if [ -f "$ST/$LANE.go" ] && [ "$ST/$LANE.go" -nt "$ST/$LANE.round" ]; then
        echo "[$(date +%F' '%T)] $LANE keeps fresh .go through the bookkeeping" >> "$LOG"
        echo "$now" > "$ST/$LANE.ended"
      else
        echo "$now" > "$ST/$LANE.ended"; rm -f "$ST/$LANE.go"
      fi
    fi
    # The coordinator reads completed round and writes state/<LANE>.go (and possibly new slot/control). Without a response to 40 min the lane continues by itself.
    if [ -f "$ST/$LANE.ended" ] && [ ! -f "$ST/$LANE.go" ]; then
      e=$(cat "$ST/$LANE.ended"); [ $((now - e)) -lt 2400 ] && continue
      echo "[$(date +%F' '%T)] $LANE no coordinator review in 40 min, continuing on its own" >> "$LOG"
    fi
    rm -f "$ST/$LANE.ended" "$ST/$LANE.go"
    [ -f "$ROOT/tasks/lanes/$LANE.md" ] || { echo "[$(date +%F' '%T)] ALERT saknar brief $LANE" >> "$LOG"; continue; }
    OUT=$ROOT/results/$LANE; mkdir -p "$OUT/night_rounds" "/mnt/games-240/research/bodytwin_solnight/$LANE" 2>/dev/null
    R=$(( $(cat "$ST/$LANE.round" 2>/dev/null || echo 0) + 1 )); echo $R > "$ST/$LANE.round"
    [ -f "$D/steer/$LANE.md" ] || echo "# Styrning $LANE" > "$D/steer/$LANE.md"
    PROMPT="You're a build lane in the night's breakthrough hunt for BodyTwin, lane $LANE, round $R. You have a mandate to perform the assignment without asking; choose reasonable defaults and report them. Read and follow tasks/build_night/COMMON.md i sin helhet (<LANE> = $LANE, <N> = $R), sedan tasks/build_night/steer/$LANE.md, tasks/lanes/$LANE.md and present in results/$LANE/. Continue the same research from the latest checkpoint and NEXT_ROUND.md. Avsluta med RESULTS.md-avsnittet, WORK_STATUS.json, NEXT_ROUND.md and results/$LANE/night_rounds/r$R.json."
    date +%s > "$ST/$LANE.started"
    systemd-run --user --unit="$U" --collect --quiet -p TimeoutStopSec=60 -p RuntimeMaxSec="$ROUND_TIMEOUT" -p MemoryMax="${LANE_MEM:-4G}" --slice=bt-solnight.slice \
      -p WorkingDirectory="$ROOT" \
      -E OMP_NUM_THREADS=2 -E OPENBLAS_NUM_THREADS=2 -E MKL_NUM_THREADS=2 -E NUMEXPR_NUM_THREADS=2 -E HOME=~ -E PATH="$PATH" \
      /bin/bash -c "exec nice -n 10 lane_runner exec -m '$MODEL' -c model_reasoning_effort='$EFFORT' -c tools.web_search=true --dangerously-bypass-approvals-and-sandbox --skip-git-repo-check -C '$ROOT' \"\$0\" </dev/null >> '$OUT/night_rounds/lane_runner_r$R.log' 2>&1" "$PROMPT" \
      && echo "[$(date +%F' '%T)] start $LANE round $R ($MODEL/$EFFORT)" >> "$LOG" \
      || { echo "[$(date +%F' '%T)] ALERT systemd-run failed $LANE" >> "$LOG"; rm -f "$ST/$LANE.started"; }
  } || true
done < "$D/slots.txt"
exit 0
