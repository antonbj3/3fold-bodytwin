#!/usr/bin/env bash
# Local third execution host; no research.slice. Admission waits count as jobs.
W=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
Q=$W/tasks/lanes/bt_queue.txt; LOG=$W/tasks/lanes/bt_queue.log
exec 9>"$W/tasks/lanes/bt_queue.local-dispatcher.lock"
flock -w 30 9 || exit 1
running(){ ps -eo args | grep -E '^[^ ]*opencode run' | grep -F -- "--title $1 " | head -1; }
nown(){
  # Count each owned worker, plus other local BT runtimes. No waiting-job gap.
  python3 - "$W" <<'PY'
from pathlib import Path
import os,re,subprocess,sys
root=Path(sys.argv[1]);owned=set()
for claim in (root/'results').glob('BT-*/.local_claim'):
 try:os.kill(int(claim.read_text().strip()),0);owned.add(claim.parent.name)
 except (OSError,ValueError):continue
args=subprocess.check_output(['ps','-eo','args'],text=True)
external=set(re.findall(r'^\S*opencode run .*?--title (BT-\S+)',args,re.M))
print(len(owned|external))
PY
}
while true; do
  MAX=$(cat "$W/tasks/lanes/bt_queue.max" 2>/dev/null || echo 12)
  AUTO_LOCAL=0
  if [ "$MAX" -gt 0 ] && [ -f external_research_path ]; then
    MAX=10000
    AUTO_LOCAL=1
  fi
  while read -r P M J; do
    [ -z "$J" ] && continue; [[ "$P" == \#* ]] && continue
    [ "$M" = swarm_worker ] && continue  # preserve paid-model pause
    [ -d "$W/results/$J" ] || continue
    [ -s "$W/results/$J/RESULTS.md" ] && continue
    [ -e "$W/results/$J/.ovh_claim" ] && continue
    [ -n "$(running "$J")" ] && continue
    if [ -e "$W/results/$J/.local_claim" ]; then
      pid=$(cat "$W/results/$J/.local_claim")
      [[ "$pid" =~ ^[0-9]+$ ]] && kill -0 "$pid" 2>/dev/null && continue
      rm -f "$W/results/$J/.local_claim"
    fi
    if [ "$(grep -cF "start $J ($P $M)" "$LOG")" -ge 3 ]; then
      python3 external_research_path "$J" || continue
    fi
    avail=$(awk '/MemAvailable/{print $2}' /proc/meminfo)
    psi_ok=$(awk -F'avg10=' '/^full/{split($2,a," ");print (a[1]<10 ? 1 : 0)}' /proc/pressure/memory)
    own=$(nown)
    [ "$own" -ge "$MAX" ] && break
    [ "$avail" -le 12582912 ] && break
    [ "$psi_ok" != 1 ] && break
    (set -o noclobber; echo "$$" > "$W/results/$J/.local_claim") 2>/dev/null || continue
    if [ -e "$W/results/$J/.ovh_claim" ] || [ -n "$(running "$J")" ]; then
      rm -f "$W/results/$J/.local_claim"; continue
    fi
    if [ "$(grep -cF "start $J ($P $M)" "$LOG")" -ge 3 ]; then
      python3 external_research_path "$J" --claim || { rm -f "$W/results/$J/.local_claim"; continue; }
    fi
    if [ "$AUTO_LOCAL" -eq 1 ]; then
      python3 external_research_path --launch "$P" "$M" "$J" >> "$W/tasks/lanes/profile${P}_$J.log" 2>&1 < /dev/null
      START_STATUS=$?
      if [ "$START_STATUS" -ne 0 ]; then
        rm -f "$W/results/$J/.local_claim"
        echo "[$(date +%T)] local admission deferred $J status=$START_STATUS" >> "$LOG"
        [ "$START_STATUS" -eq 77 ] && continue
        break
      fi
    else
      setsid nohup /bin/bash "$W/tasks/lanes/oom_slots_runtime/local_worker.sh" "$P" "$M" "$J" > "$W/tasks/lanes/profile${P}_$J.log" 2>&1 < /dev/null 9>&- &
    fi
    echo "[$(date +%T)] start $J ($P $M) own=$own avail=${avail}KiB psi<10" >> "$LOG"
    sleep 8
  done < <(python3 "$W/tasks/lanes/seed_unblock_evidence/seed_queue_order.py" "$Q")
  sleep 15
done
