#!/usr/bin/env bash
# BodyTwin model job queue (after crash 2026-09-24 00:5x). Reading bt_queue.txt: lines "PROFIL MODELL JOBB_ID".
# Starts next job when own BT-opencode-jobb < MAX, MemAvailable > 12 GB and memory PSI full avg10 < 10.
# Skips jobs that already have RESULTS.md or are already running. Runs under setsid so it survives the session.
W=$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd); RP=external_research_path
Q=$W/tasks/lanes/bt_queue.txt; MAX=$(cat $W/tasks/lanes/bt_queue.max 2>/dev/null || echo 10); LOG=$W/tasks/lanes/bt_queue.log
running(){ ps -eo args | grep -E '^[^ ]*opencode run' | grep -oE -- "--title $1( |$)" | head -1; }
nown(){ ps -eo args | grep -E '^[^ ]*opencode run' | grep -cE -- '--title BT-'; }
while true; do
  left=0
  while read -r P M J; do
    [ -z "$J" ] && continue; [[ "$P" == \#* ]] && continue
    [ -s $W/results/$J/RESULTS.md ] && continue
    [ -n "$(running $J)" ] && continue
    left=$((left+1))
    avail=$(awk '/MemAvailable/{print int($2/1048576)}' /proc/meminfo); psi=$(awk -F'avg10=' '/^full/{split($2,a," ");print int(a[1])}' /proc/pressure/memory)
    if [ "$(nown)" -lt "$MAX" ] && [ "$avail" -gt 12 ] && [ "$psi" -lt 10 ]; then
      echo "[$(date +%T)] start $J ($P $M) own=$(nown) avail=${avail}G psi=$psi" >> $LOG
      setsid nohup nice -n 19 ionice -c3 timeout 5400 python3 $RP $P $M --dir $W/results/$J --title $J "Read BRIEF.md and execute the bounded task. First check what already exists in this directory from an earlier interrupted session and build on it. Write RESULTS.md starting with the line \"$J\" when done." > $W/tasks/lanes/profile${P}_$J.log 2>&1 < /dev/null &
      sleep 20
    fi
  done < $Q
  [ $left -eq 0 ] && { echo "[$(date +%T)] queue empty" >> $LOG; sleep 120; }
  sleep 60
done
