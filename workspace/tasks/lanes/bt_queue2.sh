#!/usr/bin/env bash
# BodyTwin model job queue (after crash 2026-09-24 00:5x). Reading bt_queue.txt: lines "PROFIL MODELL JOBB_ID".
# Starts next job when own BT-opencode-jobb < MAX, MemAvailable > 12 GB and memory PSI full avg10 < 10.
# Skips jobs that already have RESULTS.md or are already running. Runs under setsid so it survives the session.
W=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd); RP=external_research_path
Q=$W/tasks/lanes/bt_queue.txt; LOG=$W/tasks/lanes/bt_queue.log
running(){ ps -eo args | grep -E '^[^ ]*opencode run' | grep -oE -- "--title $1( |$)" | head -1; }
nown(){ ps -eo args | grep -E '^[^ ]*opencode run' | grep -cE -- '--title BT-'; }
while true; do [ $(date +%s) -ge 1790603285 ] && exit 0; MAX=$(cat $W/tasks/lanes/bt_queue.max 2>/dev/null || echo 10)
  left=0
  while read -r P M J; do
    [ -z "$J" ] && continue; [[ "$P" == \#* ]] && continue
    [ -s $W/results/$J/RESULTS.md ] && continue
    [ -n "$(running $J)" ] && continue
    [ -e $W/results/$J/.ovh_claim ] && continue
    [ -n "$(find $W/results/$J -type l -print -quit 2>/dev/null)" ] && { echo "[$(date +%T)] SYMLINK-GUARD skips $J (the package contains symlinks)" >> $LOG; continue; }
    [ "$M" = "$SWARM_TOKEN" ] && [ "$(ps -eo args | grep -E '^[^ ]*opencode run' | grep -- '--title BT-' | grep -c "$SWARM_TOKEN")" -ge "$(cat $W/tasks/lanes/$SWARM_CAP_FILE 2>/dev/null || echo 6)" ] && continue   # fixed slots for The swarm   # taken by the OVH driver (ovh_agents/bt_queue_ovh.sh)
    [ "$(grep -cF "start $J ($P $M)" $LOG)" -ge 3 ] && continue   # at most 3 starts per job and account/model (quota exhausted → hangs silently)
    left=$((left+1))
    avail=$(awk '/MemAvailable/{print int($2/1048576)}' /proc/meminfo); psi=$(awk -F'avg10=' '/^full/{split($2,a," ");print int(a[1])}' /proc/pressure/memory)
    disk_kb=$(df -Pk / | awk 'NR==2 {print $4}')
    if [ "$(nown)" -lt "$MAX" ] && [ "$avail" -gt 14 ] && [ "$psi" -lt 5 ] && [ "${disk_kb:-0}" -ge 5242880 ]; then
      echo "[$(date +%T)] start $J ($P $M) own=$(nown) avail=${avail}G psi=$psi" >> $LOG
      systemd-run --user --quiet --collect --unit=bt-job-$J-$(date +%s) --slice=jobs-bodytwin.slice -p MemoryMax=1100M -p MemoryHigh=900M -p StandardOutput=append:$W/tasks/lanes/profile${P}_$J.log -p StandardError=append:$W/tasks/lanes/profile${P}_$J.log -p CPUQuota=100% -p RuntimeMaxSec=5400 --setenv=OMP_NUM_THREADS=1 --setenv=OPENBLAS_NUM_THREADS=1 --setenv=MKL_NUM_THREADS=1 nice -n 19 ionice -c3 python3 $RP $P $M --dir $W/results/$J --title $J "Read BRIEF.md and execute the bounded task. First check what already exists in this directory from an earlier interrupted session and build on it. Write RESULTS.md starting with the line \"$J\" when done." < /dev/null
      sleep 8
    fi
  done < $Q
  [ $left -eq 0 ] && { echo "[$(date +%T)] queue empty" >> $LOG; sleep 30; }
  sleep 20
done
