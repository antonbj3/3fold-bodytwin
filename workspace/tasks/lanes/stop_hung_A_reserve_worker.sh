#!/usr/bin/env bash
# Stoppar BodyTwin-reserve_worker-jobb on account A (standardprofil) who did not write any log (Hanging). Bara egna BT-titlar.
L=./tasks/lanes
for pid in $(ps -eo pid=,args= | awk '/[o]pencode run/ && /reserve_worker/ && /--title BT-/ {print $1}'); do
  env=$(tr '\0' '\n' < /proc/$pid/environ 2>/dev/null | grep -c 'SPACE_SWARM_MULTIKEY_20260923/profiles/')
  [ "$env" -gt 0 ] && continue                                   # B/C-profil → do not touch
  t=$(ps -o args= -p $pid | grep -oE -- '--title BT-[A-Za-z0-9-]+' | cut -d' ' -f2)
  log=$L/profileA_$t.log; [ -s "$log" ] && continue                # har skrivit → not hanging
  echo "stoppar $pid $t (tom logg)"; kill $pid
done
