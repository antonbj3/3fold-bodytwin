#!/usr/bin/env bash
# Stop BodyTwin-swarm_worker jobs on account A (default profile) that have not written any log (hung). Only own BT titles.
L=tasks/lanes
for pid in $(ps -eo pid=,args= | awk '/[o]pencode run/ && /swarm_worker/ && /--title BT-/ {print $1}'); do
  env=$(tr '\0' '\n' < /proc/$pid/environ 2>/dev/null | grep -c 'SPACE_SWARM_MULTIKEY_20260923/profiles/')
  [ "$env" -gt 0 ] && continue                                   # B/C profile → do not touch
  t=$(ps -o args= -p $pid | grep -oE -- '--title BT-[A-Za-z0-9-]+' | cut -d' ' -f2)
  log=$L/profileA_$t.log; [ -s "$log" ] && continue                # has written → not hung
  echo "stopping $pid $t (empty log)"; kill $pid
done
