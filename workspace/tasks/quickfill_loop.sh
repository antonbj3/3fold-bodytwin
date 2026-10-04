#!/usr/bin/env bash
# Keep at least 100 unfinished queue packets: run quickfill.py (+150) when the queue is short. Stops at 23:40.
W=; cd $W
while [ $(date +%H%M) -lt 2340 ]; do
  left=$(awk '{print $3}' tasks/lanes/bt_queue.txt | while read j; do [ -n "$j" ] && [ ! -s results/$j/RESULTS.md ] && echo $j; done | wc -l)
  if [ "$left" -lt 100 ]; then echo "[$(date +%T)] kvar $left → quickfill" >> tasks/lanes/quickfill.log; python3 tasks/quickfill.py 150 >> tasks/lanes/quickfill.log 2>&1; fi
  sleep 120
done
