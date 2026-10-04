#!/bin/bash
# Run lane_runner lanes in priority order, at most PAR simultaneously (cost/load). Log: tasks/lanes/lane_runner_queue.log
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd); PAR=${PAR:-3}; cd "$ROOT"
Q=("BT-IM1" "BT-IM2" "BT-FV2 results/FV2" "BT-FV3" "BT-H3b results/H3b" "BT-C4 results/C4" "BT-X3c results/X3c" "BT-P2b" "BT-V7" "BT-C1b" "BT-R2")
for item in "${Q[@]}"; do
  while [ "$(jobs -rp | wc -l)" -ge "$PAR" ]; do sleep 30; done
  echo "[$(date +%T)] start $item" >> tasks/lanes/codex_queue.log
  tasks/lanes/run_codex.sh $item &
  sleep 10
done
wait; echo "[$(date +%T)] queue finished" >> tasks/lanes/lane_runner_queue.log
