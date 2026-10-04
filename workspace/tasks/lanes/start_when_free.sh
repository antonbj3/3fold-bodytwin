#!/bin/bash
# Starts an extra lane_runner lane when fewer than 3 BodyTwin lane_runner lanes are running.
cd 
while [ "$(ps -eo args | grep -c '^/bin/bash tasks/lanes/run_codex.sh BT-')" -ge 3 ]; do sleep 60; done
echo "[$(date +%T)] start $*" >> tasks/lanes/codex_queue.log
exec tasks/lanes/run_codex.sh "$@"
