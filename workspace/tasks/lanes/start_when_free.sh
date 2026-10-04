#!/bin/bash
# Starts an extra lane_runner lane when fewer than 3 BodyTwin lane_runner lanes are running.
cd "$(dirname "$0")/../.." || exit 1
while [ "$(ps -eo args | grep -c '^/bin/bash tasks/lanes/run_lane_runner.sh BT-')" -ge 3 ]; do sleep 60; done
echo "[$(date +%T)] start $*" >> tasks/lanes/lane_runner_queue.log
exec tasks/lanes/run_lane_runner.sh "$@"
