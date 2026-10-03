#!/usr/bin/env bash
# Free disk created by the BodyTwin lane itself (Anton via graph session 11:05: disks full, delete only own).
# 1) shared_data: N40b's out-of-date cloud copies (cloud_base2 = aborted run, replaced by k10/fix; job_k20 = canceled K=20).
# 2) /: inputs/-kopior (datasetkopior) i AVSLUTADE BodyTwin-jobb (RESULTS.md exists) > 100 MB. Original data remain in /mnt/shared_data/datasets.
set -u
rm -rf /mnt/shared_data/bodytwin_work/N40b/cloud_base2 /mnt/shared_data/bodytwin_work/N40b/job_k20
R=./results
for d in $R/BT-*/inputs; do j=$(dirname $d); [ -s $j/RESULTS.md ] || continue
  ps -eo args | grep -qE -- "--title $(basename $j)( |$)" && continue
  s=$(du -sm $d | cut -f1); [ $s -gt 100 ] && { echo "removing $d ($s MB)"; rm -rf $d; echo "inputs removed $(date -Is) after completion to free disk; sources in BRIEF.md and original dataset" > $j/INPUTS_REMOVED.txt; }
done
