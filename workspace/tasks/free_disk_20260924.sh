#!/usr/bin/env bash
# Free disk created by the BodyTwin lane itself (Anton via graph session 11:05: disks full, delete only own).
# 1) shared_data: N40b's out-of-date cloud copies (cloud_base2 = aborted run, replaced by k10/fix; job_k20 = canceled K=20).
# 2) /: inputs/ copies (dataset copies) in COMPLETED BodyTwin jobs (RESULTS.md exists) > 100 MB. Original data remains in external_mount
set -u
rm -rf external_mount external_mount
R=results
for d in $R/BT-*/inputs; do j=$(dirname $d); [ -s $j/RESULTS.md ] || continue
  ps -eo args | grep -qE -- "--title $(basename $j)( |$)" && continue
  s=$(du -sm $d | cut -f1); [ $s -gt 100 ] && { echo "removing $d ($s MB)"; rm -rf $d; echo "inputs removed $(date -Is) after completion to free disk; sources in BRIEF.md and original dataset" > $j/INPUTS_REMOVED.txt; }
done
