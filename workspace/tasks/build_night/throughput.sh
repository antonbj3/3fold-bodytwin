#!/bin/bash
# Swarm throughput, measured with -mmin and never with -newermt.
# 3/10 18:45 (anton-5f): third time tonight that `find -name RESULTS.md -newermt '-30 minutes'`
# returned 0 while -mmin -30 returned 33 on the same tree. -newermt silently matches nothing here
# instead of erroring, so it reads as "the swarm is idle" when the swarm is at full rate. Twice I
# reported that as a finding. Empty is not absence.
cd "$(dirname "$0")/../.." || exit 1
printf '%-14s %8s %8s\n' window reports per_hour
for m in 15 30 60 120; do
  n=$(find results -name RESULTS.md -mmin -$m 2>/dev/null | wc -l)
  printf '%-14s %8d %8d\n' "${m} min" "$n" "$(( n * 60 / m ))"
done
printf '\nSol rounds in the last 2 h: %s\n' "$(find results -path '*/night_rounds/r*.json' -mmin -120 2>/dev/null | wc -l)"
printf 'queue head: %s\n' "$(head -1 tasks/lanes/bt_queue.txt 2>/dev/null)"
