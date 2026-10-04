#!/usr/bin/env bash
# 1/10 00:20 (anton-5f): the agent's event stream goes to agent.log, where runtime.run_guarded counts finished tools; previously it ended up in LOCAL_WORKER.log and the guard killed all local jobs after 6 min (exit 75).
# 2/10 01:45 (anton-5f, Dental's finding): opencode is a Bun binary that unpacks its native addon to
# $TMPDIR/.<hash>-00000000.so (5.6-13.7 MB) and never cleans up. With ~430 new ones per hour, / grew by 2.4 GB/h
# and 6962 files (39 GB) had accumulated. Each job therefore gets its OWN TMPDIR on games-240, torn down at exit.
W=
P=$1; M=$2; J=$3
CLAIM=$W/results/$J/.local_claim
printf '%s\n' "$$" > "$CLAIM"
TMPBASE=external_mount
mkdir -p "$TMPBASE" 2>/dev/null
JOBTMP=$(mktemp -d "$TMPBASE/$J.XXXXXX" 2>/dev/null) || JOBTMP=""
if [ -n "$JOBTMP" ]; then
  export TMPDIR="$JOBTMP" TMP="$JOBTMP" TEMP="$JOBTMP"
  trap 'rm -f "$CLAIM"; rm -rf "$JOBTMP"' EXIT
else
  # Falls back on the system /tmp if games-240 cannot be written to; cleanup is then handled by the sweep.
  trap 'rm -f "$CLAIM"' EXIT
fi
RP=external_research_path
nice -n 19 ionice -c3 timeout 5400 python3 "$RP" "$P" "$M" --dir "$W/results/$J" --title "$J" "Read BRIEF.md and execute the bounded task. First check what already exists in this directory from an earlier interrupted session and build on it. Write RESULTS.md starting with the line \"$J\" when done." >> "$W/results/$J/agent.log" 2>&1
status=$?
printf '%s\n' "$status" > "$W/results/$J/LOCAL_AGENT_EXIT"
printf '[%s] localfinish %s exit=%s results=%s\n' "$(date +%T)" "$J" "$status" "$([ -s "$W/results/$J/RESULTS.md" ] && echo yes || echo no)" >> "$W/tasks/lanes/bt_queue.log"
exit "$status"
