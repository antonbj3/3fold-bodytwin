#!/usr/bin/env bash
# Queue driver 2: runs Space The swarm/swarm_worker jobs on OVH (Anton 24/9 ~15:15: "yes", profile keys on OVH).
# Reads the same bt_queue.txt as the local driver. Claims a job with results/<J>/.ovh_claim (the local driver skips it),
# sends the packet to /opt/agents/jobs/<J>, starts its own systemd unit on OVH (1500M/100 % CPU), and fetches it
# back when AGENT_EXIT exists. Cap: tasks/lanes/ovh_agents/cap (default 30). Stop: systemctl --user stop bt-queue-ovh.
W=$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd); Q=$W/tasks/lanes/bt_queue.txt; D=${BT_CLOUD_DIR:-$W/tasks/lanes/ovh_agents}
# Runtime values outside the tree: keys, dashboard flag and helper scripts.
set -a; . "$HOME/.bodytwin/runners.env" 2>/dev/null || { echo "missing ~/.bodytwin/runners.env"; exit 1; }; set +a
LOG=$D/queue_ovh.log; LLOG=$W/tasks/lanes/bt_queue.log; RUN=$D/running; mkdir -p $RUN
SSHO=(-o ConnectTimeout=15 -o ServerAliveInterval=30 -i "$SSH_KEY" -o IdentitiesOnly=yes -o UserKnownHostsFile="$KNOWN_HOSTS" -o BatchMode=yes)
H=${BT_CLOUD_HOST:-ubuntu@51.77.110.4}; END=${BT_CLOUD_END:-1791042504}
# Shared with Field; the host lock makes the final check + unit start atomic.
# Keep this default equal to the unit's Environment=BT_AGENT_MEMORY_MIB, or the script
# default is a silent no-op. Measured 2026-10-03: our agents peak at 1282 MiB against this
# 1500 booking, so there is only 15 % slack here -- the real waste was elsewhere.
AGENT_MIB=${BT_AGENT_MEMORY_MIB:-1100}  # 3/10 22:00: measured 512-972 MiB actual across 19 live
# agents on the host (mean ~680), against a 1500 MiB booking. Booking above 1100 leaves the host
# with 19.5 GB free while admitting fewer jobs than it can run. The mem-kill requeue path already
# handles a job that exceeds it, so the downside of booking tighter is a retry and not a loss.
HEADROOM_MIB=${BT_HEADROOM_MIB:-4096}
SLICE_HEADROOM_MIB=${BT_SLICE_HEADROOM_MIB:-2048}
HOST_CAP=${BT_HOST_CAP:-17}
GUARD="sudo -u ubuntu python3 /opt/agents/launch_reserved.py --agent-mib $AGENT_MIB --slice-headroom-mib $SLICE_HEADROOM_MIB --host-headroom-mib $HEADROOM_MIB --host-cap $HOST_CAP"
exec 9>"$D/dispatcher.lock"
flock -n 9 || exit 0
localrun(){ ps -eo args | grep -E '^[^ ]*opencode run' | grep -oE -- "--title $1( |$)" | head -1; }
while true; do
  CAP=$(cat $D/cap 2>/dev/null || echo 30)
  # Shared cloud guard selects a measured hard budget and authorises every
  # start. A zero cap still means explicit pause; positive caps need no tuning.
  if [ -f "$AUTOMATIC_FLAG" ] && [ "$CAP" -gt 0 ]; then
    CAP=10000
  fi
  # 1) fetch finished jobs: done = AGENT_EXIT exists OR no agent-<J> unit running anymore (timeout/OOM leaves no AGENT_EXIT)
  ACTIVE=$(ssh -n "${SSHO[@]}" $H "systemctl list-units --no-legend --state=active,activating,deactivating 'agent-*' | awk '{print \$1}'" 2>/dev/null) || ACTIVE="__SSH_FAIL__"
  if [ "$ACTIVE" != "__SSH_FAIL__" ]; then
    for f in $RUN/*; do [ -e "$f" ] || continue; J=$(basename $f)
      echo "$ACTIVE" | grep -q "^agent-$J-" && continue
      if ! python3 "$COLLECT_VERIFIED" "$H:/opt/agents/jobs/$J/" "$W/results/$J/" "ssh ${SSHO[*]}" < /dev/null; then
        echo "[$(date +%T)] fetch failed; remote preserved $J" >> "$LOG"
        continue
      fi
      JMODEL=$([ -s "$f" ] && cat "$f" || echo swarm)
      rm -f $f $W/results/$J/.ovh_claim
      echo "[$(date +%T)] fetched $J exit=$(cat $W/results/$J/AGENT_EXIT 2>/dev/null || echo none-timeout) results=$([ -s $W/results/$J/RESULTS.md ] && echo yes || echo no)" >> $LOG
      EX0=$(cat $W/results/$J/AGENT_EXIT 2>/dev/null || echo "")
      # exit 75 is the provider rate limit. Record a per-model skip so the rotation stops asking for
      # a model that is not answering: measured 3/10 21:56, five consecutive BT-CONN starts exited 75
      # and every one ran on ling or swarm_worker while swarm answered 24 of 38, and throughput fell from
      # 104 to 48 reports an hour. A third model raises the ceiling only while it answers.
      if [ "$EX0" = "75" ]; then
        MOD=$(python3 -c "import json;print(json.load(open('$W/results/$J/MODEL_ROUTE.json')).get('model',''))" 2>/dev/null)
        case "$MOD" in
          *ling*)    echo $(( $(date +%s) + 1200 )) > "$D/.skip_ling"
                     echo "[$(date +%T)] rate-limit on ling, skipping it for 20 min" >> $LOG;;
          *swarm_worker*) echo $(( $(date +%s) + 1200 )) > "$D/.skip_swarm_worker"
                     echo "[$(date +%T)] rate-limit on swarm_worker, skipping it for 20 min" >> $LOG;;
        esac
      fi
      # exit 137 is SIGKILL, which for a unit with MemoryMax means the memory limit. Re-queue once.
      EX=$(cat $W/results/$J/AGENT_EXIT 2>/dev/null || echo "")
      if [ "$EX" = "137" ] && [ ! -s "$W/results/$J/RESULTS.md" ] && [ ! -f "$W/results/$J/.mem_retry_done" ]; then
        PREV=$(cat "$W/results/$J/.mem_retry" 2>/dev/null); case "$PREV" in (*[!0-9]*|"") PREV="$AGENT_MIB";; esac
        NEXT=$((PREV * 2)); [ "$NEXT" -gt "${RETRY_CAP_MIB:-8192}" ] && NEXT="${RETRY_CAP_MIB:-8192}"
        echo "$NEXT" > "$W/results/$J/.mem_retry"
        : > "$W/results/$J/.mem_retry_done"
        rm -f "$W/results/$J/AGENT_EXIT"
        echo "A $JMODEL $J" >> "$W/tasks/lanes/bt_queue.txt"
        echo "[$(date +%T)] mem-kill $J: requeued with ${NEXT}MiB" >> $LOG
      fi
      ssh -n "${SSHO[@]}" $H "sudo -u ubuntu rm -rf /opt/agents/jobs/$J"
    done
  else echo "[$(date +%T)] ssh fail at fetch" >> $LOG; fi
  # 2) start new jobs (not after 22:30)
  n=$(ls $RUN | wc -l)
  SLOTS=0
  if [ "$ACTIVE" != "__SSH_FAIL__" ]; then
    # Snapshot uses actual per-unit MemoryMax, plus non-agent slice memory.
    SLOTS=$(ssh -n "${SSHO[@]}" "$H" "$GUARD" 2>/dev/null) || SLOTS=0
    [[ "$SLOTS" =~ ^[0-9]+$ ]] || SLOTS=0
    # Measured 2026-10-03: the guard charges reclaimable page cache to the slice, so it refused all
    # starts while the host had 24 GiB available, 70 % idle CPU and agents using 0.6 GiB each of a
    # 1 GiB booking. research.slice held 8621 MiB of file cache against 6062 MiB of real agent memory.
    # Reclaiming the cache took slots 1 -> 6 without touching a single running agent. This is not a
    # policy change and it raises capacity for every session on the host: it only stops cache from
    # being counted as occupied. Rate-limited to once per RECLAIM_MIN_S.
    if [ "$SLOTS" -eq 0 ]; then
      NOW=$(date +%s); LAST=$(cat "$D/.last_reclaim" 2>/dev/null || echo 0)
      case "$LAST" in (*[!0-9]*|"") LAST=0;; esac
      if [ $((NOW-LAST)) -ge "${RECLAIM_MIN_S:-300}" ]; then
        echo "$NOW" > "$D/.last_reclaim"
        FILE_MIB=$(ssh -n "${SSHO[@]}" "$H" "awk '/^file /{printf \"%d\", \$2/1048576}' /sys/fs/cgroup/research.slice/memory.stat" 2>/dev/null)
        case "$FILE_MIB" in (*[!0-9]*|"") FILE_MIB=0;; esac
        if [ "$FILE_MIB" -gt 2048 ]; then
          ssh -n "${SSHO[@]}" "$H" "echo $((FILE_MIB-1024))M | sudo tee /sys/fs/cgroup/research.slice/memory.reclaim >/dev/null" 2>/dev/null
          echo "[$(date +%T)] reclaimed ${FILE_MIB}MiB-1024 page cache in research.slice" >> $LOG
          SLOTS=$(ssh -n "${SSHO[@]}" "$H" "$GUARD" 2>/dev/null) || SLOTS=0
          [[ "$SLOTS" =~ ^[0-9]+$ ]] || SLOTS=0
          echo "[$(date +%T)] slots after reclaim: $SLOTS" >> $LOG
        fi
      fi
    fi
  fi
  if [ $(date +%s) -lt $END ] && [ "$SLOTS" -gt 0 ]; then
    while read -r P M J; do
      [ "$n" -ge "$CAP" ] && break
      [ "$SLOTS" -le 0 ] && break
      [ -z "$J" ] && continue; [[ "$P" == \#* ]] && continue
      [ "$M" = swarm_worker ] && continue  # paid model remains paused
      # Rate limits are per model, so a third free model raises the ceiling without more memory:
      # measured 2026-10-03, 12.5 % of fetches returned the provider rate-limit exit. ling-3.1-flash
      # was smoke-tested through the real launcher the same day and answered at zero cost. The host
      # router still owns caps and backoff, so a request it cannot honour falls back to swarm by
      # itself. Every third free start asks for ling.
      if [ "$M" = swarm ]; then
        # The counter must persist across dispatcher cycles: the start loop runs in a subshell and a
        # cycle usually starts only one or two jobs, so an in-memory counter never reached three.
        LING_N=$(cat "$D/.ling_n" 2>/dev/null || echo 0)
        case "$LING_N" in (*[!0-9]*|"") LING_N=0;; esac
        LING_N=$((LING_N + 1)); echo "$LING_N" > "$D/.ling_n"
        # 3/10 21:56: ling and swarm_worker are BOTH rate-limited right now while swarm answers -- the
        # five most recent BT-CONN starts all exited 75 and all of them ran on ling or swarm_worker, and
        # throughput fell from 104 to 48 reports an hour. A third model raises the ceiling only while
        # it answers; when it does not, it burns a slot per start. So a model that returns the
        # rate-limit exit is skipped for 20 minutes, written as an epoch in .skip_<model>.
        LING_SKIP=$(cat "$D/.skip_ling" 2>/dev/null || echo 0)
        if [ $((LING_N % 3)) -eq 0 ] && [ "$(date +%s)" -ge "$LING_SKIP" ]; then M=ling; fi
      fi
      if [[ "$J" == BT-DW48-* ]] && [ ! -f "$AUTOMATIC_FLAG" ]; then
        DENT_N=$(find "$RUN" -maxdepth 1 -name 'BT-DW48-*' -type f | wc -l)
        [ "$DENT_N" -ge 12 ] && continue
      fi
      nds=$(for r in $RUN/*; do [ -e "$r" ] && { [ -s "$r" ] && cat "$r" || echo swarm_worker; }; done | grep -c swarm_worker); [ "$M" = swarm_worker ] && [ "$nds" -ge "$(cat $D/cap_swarm_worker 2>/dev/null || echo 30)" ] && continue   # fixed slots: The swarm never loses all slots
      [ -s $W/results/$J/RESULTS.md ] && continue; [ -e $W/results/$J/.ovh_claim ] && continue; [ -e $W/results/$J/.local_claim ] && continue; [ -d $W/results/$J ] || continue
      [ -n "$(localrun $J)" ] && continue
      [ "$(grep -cF "start $J ($P $M)" $LLOG)" -ge 3 ] && continue
      RATE_RETRY=0
      if [ "$(grep -cF "ovhstart $J " $LOG)" -ge 2 ]; then
        python3 "$DEFERRED_GATE" "$J" || continue
        RATE_RETRY=1
      fi
      ( set -o noclobber; : > "$W/results/$J/.ovh_claim" ) 2>/dev/null || continue
      [ -e "$W/results/$J/.local_claim" ] || [ -n "$(localrun $J)" ] && { rm -f $W/results/$J/.ovh_claim; continue; }
      LAUNCHER=/opt/agents/run_profile_ovh.py
      if [[ "$J" == BT-DW48-* || "$J" == BT-FW48-* ]] && [ -s "$W/results/$J/_run_dental.py" ]; then
        LAUNCHER=/opt/agents/jobs/$J/_run_dental.py
      fi
      if rsync -aL --rsync-path="sudo -u ubuntu rsync" -e "ssh ${SSHO[*]}" --exclude .ovh_claim $W/results/$J/ $H:/opt/agents/jobs/$J/ < /dev/null && \
      # A job killed for memory is retried once with a doubled booking. Measured 2026-10-03: zero
      # exit=137 at the current booking, and dental's peak tables put the loss at about 2.5 % of one
      # class if the booking drops to p95 -- so this path is what makes lowering the quantile safe.
      # Bounded on purpose: one retry, never above RETRY_CAP_MIB, and the file is removed on success
      # so a job cannot escalate forever.
      # Per-prefix booking. The dental session measured its own 200 most recent successful runs
      # from /opt/agents/resource_history: median 840 MiB, p95 1319, p99 2004, max 5906, with 24
      # of 200 above 1024. Our own 19 live agents measured 512-972 with a mean of 680. One number
      # cannot serve both profiles: at 1100 roughly a tenth of dental rows die once on memory and
      # each costs a rerun of a 30-60 minute job, while booking 1500 for ours wastes a third of
      # the host. Their measurement governs their rows and ours governs ours.
      case "$J" in
        BT-DW48-*|BT-DENT-*) JOB_MIB="${BT_DENTAL_MEMORY_MIB:-1500}";;
        *)                   JOB_MIB="$AGENT_MIB";;
      esac
      if [ -f "$W/results/$J/.mem_retry" ]; then
        R=$(cat "$W/results/$J/.mem_retry" 2>/dev/null)
        case "$R" in (*[!0-9]*|"") R=0;; esac
        [ "$R" -gt "$JOB_MIB" ] && JOB_MIB="$R"
      fi
         ssh -n "${SSHO[@]}" $H "${GUARD/sudo -u ubuntu/sudo} -- systemd-run --quiet --collect --slice=research.slice --unit=agent-$J-\$(date +%s) --uid=ubuntu -p MemoryMax=${JOB_MIB}M -p MemoryHigh=$((JOB_MIB * 5 / 6))M -p CPUQuota=100% -p RuntimeMaxSec=5400 --working-directory=/opt/agents/jobs/$J /bin/bash -c 'python3 $LAUNCHER $P $M --dir /opt/agents/jobs/$J --title $J \"Read BRIEF.md and execute the bounded task. First check what already exists in this directory from an earlier interrupted session and build on it. Write RESULTS.md starting with the line \\\"$J\\\" when done.\" > agent.log 2>&1 < /dev/null; echo \$? > AGENT_EXIT'"; then
        echo $M > $RUN/$J; n=$((n+1)); SLOTS=$((SLOTS-1)); echo "[$(date +%T)] ovhstart $J ($P $M) n=$n" >> $LOG
        # The shared host guard can reject a stale capacity snapshot. Charge a
        # delayed retry only after a real unit start, never for failed admission.
        if [ "$RATE_RETRY" -eq 1 ]; then
          python3 "$DEFERRED_GATE" "$J" --claim || true
        fi
      else
        rm -f "$W/results/$J/.ovh_claim"
        echo "[$(date +%T)] start failed $J" >> "$LOG"
        SLOTS=0
        break  # refresh shared capacity and fetch completed jobs on the next pass
      fi
      sleep 3
    # 3/10: dental received 6,8 % of 441 starts over 9,17 h while the queue had thousands of dental rows.
    # The cause is the SHARED ordering function, which I do not touch. The filter below sits after
    # it and guarantees a live dental row every third position. It is fail-open: on error
    # the stream is printed unchanged, so the queue cannot be starved by the filter.
    done < <(python3 "$W/tasks/lanes/seed_unblock_evidence/seed_queue_order.py" "$Q" \
             | python3 "$D/priority_interleave.py")
  fi
  [ $(date +%s) -ge $END ] && [ "$(ls $RUN | wc -l)" -eq 0 ] && { echo "[$(date +%T)] 07:30 passed and empty — exiting" >> $LOG; exit 0; }
  sleep 15
done
