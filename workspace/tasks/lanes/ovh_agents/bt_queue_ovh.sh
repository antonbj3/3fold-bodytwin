#!/usr/bin/env bash
# Queue driver 2: runs Space The swarm/reserve_worker jobs on OVH (Anton 24/9 ~15:15: "yes", profile keys on OVH).
# Reads the same bt_queue.txt as the local driver. Claims a job with results/<J>/.ovh_claim (the local driver skips it),
# sends the packet to /opt/agents/jobs/<J>, starts its own systemd unit on OVH (1500M/100 % CPU), and fetches it
# back when AGENT_EXIT exists. Cap: tasks/lanes/ovh_agents/cap (default 30). Stop: systemctl --user stop bt-queue-ovh.
W=.; Q=$W/tasks/lanes/bt_queue.txt; D=${BT_CLOUD_DIR:-$W/tasks/lanes/ovh_agents}
LOG=$D/queue_ovh.log; LLOG=$W/tasks/lanes/bt_queue.log; RUN=$D/running; mkdir -p $RUN
SSHO=(-o ConnectTimeout=15 -o ServerAliveInterval=30 -i ~/.ssh/hunt_20260923 -o IdentitiesOnly=yes -o UserKnownHostsFile=~/research/sol6_recovery_20260923/CLOUD_HUNT_20260923/infra/known_hosts -o BatchMode=yes)
H=${BT_CLOUD_HOST:-ubuntu@51.77.110.4}; END=${BT_CLOUD_END:-1791042504}
# Shared with Field; the host lock makes the final check + unit start atomic.
AGENT_MIB=${BT_AGENT_MEMORY_MIB:-1500}
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
  if [ -f ~/research/AGENT_DASHBOARD_20260930/AUTOMATIC_ENABLED.json ] && [ "$CAP" -gt 0 ]; then
    CAP=10000
  fi
  # 1) fetch finished jobs: done = AGENT_EXIT exists OR no agent-<J> unit running anymore (timeout/OOM leaves no AGENT_EXIT)
  ACTIVE=$(ssh -n "${SSHO[@]}" $H "systemctl list-units --no-legend --state=active,activating,deactivating 'agent-*' | awk '{print \$1}'" 2>/dev/null) || ACTIVE="__SSH_FAIL__"
  if [ "$ACTIVE" != "__SSH_FAIL__" ]; then
    for f in $RUN/*; do [ -e "$f" ] || continue; J=$(basename $f)
      echo "$ACTIVE" | grep -q "^agent-$J-" && continue
      if ! python3 ~/research/FREE_AUTONOMY_20260926/collect_verified.py "$H:/opt/agents/jobs/$J/" "$W/results/$J/" "ssh ${SSHO[*]}" < /dev/null; then
        echo "[$(date +%T)] fetch failed; remote preserved $J" >> "$LOG"
        continue
      fi
      rm -f $f $W/results/$J/.ovh_claim
      echo "[$(date +%T)] fetched $J exit=$(cat $W/results/$J/AGENT_EXIT 2>/dev/null || echo none-timeout) results=$([ -s $W/results/$J/RESULTS.md ] && echo yes || echo no)" >> $LOG
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
  fi
  if [ $(date +%s) -lt $END ] && [ "$SLOTS" -gt 0 ]; then
    while read -r P M J; do
      [ "$n" -ge "$CAP" ] && break
      [ "$SLOTS" -le 0 ] && break
      [ -z "$J" ] && continue; [[ "$P" == \#* ]] && continue
      [ "$M" = reserve_worker ] && continue  # paid model remains paused
      if [[ "$J" == BT-DW48-* ]] && [ ! -f ~/research/AGENT_DASHBOARD_20260930/AUTOMATIC_ENABLED.json ]; then
        DENT_N=$(find "$RUN" -maxdepth 1 -name 'BT-DW48-*' -type f | wc -l)
        [ "$DENT_N" -ge 12 ] && continue
      fi
      nds=$(for r in $RUN/*; do [ -e "$r" ] && { [ -s "$r" ] && cat "$r" || echo reserve_worker; }; done | grep -c reserve_worker); [ "$M" = reserve_worker ] && [ "$nds" -ge "$(cat $D/cap_reserve_worker 2>/dev/null || echo 30)" ] && continue   # fasta platser: The swarm tappar aldrig alla platser
      [ -s $W/results/$J/RESULTS.md ] && continue; [ -e $W/results/$J/.ovh_claim ] && continue; [ -e $W/results/$J/.local_claim ] && continue; [ -d $W/results/$J ] || continue
      [ -n "$(localrun $J)" ] && continue
      [ "$(grep -cF "start $J ($P $M)" $LLOG)" -ge 3 ] && continue
      RATE_RETRY=0
      if [ "$(grep -cF "ovhstart $J " $LOG)" -ge 2 ]; then
        python3 ~/research/FREE_AUTONOMY_20260926/SOL_SYNTHESIS_20260927/deferred_job_gate.py "$J" || continue
        RATE_RETRY=1
      fi
      ( set -o noclobber; : > "$W/results/$J/.ovh_claim" ) 2>/dev/null || continue
      [ -e "$W/results/$J/.local_claim" ] || [ -n "$(localrun $J)" ] && { rm -f $W/results/$J/.ovh_claim; continue; }
      LAUNCHER=/opt/agents/run_profile_ovh.py
      if [[ "$J" == BT-DW48-* || "$J" == BT-FW48-* ]] && [ -s "$W/results/$J/_run_dental.py" ]; then
        LAUNCHER=/opt/agents/jobs/$J/_run_dental.py
      fi
      if rsync -aL --rsync-path="sudo -u ubuntu rsync" -e "ssh ${SSHO[*]}" --exclude .ovh_claim $W/results/$J/ $H:/opt/agents/jobs/$J/ < /dev/null && \
         ssh -n "${SSHO[@]}" $H "${GUARD/sudo -u ubuntu/sudo} -- systemd-run --quiet --collect --slice=research.slice --unit=agent-$J-\$(date +%s) --uid=ubuntu -p MemoryMax=${AGENT_MIB}M -p MemoryHigh=1250M -p CPUQuota=100% -p RuntimeMaxSec=5400 --working-directory=/opt/agents/jobs/$J /bin/bash -c 'python3 $LAUNCHER $P $M --dir /opt/agents/jobs/$J --title $J \"Read BRIEF.md and execute the bounded task. First check what already exists in this directory from an earlier interrupted session and build on it. Write RESULTS.md starting with the line \\\"$J\\\" when done.\" > agent.log 2>&1 < /dev/null; echo \$? > AGENT_EXIT'"; then
        echo $M > $RUN/$J; n=$((n+1)); SLOTS=$((SLOTS-1)); echo "[$(date +%T)] ovhstart $J ($P $M) n=$n" >> $LOG
        # The shared host guard can reject a stale capacity snapshot. Charge a
        # delayed retry only after a real unit start, never for failed admission.
        if [ "$RATE_RETRY" -eq 1 ]; then
          python3 ~/research/FREE_AUTONOMY_20260926/SOL_SYNTHESIS_20260927/deferred_job_gate.py "$J" --claim || true
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
             | python3 "$D/dental_interleave.py")
  fi
  [ $(date +%s) -ge $END ] && [ "$(ls $RUN | wc -l)" -eq 0 ] && { echo "[$(date +%T)] 07:30 passed and empty — exiting" >> $LOG; exit 0; }
  sleep 15
done
