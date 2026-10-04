#!/usr/bin/env bash
# Queue driver 2: runs Space The swarm/reserve_worker jobs on OVH (Anton 24/9 ~15:15: "yes", profile keys on OVH).
# Reads the same bt_queue.txt as the local driver. Claims a job with results/<J>/.ovh_claim (the local driver skips it),
# sends the packet to /opt/agents/jobs/<J>, starts its own systemd unit on OVH (1.1G/100 % CPU), and fetches it
# back when AGENT_EXIT exists. Cap: tasks/lanes/ovh_agents/cap (default 30). Stop: systemctl --user stop bt-queue-ovh.
W=$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd); Q=$W/tasks/lanes/bt_queue.txt; D=${BT_CLOUD_DIR:-$W/tasks/lanes/ovh_agents}
LOG=$D/queue_ovh.log; LLOG=$W/tasks/lanes/bt_queue.log; RUN=$D/running; mkdir -p $RUN
SSHO=(-o ConnectTimeout=15 -o ServerAliveInterval=30 -i local_config_path/hunt_20260923 -o IdentitiesOnly=yes -o UserKnownHostsFile=external_research_path -o BatchMode=yes)
H=${BT_CLOUD_HOST:-ubuntu@51.77.110.4}; END=1790948885
localrun(){ ps -eo args | grep -E '^[^ ]*opencode run' | grep -oE -- "--title $1( |$)" | head -1; }
while true; do
  CAP=$(cat $D/cap 2>/dev/null || echo 30)
  # 1) fetch finished jobs: done = AGENT_EXIT exists OR no agent-<J> unit running anymore (timeout/OOM leaves no AGENT_EXIT)
  ACTIVE=$(ssh -n "${SSHO[@]}" $H "systemctl list-units --no-legend --state=running 'agent-*' | awk '{print \$1}'" 2>/dev/null) || ACTIVE="__SSH_FAIL__"
  if [ "$ACTIVE" != "__SSH_FAIL__" ]; then
    for f in $RUN/*; do [ -e "$f" ] || continue; J=$(basename $f)
      echo "$ACTIVE" | grep -q "^agent-$J-" && continue
      if ! python3 external_research_path "$H:/opt/agents/jobs/$J/" "$W/results/$J/" "ssh ${SSHO[*]}" < /dev/null; then
        echo "[$(date +%T)] fetch failed; remote preserved $J" >> "$LOG"
        continue
      fi
      rm -f $f $W/results/$J/.ovh_claim
      echo "[$(date +%T)] fetched $J exit=$(cat $W/results/$J/AGENT_EXIT 2>/dev/null || echo none-timeout) results=$([ -s $W/results/$J/RESULTS.md ] && echo yes || echo no)" >> $LOG
      ssh -n "${SSHO[@]}" $H "rm -rf /opt/agents/jobs/$J"
    done
  else echo "[$(date +%T)] ssh fail at fetch" >> $LOG; fi
  # 2) start new jobs (not after 22:30)
  n=$(ls $RUN | wc -l)
  SLOTS=0
  if [ "$ACTIVE" != "__SSH_FAIL__" ]; then
    TOTAL=$(printf '%s\n' "$ACTIVE" | grep -c '^agent-' || true)
    # Reserve each running agent's full systemd MemoryMax, including agents
    # started by the Field queue in the same slice. MemoryCurrent alone misses
    # agents waiting at provider admission that may later use their full limit.
    SLICE=$(ssh -n "${SSHO[@]}" "$H" "systemctl show research.slice -p MemoryMax -p MemoryCurrent" 2>/dev/null) || SLICE=
    SLICE_MAX=$(sed -n 's/^MemoryMax=//p' <<< "$SLICE")
    SLICE_CURRENT=$(sed -n 's/^MemoryCurrent=//p' <<< "$SLICE")
    MEM=$(ssh -n "${SSHO[@]}" "$H" "awk '/MemAvailable/{print int(\$2/1024)}' /proc/meminfo" 2>/dev/null) || MEM=0
    AGENT_MIB=${BT_AGENT_MEMORY_MIB:-1500}
    HEADROOM_MIB=${BT_HEADROOM_MIB:-4096}
    HOST_CAP=${BT_HOST_CAP:-42}
    if [[ "$SLICE_MAX" =~ ^[0-9]+$ && "$SLICE_CURRENT" =~ ^[0-9]+$ && "$MEM" =~ ^[0-9]+$ && "$AGENT_MIB" =~ ^[0-9]+$ && "$HEADROOM_MIB" =~ ^[0-9]+$ && "$HOST_CAP" =~ ^[0-9]+$ ]] && [ "$AGENT_MIB" -gt 0 ]; then
      AGENT_BYTES=$((AGENT_MIB*1024*1024))
      HEADROOM_BYTES=$((HEADROOM_MIB*1024*1024))
      # The reservation bound prevents idle agents from overcommitting the slice.
      # The current-use bound also covers untracked memory inside research.slice.
      SLOTS=$(( (SLICE_MAX-HEADROOM_BYTES)/AGENT_BYTES-TOTAL ))
      LIVE_SLOTS=$(( (SLICE_MAX-SLICE_CURRENT-HEADROOM_BYTES)/AGENT_BYTES ))
      HOST_SLOTS=$(( (MEM-HEADROOM_MIB)/AGENT_MIB ))
      CAP_SLOTS=$(( HOST_CAP-TOTAL ))
      [ "$SLOTS" -gt "$LIVE_SLOTS" ] && SLOTS=$LIVE_SLOTS
      [ "$SLOTS" -gt "$HOST_SLOTS" ] && SLOTS=$HOST_SLOTS
      [ "$SLOTS" -gt "$CAP_SLOTS" ] && SLOTS=$CAP_SLOTS
      [ "$SLOTS" -lt 0 ] && SLOTS=0
    fi
    DISK=$(ssh -n "${SSHO[@]}" "$H" "python3 -c \"import shutil; print(shutil.disk_usage('/opt/agents').free)\"" 2>/dev/null) || DISK=0
    [[ "$DISK" =~ ^[0-9]+$ ]] || DISK=0
    [ "$DISK" -lt 2147483648 ] && SLOTS=0
  fi
  if [ $(date +%s) -lt $END ] && [ "$SLOTS" -gt 0 ]; then
    while read -r P M J; do
      [ "$n" -ge "$CAP" ] && break
      [ "$SLOTS" -le 0 ] && break
      [ -z "$J" ] && continue; [[ "$P" == \#* ]] && continue
      if [[ "$J" == BT-DW48-* ]]; then
        DENT_N=$(find "$RUN" -maxdepth 1 -name 'BT-DW48-*' -type f | wc -l)
        [ "$DENT_N" -ge 12 ] && continue
      fi
      nds=$(for r in $RUN/*; do [ -e "$r" ] && { [ -s "$r" ] && cat "$r" || echo reserve_worker; }; done | grep -c reserve_worker); [ "$M" = reserve_worker ] && [ "$nds" -ge "$(cat $D/cap_reserve_worker 2>/dev/null || echo 30)" ] && continue   # fasta platser: The swarm tappar aldrig alla platser
      [ -s $W/results/$J/RESULTS.md ] && continue; [ -e $W/results/$J/.ovh_claim ] && continue; [ -d $W/results/$J ] || continue
      [ -n "$(localrun $J)" ] && continue
      [ "$(grep -cF "start $J ($P $M)" $LLOG)" -ge 3 ] && continue
      if [ "$(grep -cF "ovhstart $J " $LOG)" -ge 2 ]; then
        python3 external_research_path "$J" --claim || continue
      fi
      ( set -o noclobber; : > "$W/results/$J/.ovh_claim" ) 2>/dev/null || continue
      [ -n "$(localrun $J)" ] && { rm -f $W/results/$J/.ovh_claim; continue; }
      LAUNCHER=/opt/agents/run_profile_ovh.py
      if [[ "$J" == BT-DW48-* || "$J" == BT-FW48-* ]] && [ -s "$W/results/$J/_run_dental.py" ]; then
        LAUNCHER=/opt/agents/jobs/$J/_run_dental.py
      fi
      if rsync -aL -e "ssh ${SSHO[*]}" --exclude .ovh_claim $W/results/$J/ $H:/opt/agents/jobs/$J/ < /dev/null && \
         ssh -n "${SSHO[@]}" $H "sudo systemd-run --quiet --collect --slice=research.slice --unit=agent-$J-\$(date +%s) --uid=ubuntu -p MemoryMax=${AGENT_MIB}M -p MemoryHigh=1250M -p CPUQuota=100% -p RuntimeMaxSec=5400 --working-directory=/opt/agents/jobs/$J /bin/bash -c 'python3 $LAUNCHER $P $M --dir /opt/agents/jobs/$J --title $J \"Read BRIEF.md and execute the bounded task. First check what already exists in this directory from an earlier interrupted session and build on it. Write RESULTS.md starting with the line \\\"$J\\\" when done.\" > agent.log 2>&1 < /dev/null; echo \$? > AGENT_EXIT'"; then
        echo $M > $RUN/$J; n=$((n+1)); SLOTS=$((SLOTS-1)); echo "[$(date +%T)] ovhstart $J ($P $M) n=$n" >> $LOG
      else rm -f $W/results/$J/.ovh_claim; echo "[$(date +%T)] start failed $J" >> $LOG; fi
      sleep 3
    done < $Q
  fi
  [ $(date +%s) -ge $END ] && [ "$(ls $RUN | wc -l)" -eq 0 ] && { echo "[$(date +%T)] 07:30 passed and empty — exiting" >> $LOG; exit 0; }
  sleep 15
done
