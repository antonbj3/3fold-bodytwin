#!/usr/bin/env bash
# BodyTwin cloud run (OVH 16 vCPU/32 GiB, UpCloud 4 vCPU/8 GiB; OVH lease until 2026-09-24 23:00, poweroff 23:15).
# Usage:
#   tasks/cloud_run.sh submit <JOB_ID> <lokal_katalog> <ovh|upcloud> <cpus> <mem_gb> <timeout_min> -- <kommando ...>
#   tasks/cloud_run.sh status  <JOB_ID> <ovh|upcloud>
#   tasks/cloud_run.sh collect <JOB_ID> <ovh|upcloud> <lokal_utkatalog>
# The directory is sent to /opt/bt/jobs/<JOB_ID>/ (the command runs there). Python: numpy 2.2.6, scipy 1.18.1, h5py 3.16 via
# PYTHONPATH=/opt/bt/vendor (set automatically). No internet on OVH. No GPU. Unique systemd unit bt-<JOB_ID>
# with CPUQuota, MemoryMax and RuntimeMaxSec. The job MUST be complete before 07:30. Write output in the job directory.
# Receipt: tasks/cloud_receipts.jsonl
set -euo pipefail
K=external_research_path
SSHO=(-i "$HOME/.ssh/hunt_20260923" -o UserKnownHostsFile=$K -o StrictHostKeyChecking=yes -o BatchMode=yes)
REC=tasks/cloud_receipts.jsonl
host(){ case "$1" in ovh) echo ubuntu@51.77.110.4;; upcloud) echo root@95.111.211.151;; *) echo "unknown machine $1" >&2; exit 2;; esac; }
sudo_for(){ [ "$1" = ovh ] && echo sudo || echo ""; }
cmd=${1:?submit|status|collect}; shift
case "$cmd" in
submit)
  JOB=$1 DIR=$2 M=$3 CPU=$4 MEM=$5 TMIN=$6; shift 6; [ "$1" = -- ] && shift
  [[ "$JOB" =~ ^[A-Za-z0-9_.-]+$ ]] || { echo "ogiltigt JOB_ID"; exit 2; }
  H=$(host $M); S=$(sudo_for $M)
  now=$(date +%s); end=$(date -d '2026-09-25 07:30' +%s); [ $((now+TMIN*60)) -le $end ] || { echo "timeout passerar 07:30 — korta jobbet eller dela upp"; exit 3; }
  ssh "${SSHO[@]}" $H "test ! -e /opt/bt/jobs/$JOB" || { echo "JOB_ID already exists on $M"; exit 4; }
  if [ "$M" = ovh ]; then used=$(ssh "${SSHO[@]}" $H 't=0; for u in $(systemctl list-units --type=service --state=running --no-legend "bt-*" | awk "{print \$1}"); do q=$(systemctl show $u -p CPUQuotaPerSecUSec --value | tr -dc 0-9); t=$((t+q)); done; echo $t'); [ $(( ${used:-0} + CPU )) -le 12 ] || { echo "OVH quota: BodyTwin already uses ${used} vCPU out of 12 — queue (wait) or use upcloud"; exit 5; }; fi
  rsync -a --exclude __pycache__ -e "ssh ${SSHO[*]}" "$DIR"/ $H:/opt/bt/jobs/$JOB/
  q=$(printf '%q ' "$@")
  ssh "${SSHO[@]}" $H "$S systemd-run --unit=bt-$JOB --uid=\$(id -un) --working-directory=/opt/bt/jobs/$JOB \
     -p CPUQuota=$((CPU*100))% -p MemoryMax=${MEM}G -p RuntimeMaxSec=$((TMIN*60)) \
     --setenv=PYTHONPATH=/opt/bt/vendor --setenv=OMP_NUM_THREADS=$CPU --setenv=OPENBLAS_NUM_THREADS=$CPU \
     /bin/bash -c 'set -o pipefail; ( $q ) > run.log 2>&1; echo \$? > EXIT_CODE'"
  printf '{"t":"%s","job":"%s","machine":"%s","dir":"%s","cpus":%s,"mem_gb":%s,"timeout_min":%s,"cmd":%s}\n' "$(date -Is)" "$JOB" "$M" "$DIR" "$CPU" "$MEM" "$TMIN" "$(python3 -c 'import json,sys;print(json.dumps(sys.argv[1:]))' "$@")" >> $REC
  echo "started bt-$JOB on $M";;
status)
  JOB=$1 M=$2; H=$(host $M)
  ssh "${SSHO[@]}" $H "systemctl show bt-$JOB -p ActiveState -p SubState -p ExecMainStatus 2>/dev/null; cat /opt/bt/jobs/$JOB/EXIT_CODE 2>/dev/null || echo 'EXIT_CODE missing (running or died)'; tail -5 /opt/bt/jobs/$JOB/run.log 2>/dev/null; nproc; uptime";;
collect)
  JOB=$1 M=$2 OUT=$3; H=$(host $M); mkdir -p "$OUT"
  rsync -a -e "ssh ${SSHO[*]}" $H:/opt/bt/jobs/$JOB/ "$OUT"/
  (cd "$OUT" && find . -type f ! -name MANIFEST.sha256 -print0 | sort -z | xargs -0 sha256sum > MANIFEST.sha256)
  echo "collected to $OUT; EXIT_CODE=$(cat "$OUT/EXIT_CODE" 2>/dev/null || echo saknas)";;
*) echo "unknown command"; exit 2;;
esac
