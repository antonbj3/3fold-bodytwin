#!/usr/bin/env bash
# The package is uploaded ONCE to /opt/bt/shared/cxw_bundle (rsync below); each shard copies it locally.
# Run CX-WHATIF's 12 shards on OVH (2 vCPU/6 GB each, at most the BT quota of 8 vCPU), fetch results into the modal_out structure.
W=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
set -a; . "$HOME/.bodytwin/runners.env" 2>/dev/null; set +a; B="$CXW_BUNDLE_DIR"; cd $W
LOG=$W/tasks/lanes/cxw_ovh.log; T=$(date +%H%M)
jobs=(); for tag in identity vsd_z001 vsd_z009; do for s in "CCD -15 0 a" "CCD 1 15 b" "AV -15 5 a" "AV 6 25 b"; do set -- $s; jobs+=("CXW-$tag-$1-$4-$T|$tag $1 $2 $3"); done; done
rsync -a --partial -e "ssh -i $HOME/.ssh/hunt_20260923 -o IdentitiesOnly=yes -o UserKnownHostsFile=$KNOWN_HOSTS -o BatchMode=yes" $B/ rocky@54.37.4.45:/opt/bt/shared/cxw_bundle/ >> $LOG 2>&1 && echo "[$(date +%T)] bundle uploaded" >> $LOG
pending=("${jobs[@]}"); running=()
while [ ${#pending[@]} -gt 0 ] || [ ${#running[@]} -gt 0 ]; do
  if [ ${#pending[@]} -gt 0 ]; then
    j=${pending[0]}; id=${j%%|*}; args=${j#*|}
    mkdir -p /tmp/coordinator-1000/cxw_stub && if tasks/cloud_run.sh submit $id /tmp/coordinator-1000/cxw_stub ovh 2 6 200 -- bash -c "cp -a /opt/bt/shared/cxw_bundle/. . && python3 results/CX-WHATIF/cloud_shard.py $args" >> $LOG 2>&1; then running+=("$id"); pending=("${pending[@]:1}"); echo "[$(date +%T)] submit $id" >> $LOG; fi
  fi
  still=()
  for id in "${running[@]}"; do
    st=$(tasks/cloud_run.sh status $id ovh 2>/dev/null | grep -m1 ActiveState= )
    if [ "$st" = "ActiveState=active" ]; then still+=("$id"); else
      tasks/cloud_run.sh collect $id ovh $B/modal_out/$id >> $LOG 2>&1; echo "[$(date +%T)] collected $id $st" >> $LOG; fi
  done; running=("${still[@]}")
  sleep 60
done
echo "[$(date +%T)] all done" >> $LOG
