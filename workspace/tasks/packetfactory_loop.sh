#!/usr/bin/env bash
set -euo pipefail
ROOT=.
cd "$ROOT"
LOG=tasks/lanes/packetfactory.log
QUEUE=tasks/lanes/bt_queue.txt
PENDING=results/CX-PACKETFACTORY/queue_pending.txt
exec 9>results/CX-PACKETFACTORY/.queue.lock

log(){ printf '[%s] %s\n' "$(date '+%F %T %Z')" "$*" >> "$LOG"; }
while true; do
  if [[ $(date +%H%M) -ge 2340 ]]; then log 'stopp 23:40'; exit 0; fi
  unfinished=$(python3 - <<'PY'
from pathlib import Path
root=Path('results'); n=0
for line in Path('tasks/lanes/bt_queue.txt').read_text().splitlines():
    parts=line.split()
    if len(parts)==3 and not parts[0].startswith('#') and not (root/parts[2]/'RESULTS.md').is_file():n+=1
print(n)
PY
)
  if (( unfinished < 80 )); then
    log "queue remaining=$unfinished; building up to 150"
    python3 tasks/packetfactory.py generate --limit 150 >> "$LOG" 2>&1
    python3 - <<'PY'
from pathlib import Path
import json
root=Path('results')
queued={v[2] for line in Path('tasks/lanes/bt_queue.txt').read_text().splitlines() if len(v:=line.split())==3}
packs=[]
for d in root.iterdir():
    if not d.is_dir() or not d.name.startswith(('BT-AG-','BT-R-','BT-P-','BT-X-')):continue
    req=d/'inputs/REQUIREMENT.json'
    if req.is_file() and d.name not in queued and not (d/'RESULTS.md').exists():
        packs.append((d.name,json.loads(req.read_text())['family']))
priority={'AG':0,'R':1,'P':2,'X':3};packs.sort(key=lambda p:(priority[p[1]],p[0]))
lines=[];swarm=0
for ident,family in packs:
    if family=='X' or (family=='AG' and len(lines)%2==0):lines.append(f'C reserve_worker {ident}')
    else:
        lines.append(f'{"ABC"[swarm%3]} swarm {ident}');swarm+=1
Path('results/CX-PACKETFACTORY/queue_pending.txt').write_text('\n'.join(lines)+('\n' if lines else ''))
PY
    flock -x 9
    while read -r profile model ident; do
      [[ -z "${ident:-}" ]] && continue
      if ! grep -Eq "^[ABC] (swarm|reserve_worker) ${ident}$" "$QUEUE" && python3 tasks/packetfactory.py check "$ident" >/dev/null; then
        printf '%s %s %s\n' "$profile" "$model" "$ident" >> "$QUEUE"
        log "queued $ident ($profile $model)"
      fi
    done < "$PENDING"
    flock -u 9
    remaining=$(python3 - <<'PY'
import sys
sys.path.insert(0,'tasks')
import packetfactory as p
print(sum(not (p.RESULTS/s[1]).exists() for s in p.specs()))
PY
)
    if (( remaining == 0 )); then log 'all sources exhausted; stop'; exit 0; fi
  fi
  sleep 180
done
