#!/bin/bash
# One self-contained cycle: check the net, make briefs from whatever is still undetermined, queue them
# at the front. Runs unattended so it does not depend on me being awake or asked.
#
# Order matters and is the whole point. The staleness check runs FIRST because the net handed out a
# withdrawn number for hours on 2026-10-03 and every brief generated from that node carried it. A brief
# built from an unchecked edge spends free capacity on a figure we have already retracted.
set -u
cd "$(dirname "$0")/../.." || exit 1
LOG=tasks/build_night/brief_cycle.log
ts() { date '+%F %T'; }

{
  echo "[$(ts)] cycle start"

  # 1. Re-read every edge's evidence. Never let this abort the cycle: a driver that dies on a helper
  #    killed a whole night on 2026-09-26.
  if python3 tasks/assembly/net_staleness.py > /tmp/.staleness.$$ 2>&1; then
    head -1 /tmp/.staleness.$$ | sed "s/^/[$(ts)] /"
    if grep -qE '\[STALE|\[MOVED' /tmp/.staleness.$$; then
      echo "[$(ts)] STALE eller MOVED kanter — briefgenerering HOPPAS OVER THIS BIKE"
      grep -E '\[STALE|\[MOVED' /tmp/.staleness.$$ | sed "s/^/[$(ts)]   /"
      grep -E '\[STALE|\[MOVED' /tmp/.staleness.$$ >> tasks/build_night/ALERTS.log
      rm -f /tmp/.staleness.$$
      echo "[$(ts)] cycle end (blockerad)"
      exit 0
    fi
  else
    echo "[$(ts)] staleness check failed, continuing without it anyway"
    tail -3 /tmp/.staleness.$$ | sed "s/^/[$(ts)]   /"
  fi
  rm -f /tmp/.staleness.$$

  # 2. Briefs for edges that are still undetermined. Idempotent: a job with RESULTS.md is skipped, so
  #    this only ever fires on an edge that is new, or whose answer has not come back yet.
  python3 tasks/assembly/queue_creative_from_net.py 2>&1 | sed "s/^/[$(ts)] /"

  # 2b. Creative briefs from every lane's own stated obstacle, across the whole project and not just
  #     the domains I happen to be working in. A lane that gates FAIL writes down exactly what it could
  #     not get, in its own words with its own number beside it, and nothing was consuming that. 55 of
  #     the 63 creative jobs queued on 3/10 came from here, which is the only source that covers every
  #     domain automatically rather than following my attention.
  python3 tasks/assembly/queue_from_anomalies.py 2>&1 | head -2 | sed "s/^/[$(ts)] /"

  # 2c. Connection briefs are the only generator that does not run dry: one per PAIR of cells
  #     that share at least two dimensions, 257 such pairs against 120 queued, so it refills
  #     itself as the queue drains. The others are one-per-thing and are already exhausted.
  python3 tasks/assembly/queue_connection_briefs.py 2>&1 | head -1 | sed "s/^/[$(ts)] /"

  # 2d. Sync the live net into the tracked copy. Moved into tasks/build_night/sync_net.sh so it can
  #     run on every tick rather than only here; this call stays as a belt for unattended runs.
  bash tasks/build_night/sync_net.sh 2>&1 | sed "s/^/[$(ts)] /"
  # (old inline sync kept below for the record of what it reported)
  python3 -c "
import json,pathlib
live=json.load(open('CONSTRAINT_NETS.json'))['bodytwin']['tissue_constraint_net']
p=pathlib.Path('data/CONSTRAINT_NET_TISSUE.json'); d=json.loads(p.read_text())
before=len(d['bodytwin']['tissue_constraint_net']['edges'])
d['bodytwin']['tissue_constraint_net']=live
p.write_text(json.dumps(d,indent=2,ensure_ascii=False))
print(f'nat synkat: {before} -> {len(live[chr(34)+chr(34)]) if False else len(live[\"edges\"])} kanter')
" 2>&1 | sed "s/^/[$(ts)] /"

  # 2e. Harvest finished creative jobs into proposals. Reading them is not processing them: a
  #     number that stays in a report changes nothing, which is the same defect as the 131
  #     cited references nobody harvested. The gate is unit agreement AND a plausible range,
  #     because unit agreement alone attached 801 mm to axial_length from an incision job.
  python3 tasks/assembly/harvest_creative.py 2>&1 | tail -3 | sed "s/^/[$(ts)] /"

  # 3. Prune the queue. Measured 2026-10-03: 12837 of 13919 entries already had a RESULTS.md, so the
  #    queue was 92 % completed work that nothing removed, plus 35 malformed lines from a heredoc that
  #    ran away. The dispatcher re-read and skipped all of it every cycle, and it made freshly queued
  #    briefs look buried at position 915 when the real backlog was 1082 jobs.
  python3 - <<'PRUNE' 2>&1 | sed "s/^/[$(ts)] /"
import pathlib, re
CREATIVE = re.compile(r'BT-(ANOM|OBST|NET|2ND|CAP|IDEA|XSEED|CONN)')
q = pathlib.Path('tasks/lanes/bt_queue.txt')
if q.exists():
    lines = q.read_text().splitlines()
    keep, done, broken = [], 0, 0
    for l in lines:
        p = l.split()
        if len(p) < 3 or ':' in p[-1] or not re.match(r'^[A-Z]$', p[0]):
            broken += 1
            continue
        if (pathlib.Path('results') / p[-1] / 'RESULTS.md').exists():
            done += 1
            continue
        keep.append(l)
    if done or broken:
        q.write_text('\n'.join(keep) + '\n')
    # 3b. Creative families to the front, every cycle. Measured 2026-10-04 12:48: three briefs the
    #     generator had just prepended sat at position 1177, and the front of the queue was
    #     BT-DW48-AUTO again. Something between the generator and the dispatcher rotates the file, so
    #     a one-time reorder does not hold -- the ordering has to be reasserted on every cycle.
    cre = [l for l in keep if CREATIVE.search(l)]
    rest = [l for l in keep if not CREATIVE.search(l)]
    if cre and keep[:len(cre)] != cre:
        keep = cre + rest
        q.write_text('\n'.join(keep) + '\n')
        print(f'creative first: {len(cre)} of {len(keep)} (queue front was {rest[0].split()[-1] if rest else "-"})')
    print(f'queue {len(lines)} -> {len(keep)} ({done} complete, {broken} broken removed)')
PRUNE

  # 2b. Does every number an edge quotes still come out of the script that owns it? Both sides were
  #     edited repeatedly on 2026-10-04 and the first run of this gate found 16 quoted numbers across
  #     three edges that no script produced. A stale quoted number reads as evidence, so it is worse
  #     than a missing one.
  python3 tasks/assembly/decision_edge_consistency.py 2>&1 | sed "s/^/[$(ts)] /"

  # 3. Throughput, so a dead swarm is visible in the same log rather than needing a separate look.
  n=$(find results -name RESULTS.md -mmin -60 2>/dev/null | wc -l)
  echo "[$(ts)] swarm: $n reports in the last hour"
  [ "$n" -lt 5 ] && echo "[$(ts)] ALERT swarm below 5 reports/hour" | tee -a tasks/build_night/ALERTS.log

  echo "[$(ts)] cycle end"
} >> "$LOG" 2>&1
tail -8 "$LOG"
