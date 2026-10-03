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

  # 3. Throughput, so a dead swarm is visible in the same log rather than needing a separate look.
  n=$(find results -name RESULTS.md -mmin -60 2>/dev/null | wc -l)
  echo "[$(ts)] swarm: $n reports in the last hour"
  [ "$n" -lt 5 ] && echo "[$(ts)] ALERT swarm below 5 reports/hour" | tee -a tasks/build_night/ALERTS.log

  echo "[$(ts)] cycle end"
} >> "$LOG" 2>&1
tail -8 "$LOG"
