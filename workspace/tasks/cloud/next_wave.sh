#!/usr/bin/env bash
# Launch the next N unstarted, preregistered public lanes sequentially.
set -euo pipefail
cd 
count="${1:-2}"
minute="$(date -u +%M)"
if [[ "$minute" != 00 && "$minute" != 30 ]]; then
  echo 'Shared credit may be read only in the :00/:30 UTC minute' >&2
  exit 2
fi
if [[ "$count" != auto ]] && { ! [[ "$count" =~ ^[0-9]+$ ]] || (( count < 1 || count > 12 )); }; then
  echo 'usage: tasks/cloud/next_wave.sh N (1..12) or auto' >&2
  exit 2
fi
# One fresh credit check per small batch. The usage endpoint rate-limits rapid
# requests; if it does, stop closed and let the coordinator retry later.
NEXT_WAVE_REQUESTED_COUNT="$count" python - <<'PY'
import sys
import urllib.error
import datetime as dt
import json
import fcntl
import os
from pathlib import Path
sys.path.insert(0,'tasks/cloud')
from collect_cloud import api
try:
    usage=api()('/api/oauth/usage')
except urllib.error.HTTPError as error:
    raise SystemExit(f'Usage check unavailable (HTTP {error.code}); retry later')
credit=next((v for v in usage.values() if isinstance(v,dict) and v.get('limit_dollars')==250),None)
if credit is None or credit.get('locked_reason') or credit.get('remaining_dollars',0)<=0.05:
    raise SystemExit('Promotional credit unavailable or exhausted; stop')
extra=usage.get('extra_usage') or {}
spend=usage.get('spend') or {}
if extra.get('is_enabled') or spend.get('enabled') or (spend.get('used') or {}).get('amount_minor',0)>0:
    raise SystemExit('Paid usage detected or enabled; stop')
with Path('tasks/cloud/credit_log.jsonl').open('a') as out:
    fcntl.flock(out, fcntl.LOCK_EX)
    out.write(json.dumps({'time': dt.datetime.now(dt.timezone.utc).isoformat(),
                          'balance': {k: credit.get(k) for k in ('limit_dollars','used_dollars','remaining_dollars','resets_at','locked_reason')},
                          'source': 'next_wave_preflight'}) + '\n')
requested = os.environ['NEXT_WAVE_REQUESTED_COUNT']
if requested == 'auto':
    now_utc = dt.datetime.now(dt.timezone.utc)
    deadline = dt.datetime(2026, 9, 24, 21, 30, tzinfo=dt.timezone.utc)
    baseline = dt.datetime(2026, 9, 24, 16, 30, tzinfo=dt.timezone.utc)
    linear_target = 139.527881 * max(0, (deadline - now_utc).total_seconds()) / (deadline - baseline).total_seconds()
    gap = credit['remaining_dollars'] - linear_target
    proposed = max(0, min(9, 5 + round(gap / 2.5)))
    proposed = min(proposed, max(0, int(max(0, credit['remaining_dollars'] - 1) / 2)))
    Path('tasks/cloud/next_wave_count.txt').write_text(str(proposed) + '\n')
    print(f'Auto count {proposed}; linear target ${linear_target:.2f}, gap ${gap:+.2f}')
    with Path('tasks/cloud/auto_decisions.jsonl').open('a') as out:
        out.write(json.dumps({'time': now_utc.isoformat(), 'balance_usd': credit['remaining_dollars'],
                              'linear_target_usd': linear_target, 'gap_usd': gap, 'count': proposed}) + '\n')
print(f"Promotional balance: ${credit['remaining_dollars']:.2f}")
PY
if [[ "$count" == auto ]]; then
  count="$(cat tasks/cloud/next_wave_count.txt)"
fi
launched=0
# The priority file is an explicit value-qualified allowlist. Do not fall back
# to every old reserve bundle: many have only self-generated synthetic gates.
for lane in $(cat tasks/cloud/priority_order.txt); do
  bundle="tasks/cloud/lanes/$lane"
  [[ -d "$bundle" ]] || continue
  (( launched < count )) || break
  [[ ! -e "tasks/cloud/$lane.intent.json" ]] || continue
  python tasks/cloud/launch_cloud_lane.py "$lane" "$bundle"
  launched=$((launched+1))
done
echo "Launched $launched of requested $count prepared lanes"
