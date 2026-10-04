#!/usr/bin/env python3
"""Build an auditable session/credit ledger; shared balance deltas are not job costs."""
import datetime as dt
import json
from pathlib import Path
ROOT=Path(__file__).resolve().parent
OUT=ROOT.parents[1]/'results'/'CX-CLOUDWAVE2'

def stamp(x):
 return dt.datetime.fromisoformat(x.replace('Z','+00:00'))

def main():
 credits=[json.loads(x) for x in (ROOT/'credit_log.jsonl').read_text().splitlines()]
 receipts=[json.loads(x) for x in (ROOT/'receipts.jsonl').read_text().splitlines()]
 costs_file=OUT/'session_costs.json'
 costs=json.loads(costs_file.read_text()) if costs_file.exists() else {}
 rows=[]
 for receipt in receipts:
  t=stamp(receipt['time'])
  before=max((c for c in credits if stamp(c['time'])<=t),key=lambda c:stamp(c['time']),default=None)
  after=min((c for c in credits if stamp(c['time'])>=t),key=lambda c:stamp(c['time']),default=None)
  lane=receipt['lane']; dest=ROOT.parents[1]/'results'/lane
  per_session=costs.get(lane,{})
  status=json.loads((dest/'cloud_status.json').read_text()) if (dest/'cloud_status.json').exists() else {}
  shared_delta=round(before['balance']['remaining_dollars']-after['balance']['remaining_dollars'],6) if before and after else None
  rows.append({'lane':lane,'session_id':receipt.get('session_id'),'launch_time':receipt['time'],
   'requested_model':receipt.get('requested_model'),'actual_model_id':status.get('model_id'),
   'session_status':status.get('session_status','not_yet_collected'),
   'results_arrived':(dest/'RESULTS.md').is_file(),
   'balance_before_time':before['time'] if before else None,'balance_before_usd':before['balance']['remaining_dollars'] if before else None,
   'balance_after_time':after['time'] if after else None,'balance_after_usd':after['balance']['remaining_dollars'] if after else None,
   'shared_account_delta_usd':shared_delta,
   'individual_cost_usd':per_session.get('cost_usd'),
   'individual_cost_observed_at':per_session.get('observed_at'),
   'cost_note':'individual_cost_usd is API session usage at observation time; ongoing sessions can accrue more. The shared balance bracket includes Field and overlapping BodyTwin sessions and must not be summed.'})
 OUT.mkdir(exist_ok=True)
 (OUT/'session_ledger.json').write_text(json.dumps(rows,indent=2)+'\n')
 print(len(rows),'sessions',sum(r['results_arrived'] for r in rows),'RESULTS returned')
if __name__=='__main__':main()
