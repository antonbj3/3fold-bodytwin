#!/usr/bin/env python3
"""Read per-session cloud usage.cost_usd without exposing auth or event bodies."""
import json
import datetime as dt
from pathlib import Path
import sys
import urllib.error
from collect_cloud import api
ROOT=Path(__file__).resolve().parent
OUT=Path(sys.argv[1]) if len(sys.argv)>1 else ROOT.parents[1]/'results'/'CX-CLOUDWAVE2'/'session_costs.json'

def main():
 get=api()
 rows=json.loads(OUT.read_text()) if OUT.exists() else {}
 for line in (ROOT/'receipts.jsonl').read_text().splitlines():
  r=json.loads(line);session=r.get('session_id');lane=r['lane']
  if not session:continue
  try:x=get('/v1/sessions/'+session)
  except urllib.error.HTTPError as error:
   rows[lane]={'session_id':session,'error_http':error.code};continue
  usage=(x.get('external_metadata') or {}).get('usage') or {}
  rows[lane]={'session_id':session,'cost_usd':usage.get('cost_usd'),
              'observed_at':dt.datetime.now(dt.timezone.utc).isoformat(),
              'model_id':(x.get('session_context') or {}).get('model') or (x.get('external_metadata') or {}).get('last_served_model'),
              'session_status':x.get('session_status')}
 OUT.parent.mkdir(exist_ok=True)
 OUT.write_text(json.dumps(rows,indent=2)+'\n')
 print('sessions',len(rows),'with_cost',sum(isinstance(v.get('cost_usd'),(int,float)) for v in rows.values()),
       'total_cost_usd',round(sum(v.get('cost_usd') or 0 for v in rows.values()),6))
if __name__=='__main__':main()
