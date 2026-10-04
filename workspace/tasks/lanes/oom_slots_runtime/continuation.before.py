"""Fixed-code free continuation: priorities, bounded rate recovery and productive ramp."""
from pathlib import Path
import collections,fcntl,hashlib,json,os,subprocess,sys,time
import operations
R=Path(__file__).resolve().parent;A=R.parent
F=Path('../3fold-motion-engine/_private/romi_collab')
sys.path.insert(0,str(F/'lanes'))
from field_queue_store import intake,append_jobs,read_queue,enabled
def save(p,x):
 t=p.with_suffix('.next');t.write_text(json.dumps(x,ensure_ascii=False,indent=2)+'\n');t.replace(p)

def ramp(host):
 code='''from pathlib import Path
import json,time,fcntl,contextlib,io
root=Path('/opt/agents');p=root/'FREE_MODEL_POLICY.json';statefile=root/'SWARM_RAMP_STATE.json'
state=json.loads(statefile.read_text()) if statefile.exists() else {'last_rate':time.time(),'last_change':time.time()}
ns={'__name__':'snapshot'};stream=io.StringIO()
with contextlib.redirect_stdout(stream):exec(STATUS_CODE,ns)
status=json.loads(stream.getvalue());fresh=[]
events=root/'SWARM_RATE_EVENTS.jsonl'
if events.exists():
 for line in events.read_text().splitlines()[-500:]:
  try:
   row=json.loads(line)
   if row['time']>state['last_rate']:fresh.append(row)
  except (ValueError,KeyError):pass
with (root/'rate_policy.lock').open('a') as lock:
 fcntl.flock(lock,fcntl.LOCK_EX);policy=json.loads(p.read_text());cap=policy.get('bunny_active_cap',4);now=time.time()
 productive=sum(w['selected']=='swarm' and w['tool_calls']>0 and not w['network_errors'] and w['log_age_s']<600 for w in status['workers'])
 if fresh:
  cap=max(2,cap//2);state.update(last_rate=max(x['time'] for x in fresh),last_change=now)
  policy['bunny_backoff_until']=max(policy.get('bunny_backoff_until',0),now+1800)
 elif now>=policy.get('bunny_backoff_until',0) and now-state['last_change']>=180 and productive>=max(1,cap//2) and status['memory_available_mib']>=RESERVE:
  cap=min(CEILING,cap+STEP);state['last_change']=now
 policy['bunny_active_cap']=cap;p.write_text(json.dumps(policy,indent=2)+'\\n')
state.update(time=time.time(),cap=cap,productive_bunny=productive,cooldown_until=policy.get('bunny_backoff_until',0),active=status['active_workers'],models=status['by_model'],memory_mib=status['memory_available_mib'])
statefile.write_text(json.dumps(state,indent=2)+'\\n');print(json.dumps(state))
'''.replace('STATUS_CODE',repr((A/'HOME_EGRESS/worker_status.py').read_text())).replace('RESERVE',str(4096 if host=='ovh' else 1300)).replace('CEILING',str(42 if host=='ovh' else 10)).replace('STEP',str(3 if host=='ovh' else 1))
 return json.loads(operations.remote(host,code))

def priorities():
 p=R/'COVERAGE_PRIORITY';marker=p/'COMPLETE.json';target=R/'COVERAGE_INTAKE.json'
 if target.exists() or not marker.exists():return
 m=json.loads(marker.read_text());assert m['results_sha256']==hashlib.sha256((p/'RESULTS.md').read_bytes()).hexdigest()
 proposals=json.loads((p/'PRIORITIES.json').read_text())
 # Eleven exact source-bound proposals reviewed by the coordinator this turn.
 assert len(proposals)==11 and all(j['id'].startswith('FB_SOL_COVER_') for j in proposals)
 for j in proposals:
  j['body']='Use prior results as hypotheses. Build one bounded constructive test; preserve strong controls and negative outcomes. One CPU thread, <=60 s per numerical invocation; no model calls, agents, publishing or source-repo changes.\n'+j['body']
  j['origin']='SOL_COVERAGE_ROOT_REVIEWED_20260927'
 q=R/'COVERAGE_PROPOSALS.json';save(q,proposals);added=intake(q)
 save(target,{'time':time.time(),'added':added,'scientific_admission':False})
 review=A/'SOL_REVIEW_20260927/PLANNER_FEEDBACK'
 for name in ['COVERAGE.json','PRIORITIES.json','PLANNER_FEEDBACK.md']:
  (review/('ALL_WAVES_'+name)).write_bytes((p/name).read_bytes())
 latest=review/'LATEST.md';text=latest.read_text();tag='## All-wave expansion coverage, updated 2026-09-27'
 if tag not in text:latest.write_text(text+'\n\n'+tag+'\n\nRead ALL_WAVES_COVERAGE.json and ALL_WAVES_PRIORITIES.json. Distinguish cataloguing, dispositions, actual results and independent review. Recover missing high-value mechanisms and combine distant compatible operations. These priorities are planning advice, not scientific acceptance.\n\n'+(p/'PLANNER_FEEDBACK.md').read_text())

def recover_collected():
 receipt=R/'RATE_RECOVERY_LEDGER.json';state=json.loads(receipt.read_text()) if receipt.exists() else {}
 queue=read_queue();known={j['id'] for j in queue};added=[]
 # Only field jobs: BodyTwin/dental keep their own queue and graph bookkeeping.
 for mark in (F/'build').glob('FB_*/RATE_INTERRUPTION.json'):
  d=mark.parent
  if (d/'RESULTS.md').exists() or (d/'.ovh_claim').exists() or d.name in state:continue
  meta=d/'JOB.json'
  if not meta.exists():continue
  try:job=json.loads(meta.read_text());why=json.loads(mark.read_text())
  except ValueError:continue
  if time.time()<why.get('next_eligible_epoch',0):continue
  depth=job.get('rate_recovery_depth',0)
  if depth>=2:state[d.name]={'status':'needs_source_review_after_two_rate_recoveries'};continue
  jid=('FB_BPLAN_RATE_' if d.name.startswith('FB_BPLAN') else 'FB_RATE_')+hashlib.sha256(d.name.encode()).hexdigest()[:16]
  if jid in known or (F/'build'/jid).exists():state[d.name]={'job':jid,'status':'already_recorded'};continue
  j=dict(job,id=jid,model='swarm',preferred_model=job.get('preferred_model','swarm'),rate_recovery_depth=depth+1,rate_recovery_of=d.name)
  j['body']='Infrastructure continuation of '+d.name+'. Rate-limited attempts are not negative scientific evidence. Inspect supplied original sources, preserve any actual partial work, and produce the smallest decisive result. Job ID '+jid+'.\n'+job['body']
  # Retain small partial code/report files without recursively copying source trees/logs.
  packet=R/'RATE_PARTIALS'/jid;packet.mkdir(parents=True,exist_ok=True)
  for f in d.iterdir():
   if f.is_file() and f.suffix in ('.py','.md','.json','.csv','.cpp','.h') and f.name not in ('JOB.json','BRIEF.md') and f.stat().st_size<250000:(packet/f.name).write_bytes(f.read_bytes())
  j['src']=dict(job.get('src',{}),partial=str(packet));j['priority']=max(85,min(98,float(job.get('priority',50))))
  made=append_jobs([j]);known.update(made);added+=made;state[d.name]={'time':time.time(),'job':jid,'status':'queued_same_experiment'}
 save(receipt,state);return added

def main():
 if not enabled():return
 with (R/'continuation.lock').open('a') as lock:
  try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
  except BlockingIOError:return
  priorities();state={'time':time.time(),'paid_calls':0,'hosts':{},'new_rate_recoveries':recover_collected()}
  for h in operations.HOSTS:
   try:state['hosts'][h]=ramp(h)
   except Exception as e:state['hosts'][h]={'error':str(e)[:300]}
  # Home inference remains paused until a successful bounded home-route probe.
  policy=A/'MODEL_BACKOFF/POLICY.json';data=json.loads(policy.read_text())
  events=A/'MODEL_BACKOFF/SWARM_RATE_EVENTS.jsonl'
  if events.exists():
   rows=[json.loads(x) for x in events.read_text().splitlines()[-100:]]
   if rows:data['bunny_until']=max(data.get('bunny_until',0),max(x['time'] for x in rows)+1800)
  save(policy,data);save(R/'CONTINUATION_STATUS.json',state);print(json.dumps(state))

if __name__=='__main__':main()
