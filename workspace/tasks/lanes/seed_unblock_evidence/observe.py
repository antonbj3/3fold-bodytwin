from pathlib import Path
import collections,datetime,hashlib,json,re,subprocess,time
B=Path(__file__).resolve().parent
code='''from pathlib import Path
import collections,fcntl,hashlib,json,os,sys,time,subprocess
root=Path('/opt/agents');sys.path.insert(0,str(root));import slot_markers
with (root/'model_slots.lock').open('a') as lock:
 fcntl.flock(lock,fcntl.LOCK_EX);rows,cleaned=slot_markers.scan_locked(root)
 marker_count=len(list((root/'active_models').glob('*.json')))
 runtimes=[]
 for p in Path('/proc').iterdir():
  if not p.name.isdigit():continue
  try:
   args=[a.decode(errors='replace') for a in (p/'cmdline').read_bytes().split(b'\\0') if a]
   if len(args)>1 and args[0].endswith('/opencode') and args[1]=='run' and '--model' in args and '--dir' in args:
    fields=(p/'stat').read_text().rsplit(')',1)[1].split();runtimes.append({'pid':int(p.name),'ppid':int(fields[1]),'job':Path(args[args.index('--dir')+1]).name,'model':args[args.index('--model')+1]})
  except (OSError,ValueError,IndexError):continue
 rowpids={r['pid'] for r in rows}
 policy=json.loads((root/'FREE_MODEL_POLICY.json').read_text())
 out={'time':time.time(),'uid':os.geteuid(),'markers':marker_count,'rows':rows,'models':dict(collections.Counter(r['selected'] for r in rows)),'cleanup':cleaned,'runtimes':runtimes,'orphan_runtimes':[r for r in runtimes if r['ppid'] not in rowpids],'marker_setup_gaps':[r['job'] for r in rows if r['pid'] not in {p['ppid'] for p in runtimes}],'egress_sha256':hashlib.sha256((root/'EGRESS_POOL.json').read_bytes()).hexdigest(),'code_sha256':{f:hashlib.sha256((root/f).read_bytes()).hexdigest() for f in ['route_runtime.py','egress_pool.py','fallback_reserve.py','run_profile_ovh.py','slot_markers.py']},'policy':{k:policy.get(k) for k in ['swarm_worker_value_only','swarm_worker_reserve_account','swarm_worker_reserve_cap','swarm_worker_weekly_used_pct','swarm_worker_weekly_stop_pct','swarm_worker_reserve_until','swarm_worker_backoff_until','swarm_backoff_until','swarm_worker_backoff_until']}}
 print(json.dumps(out,indent=2))
'''
for name,host in [('ovh','ubuntu@51.77.110.4'),('upcloud','root@212.147.226.179')]:
 r=subprocess.run(['ssh',host,'sudo -u ubuntu python3 -'],input=code,text=True,capture_output=True,check=True);out=json.loads(r.stdout)
 assert out['markers']==len(out['rows'])
 assert out['egress_sha256']==hashlib.sha256((B/(name+'.before.EGRESS_POOL.json')).read_bytes()).hexdigest()
 for f,sha in out['code_sha256'].items():assert sha==hashlib.sha256((B/(f if f=='slot_markers.py' else name+'.'+f)).read_bytes()).hexdigest(),f
 (B/(name+'.observation.json')).write_text(r.stdout)
 print(name,'markers',out['markers'],'models',out['models'],'orphans',out['orphan_runtimes'],'setup gaps',out['marker_setup_gaps'],'policy',out['policy'])
now=time.time();before=json.loads((B/'seed.before.json').read_text());q=Path('tasks/lanes/bt_queue.txt').read_text();ids=sorted(set(re.findall(r'BT-FW48-SEED-\d+',q)));done=sorted(p.name for p in Path('results').glob('BT-FW48-SEED-*') if (p/'RESULTS.md').exists());briefs=sorted(p.stem for p in Path('tasks/free48').glob('BT-FW48-SEED-*.md'))
assert set(before['queue_seed_ids'])<=set(ids) and set(before['done_ids'])<=set(done) and set(before['brief_ids'])<=set(briefs)
after=dict(time=now,queue_seed_ids=ids,done_ids=done,brief_ids=briefs,seed253_queued='BT-FW48-SEED-253' in ids,new_done=sorted(set(done)-set(before['done_ids'])))
(B/'seed.after.json').write_text(json.dumps(after,indent=2));print('SEED queue',len(ids),'done',len(done),'new',after['new_done'])
start=json.loads((B/'measurement_time.json').read_text())['epoch'];day=datetime.datetime.fromtimestamp(start).date();metric={}
for name in ['ovh','upcloud']:
 events=[]
 # Logs only contain clock time; read backward through today's tail so older
 # days with the same HH:MM:SS cannot inflate the before/after counts.
 loglines=Path('tasks/lanes',name+'_agents','queue_ovh.log').read_text().splitlines()
 cursor=now;eventday=day
 for line in reversed(loglines):
  clock=re.match(r'\[(\d\d:\d\d:\d\d)\]',line)
  if not clock:continue
  epoch=datetime.datetime.combine(eventday,datetime.time.fromisoformat(clock[1])).timestamp()
  if epoch>cursor+60:
   eventday-=datetime.timedelta(days=1)
   epoch=datetime.datetime.combine(eventday,datetime.time.fromisoformat(clock[1])).timestamp()
  cursor=epoch
  if epoch<start-600:break
  m=re.match(r'\[(\d\d:\d\d:\d\d)\] (ovhstart|fetched) (\S+)(.*)',line)
  if m and epoch<=now:events.append({'epoch':epoch,'kind':m[2],'job':m[3],'detail':m[4]})
 events.reverse()
 def window(lo,hi):
  es=[e for e in events if lo<=e['epoch']<hi];starts=[e for e in es if e['kind']=='ovhstart'];return {'all_starts':len(starts),'seed_starts':sum('SEED-' in e['job'] for e in starts),'fetched_results':sum(e['kind']=='fetched' and 'results=yes' in e['detail'] for e in es),'fetched_exit75':sum(e['kind']=='fetched' and 'exit=75' in e['detail'] for e in es)}
 metric[name]={'before_10min':window(start-600,start),'after_so_far':window(start,min(now,start+600)),'events':events}
metrics={'measurement_start':start,'elapsed_s':now-start,'metrics':metric};(B/'metrics.json').write_text(json.dumps(metrics,indent=2));print('measurement',round(now-start),'sec', {n:{k:v for k,v in m.items() if k!='events'} for n,m in metric.items()})
