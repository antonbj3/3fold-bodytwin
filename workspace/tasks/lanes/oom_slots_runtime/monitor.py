"""15-minute kernel/cgroup evidence, sampled while all three queues run."""
from pathlib import Path
import concurrent.futures, datetime, json, os, subprocess, time
ROOT=Path('')
OUT=ROOT/'tasks/lanes/oom_slots_runtime'
SSH=['ssh','-o','ConnectTimeout=12','-i','local_config_path/hunt_20260923','-o','IdentitiesOnly=yes','-o','UserKnownHostsFile=external_research_path','-o','BatchMode=yes']
BEGAN=time.time()
SINCE=datetime.datetime.fromtimestamp(BEGAN,datetime.timezone.utc).strftime('%Y-%m-%d %H:%M:%S UTC')
REMOTE='''from pathlib import Path
import json,subprocess,time
ns={"__name__":"snapshot"};exec(Path('/opt/agents/launch_reserved.py').read_text(),ns)
s=ns['snapshot'](1500,SLICEHEADROOM,HEADROOM,CAP)
s['epoch']=time.time()
s['memory_events']=dict((k,int(v)) for k,v in (line.split() for line in Path('/sys/fs/cgroup/research.slice/memory.events').read_text().splitlines()))
journal=subprocess.run(['sudo','journalctl','-k','--since',SINCE,'--no-pager'],capture_output=True,text=True,check=True).stdout
s['kernel_oom_lines']=[l for l in journal.splitlines() if any(t in l.lower() for t in ['oom-kill:','out of memory: killed','invoked oom-killer'])]
s['research_oom_lines']=[l for l in s['kernel_oom_lines'] if 'research.slice' in l]
s['agents']=subprocess.check_output(['systemctl','list-units','--no-legend','--state=running','agent-*'],text=True).splitlines()
print(json.dumps(s))
'''
HOSTS={'ovh':('ubuntu@51.77.110.4',4096,17,2048),'upcloud':('root@212.147.226.179',1536,3,1536)}
def remote(name):
 host,head,cap,slice_head=HOSTS[name]
 code=REMOTE.replace('SLICEHEADROOM',str(slice_head)).replace('HEADROOM',str(head)).replace('CAP',str(cap)).replace('SINCE',repr(SINCE))
 r=subprocess.run(SSH+[host,'python3 -'],input=code,text=True,capture_output=True,timeout=45)
 if r.returncode:return {'error':r.stderr,'returncode':r.returncode}
 return json.loads(r.stdout)
def local():
 completed=[f'BT-FW48-SEED-{i:03}' for i in range(1,501) if (ROOT/'results'/f'BT-FW48-SEED-{i:03}'/'RESULTS.md').is_file()]
 claims=[]
 for p in (ROOT/'results').glob('BT-*/.local_claim'):
  try:os.kill(int(p.read_text().strip()),0);claims.append(p.parent.name)
  except (OSError,ValueError):pass
 args=subprocess.check_output(['ps','-eo','pid,ni,args'],text=True)
 return {'seed_completed':completed,'local_claims':claims,'memory':next(l for l in Path('/proc/meminfo').read_text().splitlines() if l.startswith('MemAvailable:')),'psi':Path('/proc/pressure/memory').read_text(),'local_bt_runtimes':[l for l in args.splitlines() if 'opencode run ' in l and '--title BT-' in l], 'services':subprocess.check_output(['systemctl','--user','is-active','bt-queue-ovh.service','bt-queue-upcloud.service','bt-queue-local.service','field-ovh.service','field-upcloud.service'],text=True).splitlines()}
(OUT/'verification_start.json').write_text(json.dumps({'epoch':BEGAN,'since':SINCE,'duration_s':930},indent=2)+'\n')
with (OUT/'verification.jsonl').open('w',buffering=1) as f:
 while True:
  row={'elapsed_s':round(time.time()-BEGAN,1),'epoch':time.time()}
  with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
   results=dict(zip(HOSTS,pool.map(remote,HOSTS)))
  row.update(results);row['local']=local();f.write(json.dumps(row)+'\n')
  print(json.dumps({'elapsed_s':row['elapsed_s'],'hosts':{k:{a:v.get(a) for a in ('total','slots','slice_current','reserved','memory_events','error')} for k,v in results.items()},'seed_done':len(row['local']['seed_completed']),'local_owned':len(row['local']['local_claims'])}),flush=True)
  if time.time()-BEGAN>=930:break
  time.sleep(30)
