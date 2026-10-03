from pathlib import Path
import json,subprocess,time
B=Path(__file__).resolve().parent
code='''from pathlib import Path
import fcntl,json,time,os
root=Path('/opt/agents');assert os.geteuid()!=0
with (root/'rate_policy.lock').open('a') as lock:
 fcntl.flock(lock,fcntl.LOCK_EX);p=root/'FREE_MODEL_POLICY.json';policy=json.loads(p.read_text());before={k:policy.get(k,0) for k in ('bunny_backoff_until','free_worker_backoff_until')};until=time.time()+120
 for k in before:policy[k]=max(before[k],until)
 tmp=p.with_suffix('.verification.tmp');tmp.write_text(json.dumps(policy,indent=2)+'\\n');tmp.replace(p)
 print(json.dumps({'time':time.time(),'before':before,'forced':{k:policy[k] for k in before},'test_until':until,'account':policy['reserve_worker_reserve_account'],'used_pct':policy['reserve_worker_weekly_used_pct'],'stop_pct':policy['reserve_worker_weekly_stop_pct'],'reserve_cap':policy['reserve_worker_reserve_cap']}))
'''
for name,host in [('ovh','ubuntu@51.77.110.4'),('upcloud','root@212.147.226.179')]:
 r=subprocess.run(['ssh',host,'sudo -u ubuntu python3 -'],input=code,text=True,capture_output=True,check=True);(B/(name+'.forced_backoff.json')).write_text(r.stdout);print(name,r.stdout)
(B/'measurement_time.json').write_text(json.dumps({'epoch':time.time()}))
