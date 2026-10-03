import json,os,runpy,sys,tempfile,time
from pathlib import Path
sys.path.insert(0,'/opt/agents')
import route_runtime
job=Path(tempfile.mkdtemp(prefix='LANE_RUNNER_reserve_worker_D_FALLBACK_',dir='/tmp'))
(job/'HIGH_VALUE').touch()
live=json.loads(Path('/opt/agents/FREE_MODEL_POLICY.json').read_text())
now=time.time()
isolated=not all(now<live.get(m+'_backoff_until',0) for m in ('swarm','free_worker'))
if isolated:
 testroot=job/'isolated_policy';testroot.mkdir()
 policy=live.copy()
 for m in ('swarm','free_worker'):policy[m+'_backoff_until']=now+120
 (testroot/'FREE_MODEL_POLICY.json').write_text(json.dumps(policy,indent=2)+'\n')
 (testroot/'active_models').symlink_to('/opt/agents/active_models',target_is_directory=True)
 route_runtime.ROOT=testroot
else:policy=live
(job/'FALLBACK_TEST_CONTEXT.json').write_text(json.dumps({'isolated_policy':isolated,'live_policy_before':live,'effective_policy':policy,'shared_live_slot_markers':True,'requested_account':'A','requested_model':'swarm','production_code':'/opt/agents/run_profile_ovh.py','global_free_backoffs_edited':False},indent=2)+'\n')
sys.argv=['/opt/agents/run_profile_ovh.py','A','swarm','--dir',str(job),'--title','LANE_RUNNER_reserve_worker_D_FALLBACK','Reply exactly reserve_worker_D_FALLBACK_OK. Do not use tools.']
log=(job/'agent.log').open('w')
old=os.dup(1);old2=os.dup(2)
os.dup2(log.fileno(),1);os.dup2(log.fileno(),2)
try:
 try:runpy.run_path('/opt/agents/run_profile_ovh.py',run_name='__main__')
 except SystemExit as e:code=e.code or 0
finally:
 sys.stdout.flush();sys.stderr.flush()
 os.dup2(old,1);os.dup2(old2,2);os.close(old);os.close(old2);log.close()
route=json.loads((job/'MODEL_ROUTE.json').read_text())
events=[json.loads(line) for line in (job/'agent.log').read_text().splitlines() if line.startswith('{')]
response=''.join(e.get('part',{}).get('text','') for e in events if e.get('type')=='text')
out={'job':str(job),'exit':code,'route':route,'response':response,'isolated_policy':isolated,'log':(job/'agent.log').read_text()}
(job/'VERIFICATION.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
assert code==0 and route['account']=='D' and route['selected']=='reserve_worker' and route['requested']=='swarm'
assert 'reserve_worker_D_FALLBACK_OK' in response
