from pathlib import Path
import subprocess,json,time
B=Path(__file__).resolve().parent
code='''from pathlib import Path
import fcntl,json,os,time,sys,hashlib
root=Path('/opt/agents');assert os.geteuid()!=0
with (root/'model_slots.lock').open('a') as slots,(root/'rate_policy.lock').open('a') as lock:
 fcntl.flock(slots,fcntl.LOCK_EX);fcntl.flock(lock,fcntl.LOCK_EX)
 p=root/'FREE_MODEL_POLICY.json';policy=json.loads(p.read_text());before=dict(policy)
 assert policy['swarm_worker_reserve_cap']==3 and policy['swarm_worker_weekly_stop_pct']==95 and policy['swarm_worker_reserve_account']=='D'
 (root/('FREE_MODEL_POLICY.before-seed-unblock-'+str(int(time.time()))+'.json')).write_bytes(p.read_bytes())
 policy['swarm_worker_value_only']=False
 policy['swarm_worker_reserve_authorization']='Anton LANE_RUNNER_SEED_UNBLOCK 2026-09-30: ordinary SEED and other queued jobs may use D reserve when no free admission (cooldown/cap/egress); 95% stop, existing stricter host cap 3 and window preserved'
 policy['allocation']='Available free model/egress first; swarm_worker reserve on unavailable free admission'
 tmp=p.with_suffix('.seed-unblock.tmp');tmp.write_text(json.dumps(policy,indent=2)+'\\n');tmp.replace(p)
 sys.path.insert(0,str(root));import slot_markers
 rows,summary=slot_markers.scan_locked(root)
 print(json.dumps({'time':time.time(),'uid':os.geteuid(),'cleanup':summary,'live_after':len(rows),'policy_changes':{k:{'before':before.get(k),'after':v} for k,v in policy.items() if before.get(k)!=v},'egress_sha256':hashlib.sha256((root/'EGRESS_POOL.json').read_bytes()).hexdigest()},indent=2))
'''
r=subprocess.run(['ssh','root@212.147.226.179','sudo -u ubuntu python3 -'],input=code,text=True,capture_output=True,check=True);(B/'upcloud.deploy.json').write_text(r.stdout);print(r.stdout)
(B/'deploy_time.json').write_text(json.dumps({'epoch':time.time()}))
