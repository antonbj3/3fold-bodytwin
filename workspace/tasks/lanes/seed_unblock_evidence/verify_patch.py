from pathlib import Path
import collections, importlib.util, json, os, signal, subprocess, sys, tempfile, time
B=Path(__file__).resolve().parent
sys.path.insert(0,str(B))
import slot_markers as slots

def load(name,path):
 spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
route=load('route_runtime',B/'ovh.route_runtime.py');egress=load('egress_pool',B/'ovh.egress_pool.py');fallback=load('fallback_reserve',B/'ovh.fallback_reserve.py')
sys.modules['egress_pool']=egress;sys.modules['fallback_reserve']=fallback
checks=[]
def check(label,ok):
 assert ok,label
 checks.append(label)
with tempfile.TemporaryDirectory() as td:
 root=Path(td);(root/'active_models').mkdir();proc=root/'proc';proc.mkdir()
 (proc/'cmdline').write_bytes(b'python3\x00/opt/agents/run_profile_ovh.py\x00--dir\x00/opt/agents/jobs/BT-FW48-SEED-999\x00');(proc/'cgroup').write_text('0::/research.slice/agent-BT-FW48-SEED-999-123.service\n')
 now=time.time();row=dict(pid=os.getpid(),started=now,account='D',selected='swarm_worker',job='BT-FW48-SEED-999',start_ticks=99,boot_id='boot')
 realidentity=slots.identity;slots.identity=lambda pid:(proc,99,'boot',now-10)
 check('matching process accepted',slots.classify(row,now)=='live')
 check('reused PID start ticks rejected',slots.classify(dict(row,start_ticks=98),now)=='reused')
 check('different boot rejected',slots.classify(dict(row,boot_id='old'),now)=='reused')
 check('legacy reused PID rejected',slots.classify(dict(row,started=now-20),now)=='reused')
 check('5700 second ceiling enforced',slots.classify(dict(row,started=now-5701),now)=='aged')
 (proc/'cmdline').write_bytes(b'python3\x00unrelated.py\x00')
 check('unrelated process rejected',slots.classify(row,now)=='inconsistent')
 (root/'active_models'/'bad.json').write_text('broken')
 slots.identity=lambda pid: (_ for _ in ()).throw(ProcessLookupError(pid))
 (root/'active_models'/'dead.json').write_text(json.dumps(row))
 with (root/'model_slots.lock').open('a') as lock:
  import fcntl;fcntl.flock(lock,fcntl.LOCK_EX);live,summary=slots.scan_locked(root)
 check('dead and malformed slots removed',not live and not list((root/'active_models').glob('*.json')) and summary=={'malformed':1,'dead':1})
 slots.identity=realidentity
 policy=dict(swarm_worker_reserve_until=now+3600,swarm_worker_reserve_account='D',swarm_worker_weekly_used_pct=0,swarm_worker_weekly_stop_pct=95,swarm_worker_reserve_cap=16)
 check('D reserve without global cooldown',fallback.choose(policy,{},now,free_unavailable=True)=='D')
 check('95 percent hard stop',fallback.choose(dict(policy,swarm_worker_weekly_used_pct=95),{},now,free_unavailable=True) is None)
 check('cannot relax 95 percent threshold',fallback.choose(dict(policy,swarm_worker_weekly_stop_pct=101,swarm_worker_weekly_used_pct=95),{},now,free_unavailable=True) is None)
 check('cap 16 enforced',fallback.choose(policy,{('D','swarm_worker'):16},now,free_unavailable=True) is None)
 check('reserve expiry enforced',fallback.choose(dict(policy,swarm_worker_reserve_until=now-1),{},now,free_unavailable=True) is None)
 check('swarm_worker cooldown enforced',fallback.choose(dict(policy,swarm_worker_backoff_until=now+10),{},now,free_unavailable=True) is None)
 egress.ROOT=root;egress.POOL=root/'EGRESS_POOL.json';egress.STATE=root/'EGRESS_POOL_STATE.json'
 pool=dict(enabled=True,egress={'direct':{'cap':2},'home_via_ssh':{'cap':2},'warp':{'cap':2}});egress.POOL.write_text(json.dumps(pool));egress.STATE.write_text(json.dumps({'backoff':{'swarm|'+n:now+60 for n in pool['egress']}}))
 (root/'active_models'/'99999999.json').write_text(json.dumps(dict(row,pid=99999999,started=time.time())))
 check('egress counter permanently cleans dead slots',egress._active_counts()=={} and not (root/'active_models'/'99999999.json').exists())
 route.ROOT=root
 (root/'active_models'/'99999999.json').write_text(json.dumps(dict(row,pid=99999999,started=time.time())))
 check('route counter permanently cleans dead slots',route._active_counts()=={} and not (root/'active_models'/'99999999.json').exists())
 check('all-IP backoff blocks free admission',not egress.has_available('swarm',[]))
 check('healthy alternate model remains free',egress.has_available('swarm_worker',[]))
 check('healthy free model wins before fallback',route.choose_free_models('swarm',{},dict(swarm_worker_fill_target=4),now,eligible=['swarm'])==['swarm'])
 # Actual reserve call for a non-high-value SEED in all-IP backoff.
 egstate=json.loads(egress.STATE.read_text());egstate['backoff'].update({'swarm_worker|'+n:now+60 for n in pool['egress']});egress.STATE.write_text(json.dumps(egstate))
 route.ROOT=root;slots.ROOT=root;directory=root/'jobs'/'BT-FW48-SEED-999';directory.mkdir(parents=True)
 (root/'FREE_MODEL_POLICY.json').write_text(json.dumps(dict(policy,swarm_worker_value_only=False)))
 account,selected,marker=route.reserve_slot('B','swarm',directory)
 check('ordinary SEED reserves D swarm_worker on egress backoff',(account,selected)==('D','swarm_worker') and marker.exists());slots.release();check('own marker released',not marker.exists())
 # Guard signal cleanup kills only its explicitly isolated child group.
 script=root/'guard_test.py';script.write_text('''import sys,time,json,os
from pathlib import Path
sys.path.insert(0,sys.argv[1])
import slot_markers as slots,importlib.util
spec=importlib.util.spec_from_file_location('route',Path(sys.argv[1])/'ovh.route_runtime.py');route=importlib.util.module_from_spec(spec);spec.loader.exec_module(route)
root=Path(sys.argv[2]);slots.ROOT=root;route.ROOT=root;slots.install_cleanup()
with (root/'model_slots.lock').open('a') as lock:
 import fcntl;fcntl.flock(lock,fcntl.LOCK_EX);slots.write_locked(root,'D','swarm_worker',root/'signal-test')
route.run_guarded([sys.executable,'-c',"import os,time;from pathlib import Path;Path('child.pid').write_text(str(os.getpid()));time.sleep(60)",'opencode-go/swarm_worker-v4.1-flash'],cwd=root,env=os.environ.copy(),runtime=root)
''')
 sentinel=subprocess.Popen([sys.executable,'-c','import time;time.sleep(60)'],start_new_session=True)
 tester=subprocess.Popen([sys.executable,str(script),str(B),str(root)])
 try:
  for _ in range(50):
   if (root/'child.pid').exists():break
   time.sleep(.1)
  childpid=int((root/'child.pid').read_text());tester.send_signal(signal.SIGTERM);tester.wait(timeout=10)
  check('SIGTERM returns 143',tester.returncode==143)
  check('SIGTERM removes own slot',not (root/'active_models'/f'{tester.pid}.json').exists())
  check('SIGTERM kills own child',not Path(f'/proc/{childpid}').exists())
  check('unrelated process group survives',sentinel.poll() is None)
 finally:
  if tester.poll() is None:tester.kill();tester.wait()
  sentinel.terminate();sentinel.wait()

 (root/'child.pid').unlink()
 alarmscript=root/'alarm_test.py'
 alarmscript.write_text(script.read_text().replace("slots.install_cleanup()", "slots.install_cleanup()\nimport signal\ndef deadline(signum,frame):raise SystemExit(124)\nsignal.signal(signal.SIGALRM,deadline);signal.setitimer(signal.ITIMER_REAL,.5)"))
 alarm=subprocess.Popen([sys.executable,str(alarmscript),str(B),str(root)])
 alarm.wait(timeout=10)
 check('deadline returns 124',alarm.returncode==124)
 check('deadline removes own slot',not (root/'active_models'/f'{alarm.pid}.json').exists())
 check('deadline kills own child',not Path('/proc/'+(root/'child.pid').read_text()).exists())
print(json.dumps({'passed':len(checks),'checks':checks},indent=2))
