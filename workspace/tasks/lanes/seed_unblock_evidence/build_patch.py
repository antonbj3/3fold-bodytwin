from pathlib import Path
B=Path(__file__).resolve().parent
helper='''"""Slot identity and cleanup. All scans/mutations hold model_slots.lock."""
from pathlib import Path
import atexit, collections, fcntl, json, os, signal, time
ROOT = Path('/opt/agents')
MAX_AGE = 5700  # systemd RuntimeMaxSec=5400 plus five minutes of cleanup margin
OWN = None

def identity(pid):
    proc = Path('/proc') / str(pid)
    fields = (proc/'stat').read_text().rsplit(')', 1)[1].split()
    if fields[0] in ('Z', 'X'): raise ProcessLookupError(pid)
    ticks = int(fields[19])
    boot = Path('/proc/sys/kernel/random/boot_id').read_text().strip()
    epoch = float(next(l.split()[1] for l in Path('/proc/stat').read_text().splitlines() if l.startswith('btime '))) + ticks/os.sysconf('SC_CLK_TCK')
    return proc, ticks, boot, epoch

def classify(row, now):
    try:
        pid = row['pid']; started = float(row['started']); job = row['job']
        if not isinstance(pid, int) or pid <= 1 or not isinstance(job, str) or not job: return 'malformed'
        if row['selected'] not in ('swarm', 'free_worker', 'reserve_worker') or row['account'] not in ('A','B','C','D'): return 'malformed'
        try: proc, ticks, boot, epoch = identity(pid)
        except (FileNotFoundError, ProcessLookupError): return 'dead'
        if not 0 <= now-started <= MAX_AGE: return 'aged'
        if epoch > started+2: return 'reused'
        if 'start_ticks' in row and (row['start_ticks'] != ticks or row.get('boot_id') != boot): return 'reused'
        args = (proc/'cmdline').read_bytes().split(b'\\0')
        # Slots belong to the profile launcher, never a later process reusing its PID.
        if not any(Path(a.decode(errors='replace')).name == 'run_profile_ovh.py' for a in args): return 'inconsistent'
        if not any(Path(a.decode(errors='replace')).name == job for a in args): return 'inconsistent'
        cg = (proc/'cgroup').read_text()
        if 'agent-' not in cg or job not in cg: return 'inconsistent'
        return 'live'
    except (FileNotFoundError, ProcessLookupError): return 'dead'
    except (OSError, ValueError, TypeError, KeyError): return 'malformed'

def scan_locked(root=None):
    root = root or ROOT
    rows = []; summary = collections.Counter(); now = time.time()
    for p in (root/'active_models').glob('*.json'):
        try: row = json.loads(p.read_text()); status = classify(row, now)
        except (OSError, ValueError): row = {}; status = 'malformed'
        summary[status] += 1
        if status == 'live': rows.append(row)
        else: p.unlink(missing_ok=True)
    return rows, dict(summary)

def live_rows(root=None):
    root = root or ROOT
    with (root/'model_slots.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        return scan_locked(root)[0]

def release():
    global OWN
    if OWN is None: return
    marker, expected = OWN
    with (ROOT/'model_slots.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        try:
            if json.loads(marker.read_text()) == expected: marker.unlink(missing_ok=True)
        except FileNotFoundError: pass
    OWN = None

def install_cleanup():
    # SIGKILL/OOM cannot run handlers; the next locked scan removes those slots.
    def terminate(signum, frame): raise SystemExit(128+signum)
    for sig in (signal.SIGTERM, signal.SIGINT, signal.SIGHUP): signal.signal(sig, terminate)

def write_locked(root, account, selected, directory, reason=None):
    global OWN
    _, ticks, boot, _ = identity(os.getpid())
    row = dict(pid=os.getpid(),account=account,selected=selected,job=directory.name,started=time.time(),start_ticks=ticks,boot_id=boot)
    if reason: row['reason'] = reason
    marker = root/'active_models'/(str(os.getpid())+'.json')
    tmp = marker.with_suffix('.tmp'); tmp.write_text(json.dumps(row)+'\\n'); tmp.replace(marker)
    OWN = marker, row
    return marker

atexit.register(release)
'''
(B/'slot_markers.py').write_text(helper)
for host in ['ovh','upcloud']:
 route=(B/(host+'.before.route_runtime.py')).read_text()
 route=route.replace("ROOT=Path('/opt/agents')", "ROOT=Path('/opt/agents')\nimport slot_markers\n\ndef _active_counts(*, locked=False):\n    rows=slot_markers.scan_locked(ROOT)[0] if locked else slot_markers.live_rows(ROOT)\n    return collections.Counter((r['account'],r['selected']) for r in rows)")
 route=route.replace("def choose_free_models(preference,counts,policy,now):", "def choose_free_models(preference,counts,policy,now,eligible=None):")
 route=route.replace("allowed=[m for m in ('swarm','free_worker') if now>=", "allowed=[m for m in ('swarm','free_worker') if (eligible is None or m in eligible) and now>=")
 route=route.replace("    began=time.time()\n    while", "    slot_markers.install_cleanup()\n    began=time.time()\n    while",1)
 start=route.index('            counts=collections.Counter()')
 end=route.index("            pref=directory/'PREFERRED_MODEL'",start)
 route=route[:start]+'            counts=_active_counts(locked=True)\n'+route[end:]
 start=route.index("            pref=directory/'PREFERRED_MODEL'")
 end=route.index('            if not models:',start)
 route=route[:start]+"""            # Atomic egress snapshot; do not acquire EGRESS_POOL.lock under slots.
            import egress_pool
            rows=slot_markers.scan_locked(ROOT)[0]
            eligible=[m for m in ('swarm','free_worker') if egress_pool.has_available(m,rows)]
            pref=directory/'PREFERRED_MODEL'
            models=choose_free_models(pref.read_text().strip() if pref.exists() else requested,counts,policy,time.time(),eligible=eligible)
"""+route[end:]
 route=route.replace("choose_go_reserve(policy,counts,time.time()) if hv else None","choose_go_reserve(policy,counts,time.time(),free_unavailable=True) if hv else None")
 start=route.index("                    marker=active/(str(os.getpid())+'.json')")
 end=route.index("                    return go_account,'reserve_worker',marker",start)
 route=route[:start]+"                    marker=slot_markers.write_locked(ROOT,go_account,'reserve_worker',directory,'free models unavailable: cooldown, capacity or egress backoff')\n"+route[end:]
 start=route.index("                marker=active/(str(os.getpid())+'.json')")
 end=route.index("                return account,selected,marker",start)
 route=route[:start]+"                marker=slot_markers.write_locked(ROOT,account,selected,directory)\n"+route[end:]
 # Preserve the entire guarded body and add unconditional own-child cleanup.
 start=route.index('    while process.poll() is None:')
 end=route.index('\nclass Progress:',start)
 body=route[start:end].rstrip()
 route=route[:start]+'    try:\n'+''.join('    '+l+'\n' for l in body.splitlines())+'    finally:\n        if process.poll() is None:stop_child(process)\n'+route[end:]
 (B/(host+'.route_runtime.py')).write_text(route)
 egress=(B/(host+'.before.egress_pool.py')).read_text()
 start=egress.index('    for p in active.glob')
 end=egress.index('        route = ROOT',start)
 egress=egress[:start]+"    import slot_markers\n    for row in slot_markers.live_rows(ROOT):\n"+egress[end:]
 # Availability does not spend an egress allocation or change pool entries.
 insert='''def has_available(model, rows):
    """Admission snapshot; no state writes and no reverse lock acquisition."""
    pool = _load(POOL, None)
    if not pool or not pool.get('enabled'): return True
    state = _normalise_state(_load(STATE, {})); now = time.time()
    live = {}
    for row in rows:
        route = _load(ROOT/'jobs'/row['job']/'EGRESS_ROUTE.json', {})
        mode = route.get('mode')
        if mode: live[mode] = live.get(mode, 0)+1
    model = _model_key(model)
    for name, entry in pool.get('egress', {}).items():
        cap = int(entry.get('cap', 0) or 0)
        if _healthy(entry, now) and state.get('backoff', {}).get(f'{model}|{name}', 0) <= now and (not cap or live.get(name, 0) < cap): return True
    return False


'''
 egress=egress.replace('def pick(model, directory=None):',insert+'def pick(model, directory=None):')
 (B/(host+'.egress_pool.py')).write_text(egress)
 fallback=(B/(host+'.before.fallback_reserve.py')).read_text()
 fallback=fallback.replace('Explicit Go fallback; only after BOTH free models enter cooldown.','Explicit Go fallback after free admission is unavailable.')
 fallback=fallback.replace('def choose(policy,counts,now):','def choose(policy,counts,now,*,free_unavailable=False):')
 fallback=fallback.replace("if not all(now<policy.get(m+'_backoff_until',0) for m in ('swarm','free_worker')):return None","if not free_unavailable and not all(now<policy.get(m+'_backoff_until',0) for m in ('swarm','free_worker')):return None")
 fallback=fallback.replace("policy.get('reserve_worker_weekly_stop_pct',101)","min(95,policy.get('reserve_worker_weekly_stop_pct',95))")
 (B/(host+'.fallback_reserve.py')).write_text(fallback)
 runner=(B/(host+'.before.run_profile_ovh.py')).read_text()
 start=runner.index("    env,runtime=environment")
 end=runner.index('    raise SystemExit(code)',start)
 runner=runner[:start]+'''    runtime=None
    try:
        env,runtime=environment(account,model,directory,(directory/'ALLOW_WEB').exists())
        executable=str(ROOT/'.opencode/bin/opencode')
        completed=run_guarded([executable,'run','--pure','--auto','--agent','build','--model',model,
            '--dir',str(directory),'--format','json','--title',args.title,args.prompt],cwd=directory,env=env,runtime=runtime)
        code=completed.returncode
    finally:
        # Own marker/store only, including setup failure, SIGTERM and timeout.
        from slot_markers import release
        release()
        if runtime is not None:shutil.rmtree(runtime,ignore_errors=True)
'''+runner[end:]
 target="    args=parser.parse_args();directory=Path(args.dir).resolve();account,selected,marker=select(args.profile,args.model,directory);model=MODELS[selected]"
 runner=runner.replace(target,"""    args=parser.parse_args();directory=Path(args.dir).resolve()
    # Packet wrapper kills this PID at 1600 s. Expire slightly earlier so the
    # isolated OpenCode group and own slot get graceful cleanup first.
    import signal
    def timeout(signum,frame):raise SystemExit(124)
    signal.signal(signal.SIGALRM,timeout);signal.setitimer(signal.ITIMER_REAL,1575)
    account,selected,marker=select(args.profile,args.model,directory);model=MODELS[selected]""")
 runner=runner.replace('        release()','        release()\n        signal.setitimer(signal.ITIMER_REAL,0)')
 (B/(host+'.run_profile_ovh.py')).write_text(runner)
for p in B.glob('*.py'): compile(p.read_text(),str(p),'exec')
print('compiled all patches')
