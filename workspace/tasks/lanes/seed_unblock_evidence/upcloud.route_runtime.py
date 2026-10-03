"""Bounded free-model concurrency and prompt termination on explicit rate limits."""
from pathlib import Path
import collections,fcntl,json,os,signal,subprocess,time

ROOT=Path('/opt/agents')
import slot_markers

def _active_counts(*, locked=False):
    rows=slot_markers.scan_locked(ROOT)[0] if locked else slot_markers.live_rows(ROOT)
    return collections.Counter((r['account'],r['selected']) for r in rows)

def _high_value(directory):
    """Paid reserve_worker is reserved for synthesis/planner jobs, not bulk AUTO churn."""
    try:
        if (directory/'HIGH_VALUE').exists():return True
        j=json.loads((directory/'JOB.json').read_text())
        if j.get('breakthrough_priority'):return True
    except (OSError,ValueError):pass
    return any(t in directory.name for t in ('PLAN','HARVEST','REVIEW','BRIDGE','CROSS','HUNT','DENSE','-NEXT-','RESERVE'))


def choose_free_models(preference,counts,policy,now,eligible=None):
    # Maximise the free tier: fill free_worker up to its fill target first, then
    # The swarm, and offer both when neither is preferred; each model is bounded by
    # its own backoff and active cap.
    counts_by_model={m:sum(n for (a,k),n in counts.items() if k==m) for m in ('swarm','free_worker')}
    allowed=[m for m in ('swarm','free_worker') if (eligible is None or m in eligible) and now>=policy.get(m+'_backoff_until',0) and counts_by_model[m]<policy.get(m+'_active_cap',1000)]
    if not allowed:return []
    if 'free_worker' in allowed and counts_by_model['free_worker']<policy.get('free_worker_fill_target',0):return ['free_worker']
    if preference in allowed:return [preference]
    return allowed

def reserve_slot(label,requested,directory):
    if requested not in ('swarm','free_worker'):return label,requested,None
    active=ROOT/'active_models';active.mkdir(exist_ok=True)
    slot_markers.install_cleanup()
    began=time.time()
    while time.time()-began<1500:
        with (ROOT/'model_slots.lock').open('a') as lock:
            fcntl.flock(lock,fcntl.LOCK_EX)
            policy=json.loads((ROOT/'FREE_MODEL_POLICY.json').read_text())
            counts=_active_counts(locked=True)
            # Atomic egress snapshot; do not acquire EGRESS_POOL.lock under slots.
            import egress_pool
            rows=slot_markers.scan_locked(ROOT)[0]
            eligible=[m for m in ('swarm','free_worker') if egress_pool.has_available(m,rows)]
            pref=directory/'PREFERRED_MODEL'
            models=choose_free_models(pref.read_text().strip() if pref.exists() else requested,counts,policy,time.time(),eligible=eligible)
            if not models:
                from fallback_reserve import choose as choose_go_reserve
                hv=(not policy.get('reserve_worker_value_only',True)) or _high_value(directory)
                go_account=choose_go_reserve(policy,counts,time.time(),free_unavailable=True) if hv else None
                if go_account is not None:
                    marker=slot_markers.write_locked(ROOT,go_account,'reserve_worker',directory,'free models unavailable: cooldown, capacity or egress backoff')
                    return go_account,'reserve_worker',marker
            accounts=[label]  # Fixed queue profile; never rotate credentials after a limit.
            if models:
                pairs=[(a,m) for a in accounts for m in models]
                account,selected=min(pairs,key=lambda am:(sum(v for (a,m),v in counts.items() if m==am[1]),counts[am],am))
                marker=slot_markers.write_locked(ROOT,account,selected,directory)
                return account,selected,marker
        # Wait without starting an OpenCode runtime or consuming a model slot.
        (directory/'ADMISSION_WAIT.json').write_text(json.dumps({'since':began,'last_checked':time.time(),'reason':'free model cooldown or concurrency ceiling'}))
        time.sleep(10)
    (directory/'RATE_INTERRUPTION.json').write_text(json.dumps({'time':time.time(),'reason':'admission wait expired','next_eligible_epoch':time.time()+900}))
    raise SystemExit(75)

def apply_route(env,routing,model_id,directory):
    mode=routing.get('model_routes',{}).get(model_id,'home_via_ssh')
    if not routing.get('enabled'):mode='direct'
    proxy=None
    # Egress pool: spread jobs of the same model over several egress IPs, each an
    # independent per-IP bucket. Falls back to the static route when disabled.
    try:
        import egress_pool
        chosen=egress_pool.pick(model_id,directory)
        if chosen:
            mode,proxy=chosen
            try:
                egress_pool.prepare_cache(mode,env,model_id)
            except Exception:
                egress_pool.record(model_id,mode)
                chosen=egress_pool.pick(model_id,directory)
                mode,proxy=chosen if chosen else ('direct',None)
    except Exception:
        proxy=None
    for k in ('HTTP_PROXY','HTTPS_PROXY','ALL_PROXY','http_proxy','https_proxy','all_proxy'):env.pop(k,None)
    env.update(NO_PROXY='localhost,127.0.0.1,::1',no_proxy='localhost,127.0.0.1,::1')
    if proxy:env.update(HTTPS_PROXY=proxy,https_proxy=proxy)
    elif mode=='home_via_ssh':env.update(HTTPS_PROXY=routing['https_proxy'],https_proxy=routing['https_proxy'])
    elif mode!='direct':mode='direct'  # unknown egress without a proxy: never crash, go direct
    (directory/'EGRESS_ROUTE.json').write_text(json.dumps({'mode':mode,'proxy':proxy,'time':time.time(),'model':model_id})+'\n')

def has_rate_error(runtime,began):
    for p in (runtime/'data/opencode/log').glob('*'):
        if p.stat().st_mtime<began-1:continue
        with p.open(errors='replace') as stream:
            stream.seek(max(0,p.stat().st_size-100000))
            for line in stream:
                if 'Rate limit exceeded' in line and ('AI_APICallError' in line or line.startswith('ERROR')):return True
    return False

def stop_child(process):
    try:os.killpg(process.pid,signal.SIGTERM)
    except ProcessLookupError:return
    try:process.wait(timeout=5)
    except subprocess.TimeoutExpired:
        try:os.killpg(process.pid,signal.SIGKILL)
        except ProcessLookupError:pass
        process.wait(timeout=5)

def run_guarded(command,*,cwd,env,runtime):
    selected='free_worker' if 'opencode/free_worker-2.5-preview-free' in command else 'swarm' if 'opencode/space-swarm-free' in command else 'reserve_worker' if 'opencode-go/reserve_worker-v4.1-flash' in command else None
    if selected is None:
        return subprocess.run(command,cwd=cwd,env=env)
    began=time.time();progress=Progress(Path(cwd)/"agent.log");process=subprocess.Popen(command,cwd=cwd,env=env,start_new_session=True)
    try:
        while process.poll() is None:
            stalled=progress.stalled(time.time()-began)
            if stalled or has_rate_error(runtime,began):
                row={'time':time.time(),'job':Path(cwd).name,'reason':'no completed tools for six minutes' if stalled else 'explicit provider rate limit','egress':(json.loads((Path(cwd)/'EGRESS_ROUTE.json').read_text()).get('mode','direct') if (Path(cwd)/'EGRESS_ROUTE.json').exists() else 'direct'),'model':selected,'next_eligible_epoch':time.time()+1800}
                metadata=Path(cwd)/'MODEL_ROUTE.json'
                if metadata.exists():row['account']=json.loads(metadata.read_text()).get('account')
                with (ROOT/(selected.upper()+'_RATE_EVENTS.jsonl')).open('a') as stream:
                    fcntl.flock(stream,fcntl.LOCK_EX);stream.write(json.dumps(row)+'\n');stream.flush()
                (Path(cwd)/'RATE_INTERRUPTION.json').write_text(json.dumps(row,indent=2)+'\n')
                suppress=False
                if selected in ('swarm','free_worker'):
                    try:
                        import egress_pool
                        if egress_pool.enabled():
                            # A per-IP limit on one egress must not pause the model; the
                            # next job is routed to a healthy egress instead.
                            suppress=not egress_pool.record(selected,row.get('egress'))
                    except Exception:
                        suppress=False
                if selected in ('swarm','free_worker','reserve_worker') and not suppress:
                    with (ROOT/'rate_policy.lock').open('a') as lock:
                        fcntl.flock(lock,fcntl.LOCK_EX)
                        policy_file=ROOT/'FREE_MODEL_POLICY.json'
                        if policy_file.exists():
                            policy=json.loads(policy_file.read_text())
                            key=selected+'_backoff_until'
                            policy[key]=max(policy.get(key,0),time.time()+1800)
                            tmp=policy_file.with_suffix('.rate.tmp');tmp.write_text(json.dumps(policy,indent=2)+'\n');tmp.replace(policy_file)
                # Only this invocation's explicitly isolated child group is stopped.
                stop_child(process)
                return subprocess.CompletedProcess(command,75)
            time.sleep(3)
        if process.returncode and progress.tools==0:
            suppress=False
            if selected in ('swarm','free_worker'):
                try:
                    import egress_pool
                    route=Path(cwd)/'EGRESS_ROUTE.json'
                    egress=json.loads(route.read_text()).get('mode','direct')
                    if egress_pool.enabled():suppress=not egress_pool.record(selected,egress)
                except Exception:
                    suppress=False
            # Proxy/tunnel startup failures isolate the failed IP just like a 429.
            # reserve_worker and the existing all-IP cooldown behavior remain unchanged.
            if not suppress:
                with (ROOT/'rate_policy.lock').open('a') as lock:
                    fcntl.flock(lock,fcntl.LOCK_EX)
                    p=ROOT/'FREE_MODEL_POLICY.json';d=json.loads(p.read_text())
                    key=selected+'_backoff_until';d[key]=max(d.get(key,0),time.time()+600)
                    tmp=p.with_suffix('.startup-error.tmp');tmp.write_text(json.dumps(d));tmp.replace(p)
            (Path(cwd)/'STARTUP_ERROR.json').write_text(json.dumps({'time':time.time(),'model':selected,'exit':process.returncode,'retry_after_seconds':600}))
        return subprocess.CompletedProcess(command,process.returncode)
    finally:
        if process.poll() is None:stop_child(process)

class Progress:
    def __init__(self,path):
        self.path=path;self.offset=path.stat().st_size if path.exists() else 0;self.tools=0
    def stalled(self,age):
        if self.path.exists():
            with self.path.open(errors='replace') as f:
                f.seek(self.offset)
                for line in f:
                    try:event=json.loads(line)
                    except ValueError:continue
                    if event.get('type')=='tool_use' and event.get('part',{}).get('state',{}).get('status')=='completed':self.tools+=1
                self.offset=f.tell()
        return self.tools==0 and age>360
