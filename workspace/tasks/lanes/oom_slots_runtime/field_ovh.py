#!/usr/bin/env python3
"Field → OVH: runs Space The swarm jobs from SWARM_QUEUE_20260924.json on OVH (51.77.110.4).\n\ndirectories. Everything else stays with the local refill.py. Unit names agent-F_<id>-<ts>, directory\n/opt/agents/jobs/F_<id>. Cap: lanes/ovh_field_cap (default 12). Stop: kill this PID."
import json, os, re, shutil, subprocess, sys, time
from pathlib import Path

def _private_input_pattern(public_pattern):
    """Combine public exclusions with required operator-maintained private terms."""
    import json
    from pathlib import Path
    import re
    config = Path.home() / ".bodytwin" / "private_source_terms.json"
    terms = json.loads(config.read_text())
    if not isinstance(terms, list) or not terms or any(not isinstance(t, str) or not t.strip() for t in terms):
        raise ValueError("A nonempty private-source exclusion list is required")
    return "(?:" + public_pattern + ")|(?:" + "|".join(re.escape(t) for t in terms) + ")"

sys.path.insert(0, '../3fold-motion-engine/_private/romi_collab/lanes')
from field_queue_store import read_queue, append_jobs, take_jobs, parent_ready, enabled, research_bridge
sys.path.insert(0, '../3fold-motion-engine/_private/romi_collab/lanes/automation_runtime')
R = '../3fold-motion-engine/_private/romi_collab'
B = R + '/build'; Q = R + '/lanes/SWARM_QUEUE_20260924.json'; LOG = R + '/lanes/field_ovh.log'
RUN = os.environ.get('FIELD_RUNNING_DIR',R + '/lanes/ovh_running'); os.makedirs(RUN, exist_ok=True)
H = os.environ.get('FIELD_CLOUD_HOST','ubuntu@51.77.110.4')
LOG = os.environ.get('FIELD_CLOUD_LOG',LOG)
CAPFILE = os.environ.get('FIELD_CAP_FILE',R+'/lanes/ovh_field_cap')
HOST_CAP=int(os.environ.get('FIELD_HOST_CAP','17'))
HEADROOM=int(os.environ.get('FIELD_HEADROOM_MIB','4096'))
SLICE_HEADROOM=int(os.environ.get('FIELD_SLICE_HEADROOM_MIB','2048'))
GUARD=f'sudo -u ubuntu python3 /opt/agents/launch_reserved.py --agent-mib 1500 --slice-headroom-mib {SLICE_HEADROOM} --host-headroom-mib {HEADROOM} --host-cap {HOST_CAP}'
SSHO = ['-o', 'ConnectTimeout=15', '-o', 'ServerAliveInterval=30', '-i', 'local_config_path/hunt_20260923', '-o', 'IdentitiesOnly=yes',
        '-o', 'UserKnownHostsFile=external_research_path', '-o', 'BatchMode=yes']
RSH = 'ssh ' + ' '.join(SSHO)
MK_PATH = '../3fold-motion-engine/_private/romi_collab/lanes/automation_runtime/mk.py'
def current_rules(last=['']):
    # Field's RULES change (1/10: TARGET_FORM, EXTERNAL_REFERENT, stub gate); read them per job, keep the last good copy.
    try:
        r = re.search(r"RULES='''(.*?)'''", open(MK_PATH).read(), re.S).group(1)
        last[0] = re.sub(r'Python: \S+', 'Python: python3 (numpy/scipy/sympy are available)', r)
    except Exception as e:
        print('rules reread failed; using last good copy', repr(e)[:120], flush=True)
    return last[0]
RULES = current_rules()
CLOUD_PLANNER_CAP = int(os.environ.get('FIELD_CLOUD_PLANNERS', '2'))
# 1/10 (anton-5f, the graph's finding via Field): the short path form bodytwin/tasks|results|notes and
# bt_memory did not match 3fold-workspaces/bodytwin. This driver sends ONLY Field jobs (F_<id>),
# so BodyTwin's own swarm jobs are not affected by the short forms being denied here.
DENY = re.compile(_private_input_pattern('3fold-workspaces/bodytwin|bodytwin/(tasks|results|notes|inputs)|bt_memory|mechanism_legacy|the collaborator(?!son)|shared_data|sdc1|L1/prep|/L1\\b|\\bL1-(operat|ruta|data)|\\bGC\\b|grand challenge|eknee|etibia|\\bJW\\d?|\\bDM\\d|\\bSC\\d|\\bPS\\d|EMG|CX-|implant'), re.I)
DENY_LANES = {'U414', 'U444', 'U447', 'U456', 'U460', 'U462', 'U470', 'U472', 'U475', 'U477', 'U478', 'U479', 'U480', 'U552', 'U554', 'U556', 'U568'}

def log(*a):
    open(LOG, 'a').write(time.strftime('%H:%M:%S ') + ' '.join(map(str, a)) + '\n')

def ssh(cmd):
    return subprocess.run(['ssh', '-n'] + SSHO + [H, cmd], capture_output=True, text=True, timeout=60)

def cloud_ok(j):
    if DENY.search(j.get('body', '') + j.get('title', '')): return False
    for p in (j.get('src') or {}).values():
        if os.path.basename(p.rstrip('/')) in DENY_LANES or DENY.search(p): return False
        if os.path.isdir(p):
            # 1/10 (anton-5f, Field's finding): the gate searched only .py/.json/.sh, so .md result text with
            # internal identifiers came along to the cloud. Now ALL text files are searched (-I skips binaries)
            # with the whole DENY pattern, not the short list.
            try:
                r = subprocess.run(['grep', '-rlIqE', DENY.pattern, p], timeout=120)
                if r.returncode == 0: return False
            except subprocess.TimeoutExpired:
                return False  # if there is no time to review the source, it is not sent
    return True

def prepare(j):
    jd = f"{B}/{j['id']}"
    storage=Path('external_mount')
    if storage.is_dir() and not os.path.lexists(jd):
        dest=storage/j['id'];dest.mkdir(exist_ok=True);Path(jd).symlink_to(dest,target_is_directory=True)
    os.makedirs(jd + '/src', exist_ok=True)
    for name, srcdir in (j.get('src') or {}).items():
        if os.path.isdir(srcdir):
            os.makedirs(f'{jd}/src/{name}', exist_ok=True)
            bridge=research_bridge()
            graph_source=bridge is not None and bridge.is_packet_source(srcdir)
            # 1/10 (anton-5f, at Field's request): a parent job used as source was copied with its own src/ (the grandparent) and binaries, 90-135 MB/job.
            parent_job=any(m in os.path.realpath(srcdir) for m in ('romi_collab/build','external_mount','external_mount'))
            subprocess.run(['rsync', '-a', '--exclude', 'opencode*.log', '--exclude', '__pycache__', '--exclude', '*.tar.zst']+
                           (['--exclude','/src','--exclude','/sims','--exclude','src.tar.zst','--max-size=20m'] if parent_job else [])+
                           ([] if graph_source else ['--exclude','work'])+
                           [srcdir.rstrip('/') + '/', f'{jd}/src/{name}/'], check=True, timeout=600)
    guard='Prior reports and planner claims are unaudited inputs. Verify their exact assumptions. Finite sampled successes or heuristic flags are not universal certificates. A statistical bound is not a deterministic force/error bound. Check units and preserve all negative controls. Use the strongest matched baseline available. If a dependency is supplied, read its actual outcome before using its claimed success.'
    principles=Path('external_research_path').read_text()
    open(jd + '/BRIEF.md', 'w').write(f"# {j['id']} — {j['title']}\n{j['body']}\n{guard}\n{current_rules()}\n\n{principles}\n")
    open(jd+'/JOB.json','w').write(json.dumps(j,ensure_ascii=False,indent=1))
    if j.get('preferred_model') in ('swarm','free_worker'):
        open(jd+'/PREFERRED_MODEL','w').write(j['preferred_model'])
    shutil.copy2(R+'/lanes/automation_runtime/execute_job.py',jd+'/_run_bunny.py')
    open(jd + '/.ovh_claim', 'w').write(time.strftime('%F %T'))
    return jd

def fetch(active):
    for f in os.listdir(RUN):
        jid = f
        if any(u.startswith(f'agent-F_{jid}-') for u in active): continue
        jd = f'{B}/{jid}'
        # Workers can create computational code under src/. Skipping that tree
        # discarded the only copy of new code before remote cleanup. Rsync is
        # incremental against the existing supplied source snapshot.
        copied=subprocess.run([sys.executable, 'external_research_path', f'{H}:/opt/agents/jobs/F_{jid}/', jd + '/', RSH], stdin=subprocess.DEVNULL)
        if copied.returncode:
            log('fetch copy failed; remote preserved',jid,copied.returncode);continue
        ok = os.path.exists(jd + '/RESULTS.md')
        log('fetched', jid, 'results' if ok else 'NO-RESULTS', open(jd + '/AGENT_EXIT').read().strip() if os.path.exists(jd + '/AGENT_EXIT') else 'no-exit')
        ssh(f'sudo -u ubuntu rm -rf /opt/agents/jobs/F_{jid}'); os.remove(f'{RUN}/{f}')
        try: os.remove(jd + '/.ovh_claim')
        except OSError: pass

def start(j, acc):
    jd = prepare(j); jid = j['id']
    try:
        r = subprocess.run(['rsync', '-a', '--rsync-path=sudo -u ubuntu rsync', '-e', RSH, '--exclude', '.ovh_claim', jd + '/', f'{H}:/opt/agents/jobs/F_{jid}/'], stdin=subprocess.DEVNULL, timeout=900)
    except subprocess.TimeoutExpired:
        log('rsync timeout 900s', jid); return False  # 1/10: an 807 MB source locked both single-threaded drivers for 47 min
    if r.returncode: log('rsync fail', jid); return False
    prompt = f'Read BRIEF.md and execute the bounded task. Write RESULTS.md with the job ID on line 1 when done.'
    cmd = (f"{GUARD.replace('sudo -u ubuntu', 'sudo')} -- systemd-run --quiet --collect --slice=research.slice --unit=agent-F_{jid}-$(date +%s) --uid=ubuntu -p MemoryMax=1500M -p MemoryHigh=1250M "
           f"-p CPUQuota=100% -p RuntimeMaxSec=5400 --working-directory=/opt/agents/jobs/F_{jid} /bin/bash -c "
           f"'python3 _run_bunny.py {acc} F_{jid} {jid} >> wrapper.log 2>&1 < /dev/null; echo $? > AGENT_EXIT'")
    r = ssh(cmd)
    if r.returncode: log('start fail', jid, r.stderr[-200:]); return False
    open(f'{RUN}/{jid}', 'w').write(acc); log('ovhstart', jid, acc); return True

n_start = 0
while True:
    try:
        if os.environ.get('FIELD_REQUIRE_READY') and ssh('test -s /opt/agents/BOOTSTRAP_OK').returncode:
            log('waiting for completed worker bootstrap');time.sleep(30);continue
        cap = int(open(CAPFILE).read().strip()) if os.path.exists(CAPFILE) else 12
        if cap > 0 and Path('external_research_path').is_file():
            cap = 10000  # Every start still passes the shared resource guard.
        r = ssh("systemctl list-units --no-legend --state=active,activating,deactivating 'agent-*' | awk '{print $1}'")
        if r.returncode == 0:
            fetch(r.stdout.split())
            free = cap - len(os.listdir(RUN))
            if free > 0 and enabled():
                disk = ssh("python3 -c \"import shutil; print(shutil.disk_usage('/opt/agents').free)\"").stdout.strip()
                if not disk.isdigit() or int(disk)<2*1024**3:
                    log('dispatch waiting for remote disk headroom',disk)
                    time.sleep(30);continue
                capacity = ssh(GUARD)
                free = min(free, int(capacity.stdout.strip())) if capacity.returncode == 0 else 0
                q = read_queue(); picked = []
                # Local planners: 189 of the last 200 finished Field dirs were rate-interrupted planners (1/10).
                cloud_planners = sum(1 for f in os.listdir(RUN) if f.startswith('FB_BPLAN'))
                for j in q:
                    if len(picked) >= free: break
                    if j['id'].startswith('FB_BPLAN') and cloud_planners >= CLOUD_PLANNER_CAP: continue
                    if os.path.exists(f"{B}/{j['id']}"): continue
                    if not parent_ready(j): continue
                    if cloud_ok(j):
                        picked.append(j)
                        if j['id'].startswith('FB_BPLAN'): cloud_planners += 1
                if picked:
                    ids = {j['id'] for j in picked}
                    for j in picked:
                        selected=take_jobs([j['id']])
                        if not selected:continue
                        j=selected[0]
                        n_start += 1
                        try: ok = start(j, 'ABC'[n_start % 3])
                        except Exception as e: ok = False; log('start error', j['id'], repr(e))
                        if not ok:
                            subprocess.run(['rm', '-rf', f"{B}/{j['id']}"]); append_jobs([j]); log('requeued', j['id'])
                        time.sleep(3)
        else:
            log('ssh fail', r.stderr[-120:])
    except Exception as e:
        log('error', repr(e))
    time.sleep(15)
