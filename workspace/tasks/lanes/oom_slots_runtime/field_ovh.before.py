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
HOST_CAP=int(os.environ.get('FIELD_HOST_CAP','44'))
HEADROOM=int(os.environ.get('FIELD_HEADROOM_MIB','4096'))
SSHO = ['-o', 'ConnectTimeout=15', '-o', 'ServerAliveInterval=30', '-i', 'local_config_path/hunt_20260923', '-o', 'IdentitiesOnly=yes',
        '-o', 'UserKnownHostsFile=external_research_path', '-o', 'BatchMode=yes']
RSH = 'ssh ' + ' '.join(SSHO)
MK = open('../3fold-motion-engine/_private/romi_collab/lanes/automation_runtime/mk.py').read()
RULES = re.search(r"RULES='''(.*?)'''", MK, re.S).group(1)
RULES = re.sub(r'Python: \S+', 'Python: python3 (numpy/scipy/sympy are available)', RULES)
DENY = re.compile(_private_input_pattern('3fold-workspaces/bodytwin|collaborator|shared_data|sdc1|L1/prep|/L1\\b|\\bL1-(operat|ruta|data)|\\bGC\\b|grand challenge|eknee|etibia|\\bJW\\d?|\\bDM\\d|\\bSC\\d|\\bPS\\d|EMG|CX-|implant'), re.I)
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
            r = subprocess.run(['grep', '-rlIqE', '--include=*.py', '--include=*.json', '--include=*.sh',
                                '3fold-workspaces/bodytwin|collaborator_WR|L1/prep|sdc1-tmp|shared_data/bodytwin', p])
            if r.returncode == 0: return False
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
            subprocess.run(['rsync', '-a', '--exclude', 'opencode*.log', '--exclude', '__pycache__', '--exclude', '*.tar.zst']+
                           ([] if graph_source else ['--exclude','work'])+
                           [srcdir.rstrip('/') + '/', f'{jd}/src/{name}/'], check=True)
    guard='Prior reports and planner claims are unaudited inputs. Verify their exact assumptions. Finite sampled successes or heuristic flags are not universal certificates. A statistical bound is not a deterministic force/error bound. Check units and preserve all negative controls. Use the strongest matched baseline available. If a dependency is supplied, read its actual outcome before using its claimed success.'
    principles=Path('external_research_path').read_text()
    open(jd + '/BRIEF.md', 'w').write(f"# {j['id']} — {j['title']}\n{j['body']}\n{guard}\n{RULES}\n\n{principles}\n")
    open(jd+'/JOB.json','w').write(json.dumps(j,ensure_ascii=False,indent=1))
    if j.get('preferred_model') in ('swarm','swarm_worker'):
        open(jd+'/PREFERRED_MODEL','w').write(j['preferred_model'])
    shutil.copy2(R+'/lanes/automation_runtime/execute_job.py',jd+'/_run_swarm.py')
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
        ssh(f'rm -rf /opt/agents/jobs/F_{jid}'); os.remove(f'{RUN}/{f}')
        try: os.remove(jd + '/.ovh_claim')
        except OSError: pass

def start(j, acc):
    jd = prepare(j); jid = j['id']
    r = subprocess.run(['rsync', '-a', '-e', RSH, '--exclude', '.ovh_claim', jd + '/', f'{H}:/opt/agents/jobs/F_{jid}/'], stdin=subprocess.DEVNULL)
    if r.returncode: log('rsync fail', jid); return False
    prompt = f'Read BRIEF.md and execute the bounded task. Write RESULTS.md with the job ID on line 1 when done.'
    cmd = (f"sudo systemd-run --quiet --collect --slice=research.slice --unit=agent-F_{jid}-$(date +%s) --uid=ubuntu -p MemoryMax=1500M -p MemoryHigh=1250M "
           f"-p CPUQuota=100% -p RuntimeMaxSec=5400 --working-directory=/opt/agents/jobs/F_{jid} /bin/bash -c "
           f"'python3 _run_swarm.py {acc} F_{jid} {jid} >> wrapper.log 2>&1 < /dev/null; echo $? > AGENT_EXIT'")
    r = ssh(cmd)
    if r.returncode: log('start fail', jid, r.stderr[-200:]); return False
    open(f'{RUN}/{jid}', 'w').write(acc); log('ovhstart', jid, acc); return True

n_start = 0
while True:
    try:
        if os.environ.get('FIELD_REQUIRE_READY') and ssh('test -s /opt/agents/BOOTSTRAP_OK').returncode:
            log('waiting for completed worker bootstrap');time.sleep(30);continue
        cap = int(open(CAPFILE).read().strip()) if os.path.exists(CAPFILE) else 12
        r = ssh("systemctl list-units --no-legend --state=running 'agent-*' | awk '{print $1}'")
        if r.returncode == 0:
            fetch(r.stdout.split())
            free = cap - len(os.listdir(RUN))
            if free > 0 and enabled():
                disk = ssh("python3 -c \"import shutil; print(shutil.disk_usage('/opt/agents').free)\"").stdout.strip()
                if not disk.isdigit() or int(disk)<2*1024**3:
                    log('dispatch waiting for remote disk headroom',disk)
                    time.sleep(30);continue
                mem = ssh("awk '/MemAvailable/{print int($2/1024)}' /proc/meminfo").stdout.strip()
                # Shared OVH host: BodyTwin and Field together must fit.
                free = min(free, max(0, HOST_CAP-len(r.stdout.split())), max(0,(int(mem or 0)-HEADROOM)//1100))
                q = read_queue(); picked = []
                for j in q:
                    if len(picked) >= free: break
                    if j['id'].startswith('FB_BPLAN'): continue
                    if os.path.exists(f"{B}/{j['id']}"): continue
                    if not parent_ready(j): continue
                    if cloud_ok(j): picked.append(j)
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
