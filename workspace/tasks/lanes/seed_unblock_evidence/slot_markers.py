"""Slot identity and cleanup. All scans/mutations hold model_slots.lock."""
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
        if row['selected'] not in ('swarm', 'swarm_worker', 'swarm_worker') or row['account'] not in ('A','B','C','D'): return 'malformed'
        try: proc, ticks, boot, epoch = identity(pid)
        except (FileNotFoundError, ProcessLookupError): return 'dead'
        if not 0 <= now-started <= MAX_AGE: return 'aged'
        if epoch > started+2: return 'reused'
        if 'start_ticks' in row and (row['start_ticks'] != ticks or row.get('boot_id') != boot): return 'reused'
        args = (proc/'cmdline').read_bytes().split(b'\0')
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
    tmp = marker.with_suffix('.tmp'); tmp.write_text(json.dumps(row)+'\n'); tmp.replace(marker)
    OWN = marker, row
    return marker

atexit.register(release)
