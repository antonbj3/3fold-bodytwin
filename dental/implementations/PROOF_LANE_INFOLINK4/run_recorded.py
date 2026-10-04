from pathlib import Path
import subprocess, sys, datetime, json, hashlib, time
H = Path(__file__).resolve().parent
P = H / 'runs'
P.mkdir(exist_ok=True)
p = P / f"run_{len(list(P.glob('run_*'))) + 1:03}"
p.mkdir()
records = []
for script in ['make_demo.py', 'interaction_demo.py', 'report.py']:
    if not (H / script).exists():
        continue
    raw = (H / script).read_bytes()
    (p / script).write_bytes(raw)
    start = time.perf_counter()
    r = subprocess.run([sys.executable, script], cwd=H, capture_output=True, text=True)
    (p / (script + '.stdout')).write_text(r.stdout)
    (p / (script + '.stderr')).write_text(r.stderr)
    records.append(dict(command=[sys.executable, script], cwd=str(H), code_sha256=hashlib.sha256(raw).hexdigest(), elapsed_seconds=time.perf_counter() - start, returncode=r.returncode, at_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
    (p / 'commands.json').write_text(json.dumps(records, indent=2) + '\n')
    print(r.stdout, end='')
    print(r.stderr, end='', file=sys.stderr)
    if r.returncode:
        sys.exit(r.returncode)
print('Run evidence:', p)
