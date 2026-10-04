"""Cold one-command reproduction into an isolated directory on the data disk."""
from dental_release.paths import expand as _release_expand
import os, shutil, subprocess, datetime
from pathlib import Path
H = Path(__file__).resolve().parents[1]
base = Path(_release_expand('@DENTAL_WORK_ROOT@/X2_occlusion_b2b'))
used = sum((p.stat().st_size for p in base.rglob('*') if p.is_file() and (not p.is_symlink())))
if used + 800000000 > 3000000000:
    raise RuntimeError('Not enough lane budget for another isolated cold replay; keep existing evidence')
out = base / 'replays' / datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%SZ')
(out / 'code').mkdir(parents=True)
(out / 'raw').mkdir()
for p in (H / 'code').glob('*.py'):
    shutil.copy2(p, out / 'code' / p.name)
for p in H.glob('PREREG_*.json'):
    shutil.copy2(p, out / p.name)
for name in ['EXTERNAL_REFERENTS.json', 'SUPPORT_INTERVALS.json', 'WITHIN_PERSON_SPREAD.json']:
    shutil.copy2(H / 'raw' / name, out / 'raw' / name)
shutil.copy2(H / 'run_all.sh', out / 'run_all.sh')
env = dict(os.environ, X2_DATA_DIR=str(out / 'arrays'))
print('Cold replay directory:', out, flush=True)
subprocess.run(['bash', str(out / 'run_all.sh')], env=env, check=True)
