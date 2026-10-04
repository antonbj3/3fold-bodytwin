"""Save prior small reports before an authorized one-command replay."""
import datetime
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
ROOT = Path(__file__).resolve().parent
used_kb = int(subprocess.check_output(['du', '-sk', str(ROOT)], text=True).split()[0])
if used_kb * 1024 > 2900000000:
    raise RuntimeError('lane approaching frozen 3 GB cap; stop before writing replay intermediates')
names = ['results.json', 'RESULTS_R1.json', 'RESULTS_R2.json', 'RESULTS_R3.json', 'RESULTS_R4.json', 'LEGACY_REPRODUCTION.json', 'BODYTWIN_HAND_CASE_COMPARISON.json', 'VERIFICATION.json']
existing = [ROOT / x for x in names if (ROOT / x).exists()]
if existing:
    stamp = datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%S.%fZ')
    dest = ROOT / 'versions/replays' / stamp
    dest.mkdir(parents=True)
    manifest = []
    for file in existing:
        shutil.copyfile(file, dest / file.name)
        manifest.append({'name': file.name, 'bytes': file.stat().st_size, 'sha256': hashlib.sha256(file.read_bytes()).hexdigest()})
    (dest / 'MANIFEST.json').write_text(json.dumps(manifest, indent=2) + '\n')
