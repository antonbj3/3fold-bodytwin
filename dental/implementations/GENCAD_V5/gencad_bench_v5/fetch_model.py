"""Fetch author code/weights only. No patient data. Pin API tree and checkpoint hashes."""
from dental_release.paths import expand as _release_expand
import requests, json, hashlib, time
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
DATA = Path(_release_expand('@DENTAL_WORK_ROOT@/PROOF_LANE_GENCAD_V5'))
DATA.mkdir(parents=True, exist_ok=True)
start = time.time()
tree = json.loads((ROOT / 'raw/ikarus1211_VISAPP_ToothCraft_tree.json').read_text())
commit = tree['sha']
dest = DATA / 'ToothCraft'
dest.mkdir(exist_ok=True)
files = []
for x in tree['tree']:
    p = x['path']
    if x['type'] == 'blob' and (p.startswith(('model/', 'configs/', 'AugmentPipeline/utils/')) or p in ['README.MD', 'requirements.txt', 'test.py', 'AugmentPipeline/generator.py']):
        q = dest / p
        q.parent.mkdir(parents=True, exist_ok=True)
        r = requests.get(f'https://raw.githubusercontent.com/ikarus1211/VISAPP_ToothCraft/{commit}/{p}', timeout=60)
        r.raise_for_status()
        q.write_bytes(r.content)
        files.append({'path': str(q), 'sha256': hashlib.sha256(r.content).hexdigest(), 'bytes': len(r.content)})
for x in json.loads((ROOT / 'raw/ToothCraft_HF_tree.json').read_text()):
    if x['path'] not in ['ToothCraft_normal_control_branch.pth', 'ToothCraft_normal_diff_branch.pth']:
        continue
    q = DATA / 'weights' / x['path']
    q.parent.mkdir(exist_ok=True)
    if not q.exists() or q.stat().st_size != x['size']:
        with requests.get('https://huggingface.co/DejvaX/ToothCraft/resolve/main/' + x['path'], stream=True, timeout=60) as r:
            r.raise_for_status()
            with q.open('wb') as f:
                for b in r.iter_content(1024 * 1024):
                    f.write(b)
    digest = hashlib.sha256(q.read_bytes()).hexdigest()
    assert digest == x['lfs']['oid']
    files.append({'path': str(q), 'sha256': digest, 'bytes': q.stat().st_size})
    print(q.name, q.stat().st_size, flush=True)
(ROOT / 'raw/EXTERNAL_MODEL_LOCK.json').write_text(json.dumps({'repository': 'https://github.com/ikarus1211/VISAPP_ToothCraft', 'commit': commit, 'weights': 'https://huggingface.co/DejvaX/ToothCraft', 'license': 'UNSPECIFIED: no licence file or model card in observed trees; private evaluation only, redistribution not granted', 'files': files, 'seconds': time.time() - start}, indent=2) + '\n')
