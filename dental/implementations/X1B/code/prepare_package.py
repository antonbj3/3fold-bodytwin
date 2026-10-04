"""One-time, bounded import of selected X1 artifacts, with unchanged origin hashes.
This preparation script records source locations; run_all does not execute it.
"""
from dental_release.paths import expand as _release_expand
from pathlib import Path
import shutil, json, hashlib, sys, platform, importlib, datetime
import numpy as np
R = Path(__file__).resolve().parents[1]
W = R.parents[1]
X = W / 'results/LANE_X1_CROWN_LOOP'
DATA = Path(_release_expand('@DENTAL_WORK_ROOT@/LANE_X1_CROWN_LOOP'))
origin = []

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def copy(src, dst):
    dst = R / dst
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)
    origin.append(dict(source=str(src), local=str(dst.relative_to(R)), sha256=sha(dst), bytes=dst.stat().st_size))

def write(p, d):
    p = R / p
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(d, indent=2) + '\n')
for rel in ['results.json', 'FROZEN_PREDICTIONS.json', 'FROZEN_PREDICTIONS.sha256', 'EXPORT_MANIFEST.json', 'README_DEMO.md', 'HANDOFF.md', 'LAB_PROTOCOL.md', 'VERIFICATION.json', 'DECOMPOSITION.json']:
    copy(X / rel, 'history/X1/' + rel)
for p in X.glob('PREREG*'):
    copy(p, 'history/X1/' + p.name)
for p in (X / 'code').glob('*.py'):
    copy(p, 'history/X1/code/' + p.name)
for p in (X / 'raw').glob('FAIL*'):
    copy(p, 'history/X1/failures/' + p.name)
for rel in ['XREVIEW.md', 'REVIEW_LANE_X1_CROWN_LOOP.json']:
    copy(W / 'results/PROOF_LANE_XREVIEW_RESTORATIVE' / rel, 'history/review/' + rel)
modules = ['design/crown_case_sts3d.py', 'design/crown_design_geometry.py', 'manufacturing/crown_fit_geometry.py', 'design/crown_design_weibull.py']
for rel in modules:
    copy(W / 'cells' / rel, 'vendor/' + Path(rel).name)
src = W / 'tasks/swarm48/sources/LIT_CROWN/crown_fracture_papers.jsonl'
copy(src, 'inputs/literature/crown_fracture_papers.jsonl')
for id in ['PMC4764450', 'PMC8558575', 'PMC6642729', 'PMC9081244', 'PMC7274823']:
    copy(W / 'results/DESIGN_crown/sources' / f'{id}.txt', f'inputs/literature/{id}.txt')
for id in ['PMC10817558', 'PMC10934854']:
    copy(W / 'data/corpus/europepmc/fulltext' / f'{id}.xml', f'inputs/literature/{id}.xml')
frozen = json.loads((X / 'FROZEN_PREDICTIONS.json').read_text())
for d in frozen['designs']:
    tid = '_'.join(d['source_id'].split('_')[:-2])
    vid = d['geometry_parameters']['id']
    des = d['design']
    copy(W / 'results/CROWN_pop/cases' / f'{tid}.npz', f'inputs/anatomy/{tid}.npz')
    wd = DATA / tid / vid
    copy(wd / 'model.npz', f'inputs/geometry/{des}_model.npz')
    a = np.load(wd / 'case_h0.12.npz')
    write(f'inputs/geometry/{des}_grid.json', dict(origin=a['origin'].tolist(), h=float(a['h']), shape=a['shape'].tolist(), z_m=float(a['z_m']), frame_R=a['R'].tolist(), frame_c=a['c'].tolist(), source=str(wd / 'case_h0.12.npz'), source_sha256=sha(wd / 'case_h0.12.npz')))
    for p in (X / 'exports' / des).iterdir():
        if p.suffix in ('.stl', '.3mf'):
            copy(p, f'exports/{des}/{p.name}')
res = json.loads((X / 'results.json').read_text())
records = []
for row in res['rows']:
    p = Path(row['source'])
    fe = json.loads(p.read_text())
    base = f"inputs/fe/{row['id']}"
    copy(p, base + '.origin.json')
    copy(fe['stress_npz'], base + '.npz')
    deck = Path(fe['stress_npz']).with_name('solve.inp')
    copy(deck, base + '.inp')
    records.append(dict(id=row['id'], material=row['material'], material_model=fe['material_model'], stress=base + '.npz', deck=base + '.inp', predictions=fe['predictions'], cases=fe['cases'], source_sha256=sha(p)))
write('inputs/fe/INDEX.json', records)
write('inputs/X1_ROWS.json', res['rows'])
write('VENDOR_MANIFEST.json', [x for x in origin if x['local'].startswith('vendor/')])
write('SOURCE_MANIFEST.json', dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), origin=origin, total_copied_bytes=sum((x['bytes'] for x in origin)), no_tree_copy=True))
runtime = dict(python_version=platform.python_version(), interpreter_locator=sys.executable, interpreter_sha256=sha(Path(sys.executable).resolve()), packages={n: importlib.import_module(n).__version__ for n in ['numpy', 'scipy', 'trimesh', 'matplotlib']})
write('RUNTIME_LOCK.json', runtime)
(R / 'requirements.lock').write_text('\n'.join((f'{k}=={v}' for (k, v) in runtime['packages'].items())) + '\n')
print('Copied MB', sum((x['bytes'] for x in origin)) / 1000000.0, 'runtime', runtime)
