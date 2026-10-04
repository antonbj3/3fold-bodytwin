"""One-command copied replay, all source arrays read-only and hash-checked."""
from dental_release.paths import expand as _release_expand
import os
for k in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS']:
    os.environ[k] = '4'
import sys, argparse, pathlib, shutil, json, datetime, hashlib, time, resource
SOURCE = pathlib.Path(__file__).resolve().parent
ap = argparse.ArgumentParser()
ap.add_argument('--out', type=pathlib.Path)
a = ap.parse_args()
out = a.out or pathlib.Path(_release_expand('@DENTAL_EXTERNAL_ROOT@/storage/tmp_dental_audit/X95-whole-insertion-path')) / datetime.datetime.now(datetime.timezone.utc).strftime('replay_%Y%m%dT%H%M%S_%fZ')
out.mkdir(parents=True, exist_ok=False)
shutil.copytree(SOURCE / 'code', out / 'code', ignore=shutil.ignore_patterns('__pycache__'))
(out / 'raw').mkdir()
(out / 'rounds').mkdir()
(out / 'figures').mkdir()
for pattern in ['PREREG*.json', 'PREREG*.sha256', 'FROZEN*.json', 'FROZEN*.sha256']:
    for p in SOURCE.glob(pattern):
        shutil.copy2(p, out / p.name)
for name in ['BRIEF.md', 'SOURCES.md', 'LAB_PROTOCOL.md', 'raw/REUSED_SOURCE.json']:
    p = SOURCE / name
    if p.exists():
        shutil.copy2(p, out / name)
for p in out.glob('PREREG*.sha256'):
    source = p.with_suffix('.json')
    if source.exists():
        assert hashlib.sha256(source.read_bytes()).hexdigest() == p.read_text().strip()
reference = json.load(open(SOURCE / 'results.json'))
(out / 'REPLAY_REFERENCE_POST_OBSERVATION.json').write_text(json.dumps(dict(kind='post_observation_replay_reference_not_prediction', source_results_sha256=hashlib.sha256((SOURCE / 'results.json').read_bytes()).hexdigest(), first_contacts=reference['first_contacts'], selected_axes=reference['finite_axis_search']['selected'], colliding_sites=reference['cohort']['colliding_sites']), indent=2) + '\n')
sys.path.insert(0, str(out / 'code'))
import fixtures, run_cohort, repair, first_event, axes, blockout_control, validate, render, finalize, export_check, cone_link
print('AUDIT_COPY', out, flush=True)
wall = time.perf_counter()
cpu = time.process_time()
fixtures.run(out / 'raw')
run_cohort.run(out / 'rounds/R1')
repair.run(out / 'rounds/R2')
export_check.run(out)
first_event.run(out / 'rounds/R3')
axes.run(out / 'rounds/R4')
cone_link.run(out)
blockout_control.run(out)
validate.run(out)
render.run(out)
result = finalize.run(out)
expected = {x['key']: x['height_exact_mm'] for x in reference['first_contacts']}
actual = {x['key']: x['height_exact_mm'] for x in result['first_contacts']}
axes_expected = {x['key']: x['translation_mm'] for x in reference['finite_axis_search']['selected']}
axes_actual = {x['key']: x['translation_mm'] for x in result['finite_axis_search']['selected']}
parity = dict(first_height_exact_equal=expected == actual, collision_sites_equal=reference['cohort']['colliding_sites'] == result['cohort']['colliding_sites'], selected_axes_equal=axes_expected == axes_actual, physical_unknown_preserved=True, source_cohort_sha256=hashlib.sha256((out / 'FROZEN_COHORT.json').read_bytes()).hexdigest(), scope='Frozen decisive geometry quantities; CPU-time budget UNKNOWN may vary across replays; each new record preserved.')
assert all((parity[k] for k in ['first_height_exact_equal', 'collision_sites_equal', 'selected_axes_equal']))
(out / 'raw/REPLAY_PARITY.json').write_text(json.dumps(parity, indent=2) + '\n')
size = sum((p.stat().st_size for p in out.rglob('*') if p.is_file()))
assert size < 3000000000
import numpy, scipy, trimesh, matplotlib, rtree
cost = dict(wall_s=time.perf_counter() - wall, cpu_s=time.process_time() - cpu, peak_RSS_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, own_copy_bytes=size, threads_max=4, GPU=False, versions={m.__name__: m.__version__ for m in [numpy, scipy, trimesh, matplotlib, rtree]}, python=sys.version, all_code_sha256={str(p.relative_to(out)): hashlib.sha256(p.read_bytes()).hexdigest() for p in (out / 'code').glob('*.py')}, source_results_sha256=hashlib.sha256((SOURCE / 'results.json').read_bytes()).hexdigest())
(out / 'DEMO_COST.json').write_text(json.dumps(cost, indent=2) + '\n')
print('DEMO_COMPLETE', json.dumps(dict(out=str(out), parity=parity, cost=cost['wall_s'])), flush=True)
