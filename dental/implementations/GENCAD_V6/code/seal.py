from dental_release.paths import expand as _release_expand
from common import *
import subprocess

def runtime(exe, names):
    code = 'import sys,json,importlib.metadata as m; print(json.dumps(dict(executable=sys.executable,python=sys.version,packages={n:m.version(n) for n in ' + repr(names) + '})))'
    return json.loads(subprocess.check_output([exe, '-B', '-c', code], text=True))
dump(ROOT / 'RUNTIME_LOCK.json', dict(geometry=runtime(_release_expand('@DENTAL_PYTHON@'), ['numpy', 'scipy', 'trimesh', 'matplotlib', 'scikit-image']), load=runtime('/usr/bin/python3', ['numpy', 'scipy', 'numba', 'llvmlite']), clean_install_verified=False))
critical = ['results.json', 'RESULTS.md', 'README_DEMO.md', 'HANDOFF.md', 'RUNTIME_LOCK.json', 'LEADERBOARD.csv', 'raw/TOOTH_LOAD_SETS.csv', 'raw/R1_ROWS.json', 'raw/R2_ROWS.json', 'raw/R3_ROWS.json', 'raw/EXTERNAL_PAIRED_OBSERVATIONS.json', 'raw/CONTACT_CONTROLS.json', 'raw/LOAD_CONTROLS.json', 'raw/OUTPUT_CONTROLS.json', 'raw/RUN_COST.json', 'FROZEN_MEASUREMENT_PREDICTIONS.json']
critical.extend((str(p.relative_to(ROOT)) for p in sorted((ROOT / 'code').glob('*.py'))))
critical.extend((str(p.relative_to(ROOT)) for p in sorted(ROOT.glob('PREREG_*.json'))))
critical.extend((str(p.relative_to(ROOT)) for p in sorted((ROOT / 'rounds').glob('*.json'))))
critical.extend(['run_all.sh', 'EXTERNAL_REFERENTS.md', 'LAB_PROTOCOL.md', 'figures/demo.png', 'figures/demo.pdf'])
critical.extend((str(p.relative_to(ROOT)) for p in sorted(ROOT.glob('DECOMPOSITION_*.json'))))
critical.extend((str(p.relative_to(ROOT)) for p in sorted((ROOT / 'raw').glob('SUFFICIENCY_*.json'))))
critical.extend((str(p.relative_to(ROOT)) for p in sorted((ROOT / 'exports').rglob('*')) if p.is_file()))
manifest = {x: dict(sha256=sha(ROOT / x), bytes=(ROOT / x).stat().st_size) for x in critical}
dump(ROOT / 'ARTIFACT_MANIFEST.json', manifest)
freeze(ROOT / 'REVIEW_SNAPSHOT.json', dict(result=read(ROOT / 'results.json'), artifact_manifest=manifest, review_state='PENDING_INDEPENDENT_REVIEW', source_locks={p.name: sha(p) for p in ROOT.glob('INPUT_LOCK_*.json')}, retrospective=True))
feedback = dict(target_id='DENT-VAL-BASELINE-COMPARISON', result_file='results/PROOF_LANE_GENCAD_V6/REVIEW_SNAPSHOT.json', sha256=sha(ROOT / 'REVIEW_SNAPSHOT.json'), review_state='PENDING_INDEPENDENT_REVIEW', status='PENDING_INDEPENDENT_REVIEW', claim_type='capability', outcome='SPATIAL_DISCRIMINATION_DELIVERED; POSITIVE_CONTACT_RESTORATION_AND_PHYSICAL_FORCE_UNVALIDATED', measured_quantity='Geometric contact polygon location/count/area; all-offset contact filtration; conditional axial tooth-load sets; published paired contact-count/relative-force observations', units='mm; mm2; mm3; per-tooth class; percentage points', uncertainty='PL gap scenarios, projected rounded mesh intervals; physical support/pose/pressure unknown; broad force hulls; unmatched population referents; floating/source enclosures missing for contact areas', population_regime='1728 roof tasks,48case clusters;243 crown requests234scored on18sites6case clusters; one published repeated-measurement subject; no matched fabricated specimen', preregistered_gate=[f'PREREG_R{i}.json' for i in range(1, 6)], baseline='Actually executed v4 geometry integrator, v5 component count, TANDLAST KKT replay and quantitative wrong-output controls; no algorithm superiority claim', negative_result=True, dispatch_scope='Independent review request; stale working-view packet currently prevents dispatch', missing_coverage='No new broad scientific node requested. Existing benchmark target remains; coordinator must refresh working-view draft index before binding.', external_referent=REFERENT, edge_contracts=read(ROOT / 'results.json')['edges'], next_operation='Matched loaded geometry/pressure/displacement acquisition; or rigorous source-supported positive-contact reconstruction with unchanged failed R5 gate')
dump(ROOT / 'GRAPH_FEEDBACK.json', feedback)
print('SEALED', sha(ROOT / 'REVIEW_SNAPSHOT.json'))
