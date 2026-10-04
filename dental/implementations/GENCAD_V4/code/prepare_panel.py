"""Author-only preparation; never invoked by the released scorer."""
from dental_release.paths import expand as _release_expand
from common import *
import shutil
V3 = ROOT.parent / 'PROOF_LANE_GENCAD_V3'
FIELD = Path(_release_expand('@DENTAL_EXTERNAL_ROOT@/engines/3fold-motion-engine/_private/romi_collab/build/SOL_FALT_GENCADDENTAL_20261003/code'))
X42 = ROOT.parent / _release_expand('X42')

def copy(src, dst):
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)

def run():
    start = time.perf_counter()
    target = Path(_release_expand('@DENTAL_WORK_ROOT@/PROOF_LANE_GENCAD_V4/payload'))
    target.mkdir(parents=True, exist_ok=True)
    if not DATA.exists():
        DATA.symlink_to(target, target_is_directory=True)
    cohort = read(V3 / 'payload/private/COHORT.json')['payload']['cases']
    selected = []
    for (split, n) in [('test', 48), ('dev', 8), ('auxiliary', 8)]:
        rows = [r for r in cohort if r['split'] == split]
        selected += sorted(rows, key=lambda r: hashlib.sha256(('PROOF_LANE-v4-panel|' + r['group']).encode()).hexdigest())[:n]
    cfg = dict(round='R2', claim_type='capability', capability='Spatial quality comparison of all available frozen participants on measured registered IOS roofs, followed by whole-crown construction', obstacle='v3 PASS fractions ignore location and shape; fields do not match physical laboratory protocols', changed_operation='Exact PL contact-mask overlap, fieldwise shape error, explicit missing physical endpoints, clustered paired Pareto comparison', consumer='Researcher selecting designs for subsequent whole-crown laboratory tests', selection=selected, selection_rule='48 test, 8 dev, 8 auxiliary, hash PROOF_LANE-v4-panel|patient_group; before new quality scoring. Dev is in-sample for X42 and never pooled with test.', participants=['population', 'parametric', 'X1B_morphology', 'X18_antagonist', 'constraint_optimizer', 'field_generator', 'x42_shape', 'x42_contact_branch'], x1b_tool_radius_mm={'easy': 0.25, 'normal': 0.5, 'hard': 0.8, 'boundary': 1.0}, metrics={'anatomy_rmse_mm': {'resolution': 'PER_TOOTH', 'comparison_tolerance': 0.01}, 'contact_symdiff_mm2': {'resolution': 'PER_SURFACE_REGION', 'comparison_tolerance': 0.05}, 'negative_gap_area_mm2': {'resolution': 'PER_SURFACE_REGION', 'comparison_tolerance': 0.05}, 'nominal_gap_mean_mm': {'resolution': 'PER_TOOTH', 'comparison_tolerance': 1e-06}, 'nominal_gap_range_mm': {'resolution': 'PER_POINT', 'comparison_tolerance': 1e-06}}, contact_band_mm=0.1, pose_scenario_mm=0.05, bootstrap_replicates=2000, bootstrap_seed=20261003, decision='Capability requires feasible constructed discriminating pairs and at least one test-family paired metric contrast whose patient bootstrap interval excludes zero. No all-endpoint or clinical claim if full crown/film/fracture unresolved.', strongest_equally_informed_control='Actual v3 constraint_optimizer and unchanged field generator; actual X1b/X18 operations; X42 frozen models trained on dev only', falsifiers=['All quality vectors collapse for all L1-PASS designs', 'Hidden reference accessible to generator', 'Wrong private geometry or postfreeze prediction accepted', 'Physical calibration silently extrapolated to generated crown'], full_cost={'preparation': 'Selected v3 public scenes and private reference bytes only; metered copy', 'fit': 'Inherited X11, v3 training templates and X42 models; costs inherited/UNKNOWN, no fresh model fit', 'discovery': 'R1 plus code and frozen panel setup', 'validation': 'Regenerate every participant, score only after freeze, input faults and copy replay', 'questions': '64 patient-case groups x 9 families x 4 difficulty levels x 8 participants', 'fallback': 'Per-task ABSTAIN/UNKNOWN/FAIL retained; no excluded failure denominators'}, external_referent={'kind': 'published_dataset', 'locator': 'https://ditto.ing.unimore.it/bite2text/', 'compared_quantity': 'Original tooth surface heights and registered antagonist geometry, geometric contact pattern at declared 0.1 mm band', 'refutes_us': True}, uncertainty='Patient bootstrap is conditional on this retrospective panel. Pose +/-0.05mm is a scenario, not measured registration uncertainty; target-domain FDI accuracy UNKNOWN.', limits={'threads': 4, 'disk_bytes': 3000000000, 'ram_gb': 3})
    freeze(ROOT / 'PREREG_R2.json', cfg)
    source = []
    lock = read(V3 / 'BENCHMARK_LOCK.json')
    expected = lock.get('files', {})
    for row in selected:
        key = row['case_key']
        tasks = read(V3 / 'payload/public/tasks' / f'{key}.json')
        files = {f'public/tasks/{key}.json', f'private/references/{key}.npz', f'private/provenance/{key}.json'}
        files.update(('public/' + t['geometry_file'] for t in tasks if t['status'] == 'READY'))
        for rel in sorted(files):
            src = V3 / 'payload' / rel
            dst = DATA / rel
            bound = expected.get('payload/' + rel)
            if not bound or sha(src) != bound['sha256']:
                raise RuntimeError('v3 manifest mismatch ' + rel)
            copy(src, dst)
            source.append(dict(source=str(src), relative='payload/' + rel, sha256=sha(dst), bytes=dst.stat().st_size))
    dump(DATA / 'private/COHORT.json', selected)
    shutil.copytree(V3 / 'payload/runtime', DATA / 'runtime', dirs_exist_ok=True)
    base = DATA / 'participant_code'
    base.mkdir(exist_ok=True)
    for name in ['task_io.py', 'fast_geometry.py', 'util.py']:
        copy(V3 / 'code' / name, base / name)
    shutil.copytree(V3 / 'code/legacy', base / 'legacy', ignore=shutil.ignore_patterns('__pycache__'), dirs_exist_ok=True)
    copy(FIELD / 'field_generator.py', base / 'field/field_generator.py')
    copy(FIELD / 'roof_manufacturing_v1.py', base / 'field/roof_manufacturing_v1.py')
    for p in (FIELD / 'field_snapshot').rglob('*'):
        if p.is_file() and '__pycache__' not in str(p):
            copy(p, base / 'field/field_snapshot' / p.relative_to(FIELD / 'field_snapshot'))
    copy(FIELD / 'native_mesh_sdf.so', base / 'field/native_mesh_sdf.so')
    copy(X42 / 'code/generator.py', base / 'x42_generator.py')
    xf = read(X42 / 'FROZEN_GENERATOR.json')
    for (rel, item) in xf['models'].items():
        if rel.split('/')[0] not in ['kernel', 'kernel_contact']:
            continue
        src = Path(xf['model_root']) / rel
        if sha(src) != item['sha256']:
            raise RuntimeError('X42 frozen model drift')
        copy(src, DATA / 'models' / rel)
    dump(DATA / 'public/PARTICIPANTS.json', {k: cfg[k] for k in ['participants', 'x1b_tool_radius_mm']})
    for (src, dst) in [(ROOT.parent / 'LANE_X13_CEMENT_GAP/FACIT.csv', 'CEMENT.csv'), (ROOT.parent / 'LANE_X36_META_REGRESSION/EXTRACTION_CROWN_CORRECTED_R2.csv', 'CROWN.csv')]:
        copy(src, DATA / 'literature' / dst)
        source.append(dict(source=str(src), relative='payload/literature/' + dst, sha256=sha(src), bytes=src.stat().st_size))
    dump(ROOT / 'SOURCE_REUSE.json', dict(files=source, field_code_sha256=sha(FIELD / 'field_generator.py'), field_operation='Unmodified code and dependency snapshot', x42_generator_freeze_sha256=sha(X42 / 'FROZEN_GENERATOR.json'), v3_lock_sha256=sha(V3 / 'BENCHMARK_LOCK.json'), copy_seconds=time.perf_counter() - start))
    state('PANEL_PREPARED', 'Frozen 64-case selection copied with source hashes', 'R1 discrimination then isolated participant generation')
if __name__ == '__main__':
    run()
