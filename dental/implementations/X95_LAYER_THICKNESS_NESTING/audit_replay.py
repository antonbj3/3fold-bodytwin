"""Relocated replay and executed-copy fault checks; read-only predecessor data."""
from dental_release.paths import expand as _release_expand
from pathlib import Path
import json, shutil, subprocess, hashlib, time, os
ROOT = Path(__file__).resolve().parent
AUDIT = Path(_release_expand('@DENTAL_EXTERNAL_ROOT@/storage/tmp_dental_audit/X95-layer-thickness-nesting'))

def science(x):
    return {'gates': x['primary_B05_gate'], 'R1': [{k: r[k] for k in ['family', 'requested_poses', 'stock_fitting_poses', 'allowed_yml_poses', 'allowed_universal0_5_poses', 'sampled_min_wall_mm']} for r in x['rounds']['R1']['rows']], 'R2': [{k: r[k] for k in ['family', 'eligible_cubes', 'source_minus_universal_removed_mm3', 'selected_guard_min_mm', 'mesh_watertight']} for r in x['rounds']['R2']['rows']], 'R3': [{k: r.get(k) for k in ['family', 'status', 'source_retained_volume_mm3', 'watertight']} for r in x['rounds']['R3']['rows']], 'R4': [{k: r[k] for k in ['family', 'source_retained_volume_mm3', 'universal0_5_retained_volume_mm3', 'source_base_xy_grid', 'source_base_z_grid', 'guard_min_used_cubes_mm', 'watertight', 'outer_identity_error_mm']} for r in x['rounds']['R4']['rows']], 'R5': [{k: r[k] for k in ['family', 'point_min_slack_mm', 'negative_point_slacks', 'fault_core_translation_min_slack_mm', 'independent_point_triangle_max_error_mm']} for r in x['rounds']['R5']['rows']], 'reference': [r['original_surface_p95_mm'] for r in x['final_validation']['rows']], 'source_table': x['rounds']['R1']['source_only_strength_control']['rows'], 'sufficiency': x['sufficiency']}

def run():
    AUDIT.mkdir(parents=True, exist_ok=True)
    dest = AUDIT / 'replay_01'
    if dest.exists():
        raise RuntimeError('Preserve existing audit copy; choose new run name')
    shutil.copytree(ROOT, dest, ignore=shutil.ignore_patterns('__pycache__'))
    start = time.perf_counter()
    with (AUDIT / 'REPLAY_COPY.log').open('w') as log:
        r = subprocess.run(['./run_all.sh'], cwd=dest, stdout=log, stderr=subprocess.STDOUT)
    assert r.returncode == 0
    original = json.loads((ROOT / 'results.json').read_text())
    copied = json.loads((dest / 'results.json').read_text())
    (a, b) = (science(original), science(copied))
    assert a == b, 'relocated scientific outputs differ'
    code = dest / 'code/pose_join.py'
    payload = code.read_bytes()
    code.write_bytes(payload + b'\n# injected copied-code mutation\n')
    fault = subprocess.run(['./run_all.sh', '--verify-only'], cwd=dest, capture_output=True, text=True)
    assert fault.returncode != 0 and 'RELEASE_HASH_MISMATCH' in fault.stderr
    code.write_bytes(payload)
    rule = dest / 'RULE_CONTRACT.json'
    payload2 = rule.read_bytes()
    mut = json.loads(payload2)
    mut['crown_rules'][0]['footnotes'] = ['*1', '*3']
    rule.write_text(json.dumps(mut, indent=2) + '\n')
    fault2 = subprocess.run(['./run_all.sh', '--verify-only'], cwd=dest, capture_output=True, text=True)
    assert fault2.returncode != 0 and 'RELEASE_HASH_MISMATCH' in fault2.stderr
    rule.write_bytes(payload2)
    restored = subprocess.run(['./run_all.sh', '--verify-only'], cwd=dest, capture_output=True, text=True)
    assert restored.returncode == 0
    result = {'audit_copy': str(dest), 'own_source_tree_copied_only': True, 'relocated_replay_exit': r.returncode, 'scientific_metrics_identical_exactly': True, 'copied_code_mutation_rejected': True, 'copied_rule_mutation_rejected': True, 'restored_copy_verify_exit': 0, 'seconds': time.perf_counter() - start, 'original_result_sha256': hashlib.sha256((ROOT / 'results.json').read_bytes()).hexdigest(), 'copy_result_sha256': hashlib.sha256((dest / 'results.json').read_bytes()).hexdigest(), 'costs_differ_and_not_compared_as_science': True, 'source_only_physical_scope': 'unchanged UNKNOWN', 'original_source_code_and_inputs_not_mutated_by_audit': True}
    (ROOT / 'raw/AUDIT_REPLAY.json').write_text(json.dumps(result, indent=2) + '\n')
    (AUDIT / 'AUDIT_REPLAY.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))
if __name__ == '__main__':
    run()
