"""Closed-form and adverse checks for the new uncertainty port."""
import copy
import hashlib
import json
import subprocess
import sys
import time
from pathlib import Path
import numpy as np
from assessment_bounds import reference_interval
from voxel_decisions import direct_cube_field
ROOT = Path(__file__).resolve().parent

def main():
    start = time.perf_counter()
    p = ROOT / 'PREREG_R2_REFERENCE_INTERVAL.json'
    if hashlib.sha256(p.read_bytes()).hexdigest() != (ROOT / 'PREREG_R2_REFERENCE_INTERVAL.sha256').read_text().strip():
        raise ValueError('R2 prereg hash drift')
    manifest = {'reference_error_bound_mm': 0.05, 'registration_error_bound_mm': 0.3, 'body_error_bound_mm': 0.0, 'coverage': 'COMPLETE_SHAPES', 'source_locators': ['OUR_ANALYTIC_FIXTURE']}
    gap = float(direct_cube_field(np.array([[0, 0, 0]]), np.array([[0, 0, 10]]), (0.3, 0.3, 0.3))[0])
    shifted = float(direct_cube_field(np.array([[0, 0, 0]]), np.array([[0, 0, 9]]), (0.3, 0.3, 0.3))[0])
    interval = reference_interval(gap, 2, manifest)
    allchecks = []

    def record(name, good, bad):
        allchecks.append({'gate': name, 'valid_pass': bool(good), 'fault_rejected': bool(not bad)})
        if not good or bad:
            raise ValueError(f'R2 gate failed: {name}')
    record('cube_closed_form', abs(gap - 2.7) < 1e-12, abs(gap + 0.05 - 2.7) < 1e-12)
    record('registration_changes_clearance', abs(gap - shifted - 0.3) < 1e-12, abs(gap - shifted + 0.1 - 0.3) < 1e-12)
    record('interval_endpoints', abs(interval['clearance_interval_mm'][0] - (gap - 0.35 - 1e-09)) < 1e-12, abs(interval['clearance_interval_mm'][0] + 0.1 - (gap - 0.35 - 1e-09)) < 1e-12)
    record('missing_is_unknown', reference_interval(gap, 2)['decision'] == 'UNKNOWN_MISSING_BOUND_INFORMATION', reference_interval(gap, 2)['decision'] == 'ABOVE')
    partial = dict(manifest, coverage='SUBSET')
    record('subset_is_unknown', reference_interval(gap, 2, partial)['decision'] == 'UNKNOWN_INCOMPLETE_BOUND_COVERAGE', reference_interval(gap, 2, partial)['decision'] == 'ABOVE')
    near = reference_interval(2.2, 2, manifest)
    record('near_threshold_abstains', near['decision'] == 'ABSTAIN', near['decision'] == 'ABOVE')
    lower = reference_interval(1.0, 2, manifest)
    record('below_threshold', lower['decision'] == 'BELOW', lower['decision'] == 'ABOVE')
    record('physical_unknown', interval['physical_status'] == 'UNKNOWN', interval['physical_status'] == 'CERTIFIED')
    for k in ['reference_error_bound_mm', 'registration_error_bound_mm', 'body_error_bound_mm']:
        removed = dict(manifest)
        removed.pop(k)
        record('missing:' + k, reference_interval(gap, 2, removed)['decision'] == 'UNKNOWN_MISSING_BOUND_INFORMATION', False)
        for value in [-1, float('nan'), True]:
            invalid = dict(manifest)
            invalid[k] = value
            rejected = False
            try:
                reference_interval(gap, 2, invalid)
            except ValueError:
                rejected = True
            record('invalid:' + k + ':' + str(value), rejected, False)
    (ROOT / 'REFERENCE_BOUND_TEMPLATE.json').write_text(json.dumps({k: None for k in manifest if k.endswith('_mm')} | {'coverage': 'UNKNOWN', 'source_locators': [], 'note': 'Nulls deliberately require independent same-site complete-shape evidence; do not insert SDs or population means'}, indent=2) + '\n')
    command = [sys.executable, str(ROOT / 'assess_model.py'), '--input', str(ROOT / 'example_input.npz'), '--threshold-mm', '2', '--reference-kind', 'annotation', '--bound-manifest', str(ROOT / 'REFERENCE_BOUND_TEMPLATE.json'), '--output', str(ROOT / 'EXAMPLE_ASSESSMENT_WITH_BOUNDS.json')]
    subprocess.run(command, check=True, cwd=ROOT)
    result = json.loads((ROOT / 'EXAMPLE_ASSESSMENT_WITH_BOUNDS.json').read_text())
    record('CLI_missing_bounds', result['bounded_reference_decision']['decision'] == 'UNKNOWN_MISSING_BOUND_INFORMATION', result['bounded_reference_decision']['decision'] == 'ABOVE')
    out = {'id': 'X41_R2_REGISTRATION_AWARE_REFERENCE_INTERVAL', 'claim_type': 'capability', 'outcome': 'CONDITIONAL_REFERENCE_INTERVAL_PORT_DELIVERED_PHYSICAL_UNKNOWN', 'checks': allchecks, 'closed_form_gap_mm': gap, 'translated_gap_mm': shifted, 'interval': interval, 'near_threshold': near, 'physical_certificates': 0, 'wall_seconds': time.perf_counter() - start, 'external_referent': {'kind': 'closed_form', 'locator': 'closed voxel cube separation; complete-shape distance triangle inequality', 'compared_quantity': 'deterministic clearance change under bounded translation; no lab data'}}
    (ROOT / 'RESULTS_R2.json').write_text(json.dumps(out, indent=2) + '\n')
    print(json.dumps({'R2': 'PASS', 'controls': len(allchecks), 'physical': 'UNKNOWN', 'wall_seconds': out['wall_seconds']}))
if __name__ == '__main__':
    main()
