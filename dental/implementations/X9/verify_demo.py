"""Consequential guards and frozen experiment integrity; failures stay visible."""
from dental_release.paths import expand as _release_expand
import json
from pathlib import Path
from ipr import sha
import numpy as np
from segmented_tooth_tool import build_protection, envelope

def run():
    r = json.loads(Path('results.json').read_text())
    required_faults = {'R1': ['population_plus_2_mm_rejected', 'control_plus_0_1_mm_rejected', 'coverage_all_rays_removed_rejected', 'praxis_equal_0_5_all_patches_rejected', 'uncalibrated_physical_safety_claim_rejected'], 'R2': ['external_plus_2_mm_rejected', 'control_plus_0_1_rejected', 'zero_contrast_profile_rejected', 'cut_plus_1_rejected_all_complete_patches'], 'R3': ['solver_plus_0_1_rejected', 'cut_plus_2_rejected_all', 'empty_tissue_mask_rejected', 'external_thickness_plus_2_rejected', 'missing_anatomical_calibration_blocks_positive_safety']}
    assert {x['round'] for x in r['rounds']} == set(required_faults) and len(r['rounds']) == len(required_faults)
    for result in r['rounds']:
        assert sha(f"PREREG_{result['round']}.json") == result['prereg_sha256']
        assert all((result.get('injected_faults', {}).get(k) is True for k in required_faults[result['round']])), result.get('injected_faults')
        for item in result.get('array_manifest', []):
            assert sha(item['path']) == item['sha256']
    tool = json.loads(Path('artifacts/tool_output.json').read_text())
    assert tool['calibrated_safe_ipr_mm'] is None
    assert tool['parity_error_mm'] <= 1e-07
    assert tool['protected_overlap_cells'] == 0
    assert tool['injected_plus_2_mm_rejected']
    assert all((p.stat().st_size <= 50000000 for p in Path('.').rglob('*') if p.is_file()))
    data = Path(_release_expand('@DENTAL_WORK_ROOT@/X9-ipr-safety'))
    assert sum((p.stat().st_size for p in data.rglob('*') if p.is_file())) < 3000000000
    h = 0.15
    axis = np.arange(-4.5, 4.5001, h)
    (x, y, z) = np.meshgrid(axis, axis, axis, indexing='ij')
    radius = np.sqrt(x * x + y * y + z * z)
    labels = np.where(radius < 4.0, 1, 0).astype(np.uint8)
    labels[radius < 1.0] = 2
    (protected, _) = build_protection(labels, np.full(3, h))
    (result, _) = envelope(labels, protected, np.full(3, h), np.full(3, axis[0]), 1, 0.0, 0.0)
    rp = max(1.0 + 0.5 + 0.3 + np.sqrt(3) * h, (4.0 + 1.0 + 0.3) / 2)
    ideal_capacity = 4.0 - rp - 0.05 - np.sqrt(3) * h / 2
    assert result['scenario_capacity_mm'] > 0.5
    assert abs(result['scenario_capacity_mm'] - ideal_capacity) <= 2 * h
    assert result['protected_overlap_cells'] == 0
    assert result['parity_error_mm'] <= 1e-07
    assert result['injected_plus_2_mm_rejected']
    Path('rounds/operator_analytic_probe.json').write_text(json.dumps({'external_referent': {'kind': 'our_own_fixture', 'locator': 'verify_demo.py concentric spheres', 'compared_quantity': 'closed-form model cut-envelope radius, mm', 'refutes_us': False}, 'used_as_anatomical_facit': False, 'ideal_capacity_mm': ideal_capacity, 'voxel_capacity_mm': result['scenario_capacity_mm'], 'frozen_numerical_tolerance_mm': 2 * h, 'injected_overcut_rejected': result['injected_plus_2_mm_rejected']}, indent=2) + '\n')
    print('Verification PASS: prereg hashes, fault rejections, tool CLI, protected tissue, disk limits')
    print('Scientific coverage/type gates may FAIL as recorded; anatomical safety remains UNKNOWN')
if __name__ == '__main__':
    run()
