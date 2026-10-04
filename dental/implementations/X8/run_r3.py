import csv, hashlib, json, time
from pathlib import Path
import numpy as np
from scipy.stats import binom, beta
from calibrate_clearance import required_blocks, tail_upper, calibrate, REQUIRED
ROOT = Path(__file__).resolve().parent

def dump(name, obj):
    (ROOT / name).write_text(json.dumps(obj, indent=2, ensure_ascii=False) + '\n')

def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()
start = time.perf_counter()
assert digest(ROOT / 'PREREG_R3.json') == (ROOT / 'PREREG_R3.sha256').read_text().strip()
pr = json.loads((ROOT / 'PREREG_R3.json').read_text())
design = pr['design']
x = design['x']
alpha = 1 - design['family_confidence']
rows = []
for m in [1, 4]:
    n = required_blocks(x, alpha, m)
    for count in design['sample_grid']:
        cov = 1 - (1 - x) ** count
        baseline = 1 - float(binom.pmf(0, count, x))
        rows.append({'candidates': m, 'n': count, 'confidence_top_x_seen': cov, 'binomial_control': baseline, 'tail_upper_at_sample_max': tail_upper(count, alpha, m), 'required_confidence': 1 - alpha / m, 'meets_1pct_tail': tail_upper(count, alpha, m) < x})
    assert tail_upper(n, alpha, m) < x and tail_upper(n - 1, alpha, m) >= x
gap = max((abs(r['confidence_top_x_seen'] - r['binomial_control']) for r in rows))
published = {'source': 'doi:10.1111/clr.13578 Table4 MANDIBLE Fu', 'n_implants': 37, 'apical_max_mm': 2.99, 'entry_max_mm': 2.94, 'angle_max_deg': 6.0, 'optimistic_iid_tail_upper_95pct': tail_upper(37, 0.05, 1), 'tail_gate_pass': tail_upper(37, 0.05, 1) < 0.01, 'effective_independent_n': 'UNKNOWN; patient clustering prevents treating 37 implants as 37 independent blocks', 'signed_loss_observations': 'NOT_AVAILABLE', 'phantom_or_clinical_transfer': 'UNKNOWN'}
base = {'independent_block_id': 'FAULT_TEST_ONLY', 'guide': 'fully_guided', 'planned_gap_mm': '2', 'measured_min_tool_gap_mm': '1.5', 'metrology_loss_upper_mm': '.1', 'source_locator': 'FAULT_TEST_NOT_DATA', 'measurement_kind': 'SYNTHETIC', 'population_transfer_verified': 'true'}

def rejected(rr):
    try:
        calibrate(rr)
    except ValueError:
        return True
    return False
dup = {**base, 'measurement_kind': 'INDEPENDENT_PHYSICAL'}
inj = {'synthetic_measurement_rejected': rejected([base]), 'duplicate_independent_block_rejected': rejected([dup, dup]), 'negative_metrology_bound_rejected': rejected([{**dup, 'metrology_loss_upper_mm': '-.1'}]), 'missing_population_transfer_rejected': rejected([{**dup, 'population_transfer_verified': 'false'}]), 'sample_size_one_too_small_rejected': all((tail_upper(required_blocks(x, alpha, m) - 1, alpha, m) >= x for m in [1, 4])), 'confidence_plus_.02_rejected_by_binomial': abs(rows[3]['confidence_top_x_seen'] + 0.02 - rows[3]['binomial_control']) > 1e-12}
assert gap <= 1e-12 and all(inj.values())
with (ROOT / 'CALIBRATION_TEMPLATE.csv').open('w') as f:
    csv.writer(f).writerow(REQUIRED)
dump('RAW_R3.json', {'construction': 'R3', 'formula_gate_pass': gap <= 1e-12, 'max_abs_binomial_error': gap, 'required_independent_blocks_one_guide': required_blocks(x, alpha, 1), 'required_independent_blocks_four_guides': required_blocks(x, alpha, 4), 'published_external_sample': published, 'published_sample_gate': 'FAIL', 'measurement_port': calibrate([]), 'controls': {'same_information_outcome': 'TIE: exact order-statistic coverage equals conventional binomial calculation', 'injected_faults': inj}, 'external_referent': pr['external_referent'], 'secondary_referent': pr['secondary_referent'], 'samples': rows, 'physical_measurements_performed': 0, 'cost': {'wall_seconds': time.perf_counter() - start, 'physical_acquisition': '299 independent blocks for one fixed guide; 437 per guide for four-way selection; not executed'}, 'next_operation': 'Correct source geometry contract: K3 samples only apical 6mm, so certify full-cylinder-to-voxel distance before phantom export'})
dump('CURRENT_WORK_STATE.json', {'lane': 'X8-guide-nerve-risk', 'status': 'R3_COMPLETE_FULL_GEOMETRY_CONSTRUCTION', 'latest_gate': 'Published sample FAIL; actual metrology UNKNOWN', 'next_operation': 'Freeze R4 full-cylinder distance and raw-label phantom export, replacing K3 partial tool surface'})
(ROOT / 'HANDOFF_R3.md').write_text('R3 mathematical calibration operator implemented; no physical measurements. 1% tail requires 299 independent blocks for one fixed guide or 437/guide for four guide types at family confidence95%. Published n=37 FAIL even under optimistic iid. Exact binomial control TIE. New source inspection: K3 nominal body surface samples apical6mm, not the whole cylinder, and EDT-minus-halfvoxel is not exact voxel-union distance. R4 must replace this geometric input before using the R2 moment bound as a whole-body statement. R2 numeric values preserved; full-body eligibility UNKNOWN.\n')
print((ROOT / 'RAW_R3.json').read_text())
