import json, hashlib, csv, math, statistics
from pathlib import Path
from datetime import datetime, timezone
R = Path(__file__).resolve().parent

def dump(p, d):
    p.write_text(json.dumps(d, ensure_ascii=False, indent=2, allow_nan=False) + '\n')
assert hashlib.sha256((R / 'PREREG_R1_IDENTIFIABILITY.json').read_bytes()).hexdigest() == (R / 'PREREG_R1_IDENTIFIABILITY.sha256').read_text().strip()
A = [0, 3, 3] * 6 + [0, 0]
B = [1, 1, 4] * 6 + [0, 0]

def summary(x):
    mean = sum(x) / len(x)
    v = sum((i * i for i in x)) / len(x) - mean * mean
    return {'mean_mm': mean, 'variance_mm2': v, 'sd_mm': math.sqrt(v)}
(sA, sB) = (summary(A), summary(B))
assert sA == sB
q = lambda a: sorted(a)[math.ceil(0.95 * len(a)) - 1]
sites = json.loads((R / 'inputs/lineage_sites.json').read_text())
groups = {}
for s in sites:
    groups.setdefault((s['image_group'], s['fdi']), []).append(s)
flips = sum((any((s['classes']['tf2'] == 'ABOVE' and s['union']['class_2mm'] == 'BELOW' for s in rr)) for rr in groups.values()))
segments = {}
for s in sites:
    d = segments.setdefault(s['segment'], {'pose_rows': 0, 'image_groups': set(), 'losses': []})
    d['pose_rows'] += 1
    d['image_groups'].add(s['image_group'])
    d['losses'].append(s['loss_upper_mm'])
for d in segments.values():
    d['image_groups'] = len(d['image_groups'])
    l = d.pop('losses')
    d['observed_loss_max_mm'] = max(l)
out = {'claim_type': 'information_link', 'physical_margin': 'UNKNOWN', 'physical_gate': 'FAIL_IDENTIFIABILITY', 'sufficiency': {'moment_test': {'A_mm': A, 'B_mm': B, 'summary_A': sA, 'summary_B': sB, 'identity_error': 0.0, 'quantile95_A_mm': q(A), 'quantile95_B_mm': q(B), 'quantile_difference_mm': q(B) - q(A), 'nominal_gap_mm': 3.5, 'breach_fraction_A': sum((v > 3.5 for v in A)) / 20, 'breach_fraction_B': sum((v > 3.5 for v in B)) / 20, 'minimum_extension': 'Signed whole-body clearance-loss distribution or externally valid upper-tail bound'}, 'direction_test': {'magnitudes_A_mm': [2.0] * 20, 'magnitudes_B_mm': [2.0] * 20, 'identity_error': 0.0, 'canal_normal_loss_A_mm': 2.0, 'canal_normal_loss_B_mm': 0.0, 'downstream_difference_mm': 2.0, 'minimum_extension': 'Canal-normal orientation plus loss at the active wall witness'}, 'external_referent': {'kind': 'our_own_fixture', 'locator': 'r1_identify.py', 'compared_quantity': 'Logical nonidentifiability, not anatomy', 'refutes_us': True}}, 'local_census': {'pose_rows': len(sites), 'unique_image_FDI_sites': len(groups), 'reversals': flips, 'reversals_per100': 100 * flips / len(groups), 'resolution': 'POPULATION', 'raw_resolution': 'PER_TOOTH', 'segments': segments}, 'next_operation': 'R2 finite conditional 90/95 union bounds using existing lineage split; explicit heldout wall gate', 'external_referent': {'kind': 'published_dataset', 'locator': str(R.parent / 'LANE_X73_CANAL_WALL_FULL/raw/LINEAGE_SITES.json'), 'compared_quantity': 'Digital mask-version cylinder clearance and decision reversals', 'refutes_us': True}}
dump(R / 'raw/R1_RESULT.json', out)
dump(R / 'CURRENT_WORK_STATE.json', {'updated_utc': datetime.now(timezone.utc).isoformat(), 'stage': 'R1_COMPLETE', 'latest_gate': 'Physical 90/95 table not identified; same moments differ by 1mm in P95 and 30 percentage points in breach fraction', 'next_operation': out['next_operation']})
(R / 'HANDOFF_R1.md').write_text("# R1 : tail and direction missing\n\nModernistic own counter-examples have exactly 0 identity error but P95 3 and 4 mm respectively and breakthrough share 0 and 0,30 respectively. This is a logical counter-test, no clinical reference. The local digital mask revisions change 12 / 741 sites ( 1,619 / 100 ). Physical 90/95 margins are UNKNOWN. The next operation will change the contract to conditional's full-body limits, with targeted audit losses and preserved wholedout-gate.\n")
print(json.dumps({'identity_error': 0.0, 'q95_A': q(A), 'q95_B': q(B), 'sites': len(groups), 'flips': flips, 'physical_margin': 'UNKNOWN'}))
