"""Keep published load cases separate. No universal aligner relaxation curve."""
from dental_release.paths import expand as _release_expand
import json, hashlib
from pathlib import Path
R = Path(__file__).resolve().parents[1]
p = Path(_release_expand('@DENTAL_INPUT_ROOT@/workspace/notes/aligner_intervals/ALIGNER_INTERVALS.json'))
src = json.load(open(p))
out = []
for (idx, kind, propagate) in [(122, 'tensile_coupon_37C_3h', True), (101, 'wet_bending_coupon_37C_24h', True), (131, 'saliva_only_no_applied_load_day1', False), (137, 'saliva_only_no_applied_load_day7', False), (138, 'dry_air_no_applied_load_day1', False), (139, 'dry_air_no_applied_load_day7', False)]:
    a = src[idx]
    lo = a['value_min']
    hi = a['value_max']
    v = a['value_central']
    factor = [v / 100, v / 100] if v is not None else [lo / 100, hi / 100]
    out.append({'load_case': kind, 'source_record_index': idx, 'doi': a['source_doi'], 'locator': a['locator'], 'condition': a['condition'], 'quantity': a['quantity'], 'retained_fraction': factor, 'source_statistic': a['statistic'], 'scale_initial_rigid_ribbon_scenario_only': propagate, 'clinical_tooth_force_propagation': 'UNKNOWN_NOT_MATCHED', 'is_relaxation_under_load': propagate, 'prohibited_operation': 'Do not average with other load cases; do not use unloaded aging as stress relaxation; no interpolation between these observation times'})
base = json.load(open(R / 'raw/PREDICTIONS_R1.json'))['rows']
samples = []
for a in base:
    if a['status'] != 'CONDITIONAL_RIBBON_MODEL' or a['activation_mm'] != 0.2:
        continue
    for r in out:
        if r['scale_initial_rigid_ribbon_scenario_only']:
            samples.append({'case': a['case'], 'jaw': a['jaw'], 'active_fdi': a['active_fdi'], 'load_case': r['load_case'], 'conditional_force_N': [a['formed']['active_force_norm_N'] * q for q in r['retained_fraction']], 'conditional_active_crown_moment_Nmm': [a['formed']['active_moment_norm_Nmm'] * q for q in r['retained_fraction']], 'physical_validation': 'UNKNOWN; assuming coupon factor uniformly scales rigid ribbon E; no patient force claim'})
(R / 'raw/RELAXATION_CASES.json').write_text(json.dumps({'claim_type': 'information_link', 'source_sha256': hashlib.sha256(p.read_bytes()).hexdigest(), 'load_cases': out, 'conditional_scenarios': samples, 'pooled_curve': None}, indent=2) + '\n')
print('Separate relaxation/aging cases:', len(out), 'conditional rows', len(samples))
