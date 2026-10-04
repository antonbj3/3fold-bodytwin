from dental_release.paths import expand as _release_expand
import csv, datetime, hashlib, json, resource, sys, time, traceback, zipfile
from pathlib import Path
import numpy as np
from full_geometry import voxel_cylinder_bracket
from risk_operator import body_bound, scenario_probability
ROOT = Path(__file__).resolve().parent
DATA = Path(_release_expand('@DENTAL_WORK_ROOT@/X8-guide-nerve-risk'))
sys.path.insert(0, str(Path(__import__('os').environ.get('DENTAL_PROJECT_ROOT', str(ROOT.parent.parent))) / 'cells/geometry'))
import tf2_io

def dump(name, obj):
    (ROOT / name).write_text(json.dumps(obj, indent=2, ensure_ascii=False) + '\n')

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()
start = time.perf_counter()
assert sha(ROOT / 'PREREG_R5.json') == (ROOT / 'PREREG_R5.sha256').read_text().strip()
previous_cost = json.loads((ROOT / 'RAW_R5.json').read_text()).get('cost', {}) if (ROOT / 'RAW_R5.json').exists() else {}
pr = json.loads((ROOT / 'PREREG_R5.json').read_text())
profiles = json.loads((ROOT / 'GUIDE_PROFILES.json').read_text())['profiles']
sites = [r for r in json.loads((DATA / 'K3_SITES.json').read_text())['sites'] if r.get('L_plan') is not None and 'd_nom_implant_body' in r]
groups = {}
for r in sites:
    groups.setdefault(r['case'], []).append(r)
dest = DATA / 'R5_FULL_GEOMETRY.jsonl'
old_by_key = {(r['case'], r['fdi']): r for r in sites}
completed = {}
failures = []
if dest.exists():
    for line in dest.read_text().splitlines():
        r = json.loads(line)
        if r['prereg_sha256'] != sha(ROOT / 'PREREG_R5.json'):
            raise ValueError('R5 resume hash mismatch')
        completed[r['case'], r['fdi']] = r
resumed_count = len(completed)
with zipfile.ZipFile(tf2_io.ZIP) as z, dest.open('a') as rawout:
    for (count, (case, rr)) in enumerate(sorted(groups.items()), 1):
        pending = [r for r in rr if (case, r['fdi']) not in completed]
        if pending:
            member = tf2_io.ROOT + '/labelsTr/' + case + '.mha'
            b = z.read(member)
            source_hash = hashlib.sha256(b).hexdigest()
            (lab, sp, hdr) = tf2_io.read_mha_bytes(b)
            sp = np.array(sp)
            centers = np.argwhere(np.isin(lab, [3, 4])) * sp
            nv = Path(__import__('os').environ['DENTAL_PROJECT_ROOT']) / 'results/NV1_canals/per_case' / f'{case}.json'
            recs = json.loads(nv.read_text())['teeth']
            for s in pending:
                try:
                    rec = recs[str(s['fdi'])]
                    E = np.array(rec['crest_entry'])
                    a = -np.array(rec['axis'])
                    a /= np.linalg.norm(a)
                    ans = voxel_cylinder_bracket(centers, sp, E, a, s['L_plan'], 2.0)
                    result = {'case': case, 'fdi': s['fdi'], 'source_member': member, 'source_member_sha256': source_hash, 'plan_file_sha256': sha(nv), 'entry_zyx_mm': E.tolist(), 'axis_zyx': a.tolist(), 'length_mm': s['L_plan'], 'radius_mm': 2.0, 'old_K3_gap_mm': s['d_nom_implant_body'], 'prereg_sha256': sha(ROOT / 'PREREG_R5.json'), **ans}
                    rawout.write(json.dumps(result) + '\n')
                    rawout.flush()
                    completed[case, s['fdi']] = result
                except Exception as exc:
                    failures.append({'case': case, 'fdi': s['fdi'], 'error': repr(exc), 'traceback': traceback.format_exc()})
            del b, lab, centers
        if count % 25 == 0 or count == len(groups):
            dump('CURRENT_WORK_STATE.json', {'lane': 'X8-guide-nerve-risk', 'status': 'R5_RUNNING', 'updated_at': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'cases_seen': count, 'total_cases': len(groups), 'completed_sites': len(completed), 'latest_gate': 'IN_PROGRESS', 'next_operation': 'Complete whole-body cohort; retained raw certificates'})
            print(f'{count}/{len(groups)} cases; {len(completed)}/{len(sites)} sites; {time.perf_counter() - start:.1f}s', flush=True)
valid = []
uncertain = []
for s in sites:
    r = completed.get((s['case'], s['fdi']))
    if r is None or r['gap_mm'] > 0.0001:
        uncertain.append({'case': s['case'], 'fdi': s['fdi'], 'reason': 'MISSING' if r is None else 'NUMERICAL_GAP', 'gap_mm': None if r is None else r['gap_mm']})
    else:
        valid.append(r)
dlo = np.array([r['lower_mm'] for r in valid])
dhi = np.array([r['upper_mm'] for r in valid])
dmid = (dlo + dhi) / 2
rows = []
sums = []
for p in profiles:
    scenario = scenario_probability(dmid, p)
    bounds = {str(B): body_bound(dlo, p, B) for B in [0, 0.3, 0.7]}
    sums.append({'guide': p['id'], 'n_complete_sites': len(valid), 'scenario_eligible_1pct': int((scenario < 0.01).sum()), 'moment_condition_eligible_1pct_B0p3': int((bounds['0.3'] < 0.01).sum()), 'scenario_gt1pct_among_full2mm_accepted': int(((dlo >= 2) & (scenario > 0.01)).sum()), 'clinical_guide_selection': 'UNKNOWN'})
    for (j, r) in enumerate(valid):
        rows.append({'case': r['case'], 'fdi': r['fdi'], 'guide': p['id'], 'old_K3_gap_mm': r['old_K3_gap_mm'], 'whole_cylinder_label_gap_lower_mm': r['lower_mm'], 'whole_cylinder_label_gap_upper_mm': r['upper_mm'], 'nominal_2mm_class': 'BELOW' if r['upper_mm'] < 2 else 'ABOVE' if r['lower_mm'] >= 2 else 'UNKNOWN', 'local_apex_scenario_p_lt1': float(scenario[j]), 'whole_body_moment_upper_B0p3': float(bounds['0.3'][j]), 'whole_body_moment_upper_B0': float(bounds['0'][j]), 'whole_body_moment_upper_B0p7': float(bounds['0.7'][j]), 'tail_population_validated': False, 'true_anatomy_validated': False, 'clinical_risk': 'UNKNOWN'})
with (ROOT / 'PER_SITE_FULL_GEOMETRY_RISK.csv').open('w') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]))
    w.writeheader()
    w.writerows(rows)
changed = int((dhi < 2).sum())
coverage = len(valid) / len(sites)
out = {'construction': 'R5', 'planned_sites': len(sites), 'cases': len(groups), 'completed_sites': len(completed), 'closed_brackets': len(valid), 'coverage': coverage, 'maximum_bracket_gap_mm': max((r['gap_mm'] for r in completed.values())), 'old_2mm_accepted': len(sites), 'full_label_gap_below2mm': changed, 'practice_class_change_fraction': changed / len(sites), 'full_label_gap_above2mm': int((dlo >= 2).sum()), 'median_old_minus_full_gap_mm': float(np.median(np.array([r['old_K3_gap_mm'] for r in valid]) - dmid)), 'gates': {'cohort_coverage': coverage >= 0.99, 'distance_bracket': all((r['gap_mm'] <= 0.0001 for r in completed.values())), 'practice_class_change': changed / len(sites) >= 0.05, 'clinical_guide_choice': 'UNKNOWN'}, 'profiles': sums, 'external_referent': pr['external_referent'], 'method_comparison': 'TIE; gain is consuming complete 3D geometry, not defeating a matched conventional solver', 'cost': {'wall_seconds': time.perf_counter() - start, 'maxrss_kib': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, 'workers': 1, 'threads': 1, 'gpu': False, 'resumed_certificate_count': resumed_count, 'cold_geometry_wall_seconds': previous_cost.get('cold_geometry_wall_seconds', previous_cost.get('wall_seconds', time.perf_counter() - start))}, 'failures': failures, 'uncertain': uncertain, 'raw_certificates': {'path': str(dest), 'bytes': dest.stat().st_size, 'sha256': sha(dest)}, 'limits': ['Actual canal boundary and completeness UNKNOWN', 'Moment population transfer and guide-tail validation UNKNOWN', 'No sensory outcomes', 'Virtual immediate tooth-axis plans', 'No patient recommendations']}
dump('RAW_R5.json', out)
frozen = {'frozen_at': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'prereg_sha256': sha(ROOT / 'PREREG_R5.json'), 'per_site_prediction_sha256': sha(ROOT / 'PER_SITE_FULL_GEOMETRY_RISK.csv'), 'physical_measurements_performed': False, 'use': 'Future independent research measurements; preserve source-exposed selection'}
if (ROOT / 'FROZEN_FULL_GEOMETRY_PREDICTIONS.json').exists():
    assert json.loads((ROOT / 'FROZEN_FULL_GEOMETRY_PREDICTIONS.json').read_text())['per_site_prediction_sha256'] == frozen['per_site_prediction_sha256']
else:
    dump('FROZEN_FULL_GEOMETRY_PREDICTIONS.json', frozen)
dump('CURRENT_WORK_STATE.json', {'lane': 'X8-guide-nerve-risk', 'status': 'R5_COMPLETE', 'latest_gate': out['gates'], 'next_operation': 'Independent metrology and signed clearance-loss acquisition; clinical outcome closure remains UNKNOWN'})
print(json.dumps({k: v for (k, v) in out.items() if k not in ['failures', 'uncertain']}, indent=2), flush=True)
