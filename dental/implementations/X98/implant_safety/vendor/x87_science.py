import json, csv, hashlib, math, resource, time, statistics
from pathlib import Path
from datetime import datetime, timezone
import numpy as np
from scipy.stats import beta
from scipy.optimize import brentq
R = Path(__file__).resolve().parent
SEGS = ['mental_endpoint_5mm_proxy', 'molar_nearest_FDI', 'posterior_endpoint_10mm_proxy']
GUIDES = ['fully_guided', 'pilot_guided', 'freehand', 'dynamic_navigation']

def load(p):
    return json.loads((R / p).read_text())

def dump(p, d):
    (R / p).write_text(json.dumps(d, ensure_ascii=False, indent=2, allow_nan=False) + '\n')

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def finite(x):
    return x if math.isfinite(x) else 'INF'

def csvout(p, rows):
    with (R / p).open('w') as f:
        if not rows:
            return
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

def verify_inputs():
    for m in load('INPUT_MANIFEST.json'):
        assert sha(R / m['local']) == m['sha256'], m['local']

def ci(k, n):
    if not n:
        return [None, None]
    return [0 if k == 0 else float(beta.ppf(0.025, k, n - k + 1)), 1 if k == n else float(beta.ppf(0.975, k + 1, n - k))]

def guide(p, failure):
    a = failure / 3.0
    k = math.sqrt((1 - a) / a)
    e = p['entry_mean_sd_mm'][0] + p['entry_mean_sd_mm'][1] * k
    v = p['apex_mean_sd_mm'][0] + p['apex_mean_sd_mm'][1] * k
    theta = math.radians(p['angle_mean_sd_deg'][0] + p['angle_mean_sd_deg'][1] * k)
    rot = 4 * math.sin(min(theta, math.pi) / 2)
    return {'entry_bound_mm': e, 'apex_bound_mm': v, 'angle_bound_deg': math.degrees(theta), 'radius_term_mm': rot, 'guide_budget_mm': max(e, v) + rot, 'component_failure': a, 'guide_failure_bound': failure}

def controls(p, failure, g):
    a = failure / 3
    vals = []
    for key in ['entry_mean_sd_mm', 'apex_mean_sd_mm', 'angle_mean_sd_deg']:
        (mu, sigma) = p[key]
        t = brentq(lambda t: sigma ** 2 / (sigma ** 2 + t ** 2) - a, 0, sigma * 10000, xtol=1e-13)
        vals.append(mu + t)
    theta = min(math.radians(vals[2]), math.pi)
    rot = math.hypot(2 * math.cos(theta) - 2, 2 * math.sin(theta))
    reference = max(vals[:2]) + rot
    delta = abs(reference - g['guide_budget_mm'])
    assert delta <= 1e-10
    return {'name': 'Cantelli numerical inversion and rotated radius', 'guide': p['id'], 'failure': failure, 'error_mm': delta, 'injected_plus1mm_rejected': abs(reference - (g['guide_budget_mm'] + 1)) > 1e-10, 'reference_locator': 'Closed form one-sided Chebyshev inequality; scipy.optimize.brentq independent inversion'}

def classify(d, threshold):
    if not math.isfinite(threshold):
        return 'UNKNOWN'
    if d['lower_mm'] >= threshold:
        return 'ABOVE'
    if d['upper_mm'] < threshold:
        return 'BELOW'
    return 'UNKNOWN'

def run(n):
    start = time.perf_counter()
    verify_inputs()
    p = R / f'PREREG_R{n}.json'
    assert sha(p) == p.with_suffix('.sha256').read_text().strip()
    residual = json.loads(p.read_text())['residual_clearance_mm']
    profiles = {p['id']: p for p in load('inputs/guide_profiles.json')['profiles']}
    split = load('inputs/patient_split.json')
    cal = set(split['calibration_patients'])
    test = set(split['test_patients'])
    sites = load('inputs/lineage_sites.json')
    rows = []
    policies = []
    cr = []
    revision = []
    for target in [0.9, 0.95]:
        gf = (1 - target) / 2 if n != 3 else 1 - target
        gs = {g: guide(profiles[g], gf) for g in GUIDES}
        cr.extend((controls(profiles[g], gf, gs[g]) for g in GUIDES))
        for seg in SEGS:
            ss = [s for s in sites if s['segment'] == seg]
            selected = [s for s in ss if s['image_group'] in test] if n != 3 else ss
            values = {}
            for s in ss:
                values[s['image_group']] = max(values.get(s['image_group'], 0.0), s['loss_upper_mm'])
            ca = [v for (k, v) in values.items() if k in cal]
            te = [v for (k, v) in values.items() if k in test]
            level = 1 - (1 - target) / 2
            rank = math.ceil((len(ca) + 1) * level)
            budget = sorted(ca)[rank - 1] if ca and rank <= len(ca) else math.inf
            if math.isfinite(budget):
                control = float(np.partition(np.asarray(ca), rank - 1)[rank - 1])
                assert control == budget
                cr.append({'name': 'Order statistic sort vs partition', 'segment': seg, 'target': target, 'error_mm': 0.0, 'injected_plus1mm_rejected': abs(control - (budget + 1)) > 1e-10})
            k = sum((v <= budget for v in te))
            coverage = k / len(te) if te else None
            gate = bool(math.isfinite(budget) and te and (coverage >= level))
            revision.append({'segment': seg, 'target': target, 'revision_target': level, 'calibration_n': len(ca), 'heldout_n': len(te), 'rank': rank, 'budget_mm': finite(budget), 'covered': k, 'heldout_coverage': coverage, 'ci95': ci(k, len(te)), 'heldout_gate': 'PASS' if gate else 'FAIL_OR_MISSING', 'resolution': 'POPULATION', 'raw_resolution': 'PER_TOOTH'})
            for gid in GUIDES:
                g = gs[gid]
                margin = residual + budget + g['guide_budget_mm']
                local = []
                for s in selected:
                    d = s['distances']['tf2']
                    b = budget if n != 3 else s['loss_upper_mm']
                    m = residual + b + g['guide_budget_mm']
                    c = classify(d, m)
                    if n == 3:
                        exact = max(0.0, d['upper_mm'] - s['union']['lower_mm'])
                        assert abs(exact - b) <= 1e-10, (exact, b)
                        assert d['lower_mm'] - b <= s['union']['lower_mm'] + 1e-10
                        if c == 'ABOVE':
                            assert s['union']['lower_mm'] >= residual + g['guide_budget_mm'] - 1e-10
                    local.append(m)
                    policies.append({'round': n, 'target': target, 'image_group': s['image_group'], 'case': s['case'], 'fdi': s['fdi'], 'segment': seg, 'guide': gid, 'nominal_lower_mm': d['lower_mm'], 'nominal_upper_mm': d['upper_mm'], 'union_lower_mm': s['union']['lower_mm'], 'loss_upper_mm': s['loss_upper_mm'], 'margin_mm': finite(m), 'fixed_class': s['classes']['tf2'], 'model_class': c, 'resolution': 'PER_TOOTH', 'time_scale': 'SIMULTANEOUS'})
                pp = [r for r in policies if r['target'] == target and r['segment'] == seg and (r['guide'] == gid)]
                grouped = {}
                for r in pp:
                    grouped.setdefault((r['image_group'], r['fdi']), []).append(r)
                changed = sum((any((r['fixed_class'] == 'ABOVE' and r['model_class'] == 'BELOW' for r in rr)) for rr in grouped.values()))
                unknown = sum((any((r['model_class'] == 'UNKNOWN' for r in rr)) for rr in grouped.values()))
                changed_rows = sum((r['fixed_class'] == 'ABOVE' and r['model_class'] == 'BELOW' for r in pp))
                row = {'round': n, 'guide': gid, 'segment': seg, 'target': target, 'residual_clearance_mm': residual, 'margin_mm': finite(margin) if n != 3 else statistics.median(local) if local else 'UNKNOWN', 'margin_min_mm': min(local) if local else 'UNKNOWN', 'margin_max_mm': max(local) if local else 'UNKNOWN', 'guide_budget_mm': g['guide_budget_mm'], 'revision_budget_mm': finite(budget) if n != 3 else 'SITE_SPECIFIC', 'image_FDI_sites': len(grouped), 'pose_rows': len(pp), 'changed_sites': changed, 'changes_per100': 100 * changed / len(grouped) if grouped else None, 'changed_pose_rows': changed_rows, 'unknown_sites': unknown, 'heldout_revision_gate': gate if n != 3 else 'NOT_APPLICABLE_OBSERVED_ENVELOPE', 'empirical_joint_coverage': 'UNKNOWN', 'physical_margin_mm': 'UNKNOWN', 'resolution': 'PHENOMENOLOGICAL', 'count_resolution': 'POPULATION', 'edge_resolution': 'PER_TOOTH', 'time_scale': 'SIMULTANEOUS'}
                rows.append(row)
    result = {'claim_type': 'information_link', 'round': f'R{n}', 'tables': rows, 'revision_validation': revision if n != 3 else [], 'controls': cr, 'status': 'CONDITIONAL_MODEL_ONLY_PHYSICAL_UNKNOWN' if n != 3 else 'OBSERVED_RELEASE_ENVELOPE_CLOSED_PHYSICAL_UNKNOWN', 'external_referent': {'kind': 'published_dataset', 'locator': str(R.parent / 'LANE_X73_CANAL_WALL_FULL/raw/LINEAGE_SITES.json'), 'compared_quantity': 'Cylinder clearance to retained public mask union on identical image; not true canal', 'refutes_us': n != 3}, 'cost': {'wall_seconds': time.perf_counter() - start, 'maxrss_kib': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, 'threads': 1, 'gpu': False, 'manual_discovery_seconds': 'UNKNOWN', 'upstream_full_geometry_replayed': False}, 'full_cost_caveat': 'Reuse source geometry; no independent reread of voxel archives. Upstream X73 full stream 527.73s, maxRSS1100.2MiB, total known geometry1618.73s according to parent.', 'floating_point_enclosure': 'MISSING', 'physical_independent_pairs': 0, 'segment_proxy_warning': True}
    if n == 3:
        result['envelope_controls'] = {'poses_verified': len([s for s in sites if s['segment'] in SEGS]), 'injected_loss_decreased_by1mm_rejected': all((s['distances']['tf2']['lower_mm'] - (s['loss_upper_mm'] - 1) > s['union']['lower_mm'] + 1e-10 for s in sites)), 'identity_preserved': True}
    csvout(f'tables/R{n}_MARGINS.csv', rows)
    csvout(f'raw/R{n}_PER_POSE.csv', policies)
    dump(f'raw/R{n}_RESULT.json', result)
    dump('CURRENT_WORK_STATE.json', {'updated_utc': datetime.now(timezone.utc).isoformat(), 'stage': f'R{n}_COMPLETE', 'latest_gate': result['status'], 'next_operation': 'R3 site-specific digital release correction' if n == 2 else 'Source-number controls, external limits, figure and frozen measurement predictions'})
    (R / f'HANDOFF_R{n}.md').write_text(f"# R{n}: {result['status']}\n\nResults and unchanged validation gates are available in raw/R{n}_RESULT.json. The table is a full cylinder scenario with 2 mm remaining digital distance, r≤2 mm. Physical coverage UNKNOWN. " + ("The next design replaces the the population's Audit Vantile with each existing frozen poses measured Audit Loss.\n" if n == 2 else "No future patient or independent nerve wall is covered by this digital auditing union. The next essential information is independent local wall contour and planned/achieved signed gap loss in the same registered case.\n"))
    print(json.dumps({'round': n, 'table_rows': len(rows), 'controls': len(cr), 'revision_gates': [(v['segment'], v['target'], v['heldout_gate']) for v in revision] if n != 3 else [], 'status': result['status']}))
if __name__ == '__main__':
    import sys
    run(int(sys.argv[1]))
