"""R3: finite solid capacity over a spatial query domain; conditional label geometry."""
import itertools
import json
from pathlib import Path
import time
import numpy as np
from measure import HERE, ROOT, DATA, Local, check_frozen, dump, sha, trace, restoration_near

def load_local(row):
    suffix = '_mirror' if row['construction'] == 'R2' else ''
    path = DATA / f"{row['case']}_{row['fdi']}{suffix}.npz"
    a = np.load(path)
    obj = Local.__new__(Local)
    obj.lab = a['lab']
    obj.lo = a['lo']
    obj.sp = a['spacing']
    obj.hi = obj.lo + obj.lab.shape
    obj.occ = {k: v.astype(np.float32) for (k, v) in {'bone': obj.lab == 2, 'sinus': np.isin(obj.lab, (5, 6)), 'restoration': np.isin(obj.lab, (8, 9, 10)), 'neighbor': (obj.lab >= 11) & (obj.lab <= 28)}.items()}
    obj.sdf = {}
    obj.signed('bone')
    obj.distance('sinus')
    if obj.distance('neighbor') is None:
        obj.sdf['neighbor_dist'] = np.full(obj.lab.shape, 99, dtype=np.float32)
    return (obj, a['anchor'], a['axis'], a['normal'], a['tangent'])

def body_points(p, axis, normal, tangent, D, L):
    rs = np.linspace(0, D / 2, int(np.ceil(D / 2 / 0.25)) + 1)
    angles = np.arange(64) * 2 * np.pi / 64
    offsets = np.concatenate([np.zeros((1, 3))] + [r * (np.cos(angles)[:, None] * normal + np.sin(angles)[:, None] * tangent) for r in rs[1:]])
    depths = np.unique(np.r_[np.linspace(0, 2, 15), np.linspace(2, L, int(np.ceil((L - 2) / 0.15)) + 1)])
    points = (p + depths[:, None, None] * axis + offsets[None]).reshape(-1, 3)
    return (points, np.repeat(depths, len(offsets)))

def candidate(local, p, axis, normal, tangent, D, L):
    (points, depth) = body_points(p, axis, normal, tangent, D, L)
    q = points / local.sp - local.lo
    inside = ((q >= 1) & (q <= np.asarray(local.lab.shape) - 2)).all()
    if not inside:
        return {'D_mm': D, 'L_mm': L, 'status': 'unknown', 'reason': 'truncated_roi', 'sample_points': len(points)}
    constraints = {'bone': ('bone_sdf', 1.0, depth >= 2 - 1e-09), 'sinus': ('sinus_dist', 1.0, np.ones(len(points), bool)), 'neighbor': ('neighbor_dist', 1.5, np.ones(len(points), bool))}
    (margins, control) = ({}, {})
    for (k, (field, required, mask)) in constraints.items():
        margins[k] = float(local.sample(field, points[mask], order=0, fill=-999).min() - required)
        control[k] = float(local.sample(field, points[mask], order=1, fill=-999).min() - required)
    numeric = float(2 * np.linalg.norm(local.sp / 2) - local.sp.min() / 2 + 0.25)
    minimum = min(margins.values())
    status = 'feasible' if minimum - numeric - 0.3 >= 0 else 'infeasible' if minimum + numeric + 0.3 < 0 else 'unknown'
    errors = {k: abs(margins[k] - control[k]) for k in margins}
    return {'D_mm': D, 'L_mm': L, 'status': status, 'margins_mm': margins, 'trilinear_control_margins_mm': control, 'cross_read_errors_mm': errors, 'control_pass': max(errors.values()) <= 0.3696153, 'numeric_bound_mm': numeric, 'physical_box_assumed_mm': 0.3, 'nominal_fits': minimum >= 0, 'max_edge_displacement_for_conditional_pass_mm': max(0, minimum - numeric), 'limiting_constraint': min(margins, key=margins.get), 'sample_points': len(points)}

def run():
    start = time.monotonic()
    check_frozen('PREREG_R3')
    check_frozen('FROZEN_PREDICTIONS_R3')
    frozen = json.loads((HERE / 'FROZEN_PREDICTIONS_R3.json').read_text())
    for (n, h) in frozen['source_raw_sha256'].items():
        assert sha(HERE / n) == h, n
    rows = [json.loads(s) for s in (HERE / 'COMBINED_SITES.jsonl').read_text().splitlines()]
    catalog = json.loads((ROOT / 'results/DESIGN_implant/catalog.json').read_text())
    pairs = sorted(set(((float(D), float(L)) for (D, L, n) in catalog['catalog_pairs'] if D in (3.5, 4.0, 4.5, 5.0) and L <= 8)), key=lambda x: (x[1], x[0]))
    raw = []
    total_queries = 0
    total_samples = 0
    control_failed = 0
    scalar_refutations = []
    for row in rows:
        out = {'case': row['case'], 'fdi': row['fdi'], 'construction': row['construction'], 'source_valid': row['valid'], 'clinical_decision': 'UNKNOWN', 'pose_domain_is_confidence_region': False, 'positions': []}
        if not row['valid']:
            out.update(status='uncertain', reason=row['reason'])
            raw.append(out)
            continue
        (local, anchor, axis, normal, tangent) = load_local(row)
        for (ob, om) in itertools.product(range(-4, 5), repeat=2):
            a = anchor + ob * normal + om * tangent
            tr = trace(local, a, axis, normal)
            (near, n) = restoration_near(local, a, axis)
            pose = {'offset_b_mm': ob, 'offset_m_mm': om, 'trace': tr, 'candidates': []}
            if not tr['valid'] or near:
                pose.update(status='unknown', reason=tr['reason'] if not tr['valid'] else 'restoration_at_query_pose')
                out['positions'].append(pose)
                continue
            p = np.asarray(tr['crest_mm'])
            for (D, L) in pairs:
                c = candidate(local, p, axis, normal, tangent, D, L)
                pose['candidates'].append(c)
                total_queries += 1
                total_samples += c['sample_points']
                control_failed += not c.get('control_pass', False)
            feasible = [c for c in pose['candidates'] if c['status'] == 'feasible']
            if feasible:
                best = min(feasible, key=lambda c: (c['L_mm'], c['D_mm']))
                pose.update(status='short_geometry', chosen_pair_mm=[best['D_mm'], best['L_mm']], max_edge_displacement_mm=best['max_edge_displacement_for_conditional_pass_mm'])
            elif all((c['status'] == 'infeasible' for c in pose['candidates'])):
                pose.update(status='no_catalog_fit_geometry')
            else:
                pose.update(status='unknown')
            if ob == om == 0:
                main = next((c for c in pose['candidates'] if c['D_mm'] == 4 and c['L_mm'] == 6))
                out['origin_main_candidate'] = main
                out['origin_status'] = pose['status']
                out['origin_chosen_pair_mm'] = pose.get('chosen_pair_mm')
                if row.get('scalar_nominal') == 'short_geometry' and (not main.get('nominal_fits', False)):
                    scalar_refutations.append({'case': row['case'], 'fdi': row['fdi'], 'main_candidate': main})
            out['positions'].append(pose)
        counts = {k: sum((p['status'] == k for p in out['positions'])) for k in ('short_geometry', 'no_catalog_fit_geometry', 'unknown')}
        out['map_counts'] = counts
        out['position_sensitive'] = counts['short_geometry'] > 0 and counts['short_geometry'] < 81
        out['status'] = 'conditional_candidate_region' if counts['short_geometry'] else 'uncertain'
        raw.append(out)
        print(row['case'], row['fdi'], 'map', counts, 'origin', out.get('origin_status', 'unknown'), flush=True)
    with open(HERE / 'RAW_R3_CAPACITY_MAPS.jsonl', 'w') as f:
        for r in raw:
            f.write(json.dumps(r, allow_nan=False) + '\n')
    summary = {'n_combined_sites': len(rows), 'n_source_valid': sum((r['valid'] for r in rows)), 'n_maps_executed': sum((r['source_valid'] for r in raw)), 'n_sites_with_conditional_candidate_region': sum((r['status'] == 'conditional_candidate_region' for r in raw)), 'n_sites_position_sensitive': sum((r.get('position_sensitive', False) for r in raw)), 'origin_counts': {k: sum((r.get('origin_status', 'unknown') == k for r in raw)) for k in ('short_geometry', 'no_catalog_fit_geometry', 'unknown')}, 'scalar_nominal_refutations': scalar_refutations, 'scalar_sufficiency_gate': 'REFUTED' if scalar_refutations else 'NOT_REFUTED_ON_THIS_SAMPLE', 'map_capability_gate': 'PASS' if any((r.get('position_sensitive', False) for r in raw)) else 'FAIL', 'paired_catalog_dimensions_mm': pairs, 'candidate_queries': total_queries, 'point_queries': total_samples, 'lookup_control_failed_candidates': control_failed, 'lookup_control_gate': 'PASS' if control_failed == 0 else 'FAIL', 'elapsed_s': time.monotonic() - start, 'clinical_and_physical_validity': 'UNKNOWN', 'external_distribution_gate': 'UNKNOWN_MINIMUM_N'}
    dump(HERE / 'SUMMARY_R3.json', summary)
    dump(HERE / 'CURRENT_WORK_STATE.json', {'lane': 'X20-short-implant-sinus', 'milestone': 'R3_completed', 'latest_gate': summary, 'next_operation': 'Verify corruptions, publish reproducible figure and minimum local measurement specification'})
    print(json.dumps({k: v for (k, v) in summary.items() if k != 'scalar_nominal_refutations'}, indent=2))
if __name__ == '__main__':
    run()
