"""Consume located pairs. Guide transfer remains a labelled population scenario."""
from dental_release.paths import expand as _release_expand
from measure import *
sys.path.insert(0, str(DENT / 'results/LANE_X8_GUIDE_NERVE_RISK'))
from full_geometry import voxel_cylinder_bracket

def classify(lo, hi, t=2):
    return 'BELOW' if hi < t else 'ABOVE' if lo >= t else 'UNKNOWN'

def cylinder_record(points, source):
    r = voxel_cylinder_bracket(points * SP, np.array([SP] * 3), np.array(source['entry_zyx_mm']), np.array(source['axis_zyx']), source['length_mm'], source['radius_mm'])
    return {k: r[k] for k in ('lower_mm', 'upper_mm', 'gap_mm', 'active_voxel_boxes')} | {'witness': r['witness']}

def main():
    start = time.perf_counter()
    patients = json.loads((ROOT / 'raw/PAIRED_PATIENTS.json').read_text())
    bycase = {r['case_tf2']: r for r in patients if r['image_identity']['exact_numeric_equal']}
    certpath = Path(_release_expand('@DENTAL_WORK_ROOT@/X8-guide-nerve-risk/R5_FULL_GEOMETRY.jsonl'))
    certificates = [json.loads(l) for l in certpath.read_text().splitlines()]
    cert_hash = sha_file(certpath)
    groups = {}
    for r in certificates:
        if r['case'] in bycase:
            groups.setdefault(r['case'], []).append(r)
    rows = []
    failed = []
    cases = []
    for (ii, (case, ss)) in enumerate(sorted(groups.items()), 1):
        p = bycase[case]
        artifact = p['point_artifact']
        assert sha_file(artifact['path']) == artifact['sha256']
        cache = ROOT / 'raw/transfer_cases' / f'{case}.json'
        if cache.exists():
            rr = json.loads(cache.read_text())
            rows += rr
            cases.append(case)
            continue
        data = np.load(artifact['path'])
        old = data['old_voxels_zyx']
        new = data['new_voxels_zyx']
        tf = data['tf2_voxels_zyx']
        rr = []
        for s in ss:
            try:
                target = {'lower_mm': s['lower_mm'], 'upper_mm': s['upper_mm'], 'gap_mm': s['gap_mm'], 'witness': s['witness']}
                oldr = target if p.get('tf2_mask_equal_maxillo') else cylinder_record(old, s)
                newr = target if p.get('tf2_mask_equal_tf1') else oldr if p['numeric_mask_equal'] else cylinder_record(new, s)
                unionlo = min(oldr['lower_mm'], newr['lower_mm'], s['lower_mm'])
                unionhi = min(oldr['upper_mm'], newr['upper_mm'], s['upper_mm'])
                witness = target['witness']['witness_zyx_mm']
                (fd, sg, info) = region_map(case, np.array([witness]) / SP)
                names = {0: 'mental_endpoint_5mm_proxy', 1: 'molar_nearest_FDI', 2: 'posterior_endpoint_10mm_proxy', -1: 'UNRESOLVED'}
                r = {'case': case, 'patient': p['patient'], 'fdi': s['fdi'], 'resolution': 'PER_TOOTH', 'time_scale': 'SIMULTANEOUS', 'revised': not p['numeric_mask_equal'], 'old': oldr, 'new': newr, 'tf2': target, 'observed_union': {'lower_mm': unionlo, 'upper_mm': unionhi, 'class_2mm': classify(unionlo, unionhi)}, 'old_class_2mm': classify(oldr['lower_mm'], oldr['upper_mm']), 'new_class_2mm': classify(newr['lower_mm'], newr['upper_mm']), 'tf2_class_2mm': classify(s['lower_mm'], s['upper_mm']), 'old_minus_new_mid_mm': 0.5 * (oldr['lower_mm'] + oldr['upper_mm'] - newr['lower_mm'] - newr['upper_mm']), 'tf2_minus_union_mid_mm': 0.5 * (s['lower_mm'] + s['upper_mm'] - unionlo - unionhi), 'canal_witness_segment': names[int(sg[0])], 'canal_witness_fdi': int(fd[0]), 'clinical_safe_unsafe': 'UNKNOWN', 'pose': {k: s[k] for k in ('entry_zyx_mm', 'axis_zyx', 'length_mm', 'radius_mm')}, 'point_artifact': artifact, 'reference_certificate_sha256': cert_hash}
                rr.append(r)
                rows.append(r)
            except Exception as e:
                failed.append({'case': case, 'fdi': s['fdi'], 'error': repr(e)})
        dump(cache, rr)
        cases.append(case)
        if ii % 5 == 0 or ii == len(groups):
            dump(ROOT / 'raw/PER_SITE_LOCATED.json', rows)
            state('R3_TRANSFERRING', 'R1_INDEPENDENCE_FAIL', 'Complete actual paired-site cylinder distance queries', sites_complete=len(rows), cases_complete=ii)
        print(f'X8 {ii}/{len(groups)} {case} sites={len(rr)}', flush=True)
    segment_names = ['mental_endpoint_5mm_proxy', 'molar_nearest_FDI', 'posterior_endpoint_10mm_proxy']
    margins = []
    x8base = json.loads((DENT / 'results/LANE_X8_GUIDE_NERVE_RISK/RAW_R2.json').read_text())['profiles']
    all_sites = {}
    with (DENT / 'results/LANE_X8_GUIDE_NERVE_RISK/PER_SITE_FULL_GEOMETRY_RISK.csv').open() as f:
        for r in csv.DictReader(f):
            all_sites[r['case'], int(r['fdi'])] = r
    all_gap = np.array([0.5 * (float(r['whole_cylinder_label_gap_lower_mm']) + float(r['whole_cylinder_label_gap_upper_mm'])) for r in all_sites.values()])
    for segment in segment_names:
        pat = [r for p in patients if not p['numeric_mask_equal'] for r in p['regional'] if r['segment'] == segment and r['n'] > 0]
        if not pat:
            continue
        b = float(np.median([r['p95_mm'] for r in pat]))
        local = [r for r in rows if r['canal_witness_segment'] == segment]
        for profile in x8base:
            base = profile['scenario_required_clearance_mm_1pct']
            required = base + b
            margins.append({'segment': segment, 'guide': profile['guide'], 'patients': len(pat), 'revision_budget_mm': b, 'guide_base_scenario_mm': base, 'scenario_margin_mm': required, 'resolution': 'PHENOMENOLOGICAL', 'input_resolution': 'PER_SURFACE_REGION', 'budget_statistic': 'median of revised-patient segment P95', 'aggregate_2581_count_pass_2mm': int((all_gap >= 2).sum()), 'aggregate_2581_count_pass_guide_base': int((all_gap >= base).sum()), 'aggregate_2581_count_pass_guide_plus_revision': int((all_gap >= required).sum()), 'aggregate_2581_guide_only_to_added_budget_changes': int(((all_gap >= base) & (all_gap < required)).sum()), 'aggregate_2581_2mm_to_combined_changes': int(((all_gap >= 2) & (all_gap < required)).sum()), 'matched_local_sites': len(local), 'matched_local_base_to_added_changes': sum((r['tf2']['lower_mm'] >= base and r['observed_union']['upper_mm'] < required for r in local)), 'physical_guarantee': False, 'clinical_validation': 'UNKNOWN', 'replacement_measurement': 'blinded repeated wall contours per patient, plus signed guide-error vectors in a matched placement series', 'limits': 'Transport to all 2581 sites is unvalidated; guide base already contains a Gaussian canal-error closure and the extra budget can double-count it'})
    dump(ROOT / 'raw/GUIDE_SEGMENT_SCENARIOS.json', margins)
    nextrows = [json.loads(l) for l in (DENT / 'results/LANE_NEXT_N_CANAL_PROXIMITY/RAW_N1.jsonl').read_text().splitlines()]
    nxt = []
    with zipfile.ZipFile(TF2) as z:
        nextgroups = {}
        for r in nextrows:
            if r['case'] in bycase:
                nextgroups.setdefault(r['case'], []).append(r)
        for (case, ccrows) in nextgroups.items():
            p = bycase[case]
            dat = np.load(p['point_artifact']['path'])
            sets = {'maxillo': dat['old_voxels_zyx'], 'tf1': dat['new_voxels_zyx'], 'tf2': dat['tf2_voxels_zyx']}
            (lab, sp, h, lh) = mha(z, f'Dataset112_ToothFairy2/labelsTr/{case}.mha')
            (cc, n) = ndi.label(lab == 10)
            for r in ccrows:
                x = np.argwhere(cc == r['component'])
                vals = {}
                witness = {}
                for (key, y) in sets.items():
                    tree = cKDTree(y)
                    (d, k) = tree.query(x, workers=1)
                    best = float(np.linalg.norm(np.maximum(np.abs(x - y[k]) - 1, 0), axis=1).min())
                    rad = best + np.sqrt(3) + 1e-12
                    inds = np.flatnonzero(d <= rad)
                    v = best
                    for j in inds:
                        jj = tree.query_ball_point(x[j], rad)
                        q = np.asarray(jj, dtype=int)
                        ds = np.linalg.norm(np.maximum(np.abs(x[j] - y[q]) - 1, 0), axis=1)
                        a = int(np.argmin(ds))
                        if ds[a] <= v:
                            v = float(ds[a])
                            witness[key] = [x[j].tolist(), y[q[a]].tolist()]
                    vals[key] = v * SP
                nxt.append({'case': case, 'component': r['component'], 'resolution': 'PER_TOOTH', 'revised': not p['numeric_mask_equal'], 'gap_mm': vals, 'observed_union_gap_mm': min(vals.values()), 'reference_next_n_tf2_gap_mm': r['mask_replay_gap_mm'], 'replay_abs_error_mm': abs(vals['tf2'] - r['mask_replay_gap_mm']), 'old_to_new_2mm_change': (vals['maxillo'] >= 2) != (vals['tf1'] >= 2), 'tf2_to_union_2mm_change': (vals['tf2'] >= 2) != (min(vals.values()) >= 2), 'witness_pair_voxels': witness, 'clinical_safe_unsafe': 'UNKNOWN'})
            print('NEXT_N', case, len(ccrows), flush=True)
    dump(ROOT / 'raw/NEXT_N_LOCATED.json', nxt)
    out = {'claim_type': 'information_link', 'x8_total_sites': len(certificates), 'paired_sites': len(rows), 'paired_cases': len(groups), 'paired_revised_sites': sum((r['revised'] for r in rows)), 'old_to_new_2mm_changes': sum((r['old_class_2mm'] != r['new_class_2mm'] for r in rows)), 'tf2_to_observed_union_2mm_changes': sum((r['tf2_class_2mm'] != r['observed_union']['class_2mm'] for r in rows)), 'max_located_old_minus_new_mm': max((r['old_minus_new_mid_mm'] for r in rows), default=None), 'mean_located_abs_change_mm': float(np.mean([abs(r['old_minus_new_mid_mm']) for r in rows])) if rows else None, 'gates': {'actual_paired_2mm_reversal': any((r['old_class_2mm'] != r['new_class_2mm'] for r in rows)), 'all_distance_brackets_closed': all((max(r['old']['gap_mm'], r['new']['gap_mm']) <= 0.0001 for r in rows)), 'clinical_margin_calibrated': False, 'independent_annotator_floor': False}, 'dropout': {'x8_unmatched_sites': len(certificates) - len(rows), 'x8_unmatched_fraction': 1 - len(rows) / len(certificates), 'reason': 'No original dense pair or no verified image identity; no population budget imputation in located analysis', 'numerical_failures': failed, 'next_n_total': len(nextrows), 'next_n_paired': len(nxt), 'next_n_unmatched_fraction': 1 - len(nxt) / len(nextrows)}, 'next_n': {'paired_components': len(nxt), 'old_to_new_2mm_changes': sum((r['old_to_new_2mm_change'] for r in nxt)), 'tf2_to_union_2mm_changes': sum((r['tf2_to_union_2mm_change'] for r in nxt)), 'max_replay_error_mm': max((r['replay_abs_error_mm'] for r in nxt), default=None)}, 'independent_model_test': {'status': 'UNKNOWN_NO_PREDICTION_VOLUME', 'reason': 'X8 manual label consumer; X30 concerns tooth root pulp canals; neither provides independent inferior-alveolar-canal predictions for these cases'}, 'cost': {'wall_seconds': time.perf_counter() - start, 'maxrss_kib': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, 'fit_seconds': 0, 'threads': 1, 'questions': 0}, 'external_referent': {'kind': 'published_dataset', 'locator': str(TF1) + '; ' + str(MAX) + '; ' + str(TF2), 'compared_quantity': 'same-patient canal occupancy and its geometric implant clearance', 'refutes_us': False}}
    dump(ROOT / 'raw/PER_SITE_LOCATED.json', rows)
    dump(ROOT / 'raw/R3_RESULT.json', out)
    state('R3_COMPLETE', out['gates'], 'Run falsifying controls, identical-summary counterexample and one-command report')
    print(json.dumps(out, indent=2), flush=True)
if __name__ == '__main__':
    main()
