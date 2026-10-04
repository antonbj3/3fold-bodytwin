from common import *
import math, resource
from scipy.stats import beta
SEG = ['mental_endpoint_5mm_proxy', 'molar_nearest_FDI', 'posterior_endpoint_10mm_proxy']

def number(x):
    return float('inf') if x == 'INF' else float(x)

def finite(x):
    return float(x) if math.isfinite(x) else 'INF'

def quantile95(a):
    n = len(a)
    rank = math.ceil((n + 1) * 0.95)
    return (sorted(a)[rank - 1] if rank <= n else float('inf'), rank)

def ci(k, n):
    if not n:
        return [None, None]
    return [0.0 if not k else float(beta.ppf(0.025, k, n - k + 1)), 1.0 if k == n else float(beta.ppf(0.975, k + 1, n - k))]

def coverage(a, q):
    k = sum((x <= q for x in a))
    n = len(a)
    return dict(covered_patients=k, test_patients=n, coverage=k / n if n else None, exact_binomial95_ci=ci(k, n), interpretation='Descriptive patient coverage; CI assumes independent patient draws; exchangeability unverified', vacuous=not math.isfinite(q))

def guide_budget(p):
    a = 0.05 / 3
    k = math.sqrt((1 - a) / a)
    entry = p['entry_mean_sd_mm'][0] + p['entry_mean_sd_mm'][1] * k
    apex = p['apex_mean_sd_mm'][0] + p['apex_mean_sd_mm'][1] * k
    theta = math.radians(p['angle_mean_sd_deg'][0] + p['angle_mean_sd_deg'][1] * k)
    angular = 4 * math.sin(min(theta, math.pi) / 2)
    return dict(guide=p['id'], entry_bound_mm=entry, apex_bound_mm=apex, angle_bound_deg=math.degrees(theta), radius_mm=2.0, angular_bound_mm=angular, hausdorff_budget_mm=max(entry, apex) + angular, allocated_failure_probability_per_component=a, guide_conditional_coverage=0.95, resolution='PHENOMENOLOGICAL', source=p['source'], sample_n_implants=p['n_implants'], sample_to_population_moments='UNKNOWN; conditional closure only', replacement_measurement='Same intended/achieved whole implant, independent signed entry/apex/angular vectors by guide and matched arch/segment; confidence bounds on population moments')

def main():
    check_frozen()
    start = time.perf_counter()
    patients = json.loads((ROOT / 'raw/PATIENTS.json').read_text())
    sites = json.loads((ROOT / 'raw/LINEAGE_SITES.json').read_text())
    lineage = json.loads((ROOT / 'raw/IMAGE_GROUPS.json').read_text())
    canonical = lineage['canonical']
    ids = sorted(set(canonical.values()), key=lambda k: int(k[1:]))
    ordered = sorted(ids, key=lambda k: hashlib.sha256(('X73|' + k).encode()).hexdigest())
    cut = math.floor(0.6 * len(ordered))
    cal = set(ordered[:cut])
    test = set(ordered[cut:])
    assert not cal & test
    for r in sites:
        r['patient'] = canonical[r['patient']]
    split = dict(rule='R4 canonical exact-image-hash group; sha256("X73|"+canonical) ordered, floor(0.6*Ngroups) calibration; all aliases together', calibration_patients=sorted(cal), test_patients=sorted(test), patient_overlap=0, calibration_n=len(cal), test_n=len(test), dataset_id_to_image_group=canonical, calibration_dataset_ids=sorted((p for (p, g) in canonical.items() if g in cal)), test_dataset_ids=sorted((p for (p, g) in canonical.items() if g in test)))
    dump(ROOT / 'raw/PATIENT_SPLIT.json', split)
    radial = []
    query = []
    score_rows = []
    for seg in SEG:
        rr = {p['patient']: next((r for r in p['regions'] if r['segment'] == seg), None) for p in patients}
        scores = {p: max((number(r['strict_score_mm']) if r and r['all_rays'] else float('inf') for (alias, r) in rr.items() if canonical[alias] == p)) for p in ids}
        (q, rank) = quantile95([scores[p] for p in cal])
        eval = coverage([scores[p] for p in test], q)
        group_p95 = {p: max((r['p95_mm'] for (alias, r) in rr.items() if canonical[alias] == p and r and r['n'])) for p in ids if any((canonical[alias] == p and r and r['n'] for (alias, r) in rr.items()))}
        radial.append(dict(segment=seg, resolution='POPULATION', raw_resolution='PER_POINT', patient_median_of_supported_p95_mm=float(np.median(list(group_p95.values()))) if group_p95 else None, duplicate_alias_summary='maximum supported P95 within exact-image group before patient median', supported_patient_count=len(group_p95), strict_calibration_n=len(cal), strict_test_n=len(test), strict_calibration_finite=sum((math.isfinite(scores[p]) for p in cal)), strict_test_finite=sum((math.isfinite(scores[p]) for p in test)), strict_quantile95_mm=finite(q), rank=rank, evaluation=eval, gate_finite_95percent_margin=math.isfinite(q) and eval['coverage'] >= 0.95))
        grouped = {}
        for s in sites:
            if s['segment'] == seg:
                grouped.setdefault(s['patient'], []).append(s)
        loss = {p: max((r['loss_upper_mm'] for r in ss)) for (p, ss) in grouped.items()}
        ca = [loss[p] for p in cal if p in loss]
        te = [loss[p] for p in test if p in loss]
        (qq, qr) = quantile95(ca)
        qe = coverage(te, qq)
        query.append(dict(segment=seg, resolution='POPULATION', raw_resolution='PER_TOOTH', calibration_patients=len(ca), heldout_patients=len(te), calibration_sites=sum((len(grouped[p]) for p in cal if p in grouped)), heldout_sites=sum((len(grouped[p]) for p in test if p in grouped)), wall_loss_budget95_mm=finite(qq), rank=qr, evaluation=qe, finite_empirical_95_gate=math.isfinite(qq) and qe['coverage'] is not None and (qe['coverage'] >= 0.95), coverage_scope='Patient maximum clearance loss among frozen X8 poses whose TF2 nearest-wall witness maps to this segment; no claim for all other poses or true anatomical walls'))
        for p in ids:
            score_rows.append(dict(patient=p, segment=seg, split='CALIBRATION' if p in cal else 'HELDOUT', strict_radial_score_mm=finite(scores[p]), query_loss_score_mm=loss.get(p), query_sites=len(grouped.get(p, [])), resolution='PER_TOOTH'))
    guides = [guide_budget(p) for p in json.loads((X8 / 'GUIDE_PROFILES.json').read_text())['profiles']]
    margins = []
    policies = []
    for q in query:
        seg = q['segment']
        b = number(q['wall_loss_budget95_mm'])
        ss = [s for s in sites if s['patient'] in test and s['segment'] == seg]
        for g in guides:
            m = 2 + b + g['hausdorff_budget_mm']
            accepted = [s for s in ss if s['classes']['tf2'] == 'ABOVE']
            changes = sum((s['distances']['tf2']['upper_mm'] < m for s in accepted)) if math.isfinite(m) else len(accepted)
            margins.append(dict(segment=seg, guide=g['guide'], wall_loss_budget95_mm=finite(b), guide_budget95_mm=g['hausdorff_budget_mm'], conditional_margin_for_retaining_2mm_mm=finite(m), heldout_sites=len(ss), fixed_2mm_accepted=len(accepted), fixed_2mm_to_conditional_policy_changes=changes, changes_per100_paired_sites=100 * changes / len(ss) if ss else None, combined_target_coverage=0.9, empirical_wall_heldout_coverage=q['evaluation']['coverage'], coverage_claim='CONDITIONAL only: patient exchangeability + exact population guide moments + r<=2mm. Physical independent-annotation coverage UNKNOWN.', resolution='PHENOMENOLOGICAL', edge_resolution='PER_TOOTH', time_scale='SIMULTANEOUS', finite=math.isfinite(m)))
            for s in ss:
                d = s['distances']['tf2']
                c = classify(d['lower_mm'], d['upper_mm'], m) if math.isfinite(m) else 'UNKNOWN_INFINITE_MARGIN'
                policies.append(dict(patient=s['patient'], fdi=s['fdi'], segment=seg, guide=g['guide'], nominal_lower_mm=d['lower_mm'], nominal_upper_mm=d['upper_mm'], conditional_margin_mm=finite(m), fixed_class=s['classes']['tf2'], conditional_class=c, resolution='PER_TOOTH', time_scale='SIMULTANEOUS'))
    direct = []
    for (label, ss) in [('ALL_PAIRED', sites), ('CALIBRATION', [s for s in sites if s['patient'] in cal]), ('HELDOUT', [s for s in sites if s['patient'] in test])]:
        for seg in ['ALL'] + SEG + ['UNRESOLVED']:
            dd = ss if seg == 'ALL' else [s for s in ss if s['segment'] == seg]
            changes = [s for s in dd if s['classes']['tf2'] == 'ABOVE' and s['union']['class_2mm'] == 'BELOW']
            unknown = sum((s['classes']['tf2'] == 'UNKNOWN' or s['union']['class_2mm'] == 'UNKNOWN' for s in dd))
            direct.append(dict(split=label, segment=seg, sites=len(dd), patients=len(set((s['patient'] for s in dd))), fixed_2mm_to_observed_union_reversals=len(changes), reversals_per100_paired_sites=100 * len(changes) / len(dd) if dd else None, ambiguous=unknown, resolution='POPULATION', raw_resolution='PER_TOOTH', time_scale='SIMULTANEOUS', not_clinical_prevalence=True))
    oldnew = [s for s in sites if 'maxillo' in s['classes'] and 'tf1' in s['classes']]
    unique = []
    for (label, ss) in [('ALL_PAIRED', sites), ('HELDOUT', [s for s in sites if s['patient'] in test])]:
        gg = {}
        for s in ss:
            gg.setdefault((s['patient'], s['fdi']), []).append(s)
        flips = sum((any((s['classes']['tf2'] == 'ABOVE' and s['union']['class_2mm'] == 'BELOW' for s in rr)) for rr in gg.values()))
        unique.append(dict(split=label, unique_lineage_FDI_sites=len(gg), annotation_pose_rows=len(ss), affected_sites_any_frozen_pose=flips, affected_per100_unique_sites=100 * flips / len(gg) if gg else None, resolution='POPULATION', consumer_resolution='PER_TOOTH', aggregation='Same image group and FDI counted once; if any retained frozen pose reverses, consumer site affected; pose alternatives remain in raw outputs', not_clinical_prevalence=True))
    reversals = [dict(patient=s['patient'], case=s['case'], fdi=s['fdi'], segment=s['segment'], tf2_lower_mm=s['distances']['tf2']['lower_mm'], union_upper_mm=s['union']['upper_mm'], loss_upper_mm=s['loss_upper_mm'], split='CALIBRATION' if s['patient'] in cal else 'HELDOUT', resolution='PER_TOOTH') for s in sites if s['classes']['tf2'] == 'ABOVE' and s['union']['class_2mm'] == 'BELOW']
    result = dict(claim_type='information_link', construction='R4_IMAGE_LINEAGE', image_lineage=lineage, radial_segments=radial, query_segments=query, guide_budgets=guides, guide_segment_margins=margins, direct_decisions=direct, old_new=dict(sites=len(oldnew), reversals=sum((s['classes']['maxillo'] != s['classes']['tf1'] for s in oldnew))), reversal_witnesses=reversals, patient_split=split, dropout=dict(heldout_sites_total=sum((s['patient'] in test for s in sites)), heldout_segment_unresolved=sum((s['patient'] in test and s['segment'] == 'UNRESOLVED' for s in sites)), segment_missing_patients=[dict(segment=q['segment'], no_query_cal=len(cal) - q['calibration_patients'], no_query_test=len(test) - q['heldout_patients']) for q in query]), physical_independent_margin='UNKNOWN', cost=dict(wall_seconds=time.perf_counter() - start, maxrss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, fit_included=True, questions=0), formal_floating_point_enclosure='MISSING; no affine sensitivity claim')
    dump(ROOT / 'raw/CALIBRATION_RESULT.json', result)
    writecsv(ROOT / 'raw/PATIENT_SCORES.csv', score_rows)
    writecsv(ROOT / 'raw/DECISION_REVERSALS.csv', reversals)
    writecsv(ROOT / 'raw/GUIDE_SEGMENT_MARGINS.csv', margins)
    writecsv(ROOT / 'raw/HELDOUT_POLICIES.csv', policies)
    result['unique_tooth_consumers'] = unique
    dump(ROOT / 'raw/CALIBRATION_RESULT.json', result)
    (ROOT / 'HANDOFF_R2.md').write_text("# R2 : coverage must follow the wall dimension\n\nAll candidate patients and support -/ directional states are present in raw/PATIENTS.json and hashade puncture arrays. Frozen strict patient maximum gate retains each missing beam as INF; its outcome is found in radial_segments in raw/CALIBRATION_RESULT.json.\n\nR3 changes the operation: the targeted distance loss of the issue to the Union of actually published masks instead of a placeless wall - P95 . The guide coupling replaces no unknown anatomy errors.\n")
    state('R3_CALIBRATION_COMPLETE', dict(finite_strict_radial_segments=sum((r['gate_finite_95percent_margin'] for r in radial)), finite_query_budgets=sum((q['wall_loss_budget95_mm'] != 'INF' for q in query)), independent_margin=False), 'Run independent controls and injected corruption probes')
    print(json.dumps(dict(radial=radial, query=query, direct=[r for r in direct if r['segment'] == 'ALL']), indent=2))
if __name__ == '__main__':
    main()
