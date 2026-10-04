"""Frozen lane measurements; local derived outputs only."""
from dental_release.paths import expand as _release_expand
import csv, json, time, resource
from collections import Counter
from pathlib import Path
import numpy as np
from model import *

def csv_write(path, rows):
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

def public_row(case, jaw, key, answer):
    return dict(case=case, jaw=jaw, tooth_or_region=key, classification=answer['classification'], force_lower_N=float(answer['force_hull_N'][0]), force_upper_N=float(answer['force_hull_N'][1]), share_lower_pp=float(answer['share_hull_pp'][0]), share_upper_pp=float(answer['share_hull_pp'][1]), lower_is_attained=answer['infimum_attained'], necessity_margin_exact_mm=str(answer['necessity_margin_mm']), license=answer['license'], physical_classification='UNCERTAIN')

def references():
    source = json.load(open(X2 / 'raw/EXTERNAL_REFERENTS.json'))
    (ferr, hatt) = source['studies'][:2]
    upper = {int(t): Q(str(v)) for (t, v) in zip(ferr['fdi'], ferr['mean_percent'])}
    lower = {10 * q + i + 1: Q(str(v)) for q in (3, 4) for (i, v) in enumerate(hatt['single_side_mean_percent'])}
    return (dict(upper=upper, lower=lower), source)

def external_metrics(results):
    (refs, src) = references()
    out = []
    for (jaw, ref) in refs.items():
        for level in ('tooth', 'region8', 'region4'):
            if level == 'tooth':
                target = {str(t): v for (t, v) in ref.items()}
                intervals = [{str(t): a['teeth'][str(t)]['share_hull_pp'] for t in ref} for a in results]
            elif level == 'region8':
                target = {r: sum((v for (t, v) in ref.items() if region(t) == r)) for r in REGIONS}
                intervals = [{r: a['regions'][jaw][r]['share_hull_pp'] for r in REGIONS} for a in results]
            else:
                target = {r: sum((v for (t, v) in ref.items() if region(t).startswith(r + '_'))) for r in ('I', 'C', 'P', 'M')}
                intervals = []
                for a in results:
                    (e, _) = edges_from_x21(a['case'])
                    intervals.append({r: classify(e, [t for t in FDI if (t < 30) == (jaw == 'upper') and region(t).startswith(r + '_')])['share_hull_pp'] for r in target})
            keys = list(target)
            included = [Q(b[k][0]) <= target[k] <= Q(b[k][1]) for b in intervals for k in keys]
            overlap = []
            for b in intervals:
                for k in keys:
                    nt = 1 if level == 'tooth' else sum((1 for t in ref if (region(t) == k if level == 'region8' else region(t).startswith(k + '_'))))
                    tol = Q(nt) if jaw == 'lower' else Q(0)
                    overlap.append(Q(b[k][0]) <= target[k] + tol and Q(b[k][1]) >= target[k] - tol)
            widths = [float(Q(b[k][1]) - Q(b[k][0])) for b in intervals for k in keys]
            means = {k: [sum((Q(b[k][side]) for b in intervals)) / len(intervals) for side in (0, 1)] for k in keys}
            mean_hits = {k: means[k][0] <= target[k] <= means[k][1] for k in keys}
            out.append(dict(jaw=jaw, level=level, n_cases=len(results), n_comparisons=len(included), marginal_covered=sum(included), marginal_coverage=float(np.mean(included)), marginal_digitization_overlap=float(np.mean(overlap)), mean_width_pp=float(np.mean(widths)), median_width_pp=float(np.median(widths)), reference_pp=target, reference_sum_pp=sum(target.values()), cohort_mean_intervals_pp=means, cohort_mean_coverage=sum(mean_hits.values()) / len(mean_hits), cohort_mean_width_pp=float(np.mean([float(v[1] - v[0]) for v in means.values()])), full_vector_normalization='INCOMPATIBLE_93_PERCENT' if jaw == 'upper' else '100_PERCENT; inclusion in marginal hull is not joint mechanical feasibility', matched_subject_validation='NOT_AVAILABLE', individual_transfer='UNKNOWN'))
    rows = [json.load(open(p)) for p in sorted((X2 / 'raw/cases_3d').glob('*.json'))]
    practice = []
    for (jaw, offset) in [('upper', 0), ('lower', 8)]:
        ref = np.array([float(sum((v for (t, v) in refs[jaw].items() if region(t).startswith(r + '_')))) for r in ('I', 'C', 'P', 'M')])
        mech = []
        area = []
        for r in rows:
            m = r['metrics'][1]
            if m['mechanics'] is None:
                continue

            def aggregate(v):
                v = np.asarray(v)[offset:offset + 8]
                return v[:4] + v[4:]
            mech.append(aggregate(m['mechanics']['shares_pp']))
            area.append(aggregate(m['area_shares_pp']))
        mm = np.mean(mech, axis=0)
        aa = np.mean(area, axis=0)
        practice.append(dict(jaw=jaw, n=len(mech), reference_pp=ref.tolist(), point_pp=mm.tolist(), area_pp=aa.tolist(), point_MAE_pp=float(np.abs(mm - ref).mean()), area_MAE_pp=float(np.abs(aa - ref).mean()), point_exact_coverage=float(np.mean(mm == ref)), area_exact_coverage=float(np.mean(aa == ref)), point_interval_width_pp=0, area_interval_width_pp=0, notes='Original X2 Bits2Bites cohort; no X21 subject linkage. Singleton intervals have no robustness guarantee.'))
    return dict(descriptive_external_test=out, practice_X2_recomputed=practice, calibration_to_facit=False, patient_validation='UNKNOWN', dominance_gate=dict(selector_built=False, reason='capability, no competing sound routes', conservation_hull_coverage_ceiling=1.0, conservation_hull_width_pp=100, point_area_with_common_conservation_hull='IDENTICAL_CEILING; no selector can narrow epistemic bounds'), source_sha256=sha(X2 / 'raw/EXTERNAL_REFERENTS.json'))

def x54_edges(case, tag):
    path = X54 / 'raw' / f'geometry_{case:03d}_{tag}.json'
    r = json.load(open(path))
    d = {}
    for p in r['patches']:
        key = (p['upper_fdi'], p['lower_fdi'])
        v = Q(p['vertical_gap_mm'])
        d[key] = min(d.get(key, v), v)
    return ([(u, l, g, g) for ((u, l), g) in sorted(d.items())], path)

def resolution(results, panel):
    changed = []
    comparisons = []
    for case in panel:
        (old, p) = edges_from_x21(case, False)
        (new, q) = edges_from_x21(case)
        for t in FDI:
            (a, b) = (classify(old, [t]), classify(new, [t]))
            if a['classification'] != b['classification']:
                changed.append(dict(case=case, fdi=t, R3=a['classification'], R4=b['classification']))
    for case in (1, 2, 3):
        models = {tag: x54_edges(case, tag) for tag in ('h02', 'h01', 'continuous', 'continuous_union')}
        classes = {tag: {t: classify(e, [t])['classification'] for t in FDI} for (tag, (e, p)) in models.items()}
        changes = {f'{a}_to_{b}': [t for t in FDI if classes[a][t] != classes[b][t]] for (a, b) in [('h02', 'h01'), ('h02', 'continuous'), ('continuous', 'continuous_union')]}
        comparisons.append(dict(case=case, classes=classes, class_flips=changes, pair_counts={tag: len(e) for (tag, (e, p)) in models.items()}, source_sha256={str(p): sha(p) for (e, p) in models.values()}, scope='X54 gaps projected onto the NEW AXIAL NORMAL MODEL; not X54 Coulomb output or X21 cases'))
    return dict(X21_continuous_gap_h_independent=True, X21_R3_to_R4_class_changes=changed, X21_refinement_is_more_information_not_h_convergence=True, X54_sampled_transfer=comparisons, X54_h_gate='PASS' if all((not c['class_flips']['h02_to_h01'] for c in comparisons)) else 'FAIL', limitations='X21 h-independence applies only to fixed continuous declared pair list; does not certify contact completeness.')

def main():
    stopped()
    start = time.perf_counter()
    cpu = time.process_time()
    assert sha(ROOT / 'PREREG.md') == (ROOT / 'PREREG.sha256').read_text().split()[0]
    manifest = json.load(open(X21 / 'raw/INPUT_MANIFEST.json'))
    results = []
    counts = Counter()
    toothrows = []
    regionrows = []
    sources = {str(p): sha(p) for p in [X21 / 'RESULTS.md', X21 / 'raw/INPUT_MANIFEST.json', X2 / 'RESULTS.md', X2 / 'raw/EXTERNAL_REFERENTS.json', X54 / 'RESULTS.md', SUPPORT]}
    for case in manifest['cases']:
        stopped()
        (e, p) = edges_from_x21(case)
        sourcecase = X21 / 'raw/cases' / (case + '.json')
        sources[str(p)] = sha(p)
        sources[str(sourcecase)] = sha(sourcecase)
        s = json.load(open(sourcecase))
        teeth = {str(t): classify(e, [t]) for t in FDI}
        regions = {jaw: {r: classify(e, [t for t in FDI if (t < 30) == (jaw == 'upper') and region(t) == r]) for r in REGIONS} for jaw in ['upper', 'lower']}
        a = dict(case=case, source_gap_path=str(p), source_gap_sha256=sha(p), n_pairs=len(e), gap_bound_version='R4' if 'refined' in str(p) else 'R3', maximum_gap_width_mm=float(max((x[3] - x[2] for x in e))), pose_overlap_below_minus03=any((x[3] < Q(-0.03) for x in e)), minimum_labelled_face_fraction=min((s['stats'][j]['labelled_face_fraction'] for j in ['upper', 'lower'])), teeth=teeth, regions=regions, physical_status='UNKNOWN; conditional graph only')
        results.append(a)
        for t in FDI:
            counts[teeth[str(t)]['classification']] += 1
            toothrows.append(public_row(case, 'upper' if t < 30 else 'lower', t, teeth[str(t)]))
        for jaw in regions:
            for (r, b) in regions[jaw].items():
                regionrows.append(public_row(case, jaw, r, b))
    query_s = time.perf_counter() - start
    csv_write(ROOT / 'raw/TOOTH_CLASSES.csv', toothrows)
    csv_write(ROOT / 'raw/REGION_INTERVALS.csv', regionrows)
    write(ROOT / 'raw/CASE_CERTIFICATES.json', results)
    srows = json.load(open(SUPPORT))
    unknown = [r['id'] for r in srows if r['evidence'] == 'UNKNOWN']
    write(ROOT / 'raw/SUPPORT_AUDIT.json', dict(source_path=str(SUPPORT), source_sha256=sha(SUPPORT), source_count=len(srows), unknown_count=len(unknown), unknown_ids=unknown, lane_declared_unknown_count=8, source_explicit_unknown_count=len(unknown), unknown_count_contract='DISCREPANCY: frozen lane says 8; local source has 7 explicitly UNKNOWN records', unidentified_eighth_record='UNKNOWN; no invented source ID or stiffness', preserved_source_ids=[r['id'] for r in srows], per_tooth_axial_interval=dict(unit='N/mm', lower=0, lower_open=True, upper=None, upper_kind='UNBOUNDED', status='UNKNOWN_POSITIVE_MODEL_ASSUMPTION'), incisor_prior=dict(interval_N_per_mm=[500, 1130], status='DERIVED_HYPOTHESIS_NOT_USED_AS_AXIAL_CERTIFICATE'), preserved=True, clinical_bounded_axial_stiffness_available=False))
    write(ROOT / 'raw/INPUT_CONTRACT_DISCREPANCIES.json', dict(dataset_link=dict(lane_impression='X21 to X2 force reference', observed='Bite2Text vs Bits2Bites', matched_case_count=0), support_count=dict(lane_and_frozen_prereg=8, source_explicit_unknown_count=len(unknown), source_total=len(srows), source_sha256=sha(SUPPORT), missing_record_id=None, status='UNRESOLVED_SOURCE_COUNT; all axial inputs remain unknown'), scientific_policy='No frozen criterion changed; preserve discrepant expectations and actual source evidence.'))
    st = time.perf_counter()
    ext = external_metrics(results)
    external_s = time.perf_counter() - st
    write(ROOT / 'raw/EXTERNAL_COMPARISON.json', ext)
    st = time.perf_counter()
    res = resolution(results, manifest['numerical_panel'])
    resolution_s = time.perf_counter() - st
    write(ROOT / 'raw/RESOLUTION.json', res)
    for c in res['X54_sampled_transfer']:
        sources.update(c['source_sha256'])
    for p in sorted((X2 / 'raw/cases_3d').glob('*.json')):
        sources[str(p)] = sha(p)
    for case in manifest['numerical_panel']:
        p = X21 / 'raw/certificates' / (case + '.json')
        sources[str(p)] = sha(p)
    sources.update({str(p): sha(p) for p in (ROOT / 'engine_snapshot').rglob('*.py')})
    write(ROOT / 'raw/SOURCE_MANIFEST.json', dict(sha256=sources, engine_commit='a357dcc', source_scope='hashes/derived quantities only; no mesh or raw array redistribution', study_links=dict(X21='Bite2Text', X2='Bits2Bites', X54='Bits2Bites', matched_X21_X2_cases=0), licenses=dict(Bits2Bites='CC BY-NC-SA 4.0, local read only', Bite2Text='UNKNOWN exact terms, local read only', source_label_provenance='Teeth3DS local binding UNKNOWN; upstream CC BY-NC-ND')))
    with open(ROOT / 'CASE_TABLES.md', 'w') as f:
        f.write('# One table per X21 case\n\nAll classes apply to the axial model and declared the pair list.\nPhysical patient class is UNCERTAIN for all. A FDI without a candidate is\nNever loaded in this model; physical absence or zeroest is unknown .\nShares are closed infimum/superremum casings. MUST with lower zero can\nhave strictly positive force everywhere without positive common minimum.\n\n')
        labels = {'MUST': 'must:', 'CAN': 'may', 'NEVER': 'never'}
        for a in results:
            f.write(f"## {a['case']} ({a['gap_bound_version']}, {a['n_pairs']} par)\n\n")
            f.write('| FDI | Klass i modellen | Force N | Andel % | Physical Class |\n|---:|---|---|---|---|\n')
            for (t, b) in a['teeth'].items():
                fn = '–'.join((str(float(Q(x))).removesuffix('.0') for x in b['force_hull_N']))
                pp = '–'.join((str(float(Q(x))).removesuffix('.0') for x in b['share_hull_pp']))
                f.write(f"| {t} | {labels[b['classification']]} | {fn} | {pp} | UNCERTAIN |\n")
            f.write('\n')
    summary = dict(lane=_release_expand('FALT_TANDLAST'), claim_type='capability', n_cases=len(results), n_tooth_slots=len(toothrows), n_region_intervals=len(regionrows), counts=dict(counts), cases_with_must=sum((any((b['classification'] == 'MUST' for b in a['teeth'].values())) for a in results)), cases_without_must=sum((not any((b['classification'] == 'MUST' for b in a['teeth'].values())) for a in results)), physical_classifications_certified=0, external_validation='UNKNOWN_UNMATCHED', prereg_sha256=sha(ROOT / 'PREREG.md'), resolution_gate=res['X54_h_gate'], wall_s=time.perf_counter() - start, cpu_s=time.process_time() - cpu, classification_plus_source_audit_s=query_s, external_s=external_s, resolution_s=resolution_s, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024)
    write(ROOT / 'raw/SUMMARY.json', summary)
    write(ROOT / 'COSTS.json', dict(measured_run=summary, pilot=json.load(open(ROOT / 'raw/PILOT_VS_CONTROL.json'))['point']['total_s'], setup='engine archive/read and human discovery not separately timed', fit='NOT_RUN', discovery_human_s='UNKNOWN', upstream_segmentation_geometry_cost='UNKNOWN_IN_THIS_LANE; source reports retain their costs', physical_validation='NOT_RUN', update='NOT_RUN', future_queries='NOT_RUN', fallback='UNCERTAIN, no physical solver run', acceleration='CPU only, 4-thread caps; no cloud/GPU', speedup_claim=None))
    print(json.dumps(summary, indent=2))
if __name__ == '__main__':
    main()
