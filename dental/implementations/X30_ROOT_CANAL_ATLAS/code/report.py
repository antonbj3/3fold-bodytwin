import csv, datetime, hashlib, json, collections
from pathlib import Path
import numpy as np
from atlas import ROOT, DATA, X12, dump, csvwrite, sha

def read(n):
    return list(csv.DictReader((ROOT / n).open()))

def number(v):
    return float(v) if v not in ['', None] else None

def clustered_ci(rows, field):
    cases = sorted({r['case'] for r in rows})
    rng = np.random.default_rng(30)
    sums = np.array([sum((number(r[field]) >= 2 for r in rows if r['case'] == c)) for c in cases])
    ns = np.array([sum((r['case'] == c for r in rows)) for c in cases])
    idx = rng.integers(0, len(cases), size=(1000, len(cases)))
    values = sums[idx].sum(1) / ns[idx].sum(1)
    return np.quantile(values, [0.025, 0.975]).tolist()

def summarize(rows, name):
    out = {'group': name, 'resolution_level': 'POPULATION', 'n_teeth': len(rows), 'n_CT_cases': len({r['case'] for r in rows}), 'nominal_section_counts': dict(collections.Counter((r['persistent_section_count'] for r in rows))), 'annotation_second_path_fraction': float(np.mean([int(r['persistent_section_count']) >= 2 for r in rows])), 'sampling_cluster_bootstrap95': clustered_ci(rows, 'persistent_section_count')}
    for field in ['max_annotation_path_mm', 'min_window_radius_mm', 'section_aspect_median', 'section_area_median_mm2', 'rootward_clearance_p05_mm']:
        x = [number(r[field]) for r in rows if number(r[field]) is not None]
        out[field + '_P5_P50_P95'] = np.quantile(x, [0.05, 0.5, 0.95]).tolist() if x else None
        out[field + '_observed_n'] = len(x)
    return out

def main():
    rows = read('raw/R1_teeth.csv')
    r2 = json.loads((ROOT / 'raw/R2_result.json').read_text())
    r3 = json.loads((ROOT / 'raw/R3_result.json').read_text())
    controls = json.loads((ROOT / 'raw/CONTROLS.json').read_text())
    r4 = json.loads((ROOT / 'raw/R4_result.json').read_text())
    fdi = [summarize([r for r in rows if int(r['fdi']) == f], str(f)) for f in range(31, 49) if any((int(r['fdi']) == f for r in rows))]
    types = [summarize([r for r in rows if r['tooth_type'] == t], t) for t in dict.fromkeys((r['tooth_type'] for r in rows))]
    dump(ROOT / 'population_by_FDI.json', fdi)
    dump(ROOT / 'population_by_type.json', types)
    facit = json.loads((ROOT / 'FACIT.json').read_text())
    expected = facit['entries'][0]['values']
    compar = []
    for s in fdi:
        if s['group'] not in expected:
            continue
        ref = expected[s['group']]
        p = ref['k'] / ref['n']
        delta = s['annotation_second_path_fraction'] - p
        compar.append({'fdi': int(s['group']), 'n': s['n_teeth'], 'annotation_second_path_fraction': s['annotation_second_path_fraction'], 'sampling_cluster_bootstrap95': s['sampling_cluster_bootstrap95'], 'external_second_canal_fraction': p, 'difference_percentage_points': 100 * delta, 'absolute15pp_gate_pass': abs(delta) <= 0.15, 'locator': 'https://doi.org/10.12659/PJR.901840', 'table_locator': 'Tables1-4', 'resolution_level': 'POPULATION', 'operator_matched': False, 'population_matched': False, 'meaning': 'Population plausibility stress test; failure rejects assuming these labels yield a calibrated canal atlas, not cohort equality or each tooth segmentation'})
    dump(ROOT / 'raw/R1_external_comparison.json', compar)
    cal = json.loads((X12 / 'raw/R8_calibration.json').read_text())
    ex = json.loads((X12 / 'raw/R8_excluded.json').read_text())
    orig = list(csv.DictReader((X12 / 'raw/R8_teeth.csv').open()))
    sources_manifest = json.loads((ROOT / 'sources/UPSTREAM_MANIFEST.json').read_text())
    for item in sources_manifest:
        if sha(item['path']) != item['sha256']:
            raise RuntimeError('Upstream SHA changed')
    artifacts = []
    for p in sorted(DATA.glob('R1_*.jsonl')):
        artifacts.append({'path': str(p), 'bytes': p.stat().st_size, 'sha256': sha(p)})
    for p in sorted((ROOT / 'raw').glob('R*.csv')):
        artifacts.append({'path': str(p), 'bytes': p.stat().st_size, 'sha256': sha(p)})
    for p in sorted((ROOT / 'raw').glob('*result.json')):
        artifacts.append({'path': str(p), 'bytes': p.stat().st_size, 'sha256': sha(p)})
    for p in [ROOT / 'raw/DEMO_P1_FDI36.json', ROOT / 'raw/DEMO_P1_FDI36_points.csv', ROOT / 'raw/R1_manifest.json', ROOT / 'raw/R4_CT_review_manifest.json']:
        artifacts.append({'path': str(p), 'bytes': p.stat().st_size, 'sha256': sha(p)})
    rejected_cases = [r for r in cal if not r['accepted']]
    result = {'lane': 'X30-root-canal', 'claim_type': 'capability', 'timestamp_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'outcome': 'ANNOTATION_ATLAS_DELIVERED; ANATOMICAL_CALIBRATION_REJECTED; CLINICAL_QUANTITIES_UNKNOWN', 'capability_delivered': 'Per-FDI annotation section graph, branch window geometry, shape, clearance and explicit measurement uncertainty/debt', 'external_referent': {'kind': 'independent_measurement', 'locator': 'https://doi.org/10.12659/PJR.901840;Tables1-4', 'compared_quantity': 'prevalence of second canal per FDI31/32/41/42 vs persistent rootward section multiplicity', 'refutes_us': any((not r['absolute15pp_gate_pass'] for r in compar)), 'refuted_claim': 'These annotations can be assumed to yield a calibrated canal-population atlas without independent landmarks and topology validation', 'qualification': 'Different cohorts and operators; neither patient-level false-negative rate nor true Pulpy population prevalence is identified'}, 'coverage': {'advertised_CT': 423, 'exact_CT_pairs': 397, 'semantically_accepted_CT': 283, 'nontruncated_CT': len({r['case'] for r in rows}), 'nontruncated_teeth': len(rows), 'upstream_measured_teeth': len(orig), 'measurement_coverage_of_nontruncated_rows': len(rows) / sum((r['domain_truncated'] == 'False' for r in orig)), 'case_rejection_fraction': len(rejected_cases) / len(cal), 'case_rejection_reasons': dict(collections.Counter((r.get('reason', r.get('calibration', {}).get('reason', 'UNKNOWN')) for r in rejected_cases))), 'tooth_rejection_in_accepted_cases_n': len(ex), 'tooth_rejection_reasons': dict(collections.Counter((r['reason'] for r in ex))), 'edge_truncated_teeth_n': sum((r['domain_truncated'] == 'True' for r in orig)), 'fraction_edge_truncated_of_measured_teeth': sum((r['domain_truncated'] == 'True' for r in orig)) / len(orig), 'physical_tooth_population_denominator': 'UNKNOWN: absent tooth vs absent label cannot be distinguished'}, 'population_by_FDI': fdi, 'population_by_type': types, 'external_population_comparison': compar, 'R2': r2, 'R3': r3, 'R4': r4, 'controls': controls, 'quantity_contract': {'counts': 'Digitized persistent section counts only; clinical canal count null', 'length': 'rootward10%-75% section-graph path, not full canal nor clinical working length', 'angle_radius': 'Three-point3/6/9mm-window turn and circumradius; Schneider/Pruett anatomical measurements null', 'precision_budget': 'Sufficient point-error budget for the conservative interval; not a proven necessary accuracy, scanner voxel size, or physical impossibility result', 'shape': 'axial moment-equivalent ellipse, not tangent-normal histological section', 'clearance': 'two voxel-center boundaries on whole tooth; no DEJ/emalj/dentin split; voxelunion halfwidthsqrt(3)*0.3mm', 'IAN': 'rootward section-graph endpoint to TF2IAN boundary; endpoint at10%cut is not anatomical apex', 'sinus': 'null; Pulpy mandibular only', 'resolution': 'point distances PER_POINT; branch/sections PER_SURFACE_REGION; tooth summaries PER_TOOTH; distributions POPULATION'}, 'missing_quantities': {'working_length': 'needs coronal reference,orifice,apical constriction', 'Schneider_Pruett_landmarks': 'needs expert curve-start/end', 'anatomical_apex_distances': 'needs root apex identification', 'microanatomical_canal_count': 'needs matched same-specimen high-resolution/expert topology', 'clinical_fracture_probability': 'needs alloy,instrument,temperature,load,history and calibrated fatigue law'}, 'source_selection': {'accepted': len(facit['entries']), 'rejected': len(facit['rejected_sources']), 'rejection_fraction': len(facit['rejected_sources']) / (len(facit['entries']) + len(facit['rejected_sources'])), 'reasons': facit['rejected_sources']}, 'cost': {'R1': json.loads((ROOT / 'raw/R1_cost.json').read_text()), 'R2': {k: r2[k] for k in ['wall_s', 'cpu_s', 'peak_RSS_MiB']}, 'R3': {k: r3[k] for k in ['wall_s', 'cpu_s', 'peak_RSS_MiB']}, 'R4': {k: r4[k] for k in ['wall_s', 'cpu_s', 'peak_RSS_MiB']}, 'discovery_preparation_agent_tokens': 'UNKNOWN', 'data_acquisition': 'Existing local datasets; no dataset download', 'validation_fallback': 'X12 costs charged; labelCRC/SHA+distance checks+morphology+interval sampling+injected errors; blind review not yet performed'}, 'phenomenological_debt': [{'quantity': '0.3mm boundary edit scenario', 'resolution_level': 'PHENOMENOLOGICAL', 'replace_with': 'signed mask-boundary errors from expert vs matched microCT'}, {'quantity': 'centroid point epsilon', 'resolution_level': 'PHENOMENOLOGICAL', 'replace_with': 'independent centerline localization error per branch/window'}, {'quantity': 'fatigue risk', 'resolution_level': 'PHENOMENOLOGICAL', 'replace_with': 'instrument-specific same-geometry cyclic-fatigue experiments'}], 'edges': [{'from': 'Pulpy+TF2 exact CT mask frame', 'to': 'section graph and local clearance', 'resolution_level': 'PER_POINT', 'timescale': 'SIMULTANEOUS', 'status': 'PENDING_INDEPENDENT_REVIEW'}, {'from': 'branch window radius/turn/uncertainty', 'to': 'NiTi fatigue geometry consumer', 'resolution_level': 'PER_SURFACE_REGION', 'timescale': 'SIMULTANEOUS', 'status': 'QUALITATIVE_INFORMATION_LINK_ONLY'}, {'from': 'instrument cyclic damage state', 'to': 'next-use instrument state', 'resolution_level': 'PER_SURFACE_REGION', 'timescale': 'HANDOVER', 'status': 'UNKNOWN_NOT_IMPLEMENTED'}], 'numerical_artifact_manifest': artifacts, 'scientific_admission': False, 'review_state': 'PENDING_INDEPENDENT_REVIEW'}
    for item in json.loads((ROOT / 'raw/R4_CT_review_manifest.json').read_text()):
        result['numerical_artifact_manifest'].append({'path': item['path'], 'bytes': item['bytes'], 'sha256': item['sha256']})
    slots = len(orig) + len(ex)
    assert slots == 16 * sum((bool(r['accepted']) for r in cal))
    result['coverage'].update(candidate_FDI_slots_in_accepted_283_cases=slots, tooth_rejection_fraction_of_candidate_FDI_slots=len(ex) / slots, FDI_slot_denominator_note='includes missing-mask slots,not physical tooth census')
    dump(ROOT / 'results.json', result)
    (ROOT / 'results.sha256').write_text(sha(ROOT / 'results.json') + '  results.json\n')
    print(result['outcome'])
    print('teeth', len(rows), 'cases', result['coverage']['nontruncated_CT'])
    print(compar)
if __name__ == '__main__':
    main()
