"""Replay frozen decision evidence without running or mutating original lanes."""
import os
for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[name] = '1'
from pathlib import Path
import collections
import csv
import datetime
import hashlib
import itertools
import json
import math
import platform
import resource
import subprocess
import sys
import time
import xml.etree.ElementTree as ET
import numpy as np
import scipy
from scipy import ndimage as ndi
from scipy.spatial import cKDTree
from voxel_decisions import local_cube_field, direct_cube_field, hd95_points
ROOT = Path(__file__).resolve().parent
CHECKS = []
LEDGER = []

def read(p):
    return json.loads((ROOT / p).read_text())

def save(p, d):
    (ROOT / p).parent.mkdir(parents=True, exist_ok=True)
    (ROOT / p).write_text(json.dumps(d, indent=2, ensure_ascii=False, allow_nan=False) + '\n')

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def check(name, observed, expected, tolerance=0, bad=None):
    """Each claim control must reject an explicitly injected wrong output."""
    if isinstance(expected, (int, float)) and (not isinstance(expected, bool)):
        predicate = lambda x: isinstance(x, (int, float)) and math.isfinite(x) and (abs(x - expected) <= tolerance)
        if bad is None:
            bad = expected + max(1.0, 100 * tolerance)
    else:
        predicate = lambda x: x == expected
        if bad is None:
            bad = not expected if isinstance(expected, bool) else 'INJECTED_WRONG_VALUE'
    row = {'gate': name, 'observed': observed, 'expected': expected, 'tolerance': tolerance, 'pass': bool(predicate(observed)), 'injected_value': bad, 'injected_rejected': bool(not predicate(bad))}
    CHECKS.append(row)
    if not row['pass'] or not row['injected_rejected']:
        save('raw/FAILED_CONTROLS.json', CHECKS)
        raise ValueError(f'Gate failed: {row}')

def claim(name, value, locator, level, unit, kind='DERIVED'):
    LEDGER.append({'claim': name, 'value': value, 'source_locator': locator, 'resolution_level': level, 'unit': unit, 'evidence_type': kind, 'time_scale': 'SIMULTANEOUS'})

def table(name, rows):
    p = ROOT / 'tables' / name
    p.parent.mkdir(exist_ok=True)
    with p.with_suffix('.csv').open('w') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    keys = list(rows[0])
    lines = ['| ' + ' | '.join(keys) + ' |', '| ' + ' | '.join(['---'] * len(keys)) + ' |']
    lines.extend(('| ' + ' | '.join((str(r[k]) for k in keys)) + ' |' for r in rows))
    p.with_suffix('.md').write_text('\n'.join(lines) + '\n')

def attached_to_source(source, added):
    """Every connected added component must reach an ORIGINAL source voxel."""
    original = {tuple(x) for x in source}
    remaining = {tuple(x) for x in added}
    if len(remaining) != len(added) or remaining & original:
        return False
    steps = [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)]
    queue = [x for x in remaining if any((tuple((x[i] + s[i] for i in range(3))) in original for s in steps))]
    seen = set(queue)
    while queue:
        x = queue.pop()
        for s in steps:
            y = tuple((x[i] + s[i] for i in range(3)))
            if y in remaining and y not in seen:
                seen.add(y)
                queue.append(y)
    return len(seen) == len(remaining)

def surface(points):
    lo = points.min(0) - 1
    hi = points.max(0) + 2
    mask = np.zeros(tuple(hi - lo), bool)
    mask[tuple((points - lo).T)] = True
    out = np.argwhere(mask & ~ndi.binary_erosion(mask)) + lo
    return out

def independent_minimum_additions(body, hazard, spacing, threshold):
    """Enumerate exact danger cubes; nearest L1 grid distance is shortest path.

    Independent of inherited taxicab distance-transform and binary-dilation
    implementations. No forbidden-tissue or image-boundary constraint is assumed.
    """
    radius = np.ceil(threshold / spacing).astype(int) + 1
    offsets = np.array(list(itertools.product(*(range(-r, r + 1) for r in radius))), int)
    gaps = np.linalg.norm(np.maximum(abs(offsets) * spacing - spacing, 0), axis=1)
    offsets = offsets[gaps < threshold]
    targetset = set()
    for start in range(0, len(body), 32):
        part = (body[start:start + 32, None, :] + offsets[None, :, :]).reshape(-1, 3)
        targetset.update(map(tuple, np.unique(part, axis=0)))
    targets = np.array(list(targetset), int)
    distances = cKDTree(hazard).query(targets, p=1, workers=1)[0]
    return (int(distances.min()), len(targets))

def bootstrap(endpoint, seed, draws):
    """Independent ratio-of-cluster-sums replay with inherited RNG schedule."""
    values = np.array(endpoint['values'], float).reshape(-1)
    ids = np.array(endpoint['cluster_ids'])
    (u, inv) = np.unique(ids, return_inverse=True)
    size = np.bincount(inv)
    hits = np.bincount(inv, weights=values)
    metricseed = seed + int(hashlib.sha256((endpoint['demo'] + endpoint['metric']).encode()).hexdigest()[:8], 16)
    rng = np.random.default_rng(metricseed)
    samples = []
    for start in range(0, draws, 100):
        b = min(100, draws - start)
        idx = rng.integers(len(u), size=(b, len(u)))
        samples.extend((hits[idx].sum(1) / size[idx].sum(1)).tolist())
        rng.integers(len(ids), size=(b, len(ids)))
    return (float(values.mean()), np.quantile(samples, [0.025, 0.975]).tolist(), len(u))

def main():
    start = time.perf_counter()
    cpu = time.process_time()
    for d in ['raw', 'figures', 'tables', 'logs']:
        (ROOT / d).mkdir(exist_ok=True)
    prereg = read('PREREG_R1_MANUSCRIPT.json')
    metric = prereg['metrics']
    for name in ['PREREG_R1_MANUSCRIPT', 'FROZEN_PREDICTIONS']:
        check(name + '_hash', sha(ROOT / (name + '.json')), (ROOT / (name + '.sha256')).read_text().strip())
    frozen = read('FROZEN_PREDICTIONS.json')
    check('source_manifest_hash', sha(ROOT / 'SOURCE_MANIFEST.json'), frozen['source_manifest_sha256'])
    for entry in read('SOURCE_MANIFEST.json')['files'] + read('SOURCE_MANIFEST_REFERENCES.json'):
        check('input_hash:' + entry['local'], sha(ROOT / entry['local']), entry['sha256'])
    X5 = 'inputs/LANE_X5_DECISION_SEG_ERROR/'
    X26 = 'inputs/LANE_X26_DECIDABILITY/'
    X31 = 'inputs/LANE_X31_CLINICAL_ANSWERS/'
    X35 = 'inputs/LANE_X35_PACKAGE_STATISTICS/'
    x5 = read(X5 + 'results.json')
    x26 = read(X26 + 'results.json')
    x31 = read(X31 + 'results.json')
    stats = read(X35 + 'results.json')
    raw = read(X5 + 'R3_CONNECTED_COUNTEREXAMPLES.json')
    probes = []
    tab = []
    check('counterexample_count', len(raw), frozen['predictions']['connected_probe_count'])
    for (i, r) in enumerate(raw):
        p = ROOT / 'inputs/geometry' / f"{r['case']}_cc{r['component']}_local.npz"
        with np.load(p, allow_pickle=False) as data:
            body = data['body_voxel']
            hazard = data['hazard_solid_voxel']
            s = data['spacing_zyx_mm']
            old_surface = data['hazard_surface_voxel']
        near = np.array(r['adverse_added_voxel'], int)
        far = np.array(r['safe_added_voxel'], int)
        refgap = float(direct_cube_field(body, hazard, s).min())
        ngap = min(refgap, float(direct_cube_field(body, near, s).min()))
        fgap = min(refgap, float(direct_cube_field(body, far, s).min()))
        near_hd = hd95_points(old_surface, surface(np.vstack([hazard, near])), s)
        far_hd = hd95_points(old_surface, surface(np.vstack([hazard, far])), s)
        dice_near = 2 * len(hazard) / (2 * len(hazard) + len(near))
        dice_far = 2 * len(hazard) / (2 * len(hazard) + len(far))
        (k, targetn) = independent_minimum_additions(body, hazard, s, 2)
        for (key, actual) in [('reference_gap_mm', refgap), ('adverse_gap_mm', ngap), ('safe_gap_mm', fgap)]:
            check(r['case'] + ':' + key, actual, r[key], metric['cube_distance_replay_tolerance_mm'])
            claim(r['case'] + ':' + key, r[key], X5 + 'R3_CONNECTED_COUNTEREXAMPLES.json#/' + str(i) + '/' + key, 'PER_SURFACE_REGION', 'mm', 'COMPUTED_ON_ANNOTATIONS')
        for (label, val) in [('Dice_adverse', dice_near), ('Dice_safe', dice_far)]:
            check(r['case'] + ':' + label, val, r[label], metric['Dice_pair_tolerance'])
            claim(r['case'] + ':' + label, r[label], X5 + 'R3_CONNECTED_COUNTEREXAMPLES.json#/' + str(i) + '/' + label, 'PER_ARCH', 'fraction', 'INJECTED_ERROR')
        for (label, val) in [('HD95_adverse_mm', near_hd), ('HD95_safe_mm', far_hd)]:
            check(r['case'] + ':' + label, val, r[label], metric['HD95_tolerance_mm'])
        check(r['case'] + ':different_decisions', ngap < 2 <= fgap, True)
        check(r['case'] + ':minimality', k, r['minimum_connected_added_voxels'])
        check(r['case'] + ':minimum_prediction', k, frozen['predictions']['minimum_added_voxels'][i])
        check(r['case'] + ':adverse_source_attachment', attached_to_source(hazard, near), True)
        check(r['case'] + ':control_source_attachment', attached_to_source(hazard, far), True)
        disconnected = near.copy()
        disconnected[-1] += 10000
        check(r['case'] + ':topology_fault_rejected', not attached_to_source(hazard, disconnected), True)
        volume = float(k * np.prod(s))
        check(r['case'] + ':volume', volume, r['minimum_added_volume_mm3'], 1e-12)
        claim(r['case'] + ':minimum_additions', k, X5 + 'R3_CONNECTED_COUNTEREXAMPLES.json#/' + str(i) + '/minimum_connected_added_voxels', 'PER_SURFACE_REGION', 'voxels', 'INJECTED_ERROR')
        out = {'case': r['case'], 'component': r['component'], 'reference_gap_mm': refgap, 'adverse_gap_mm': ngap, 'control_gap_mm': fgap, 'minimum_added_voxels': k, 'added_volume_mm3': volume, 'Dice': dice_near, 'HD95_mm': near_hd, 'target_cubes': targetn, 'reference_kind': 'annotation', 'physical_status': 'UNKNOWN'}
        probes.append(out)
        tab.append({'Case/component': r['case'].replace('ToothFairy2', '') + '/cc' + str(r['component']), 'Reference gap (mm)': f'{refgap:.3f}', 'Added voxels': k, 'Added volume (mm³)': f'{volume:.3f}', 'Adverse gap (mm)': f'{ngap:.3f}', 'Control gap (mm)': f'{fgap:.3f}', 'Dice both (%)': f'{100 * dice_near:.4f}', 'HD95 both (mm)': f'{near_hd:.1f}'})
        print(f"Replay {r['case']}: k={k}, gap={ngap:.6f} mm, source-attached controls PASS", flush=True)
        if i == 0:
            np.savez_compressed(ROOT / 'example_input.npz', body_voxel=body, reference_hazard_voxel=hazard, prediction_hazard_voxel=np.vstack([hazard, near]), spacing_zyx_mm=s)
    table('TABLE_1_CONNECTED_ERRORS', tab)
    save('raw/GEOMETRY_REPLAY.json', probes)
    canals = read(X26 + 'RAW_R2_CANAL.json')
    field = read(X26 + 'inputs/GEOM_field.json')['canal']
    floor = abs(field['all']['bias']) + 3 * (field['all']['sigma_w'] + field['tau'])
    expected_by_policy = collections.defaultdict(collections.Counter)
    current = collections.Counter()
    class_by_query = []
    with (ROOT / X26 / 'inputs/X5_margins.csv').open() as f:
        marginrows = list(csv.DictReader(f))
    for r in canals:
        m = r['margin']
        C = math.sqrt(3)
        limit = (abs(m) - floor) / C if abs(m) > floor else None
        klass = 1 if limit is not None else 2
        check(f"query:{r['record_index']}:class", klass, r['class'])
        if limit is not None:
            check(f"query:{r['record_index']}:h_limit", limit, r['pitch_limit_strict'], metric['endpoint_arithmetic_tolerance_mm'])
        orig = marginrows[r['record_index']]
        check(f"query:{r['record_index']}:raw_margin", float(orig['signed_margin_mm']), m, metric['endpoint_arithmetic_tolerance_mm'])
        check(f"query:{r['record_index']}:physical_status", r['physical_status'], 'UNKNOWN')
        expected_by_policy[r['selection']][klass] += 1
        current[r['conditional_status']] += 1
        class_by_query.append(klass)
    counts = collections.Counter(class_by_query)
    for (k, n) in counts.items():
        check('class_count:' + str(k), n, frozen['predictions']['conditional_class_counts'][str(k)])
    check('canal_scans', len({r['case'] for r in canals}), x26['canal_conditional']['represented_scans'])
    check('conditional_floor', floor, x26['canal_conditional']['floor_mm'], 1e-10)
    claim('conditional_floor', floor, X26 + 'results.json#/canal_conditional/floor_mm', 'PHENOMENOLOGICAL', 'mm', 'SCENARIO_CLOSURE')
    endpoints = read(X35 + 'raw/ENDPOINT_RECORDS.json')
    statsrows = stats['metrics']
    intervals = []
    guide_raw = list(csv.DictReader((ROOT / X31 / 'raw/GUIDE_PER_SITE.csv').open()))
    for ep in endpoints:
        if not (ep['demo'] == 'X26' or (ep['demo'] == 'X31' and 'apex-scenario' in ep['metric'])):
            continue
        row = next((v for v in statsrows if v['demo'] == ep['demo'] and v['metric'] == ep['metric']))
        if ep['demo'] == 'X26':
            recomputed_values = [float(k == 1) for k in class_by_query]
            ids = [r['case'] for r in canals]
        else:
            guide = ep['metric'].split()[0]
            rs = [r for r in guide_raw if r['guide'] == guide]
            recomputed_values = [float(r['two_mm_accept'] != r['local_apex_scenario_accept']) for r in rs]
            ids = [r['case'] for r in rs]
        check('endpoint_raw_values:' + ep['metric'], np.array(ep['values']).reshape(-1).tolist(), recomputed_values)
        check('endpoint_raw_clusters:' + ep['metric'], ep['cluster_ids'], ids)
        (point, ci, k) = bootstrap(ep, stats['seed'], stats['bootstrap_draws'])
        check('endpoint_point:' + ep['metric'], point, row['point'], metric['fraction_point_tolerance'])
        for j in range(2):
            check(f'endpoint_CI{j}:' + ep['metric'], ci[j], row['cluster_ci95'][j], metric['reported_interval_tolerance'])
        check('endpoint_cluster_count:' + ep['metric'], k, row['k_clusters'])
        intervals.append({'demo': ep['demo'], 'metric': ep['metric'], 'point': point, 'ci95': ci, 'n_records': len(ep['values']), 'clusters': k, 'scope': 'CONDITIONAL_DESCRIPTIVE; not anatomical/clinical uncertainty'})
        claim(ep['metric'], {'point': point, 'ci95': ci}, X35 + 'results.json#/metrics/' + str(statsrows.index(row)), 'POPULATION', 'fraction', 'CONDITIONAL_CASE_BOOTSTRAP')
    ci = next((r for r in intervals if r['demo'] == 'X26'))['ci95']
    n = len(canals)
    table('TABLE_2_DECIDABILITY', [{'Policy': key, 'Queries': sum(count.values()), 'Class 1': count[1], 'Class 2': count[2], 'Class 1 (%)': f'{100 * count[1] / sum(count.values()):.2f}', '95% case interval (%)': 'Not separately reported'} for (key, count) in expected_by_policy.items()] + [{'Policy': 'Combined dependent queries', 'Queries': n, 'Class 1': counts[1], 'Class 2': counts[2], 'Class 1 (%)': f'{100 * counts[1] / n:.2f}', '95% case interval (%)': f'{100 * ci[0]:.2f}–{100 * ci[1]:.2f}'}])
    table('TABLE_3_GUIDE_SCENARIOS', [{'Guide': guide, 'Apex scenario gap (mm)': f"{r['scenario_apex_required_gap_mm']:.3f}", 'Apex B sensitivity (mm)': f"{r['scenario_boundary_sensitivity_range_mm'][0]:.3f}–{r['scenario_boundary_sensitivity_range_mm'][1]:.3f}", 'Whole-body sufficient B=0.3 gap (mm)': f"{r['whole_body_sufficient_gap_B0p3_mm']:.3f}", 'Scenario decision changes (%)': f"{100 * next((x['point'] for x in intervals if x['metric'] == guide + ' apex-scenario decision change')):.2f}", '95% case interval (%)': '–'.join((f'{100 * v:.2f}' for v in next((x['ci95'] for x in intervals if x['metric'] == guide + ' apex-scenario decision change')))), 'Physical margin': 'UNKNOWN'} for (guide, r) in x31['guide_parameters'].items()])
    for (guide, r) in x31['guide_parameters'].items():
        claim(guide + ':apex_scenario_gap', r['scenario_apex_required_gap_mm'], X31 + 'results.json#/guide_parameters/' + guide + '/scenario_apex_required_gap_mm', 'PHENOMENOLOGICAL', 'mm', 'SCENARIO_CLOSURE')
        claim(guide + ':whole_body_bound', r['whole_body_sufficient_gap_B0p3_mm'], X31 + 'results.json#/guide_parameters/' + guide + '/whole_body_sufficient_gap_B0p3_mm', 'PHENOMENOLOGICAL', 'mm', 'SCENARIO_CLOSURE')
    table('TABLE_4_IDENTIFIABILITY', [{'Guide': r['guide'], 'Same radial mean±SD (mm)': f"{r['same_radial_mean_sd_mm'][0]:.2f}±{r['same_radial_mean_sd_mm'][1]:.2f}", 'Toward-plane sufficient scenario gap (mm)': f"{r['toward_plane_required_apex_gap_mm']:.3f}", 'Away-plane sufficient scenario gap (mm)': f"{r['away_plane_sufficient_apex_gap_mm']:.3f}", 'Level': 'PHENOMENOLOGICAL'} for r in x31['rounds']['R1']['orientation_counterexample']])
    save('raw/INTERVAL_REPLAY.json', intervals)
    root = ET.parse(ROOT / 'references/PMC9633839.xml').getroot()
    matches = [' '.join(e.itertext()) for e in root.iter('p') if '0.77' in ''.join(e.itertext()) and 'highest' in ''.join(e.itertext()) and ('median' in ''.join(e.itertext()))]
    check('external_SMCD_locator_found', bool(matches), True)
    save('EXTERNAL_FACIT_CHECK.json', {'doi': '10.1038/s41598-022-20605-w', 'local': 'references/PMC9633839.xml', 'locator': 'Results, highest-variability paragraph; Figure 4', 'reported_median_mm': 0.77, 'reported_quantity': 'median of highest interobserver symmetric MEAN curve distances', 'required_quantity': 'same-site directed maximum signed anatomical boundary displacement', 'transport_as_local_bound': 'REJECTED_METRIC_AND_COHORT_MISMATCH', 'matched_physical_facit_n': 0})
    for ref in read('REFERENCES_FINAL.json')['references']:
        check('DOI_identity:' + ref['doi'], ref['identity_verified'], True)
        check('DOI_metadata_hash:' + ref['doi'], sha(ROOT / ref['metadata_file']), ref['sha256'])
    command = [sys.executable, str(ROOT / 'assess_model.py'), '--input', str(ROOT / 'example_input.npz'), '--threshold-mm', '2', '--reference-kind', 'annotation', '--output', str(ROOT / 'EXAMPLE_ASSESSMENT.json')]
    subprocess.run(command, check=True, cwd=ROOT)
    assessed = read('EXAMPLE_ASSESSMENT.json')
    check('CLI_gap', assessed['prediction']['gap_mm'], probes[0]['adverse_gap_mm'], 1e-09)
    check('CLI_flip', assessed['decision_flip'], True)
    check('CLI_physical_status', assessed['physical_status'], 'UNKNOWN_UNLESS_REFERENCE_INDEPENDENTLY_VALIDATED')
    save('NUMERIC_CLAIMS.json', LEDGER)
    save('raw/CONTROLS.json', CHECKS)
    env = {'python': sys.version, 'numpy': np.__version__, 'scipy': scipy.__version__, 'platform': platform.platform(), 'threads': 1, 'gpu': False}
    save('ENVIRONMENT.json', env)
    result = {'lane': 'X41-manuscript-decision', 'claim_type': 'capability', 'review_state': 'PENDING_INDEPENDENT_REVIEW', 'scientific_admission': False, 'outcome': 'REPRODUCIBLE_DECISION_UNIT_MANUSCRIPT_PHYSICAL_REFERENCE_UNKNOWN', 'archive': {'sites': x5['rounds'][0]['counts']['source_sites'], 'represented_scans': x5['scan_count_archived'], 'source_candidates': x5['rounds'][0]['counts']['source_candidates']}, 'counterexamples': probes, 'conditional_canal': {'queries': n, 'scans': len({r['case'] for r in canals}), 'class1': counts[1], 'class2': counts[2], 'floor_mm': floor, 'current_status': dict(current), 'point_class1': counts[1] / n, 'case_ci95_class1': ci, 'per_policy': {k: dict(v) for (k, v) in expected_by_policy.items()}, 'parameter_level': 'PHENOMENOLOGICAL', 'query_level': 'PER_SURFACE_REGION', 'physical_certificates': 0}, 'guide_scenarios': x31['guide_parameters'], 'identifiability': x31['rounds']['R1']['orientation_counterexample'], 'intervals': intervals, 'bootstrap': {'draws': stats['bootstrap_draws'], 'seed': stats['seed'], 'definition': 'ratio of cluster hit sums to cluster record sums; whole CT cases resampled with dependent policies retained'}, 'attrition': {'X5': {'source_sites': 2342, 'represented_sites': 2342, 'scope': 'archive snapshot, not all public scans'}, 'X26': x26['canal_conditional']['dropout'], 'X31': x31['attrition']['guide'], 'matched_physical_reference': {'attempted_published_metric': 1, 'accepted': 0, 'rejected': 1, 'fraction': 1.0, 'reason': 'mean curve distance cannot bound local directed maximum'}}, 'external_referent': read('EXTERNAL_FACIT_CHECK.json'), 'controls': {'declared': len(CHECKS), 'pass': sum((x['pass'] for x in CHECKS)), 'fault_rejections': sum((x['injected_rejected'] for x in CHECKS)), 'interpretation': 'computational consistency, not empirical accuracy'}, 'source_review': read('inputs/PROOF_LANE_XREVIEW_CLINICAL/REVIEW_LANE_X5_DECISION_SEG_ERROR.json')['class'], 'source_review_limits': 'Corrections were not applied to original; stronger source-attachment check executed here. No own independent scientific review.', 'cost': {'replay_wall_seconds_before_figures': time.perf_counter() - start, 'cpu_seconds_before_figures': time.process_time() - cpu, 'peak_rss_kib_before_figures': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, 'threads': 1, 'gpu': False, 'fit': 0, 'source_selection_editorial_time': 'UNKNOWN', 'inherited_X5': [r['cost'] for r in x5['rounds']], 'inherited_X26': x26['full_cost'], 'new_physical_measurements': 0, 'total_project_cost': 'UNKNOWN'}, 'protocol_sha256': sha(ROOT / 'PREREG_R1_MANUSCRIPT.json'), 'frozen_predictions_sha256': sha(ROOT / 'FROZEN_PREDICTIONS.json'), 'next_construction': 'Freeze actual model output, acquire independently registered same-site anatomical surface, then score false geometric passes with registration-aware abstention'}
    save('results.json', result)
    figure_env = dict(os.environ, PYTHONNOUSERSITE='1')
    subprocess.run([sys.executable, str(ROOT / 'make_figures.py')], check=True, cwd=ROOT, env=figure_env)
    save('LAST_RUN_RECEIPT.json', {'utc': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'wall_seconds_including_figures': time.perf_counter() - start, 'cpu_seconds_parent': time.process_time() - cpu, 'peak_parent_rss_kib': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, 'success': True})
    save('CURRENT_WORK_STATE.json', {'lane': 'X41-manuscript-decision', 'status': 'R1_COMPUTATIONAL_REPLAY_COMPLETE', 'latest_gate': 'All frozen numeric, input, topology and interval controls reject injected faults; physical facit UNKNOWN', 'next_operation': 'Reader manuscript, adapted CLAIM checklist, held-out physical measurement protocol', 'review_state': 'PENDING_INDEPENDENT_REVIEW'})
    print(json.dumps({'status': 'PASS_COMPUTATIONAL_REPLAY', 'controls': len(CHECKS), 'physical': 'UNKNOWN', 'wall_seconds': time.perf_counter() - start}), flush=True)
if __name__ == '__main__':
    main()
