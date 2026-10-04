"""External measurement tracks. No generated part inherits an unmatched physical score."""
from dental_release.paths import expand as _release_expand
import sys, csv, time
from pathlib import Path
import numpy as np
from scipy.ndimage import distance_transform_edt, binary_erosion
from .common import *
sys.path.insert(0, str(ROOT / 'gencad_bench_v2/vendor/v1'))
from gencad_bench.calibration import fit
X13 = DENTAL / 'results/LANE_X13_CEMENT_GAP'
X14 = DENTAL / 'results/LANE_X14_SINTER_POSITION'
X12 = DENTAL / 'results/LANE_X12_PULPY3D'

def predict():
    start = time.perf_counter()
    rows = read(ROOT / 'data/published_measurements.json')['rows']
    frac = []
    for held in sorted({r['study'] for r in rows}):
        train = [r for r in rows if r['study'] != held]
        m = fit(train)
        for r in rows:
            if r['study'] == held:
                frac.append(dict(id=r['row_id'], heldout=held, training_studies=sorted({a['study'] for a in train}), candidate=m['model'], control=m['control']))
    cement = list(csv.DictReader((X13 / 'measurements.csv').open()))
    cp = []
    for row in cement:
        x = row['internal_spacer_um']
        if not x:
            cp.append(dict(id=row['row_id'], status='UNKNOWN_MISSING_CAD_SPACER'))
            continue
        train = [r for r in cement if r['source_family'] != row['source_family'] and r['region'] == row['region'] and r['internal_spacer_um']]
        if not train:
            cp.append(dict(id=row['row_id'], status='UNKNOWN_NO_TRAINING_REGION'))
            continue
        residual = np.array([float(r['measured_mean_um']) - float(r['internal_spacer_um']) for r in train])
        x = float(x)
        family_res = [np.median([float(r['measured_mean_um']) - float(r['internal_spacer_um']) for r in train if r['source_family'] == s]) for s in sorted({r['source_family'] for r in train})]
        cp.append(dict(id=row['row_id'], status='PREDICTED', candidate_um=max(0, x + float(np.median(family_res))), practice_um=x, lower_um=max(0, x + float(np.quantile(residual, 0.05))), upper_um=max(0, x + float(np.quantile(residual, 0.95))), training_families=sorted({r['source_family'] for r in train}), heldout_family=row['source_family']))
    angle = read(X14 / 'FROZEN_PREDICTIONS.json')['predictions']
    shrink = read(X14 / 'raw/shrinkage_cells.json')
    sp = []
    for (i, r) in enumerate(shrink):
        train = [q for q in shrink if q['position'] != r['position'] and q['material'] == r['material'] and (q['method'] == r['method'])]
        same_axis = [q for q in train if q['axis'] == r['axis']]
        sp.append(dict(id=i, candidate_pct=float(np.mean([q['mean_shrinkage_pct'] for q in same_axis])) if same_axis else None, practice_pct=float(np.mean([q['mean_shrinkage_pct'] for q in train])) if train else None, heldout_position=r['position']))
    inputs = [ROOT / 'data/published_measurements.json', X13 / 'measurements.csv', X14 / 'FROZEN_PREDICTIONS.json', X14 / 'raw/R1_test.json', X14 / 'raw/shrinkage_cells.json', X12 / 'raw/R8_teeth.csv', X12 / 'raw/STS_COMPARISON.json', DENTAL / 'tasks/swarm48/sources/LIT_CROWN/crown_fracture_papers.jsonl']
    payload = dict(claim_type='capability', fracture=frac, cement=cp, sinter_angle=angle, sinter_shrinkage=sp, source_hashes={str(p): sha(p) for p in inputs}, prospective=False, prior_exposure='Original lanes and this author have seen literature; retrospective folds, computational freeze only', physical_scores_for_generated_parts='UNKNOWN_SETUP_MISMATCH')
    freeze(ROOT / 'FROZEN_MEASUREMENT_PREDICTIONS.json', payload)
    dump(ROOT / 'raw/measurement_fit_cost.json', dict(seconds=time.perf_counter() - start))
    return payload

def scalar_error(pred, truth):
    p = np.asarray(pred, float)
    y = np.asarray(truth, float)
    if p.shape != y.shape or not p.size or (not np.isfinite(p).all()) or (not np.isfinite(y).all()):
        raise ValueError('invalid measurement comparison')
    return dict(mae=float(np.mean(np.abs(p - y))), rmse=float(np.sqrt(np.mean((p - y) ** 2))))

def score():
    frozen = read(ROOT / 'FROZEN_MEASUREMENT_PREDICTIONS.json')
    p = frozen['payload']
    if digest(p) != frozen['payload_sha256']:
        raise ValueError('measurement prediction hash drift')
    for (path, h) in p['source_hashes'].items():
        if sha(path) != h:
            raise ValueError('measurement source drift ' + path)
    out = dict(claim_type='capability', generated_part_calibration='UNKNOWN_SETUP_MISMATCH', prospective=False)
    truth = {r['row_id']: r for r in read(ROOT / 'data/published_measurements.json')['rows']}
    fr = []
    for r in p['fracture']:
        t = truth[r['id']]
        if r['heldout'] in r['training_studies']:
            raise ValueError('study leakage')
        for arm in ['candidate', 'control']:
            a = r[arm]
            (lo, hi) = (a['lower_N'], a['upper_N'])
            y = t['mean']
            score = hi - lo + 40 * max(0, lo - y, y - hi)
            fr.append(dict(id=r['id'], study=r['heldout'], arm=arm, prediction_N=a['prediction_N'], measured_N=y, reported_sd_N=t['sd'], interval_score_N=score, covered=lo <= y <= hi, locator=t['locator'], relative_error=abs(a['prediction_N'] - y) / y))
    out['fracture'] = dict(groups=len(truth), studies=len({t['study'] for t in truth.values()}), summary={arm: dict(study_balanced_interval_score_N=float(np.mean([np.mean([x['interval_score_N'] for x in fr if x['arm'] == arm and x['study'] == s]) for s in sorted({x['study'] for x in fr})])), study_balanced_coverage=float(np.mean([np.mean([x['covered'] for x in fr if x['arm'] == arm and x['study'] == s]) for s in sorted({x['study'] for x in fr})])), all_relative_errors_le_25pct=all((x['relative_error'] <= 0.25 for x in fr if x['arm'] == arm))) for arm in ['candidate', 'control']}, external_referent=dict(kind='independent_measurement', locator='https://doi.org/10.3390/ma17020365; https://doi.org/10.1186/s12903-023-03305-5; https://pmc.ncbi.nlm.nih.gov/articles/PMC10413631/', compared_quantity='published group mean crown fracture force in N', refutes_us=True))
    out['fracture']['external_referent']['locator'] = '; '.join(sorted({'https://doi.org/' + t['doi'] for t in truth.values()}))
    ct = {r['row_id']: r for r in csv.DictReader((X13 / 'measurements.csv').open())}
    cr = []
    for r in p['cement']:
        t = ct[r['id']]
        if r['status'] != 'PREDICTED':
            cr.append(dict(id=r['id'], status=r['status']))
            continue
        if r['heldout_family'] in r['training_families']:
            raise ValueError('cement source-family leakage')
        cr.append(dict(id=r['id'], status='SCORED', region=t['region'], study=t['study'], measured_um=float(t['measured_mean_um']), candidate_um=r['candidate_um'], practice_um=r['practice_um'], covered=r['lower_um'] <= float(t['measured_mean_um']) <= r['upper_um'], locator='https://doi.org/' + t['doi'] + ' ' + t['source_table'] + ' row ' + t['source_row']))
    summary = {}
    for region in sorted({x['region'] for x in cr if x['status'] == 'SCORED'}):
        a = [x for x in cr if x.get('region') == region]
        summary[region] = dict(cells=len(a), **{arm: scalar_error([x[arm + '_um'] for x in a], [x['measured_um'] for x in a]) for arm in ['candidate', 'practice']})
    out['cement'] = dict(cells=len(ct), scored=sum((x['status'] == 'SCORED' for x in cr)), summary=summary, gate_all_primary_MAE_le_20um=all((summary.get(k, {}).get('candidate', {}).get('mae', np.inf) <= 20 for k in ['marginal', 'axial', 'occlusal'])), external_referent=dict(kind='independent_measurement', locator=str(X13 / 'FACIT.csv'), compared_quantity='148 published regional gap group-mean cells, micrometres; dry/seated methods kept by row', refutes_us=True))
    at = {(r['material'], r['area']): r for r in read(X14 / 'raw/R1_test.json')}
    ar = []
    for r in p['sinter_angle']:
        t = at[r['material'], r['area']]
        ar.append(dict(material=r['material'], area=r['area'], candidate_deg=r['candidate_deg'], practice_deg=r['isotropic_control_deg'], measured_deg=t['delta_deg'], locator=t['locator']))
    angles = {a: scalar_error([x[a + '_deg'] for x in ar], [x['measured_deg'] for x in ar]) for a in ['candidate', 'practice']}
    st = read(X14 / 'raw/shrinkage_cells.json')
    sr = []
    for r in p['sinter_shrinkage']:
        t = st[r['id']]
        if r['candidate_pct'] is not None:
            sr.append(dict(**r, measured_pct=t['mean_shrinkage_pct'], locator=t['locator']))
    out['sinter'] = dict(angle_cells=len(ar), angles=angles, angle_gate_RMSE_le_point2deg=angles['candidate']['rmse'] <= 0.2, shrinkage_cells=len(st), shrinkage={a: scalar_error([x[a + '_pct'] for x in sr], [x['measured_pct'] for x in sr]) for a in ['candidate', 'practice']}, external_referent=dict(kind='independent_measurement', locator='https://doi.org/10.3390/ma18184234 Tables 3–4; https://doi.org/10.3390/ma18143217 Table 2', compared_quantity='angle change degrees; dimensional shrinkage percentage, distinct quantities', refutes_us=True), generated_crown_transfer='UNKNOWN')
    atlas = list(csv.DictReader((X12 / 'raw/R8_teeth.csv').open()))
    valid = [r for r in atlas if r['domain_truncated'] == 'False']
    tissues = []
    for r in valid:
        d = float(r['horn_boundary_min_mm'])
        eps = float(r['digital_two_surface_radius_mm'])
        rows = []
        for reduction in [0.5, 1.0, 1.5]:
            lo = max(0, d - eps) - reduction
            hi = d + eps - reduction
            decision = 'PASS' if lo >= 0.5 else 'FAIL' if hi < 0.5 else 'UNKNOWN'
            rows.append(dict(reduction_mm=reduction, remaining_lower_mm=lo, remaining_upper_mm=hi, status=decision))
        tissues.append(dict(case_hash=digest(r['case'])[:12], fdi=r['fdi'], queries=rows))
    probes = []
    for name in ['R8_P1_36_paired.npz', 'R8_P1_46_paired.npz', 'R8_P4_36_paired.npz']:
        path = Path(_release_expand('@DENTAL_WORK_ROOT@/X12')) / name
        a = np.load(path, allow_pickle=False)
        tooth = a['tooth']
        pulp = a['pulp']
        spacing = a['spacing']
        boundary = tooth & ~binary_erosion(tooth)
        pb = pulp & ~binary_erosion(pulp)
        dist = distance_transform_edt(~boundary, sampling=spacing)
        value = float(dist[pb].min())
        probes.append(dict(file=str(path), sha256=sha(path), global_pulp_boundary_min_mm=value, pulp_inside_tooth=bool(np.all(tooth[pulp])), comparison='global boundary minimum is not the atlas horn-only operator'))
    out['tissue'] = dict(atlas_rows=len(atlas), untruncated_rows=len(valid), cases=len({r['case'] for r in valid}), queries=len(tissues) * 3, IOS_patient_join='UNKNOWN_NO_SHARED_PATIENTS', anatomical_annotation_error='UNKNOWN; digital radius is not physical error guarantee', sts=read(X12 / 'raw/STS_COMPARISON.json'), array_probes=probes, external_referent=dict(kind='published_dataset', locator='https://doi.org/10.1007/978-3-031-72111-3_2; ' + str(X12 / 'raw/R8_teeth.csv'), compared_quantity='pulp/whole-tooth label distances in independent CT cohort, not IOS-specific safety', refutes_us=True))
    dump(ROOT / 'raw/MEASUREMENT_ROWS.json', dict(fracture=fr, cement=cr, sinter_angle=ar, sinter_shrinkage=sr, tissue=tissues))
    dump(ROOT / 'raw/MEASUREMENT_SCORES.json', out)
    return out
