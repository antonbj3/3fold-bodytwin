"""Consume a cached arch pair with tooth evidence; research flags, not clinical advice."""
from dental_release.paths import expand as _release_expand
from common import *
import argparse, pickle, time, resource
import numpy as np

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--case', default=_release_expand('@DENTAL_CASE_ID@'))
    ap.add_argument('--round', default='R3', choices=['R2', 'R3'])
    ap.add_argument('--specificity', type=float, default=0.95, choices=[0.9, 0.95])
    ap.add_argument('--out', default='raw/DEMO_CASE.json')
    a = ap.parse_args()
    st = time.perf_counter()
    frozen = read(P / f'FROZEN_PREDICTIONS_{a.round}.json')
    assert sha(frozen['model_path']) == frozen['model_sha256']
    bundle = pickle.loads(Path(frozen['model_path']).read_bytes())
    geo = next((json.loads(s) for s in (X7 / 'raw/geometry.jsonl').read_text().splitlines() if json.loads(s)['case_id'] == a.case))
    row = read(X21 / 'raw/cases' / (a.case + '.json'))
    op = read(X21 / 'raw/OPPOSITION_FEATURES_R2.json')[a.case]
    feat = {'arch__' + k: v for (k, v) in geo['features'].items()}
    feat.update({'pair__' + k: v for (k, v) in row['features'].items()})
    feat.update({'opposition__' + k: v for (k, v) in op['features'].items()})
    X = np.array([[feat.get(k) for k in bundle['columns']]], float)
    th = read(P / f'raw/THRESHOLDS_{a.round}.json')['candidate']
    result = read(P / f'raw/RESULTS_{a.round}.json')
    out = {}
    for f in FINDINGS:
        obj = bundle['models'].get(f)
        p = float(obj['candidate'].predict_proba(X)[0, 1]) if obj else None
        t = th[f][str(a.specificity)]['threshold']
        agg = result['findings'][f]['methods']['candidate'].get(str(a.specificity), {})
        out[f] = dict(score=p, threshold=t, screen_positive=bool(p >= t) if p is not None and t is not None else None, score_is='Empirical model score, not calibrated clinical disease probability', test_specificity=agg.get('fixed_threshold_test', {}).get('specificity'), test_sensitivity=agg.get('fixed_threshold_test', {}).get('sensitivity'), research_gate=agg.get('primary_gate'), noise_floor='UNKNOWN')
    if (P / 'FROZEN_PREDICTIONS_R4.json').exists():
        fr = read(P / 'FROZEN_PREDICTIONS_R4.json')
        assert sha(fr['model_path']) == fr['model_sha256']
        mdl = pickle.loads(Path(fr['model_path']).read_bytes())
        pp = float(mdl['model'].predict_proba(np.array([[feat.get(k) for k in mdl['columns']]], float))[0, 1])
        tr = read(P / 'raw/THRESHOLDS_R4.json')['candidate'][str(a.specificity)]['threshold']
        rr = read(P / 'raw/RESULTS_R4.json')['methods']['candidate'][str(a.specificity)]
        out['any_open_bite'] = dict(score=pp, threshold=tr, screen_positive=pp >= tr, scope='Any positively reported anterior/lateral/unspecified open-bite region; R4 scope extension', test_specificity=rr['fixed_threshold_test']['specificity'], test_sensitivity=rr['fixed_threshold_test']['sensitivity'], research_gate=rr['primary_gate'], noise_floor='UNKNOWN')
        out['open_bite']['scope'] = 'Anterior only, original R1-R3 endpoint retained'
    evidence = [dict(upper_fdi=r['upper_fdi'], lower_fdi=r['lower_fdi'], projected_gap_mm=r['minimum_projected_gap_mm'], declared_pose_scenario_mm=r['conditional_vertical_offset_interval_mm'], scenario_is='Declared +/-0.05 mm offset, not measured physical uncertainty enclosure', source_triangles={k: r['witness'][k] for k in ['upper_source_face', 'lower_source_face']}, source_status='Predicted target FDI and conditional projected geometry; physical contact UNKNOWN') for r in row['pairs']]
    write(P / a.out, dict(case_id=a.case, round=a.round, target_specificity=a.specificity, claim_type='information_link', resolution='PER_ARCH', input_resolution='PER_TOOTH', timescale='SIMULTANEOUS', findings=out, tooth_pair_evidence=evidence, opposition_landmark_candidates=op.get('teeth'), source_zip_members=row['source_zip_members'], model_sha256=frozen['model_sha256'], threshold_sha256=sha(P / f'raw/THRESHOLDS_{a.round}.json'), patient_uncertainty='No calibrated case-specific confidence, anatomical landmark validation or pose-repeat error budget.', wall_s=time.perf_counter() - st, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024))
    print('Written', a.out, 'case', a.case, flush=True)
if __name__ == '__main__':
    main()
