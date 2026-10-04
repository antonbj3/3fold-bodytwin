import math, json, time
from analyze_r1 import R, read, dump, state, slope
import numpy as np

def bracket(rs, lo, hi):
    w = math.log(rs['t'] / lo['t']) / math.log(hi['t'] / lo['t'])
    pred = lo['mean'] ** (1 - w) * hi['mean'] ** w
    A = np.array([[1.0, math.log(lo['t'])], [1.0, math.log(hi['t'])]])
    b = np.linalg.solve(A, np.log([lo['mean'], hi['mean']]))
    control = math.exp(b[0] + b[1] * math.log(rs['t']))
    return dict(w=w, predicted_N=pred, observed_N=rs['mean'], log_error=math.log(rs['mean'] / pred), log_measurement_variance=rs['v'] + (1 - w) ** 2 * lo['v'] + w * w * hi['v'], equally_informed_control_N=control, control_difference_N=abs(control - pred), local_exponent=math.log(hi['mean'] / lo['mean']) / math.log(hi['t'] / lo['t']), practice_1_398_N=lo['mean'] * (rs['t'] / lo['t']) ** 1.398, practice_1_5_N=lo['mean'] * (rs['t'] / lo['t']) ** 1.5, derived_n2_N=lo['mean'] * (rs['t'] / lo['t']) ** 2)

def run(local=False):
    t0 = time.perf_counter()
    rs = read('raw/CURATED_ROWS.json')
    blocks = {}
    for r in rs:
        if r['eligible']:
            blocks.setdefault((r['study'], r['protocol'], r['product']), []).append(r)
    prereg = 'PREREG_R3_LOCAL_BRACKET.json' if local else 'PREREG_R2_BRACKET_CALIBRATION.json'
    h = read(prereg)['metrics']['operational_halfwidth_log']
    out = []
    rejections = []
    for ((st, prot, prod), rows) in sorted(blocks.items()):
        rr = sorted(rows, key=lambda r: r['t'])
        if len(set((r['t'] for r in rr))) < 3:
            rejections.append(dict(study=st, protocol=prot, product=prod, reason='fewer than3 distinct thicknesses; no interior validation', rows=len(rr)))
            continue
        if local and st == 'PMC8558575':
            rr = [r for r in rr if 0.8 <= r['t'] <= 1.5]
        (lo, hi) = (rr[0], rr[-1])
        pts = []
        for target in rr[1:-1]:
            q = bracket(target, lo, hi)
            q.update(t_mm=target['t'], locator='doi:' + target['doi'] + ' ' + target['locator'], covered=abs(q['log_error']) <= h, window_N=[q['predicted_N'] * math.exp(-h), q['predicted_N'] * math.exp(h)], resolution=target['resolution'], quantity_kind=target['origin'])
            pts.append(q)
        out.append(dict(held_out_study=st, protocol=prot, product=prod, calibration_endpoints_t_mm=[lo['t'], hi['t']], calibration_endpoint_cost=2, endpoint_means_N=[lo['mean'], hi['mean']], train_studies=sorted({r['study'] for r in rs if r['study'] != st}), learning_from_train='NONE: h fixed in prereg; fixed affine relation; study held-out outcomes never set tolerance', target_points=pts, mean_squared_log_error=np.mean([p['log_error'] ** 2 for p in pts]).item(), coverage=np.mean([p['covered'] for p in pts]).item()))
    studyids = sorted({f['held_out_study'] for f in out})
    stats = {}
    for key in ['candidate', 'practice_1_398_N', 'practice_1_5_N', 'derived_n2_N']:
        st_mses = []
        st_cov = []
        for st in studyids:
            ps = [p for f in out if f['held_out_study'] == st for p in f['target_points']]
            es = [p['log_error'] if key == 'candidate' else math.log(p['observed_N'] / p[key]) for p in ps]
            st_mses.append(np.mean(np.square(es)).item())
            st_cov.append(np.mean([abs(e) <= h for e in es]).item())
        stats[key] = dict(study_balanced_logRMSE=math.sqrt(np.mean(st_mses)), study_balanced_coverage=float(np.mean(st_cov)), calibration_endpoint_groups_per_block=2 if key == 'candidate' else 1)
    crowns = {f['held_out_study'] for f in out if f['protocol'].startswith('crown')}
    count = sum((len(f['target_points']) for f in out))
    rep = len(crowns) >= 2
    gates = dict(logRMSE=stats['candidate']['study_balanced_logRMSE'] <= 0.2, coverage=stats['candidate']['study_balanced_coverage'] >= 0.9, operational_width=math.exp(2 * h) <= 1.5, min_two_crown_studies=rep, wrong_force_rejected=False)
    if local:
        gates['exact_3Y_protocol_replication'] = False
    p = out[0]['target_points'][0]
    wrongforce = p['predicted_N'] * 1.6
    rejected = abs(math.log(wrongforce / p['predicted_N'])) > h
    gates['wrong_force_rejected'] = rejected
    results = dict(claim_type='information_link', round='R3' if local else 'R2', folds=out, stats=stats, gates=gates, independent_studies=len(studyids), crown_studies_with_interior_targets=len(crowns), target_points=count, dropout=dict(rejected_blocks=len(rejections), total_blocks=len(blocks), reasons=rejections), fault=dict(multiplier=1.6, log_error=math.log(1.6), rejected=rejected), largest_control_difference_N=max((p['control_difference_N'] for f in out for p in f['target_points'])), same_information_control='ordinary endpoint affine log-load fit; mathematically equivalent', outcome='OPERATIONAL_TEST_PASS' if all(gates.values()) else 'LIMITED_OR_FAILED', elapsed_seconds=time.perf_counter() - t0)
    results.pop('elapsed_seconds')
    dump('raw/R3_RESULTS.json' if local else 'raw/R2_RESULTS.json', results)
    state('R3_LOCAL_BRACKET_DECIDED' if local else 'R2_BRACKET_DECIDED', gates, 'Freeze lab coefficient formulas, matched protocols and fault-rejecting scorer' if local else 'R3: constrain the 3Y bracket to actual D1/M2 lab designs; same h and gates; preserve full-range failure')
    (R / ('HANDOFF_R3.md' if local else 'HANDOFF_R2.md')).write_text(json.dumps(dict(stats=stats, gates=gates), indent=2) + '\nSame-information endpoint control agrees. These are operational windows, never population prediction intervals. Full study-out excludes target points but uses two explicit same-study endpoint calibration groups. Next construction specified in CURRENT_WORK_STATE.\n')
    print(json.dumps(dict(round=results['round'], stats=stats, gates=gates, crown_studies=len(crowns), target_points=count), indent=2))
    return results
if __name__ == '__main__':
    run('--local' in __import__('sys').argv)
