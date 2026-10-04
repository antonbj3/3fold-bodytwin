from common import *
from parameter_bounds import Response
from scipy.optimize import linprog
import time

def run():
    start = time.perf_counter()
    p = read(ROOT / 'PREREG_C.json')
    z = np.load(DATA / 'B_FIELDS.npz')
    rows = []
    worst = None
    for support in p['scenarios']['supports']:
        ref = Response(z['reference__' + support + '__stress_basis'], z['reference__' + support + '__vol'])
        opt = Response(z['optimized__' + support + '__stress_basis'], z['optimized__' + support + '__vol'])
        for angle in p['scenarios']['angles_deg']:
            s0 = ref.principal(angle)
            s1 = opt.principal(angle)
            for m in p['scenarios']['m_grid']:
                scale = max(float(s0.max()), float(s1.max()))
                a = ref.v * (s0 / scale) ** m
                b = opt.v * (s1 / scale) ** m
                norm = float(a.sum())
                d = (b - a) / norm
                P = float(d[d > 0].sum())
                N = float(d[d < 0].sum())
                k = max(1.0, -N / P) if P > 0 else None
                if k is None:
                    rows.append(dict(support=support, angle_deg=angle, m=m, kappa=None, status='No adverse elements'))
                    continue
                w = np.where(d > 0, k, 1.0)
                res = float(w @ d)
                lp = linprog(-d, bounds=[(1, k)] * len(d), method='highs', options=dict(dual_feasibility_tolerance=1e-09, primal_feasibility_tolerance=1e-09))
                lpres = float(-lp.fun)
                err = abs(lpres - res)
                below = max(1, 1 + (k - 1) * 0.99)
                above = 1 + (k - 1) * 1.01
                rb = float(np.where(d > 0, below, 1.0) @ d)
                ra = float(np.where(d > 0, above, 1.0) @ d)
                uniform = float(d.sum())
                wrong_rejected = uniform < 0 and ra > 0 if k > 1 else True
                row = dict(support=support, angle_deg=angle, m=m, kappa=k, homogeneous_force_ratio=float((a.sum() / b.sum()) ** (1 / m)), threshold_residual=res, lp_residual=lpres, lp_error=err, below_difference=rb, above_difference=ra, uniform_difference=uniform, wrong_uniform_witness_rejected=wrong_rejected, adverse_element_count=int((d > 0).sum()), element_count=len(d), adverse_volume_fraction=float(opt.v[d > 0].sum() / opt.v.sum()), threshold_pass=bool(k >= 2), control_pass=bool(lp.success and err <= 1e-09 and (abs(res) <= 1e-10) and (rb <= 0) and (ra >= 0) and wrong_rejected), resolution='PER_TOOTH; kappa PHENOMENOLOGICAL derived requirement')
                rows.append(row)
                if worst is None or k < worst['kappa']:
                    worst = row.copy()
                    np.savez_compressed(DATA / 'C_ADVERSE_REGION.npz', vertices=z['vertices_optimized'], tetra=z['tetra'], adverse_mask=d > 0, normalized_signed_difference=d, witness_weights=w)
    out = dict(claim_type='capability', round='C', rows=rows, worst=worst, all_controls_pass=all((x.get('control_pass', False) for x in rows)), seconds=time.perf_counter() - start, lp_queries=len(rows), external_referent=p['external_referent'], physical_density_measurement='UNKNOWN', dropout=dict(total=len(rows), excluded=0), region_path=DATA / 'C_ADVERSE_REGION.npz', region_sha256=sha(DATA / 'C_ADVERSE_REGION.npz'))
    save(ROOT / 'raw/C.json', out)
    freeze(ROOT / 'FROZEN_PREDICTIONS_C.json', dict(prereg_sha256=sha(ROOT / 'PREREG_C.json'), field_sha256=sha(DATA / 'B_FIELDS.npz'), worst_scenario=worst, region_sha256=out['region_sha256'], measurement_status='NOT_RUN', prediction_scope='Distribution ranking under stipulated spatial flaw field, not individual fracture order'))
    state('C_COMPLETE', dict(kappa=worst['kappa'], controls=out['all_controls_pass']), 'Package field validation, lab falsifier and reproduction; physical feasibility remains unresolved')
    print(json.dumps(clean(dict(worst=worst, all_controls_pass=out['all_controls_pass'], seconds=out['seconds']))))
if __name__ == '__main__':
    run()
