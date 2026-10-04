from dental_release.paths import expand as _release_expand
from common import *
import numpy as np, copy, time

def main():
    start = time.perf_counter()
    freeze(ROOT / 'PREREG_SUFFICIENCY.json', dict(frozen_utc=now(), claim_type='capability', capability='Test whether a summary is sufficient for its downstream consumer', metric='summary identity0 exactly; report downstream difference, smallest extension', tolerance=0, controls='copied exact geometry and support multiset; bootstrap sample means; same X49 nominal film maximum endpoint', falsifier='any summary identity nonzero', cost='wall time measured; original preparation reused; no fit,questions or acquisition', external_referent=dict(kind='our_own_fixture', locator='algebraic counterexamples plus existing independently sourced geometry', compared_quantity='summary sufficiency, not physical patient validation', refutes_us=True)))
    p = BASE / _release_expand('X54')
    sys.path.insert(0, str(p / 'code'))
    for q in [_release_expand('@DENTAL_EXTERNAL_ROOT@/projects/3fold_release_motion_staging_20260922/reproducibility/p1/data/ncp/loader.py'), _release_expand('@DENTAL_EXTERNAL_ROOT@/projects/3fold_release_motion_staging_20260922/reproducibility/p2/solver/ncp_ref.py')]:
        record(q)
    cm = module('contact_model', p / 'code/contact_model.py')
    r2 = module('x72_live_r2', p / 'code/run_r2.py')
    g = read(p / 'raw/geometry_001_h02.json')
    gg = copy.deepcopy(g)
    (a, _) = r2.force_control(g, support_scale={('upper', 16): 0.5, ('upper', 26): 2.0})
    (b, _) = r2.force_control(gg, support_scale={('upper', 16): 2.0, ('upper', 26): 0.5})
    ident = max((abs(x['gap_mm'] - y['gap_mm']) for (x, y) in zip(g['patches'], gg['patches'])))
    aid = max((abs(x['area_mm2'] - y['area_mm2']) for (x, y) in zip(g['patches'], gg['patches'])))
    diff = max((abs(a['force_shares_pp'][j][k] - b['force_shares_pp'][j][k]) for j in ['upper', 'lower'] for k in a['force_shares_pp'][j]))
    assert ident == aid == 0 and a['numerics']['all_pass'] and b['numerics']['all_pass']
    x54 = dict(summary='all regional rigid gaps and areas plus support reciprocal multiset', gap_identity_error_mm=ident, area_identity_error_mm2=aid, support_reciprocal_sum_identity_error=0, downstream_share_difference_pp=diff, argmax_before=a['argmax'], argmax_after=b['argmax'], force_target_identity_N=0, root_output_error_N=[a['force_control']['force_error_N'], b['force_control']['force_error_N']], minimal_extension='named tooth compliance plus load/preload', rigorous_NCP_enclosure='MISSING; this live output is conditional numerical simulation')
    A = np.zeros(16)
    B = np.tile([-1.0, 1.0], 8)
    assert A.mean() == B.mean() == 0
    rng = np.random.default_rng(72)
    indices = rng.integers(0, 16, (2000, 16))
    ba = A[indices].mean(1)
    bb = B[indices].mean(1)
    x2 = dict(summary='population sample mean', identity_error=0, bootstrap_width_A=float(np.quantile(ba, 0.975) - np.quantile(ba, 0.025)), bootstrap_width_B=float(np.quantile(bb, 0.975) - np.quantile(bb, 0.025)), minimal_extension='case distribution/design and dependence; mean alone is insufficient')
    x49 = read(BASE / 'LANE_X49_DESIGN_GATE/rounds/R3.json')
    run = Path(x49['run_directory'])
    a = read(run / 'valid.json')
    b = read(run / 'near_limit.json')
    va = a['rules']['film_min']
    vb = b['rules']['film_min']
    v = va['nominal_interval_mm'][1]
    u = vb['nominal_interval_mm'][1]
    assert v == u
    wall = dict(summary='nominal upper endpoint of film-minimum enclosure', identity_error_mm=v - u, nominal_endpoint_mm=v, before=va['status'], after=vb['status'], robust_intervals_mm=[va['robust_interval_mm'], vb['robust_interval_mm']], geometry_sha_identity=all((a['inputs'][k]['sha256'] == b['inputs'][k]['sha256'] for k in a['inputs'])), minimal_extension='independent calibrated full-region error bounds and signed topology')
    write(ROOT / 'raw/SUFFICIENCY.json', dict(X54=x54, X2=x2, X49=wall, cost_wall_s=time.perf_counter() - start, physical_external_validation='NOT_CLAIMED; exact-summary tests are counterexamples, not empirical facit'))
    print('SUFFICIENCY', diff, wall['before'], wall['after'])
if __name__ == '__main__':
    main()
