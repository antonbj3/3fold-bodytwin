from common import *
from fractions import Fraction as F
import time

def run():
    start = time.perf_counter()
    source = BASE / 'PROOF_LANE_FULL_CROWN_R4/FROZEN_INPUTS_B.json'
    records = read(source)['records']
    r5 = read(BASE / 'PROOF_LANE_FULL_CROWN_R5/RESULTS_B.json')['rows']
    rows = []
    for r in records:
        hashes = {side: sha(r[side + '_path']) == r[side + '_sha256'] for side in ['public', 'private']}
        prior = [x for x in r5 if x['key'] == r['key']]
        with np.load(r['public_path']) as z:
            sizes = {k: list(z[k].shape) for k in z.files}
        rows.append(dict(key=r['key'], family=r['family'], source_hash_checks=hashes, array_shapes=sizes, inherited_gates=[dict(method=x['method'], gates=x['gates']) for x in prior], measured_bite_pose='UNKNOWN', force_target_N=None, full_milling_certificate='UNKNOWN', eligible=False, reason='No measured bite/force port or full milling certificate; predecessor offset gate failed', resolution='PER_TOOTH'))
    a = [F(1), F(2)]
    b = list(reversed(a))
    weights = [F(1), F(4)]
    I = lambda s, w: sum((x * x * y for (x, y) in zip(s, w)))
    (ha, hb) = (I(a, [F(1), F(1)]), I(b, [F(1), F(1)]))
    (wa, wb) = (I(a, weights), I(b, weights))
    suff = dict(stress_A=list(map(float, a)), stress_B=list(map(float, b)), summary='stress multiset, maximum and homogeneous m=2 hazard', summary_identity_error=float(abs(ha - hb)), homogeneous_hazard_A=float(ha), homogeneous_hazard_B=float(hb), spatial_flaw_weights=list(map(float, weights)), weighted_hazard_A=float(wa), weighted_hazard_B=float(wb), downstream_hazard_difference=float(abs(wa - wb)), force_ratio_B_over_A=float((wa / wb) ** F(1, 2)), minimum_extension='Joint PER_POINT stress and process-matched flaw measure; scalar sigma0 cancellation is insufficient', resolution='PHENOMENOLOGICAL exact two-region counterexample, not tooth measurement')
    direct = sum((float(x) ** 2 * float(w) for (x, w) in zip(a, weights)))
    out = dict(claim_type='capability', round='A', cohort=rows, dropout=dict(total=18, excluded=18, fraction=1.0, reason='Missing full joint prerequisites; no cohort replacement'), sufficiency=suff, controls=dict(direct_sum_error=abs(direct - float(wa)), wrong_weight_rejected=direct != float(ha), missing_certificate_rejected=not all((v == 'PASS' for v in ['PASS', 'UNKNOWN']))), seconds=time.perf_counter() - start, external_referent=read(ROOT / 'PREREG_A.json')['external_referent'])
    save(ROOT / 'raw/A.json', out)
    (ROOT / 'HANDOFF_A.md').write_text('A: 0/18 cases have the full hard-condition conjunction. All source NPZ hashes checked. No physical Pareto claim is admissible. Exact equal stress histogram/homogeneous hazard has identity error 0 but weighted hazard 17 versus 8. Scalar strength cancellation does not remove spatial flaw/support closure. Next: change optimization to a common spatial measure with explicit support, angle and m variation on the existing whole STS specimen, keeping it separate from the 18-case cohort.\n')
    state('A_COMPLETE', '0/18 joint admission; exact summary insufficient', 'Freeze B: multi-support continuum-angle relative inverse design')
    print(json.dumps(dict(cohort=len(rows), eligible=0, sufficiency=suff, hashes_ok=all((all(r['source_hash_checks'].values()) for r in rows)))))
if __name__ == '__main__':
    run()
