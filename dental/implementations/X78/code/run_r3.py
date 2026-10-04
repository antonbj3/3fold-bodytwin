from common import *
from tissue import exact_box_query, translated_interval, rigid_displacement_enclosure
import math, time
import numpy as np
from fractions import Fraction as F

def main():
    r = read(ROOT / 'raw/R2_RESULTS.json')
    a = np.load(r['artifact']['path'])
    freeze('PREREG_R3_DIRECTED_STATE.json', dict(frozen_utc=now(), claim_type='capability', capability='Propagate a shared finite registration state into every fixed-frame tissue query and show whether scalar error magnitude is sufficient', obstacle='Scalar fit residual/error magnitude loses translation direction, shared cancellation, gray-value location and tissue dependence', changed_operation='Apply exactly +/- one native voxel along each axis to R2 queries. Carry one shared pose; query against fixed pulp/canal voxels. Exact integer geometry and global Lipschitz enclosure replace affine sensitivity.', consumer='Same-patient preparation/implant tissue constraint consumer', metrics=dict(summary_squared_norm_identity_error=0, required_downstream_difference_mm=0.01, enclosure_violations=0, common_rigid_invariance_error=0, injected_false_acceptances=0), tolerance=0, decision='Scalar insufficient only if exact norm identity holds and at least one actual-label downstream distance differs >=.01mm. All finite perturbations must lie in precomputed Lipschitz intervals. Physical errors remain UNKNOWN.', strongest_equal_information_control='Independent SciPy native-box projection from R2 plus direct exact common-motion cancellation; no algorithm claim', falsifier='Any admissible perturbation escapes bound, common transform is charged twice, or measured registration claim is inferred from scenario', full_cost=dict(preparation='Inherit R2 full source read and failed attempts', fit='No fit; conditional exact .3mm translation scenario', discovery='All 3 axes both signs; no best-only reporting', validation='Actual rigid invariance, inequality and injected wrong values', questions=0, fallback='Directed local error measurement port; physical source epsilon still UNKNOWN'), external_referent=r['external_referent'], resolution='PER_POINT', time_scale='SIMULTANEOUS', translation_norm_mm='3/10', scenario_status='PHENOMENOLOGICAL; replace with independent directed target registration-error observations'))
    pred = []
    for (family, rows) in [('pulp', r['pulp_queries']), ('canal', r['canal_queries'])]:
        for row in rows:
            base = row['distance_interval_mm'] if family == 'pulp' else row['union_distance_interval_mm']
            pred.append(dict(family=family, query=row['query'], nominal=base, predicted_interval_mm=translated_interval(base, F(3, 10), F(0), F(0))))
    freeze('FROZEN_PREDICTIONS_R3.json', dict(frozen_utc=now(), prereg_sha256=sha(ROOT / 'PREREG_R3_DIRECTED_STATE.json'), parent_result_sha256=sha(ROOT / 'raw/R2_RESULTS.json'), predictions=pred, physical_status='Not a measurement. Tissue and true registration-error terms UNKNOWN. Zero terms only mean exact annotation model sensitivity.'))
    prior = read(ROOT / 'FROZEN_PREDICTIONS_R3.json')
    assert prior['parent_result_sha256'] == sha(ROOT / 'raw/R2_RESULTS.json') or read(ROOT / 'FROZEN_PREDICTIONS_R3.json')['predictions'] == pred
    started = time.perf_counter()
    rows = []
    fails = []
    invariants = []
    for (family, qs) in [('pulp', a['prep_query_ticks']), ('canal', a['bone_profile_indices'] * 2)]:
        sets = {'pulp': a['pulp_centres_indices']} if family == 'pulp' else {k: a[k + '_canal_indices'] for k in ('tf1', 'maxillo', 'tf2')}
        for (i, q) in enumerate(qs):
            predicted = next((p['predicted_interval_mm'] for p in pred if p['family'] == family and p['query'] == i))
            for axis in range(3):
                e = np.eye(3, dtype=int)[axis]
                states = []
                norm_squared = []
                for sign in (-1, 1):
                    shift = sign * 2 * e
                    norm_squared.append(F(int(shift @ shift)) * F(3, 20) ** 2)
                    candidates = [exact_box_query(q + sign * 2 * e, c) for c in sets.values()]
                    distance = [min((d['distance_interval_mm'][k] for d in candidates)) for k in (0, 1)]
                    valid = distance[0] >= predicted[0] and distance[1] <= predicted[1]
                    if not valid:
                        fails.append(dict(family=family, query=i, axis=axis, sign=sign))
                    states.append(distance)
                    for c in sets.values():
                        d0 = exact_box_query(q, c)
                        d1 = exact_box_query(q + sign * 2 * e, c + sign * e)
                        invariants.append(d0['squared_mm2'] == d1['squared_mm2'])
                delta_exact = max(F(0), F.from_float(states[1][0]) - F.from_float(states[0][1]), F.from_float(states[0][0]) - F.from_float(states[1][1]))
                delta_lo = float(delta_exact)
                while F.from_float(delta_lo) > delta_exact:
                    delta_lo = math.nextafter(delta_lo, -math.inf)
                summary_error = abs(norm_squared[0] - norm_squared[1])
                assert summary_error == 0
                rows.append(dict(family=family, query=i, axis=axis, summary_squared_norm_mm2=[str(n) for n in norm_squared], summary_identity_error=float(summary_error), negative_distance=states[0], positive_distance=states[1], downstream_difference_lower_mm=delta_lo, predicted_interval_mm=predicted, resolution='PER_POINT', time_scale='SIMULTANEOUS'))
    assert not fails and all(invariants)
    biggest = max(rows, key=lambda r: r['downstream_difference_lower_mm'])
    assert biggest['downstream_difference_lower_mm'] >= 0.01
    faults = []
    for p in pred:
        bad = math.nextafter(p['predicted_interval_mm'][1] + 0.01, math.inf)
        faults.append(dict(family=p['family'], query=p['query'], wrong_distance_mm=bad, rejected=bad > p['predicted_interval_mm'][1]))
    assert all((f['rejected'] for f in faults))
    pivot = a['prep_query_ticks'][0]
    rotation_bounds = [dict(query=i, enclosure_mm=rigid_displacement_enclosure(q, pivot, F(3, 10), F(1, 100)), status='CONDITIONAL_MODEL_ONLY', physical_bound_mm=None, formula='tau+||q-c|| min(theta,2); finite rotation theorem', resolution='PER_POINT') for (i, q) in enumerate(a['prep_query_ticks'])]
    out = dict(claim_type='capability', status='SCALAR_REGISTRATION_MAGNITUDE_INSUFFICIENT_DIRECTED_STATE_REQUIRED', external_referent=r['external_referent'], summary_identity_error=0, comparisons=rows, strongest_witness=biggest, common_translation_checks=len(invariants), common_rigid_invariance_error=0, enclosure_violations=len(fails), faults=faults, global_finite_rotation_enclosures=rotation_bounds, minimum_extension='One directed shared 3-vector suffices for tested pure-translation family; general rigid pose requires SE(3), actual directed spatial error and anatomy error. Full physical patient bound UNKNOWN.', registration_model='PHENOMENOLOGICAL scenario, not independent measurement', wall_seconds=time.perf_counter() - started)
    dump(ROOT / 'raw/R3_RESULTS.json', out)
    (ROOT / 'HANDOFF_R3.md').write_text('# R3 handoff\n\nTwo opposite translations have exact squared error magnitude 9/100 mm² and identity error zero. Actual fixed external P1 masks give different downstream distances. The smallest added state for this tested family is the directed common translation vector. All finite translations fit the frozen Lipschitz enclosure; common query+tissue motion cancels exactly. The scenario is not measured registration error. Next construction needs independently directed same-Demo1 tissue/registration observations.\n')
    state('R3_DECIDED', 'Exact summary counterexample and global translation enclosures PASS', 'Package one-command replay and a source-bound prospective same-Demo1 observation port')
    print(out['status'], biggest['downstream_difference_lower_mm'])
if __name__ == '__main__':
    main()
