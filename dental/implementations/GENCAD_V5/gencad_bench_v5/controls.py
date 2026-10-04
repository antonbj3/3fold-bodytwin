from common import *
import copy, math
from functional import functional_map, validate_function, quality
from physics import calibrate_candidate, PROTOCOL, ctrl, quantile
from boundary import verify_witness
from lab import icc
from contracts import observation_gate, validate_external, reconstruction_tiers

def run():
    tests = []
    suff = []

    def check(name, normal, bad):
        tests.append(dict(name=name, normal_accepted=bool(normal), injected_fault_rejected=not bool(bad)))
    for tolerance in [0.35, 0.5, 0.75, 1.0, 1.5]:
        check('reconstruction_tier_' + str(tolerance), reconstruction_tiers(tolerance)[str(tolerance)], reconstruction_tiers(tolerance + 0.001)[str(tolerance)])
    xy = np.array([[0.0, 0], [1, 0], [1, 1], [0, 1], [2, 0], [3, 0], [3, 1], [2, 1]])
    faces = np.array([[0, 1, 2], [0, 2, 3], [4, 5, 6], [4, 6, 7]])
    a = np.r_[np.full(4, 0.05), np.full(4, 0.2)]
    b = np.r_[np.full(4, 0.2), np.full(4, 0.05)]
    ca = functional_map(xy, faces, a)
    cb = functional_map(xy, faces, b)
    sa = np.array([ca['contact_area_mm2'], ca['contact_regions']])
    sb = np.array([cb['contact_area_mm2'], cb['contact_regions']])
    diff = quality.contact_map(xy, faces, a, b)['contact_symdiff_mm2']
    suff.append(dict(name='contact_area_and_count', summary_A=sa, summary_B=sb, identity_error=float(np.max(abs(sa - sb))), downstream_difference=diff, downstream_units='mm2 contact-mask symmetric difference', minimal_extension='For fixed reference-mask symmetric difference: overlap area with that reference (one scalar; one bit distinguishes this two-state fixture). Full source-addressed mask retained for arbitrary later spatial queries', resolution='PER_SURFACE_REGION', external_referent={'kind': 'our_own_fixture', 'locator': 'controls.py square pair', 'compared_quantity': 'counterexample only, no physical validity', 'refutes_us': True}))
    missing = functional_map(xy, faces, np.full(8, np.nan))
    bad = dict(missing, contact_regions=0, contact_area_mm2=0, penetration_area_mm2=0)
    check('missing_contact_support', validate_function(missing), validate_function(bad))
    bad = dict(ca, contact_area_mm2=1000)
    check('contact_area_within_support', validate_function(ca), validate_function(bad))
    check('spatial_contact_not_just_area', observation_gate(0, quality.contact_map(xy, faces, a, a)['contact_symdiff_mm2'], 1e-10), observation_gate(0, diff, 1e-10))
    va = np.array([1.0, 1.0, 0.0, 0.0])
    vb = np.array([1.0, 0.0, 1.0, 0.0])
    suff.append(dict(name='pass_rate', summary_A=[va.mean()], summary_B=[vb.mean()], identity_error=float(abs(va.mean() - vb.mean())), downstream_difference=float(np.mean(va != vb)), downstream_units='fraction different case decisions', minimal_extension='For disagreement with a fixed comparator: joint pass fraction (one scalar). Retain paired case vectors for cluster uncertainty and other comparisons', resolution='POPULATION', external_referent={'kind': 'our_own_fixture', 'locator': 'controls.py paired binary vectors', 'compared_quantity': 'counterexample only', 'refutes_us': True}))
    base = {k: 'matched' for k in PROTOCOL}
    base.update(quantity='fracture_force05', units='N', resolution='PER_TOOTH')
    cal = dict(base, characteristic_load_N=1000.0, Weibull_m=5.0)
    good = calibrate_candidate(base, cal)
    for key in PROTOCOL:
        broken = dict(base)
        broken[key] = 'WRONG'
        check('force_protocol_' + key, good['status'] == 'MATCHED_MODEL_QUANTILE', calibrate_candidate(broken, cal)['status'] == 'MATCHED_MODEL_QUANTILE')
    base.update(quantity='seated_film', units='um', resolution='PER_SURFACE_REGION')
    cal = dict(base, measurement_state='SEATED_CEMENTED', observed_um=70.0)
    good = calibrate_candidate(base, cal)
    broken = dict(cal, measurement_state='DRY_FIT')
    check('dry_gap_is_not_seated_film', good['status'] == 'MATCHED_OBSERVATION', calibrate_candidate(base, broken)['status'] == 'MATCHED_OBSERVATION')
    for r in read(V4 / 'FROZEN_LITERATURE_PREDICTIONS.json')['rows']:
        bad = dict(r, observed=r['predicted'] + 100 if r['kind'] == 'cement' else 1.6 * r['predicted'])
        tests.append(dict(name='reviewed_observation_' + r['study'] + '_' + str(r['x']), normal_accepted=ctrl.observed_gate(r), source_failure_retained=not ctrl.observed_gate(r), injected_fault_rejected=not ctrl.observed_gate(bad)))
    check('wrong_weibull_tail', observation_gate(0.05, 1 - math.exp(-(quantile(1000, 5) / 1000) ** 5), 1e-12), observation_gate(0.05, 1 - math.exp(-(1000 * (-math.log(0.05)) ** 0.2 / 1000) ** 5), 1e-12))
    boundary = read(ROOT / 'raw/R1B_BOUNDARY.json')
    for r in boundary['rows']:
        if not r['exact_core_feasible']:
            check('exact_witness_' + r['task_id'], verify_witness(r['witness_coefficients'], r['witness_lower'], r['witness_rhs']), verify_witness(r['witness_coefficients'], r['witness_lower'], '100000000000'))
    task = dict(task_id='case1', frame='native1')
    result = dict(task, units='mm', points=np.zeros((8, 3)))
    for (field, value) in [('task_id', 'wrong'), ('frame', 'wrong'), ('units', 'm'), ('points', np.ones((8, 3)) * 1000), ('points', np.full((8, 3), np.nan))]:
        check('external_' + field, validate_external(task, result), validate_external(task, dict(result, **{field: value})))
    ratings = np.array([[1, 1], [2, 2], [3, 3], [4, 4]], float)
    biased = ratings.copy()
    biased[:, 1] += 1
    normal = icc(ratings)
    bad = icc(biased)
    check('absolute_agreement_rejects_rater_bias', observation_gate(1, normal['ICC_2_1'], 1e-10), observation_gate(1, bad['ICC_2_1'], 1e-10))
    try:
        icc(np.full((4, 2), np.nan))
        reject = False
    except ValueError:
        reject = True
    tests.append(dict(name='ICC_missing_ratings', normal_accepted=True, injected_fault_rejected=reject))
    for s in suff:
        normal = np.array_equal(s['summary_A'], s['summary_B']) and s['identity_error'] == 0 and (s['downstream_difference'] > 0)
        changed = np.asarray(s['summary_B']).copy()
        changed[0] += 1e-06
        check('summary_identity_' + s['name'], normal, np.array_equal(s['summary_A'], changed))
    out = dict(tests=tests, sufficiency=suff, all_injected_faults_rejected=all((t['injected_fault_rejected'] for t in tests)), normal_source_gate_failures=sum((t.get('source_failure_retained', False) for t in tests)), normal_implementation_controls_pass=all((t['normal_accepted'] for t in tests if 'source_failure_retained' not in t)), rigorous_enclosure='Exact rational linear feasibility only; no continuous surface/elastic/physical enclosure')
    dump(ROOT / 'raw/CONTROLS.json', out)
    assert out['all_injected_faults_rejected'] and out['normal_implementation_controls_pass']
    return out
if __name__ == '__main__':
    run()
