"""Adversarial contract regressions: minimum, arithmetic, source scope and OED."""
from copy import deepcopy
import json
from itertools import combinations
from field_engine.experimental.farkas_witness import farkas_witness, verify_certificate, adapt_claim_federation, adapt_bodytwin, adapt_g3_pair, load_exact_json
from copy import deepcopy
import pytest

@pytest.fixture(autouse=True)
def require_live_bidirectional_certificate_pipeline():
    m = cases()['interval_conflict']
    r = farkas_witness(m)
    assert r['status'] == 'INCONSISTENT' and verify_certificate(m, r['certificate'])
    bad = deepcopy(r['certificate'])
    bad['decision']['witness'] = ['0'] * 4
    assert not verify_certificate(m, bad)
    m = cases()['consistent']
    r = farkas_witness(m)
    assert r['status'] == 'CONSISTENT' and verify_certificate(m, r['certificate'])
S = {'population': 'declared mathematical model', 'location': 'test domain', 'window': 'steady'}

def variable(i):
    return {'id': i, 'quantity': i, 'unit': '1'}

def row(coeff, rhs, relation='<='):
    return {'coefficients': coeff, 'rhs': rhs, 'relation': relation, 'unit': '1'}

def claim(i, rows):
    return {'id': i, 'source': {'id': 'declared:' + i}, 'scope': deepcopy(S), 'constraints': rows}

def model(variables, claims, domain='inequality'):
    return {'schema': 'farkas-claims/v1', 'scope': deepcopy(S), 'variables': [variable(i) for i in variables], 'claims': claims, 'domain': domain}

def interval(i, lo, hi, v='x', unit='1'):
    return {'id': i, 'source': {'id': 'declared:' + i}, 'scope': deepcopy(S), 'variable': v, 'unit': unit, 'interval': [lo, hi]}

def measurement(i='sensor', error='1/10', cost=1, v='x'):
    return {'id': i, 'quantity': v, 'coefficients': {v: 1}, 'unit': '1', 'where': {'location': 'test domain', 'window': 'steady'}, 'scope': deepcopy(S), 'error_bound': error, 'resolution': '1/10', 'cost': cost}

def cases():
    a = model(['x', 'y'], [interval('low', 0, 1), interval('high', 3, 4), interval('irrelevant', 0, 10, 'y')])
    a['measurements'] = [measurement()]
    b = model(['x', 'y'], [claim('xpos', [row({'x': -1}, 0)]), claim('ypos', [row({'y': -1}, 0)]), claim('sum', [row({'x': 1, 'y': 1}, -1)])])
    c = model(['x', 'y', 'z'], [interval('pairA', 0, 0, 'z'), interval('pairB', 1, 1, 'z'), *deepcopy(b['claims'])])
    d = model(['x'], [interval('one', 0, 2), interval('two', 1, 3)])
    e = model(['t', 'u', 'v'], [claim('normal', [row({'t': 1}, 1, '=')]), claim('tangent', [row({'u': 1}, 2, '=')])], 'soc3')
    f = model(['t', 'u', 'v'], [claim('face', [row({'t': 1, 'u': 1}, 0, '=')]), claim('offface', [row({'v': 1}, 2, '=')])], 'soc3')
    return {'interval_conflict': a, 'joint_triple': b, 'greedy_trap': c, 'consistent': d, 'strong_soc': e, 'weak_soc': f}

def test_minimum_certificates_and_joint_conflict():
    for (name, minimum) in [('interval_conflict', 2), ('joint_triple', 3), ('greedy_trap', 2)]:
        m = cases()[name]
        r = farkas_witness(m)
        assert r['status'] == 'INCONSISTENT'
        assert len(r['minimal_subset']) == minimum
        assert r['minimality'] == 'CARDINALITY_MINIMUM'
        assert verify_certificate(m, r['certificate'])
        assert r['dual_weights']

def test_greedy_deletion_is_not_cardinality_minimum():
    m = cases()['greedy_trap']
    remaining = deepcopy(m)
    for c in m['claims']:
        trial = deepcopy(remaining)
        trial['claims'] = [x for x in trial['claims'] if x['id'] != c['id']]
        if trial['claims'] and farkas_witness(trial)['status'] == 'INCONSISTENT':
            remaining = trial
    assert len(remaining['claims']) == 3
    assert len(farkas_witness(m)['minimal_subset']) == 2

def test_consistent_is_exact_primal_not_solver_status():
    m = cases()['consistent']
    r = farkas_witness(m)
    assert r['status'] == 'CONSISTENT' and verify_certificate(m, r['certificate'])
    c = deepcopy(r['certificate'])
    c['decision']['witness'] = ['4']
    assert not verify_certificate(m, c)
    c = deepcopy(r['certificate'])
    c['decision']['status'] = 'arbitrary'
    assert not verify_certificate(m, c)

def test_cardinality_proof_missing_or_forged_subset_is_rejected():
    m = cases()['joint_triple']
    r = farkas_witness(m)
    c = deepcopy(r['certificate'])
    c['minimum_proof']['feasible_k_minus_1'].pop()
    assert not verify_certificate(m, c)
    c = deepcopy(r['certificate'])
    c['minimum_proof']['cardinality'] = 2
    assert not verify_certificate(m, c)
    c = deepcopy(r['certificate'])
    c['minimum_proof']['feasible_k_minus_1'][0]['witness'] = ['-999', '-999']
    assert not verify_certificate(m, c)

def test_negative_weights_nonzero_residual_and_zero_margin_are_rejected():
    m = cases()['interval_conflict']
    r = farkas_witness(m)
    for w in (['-1', '0', '0', '0'], ['1', '0', '0', '0'], ['0', '0', '0', '0']):
        c = deepcopy(r['certificate'])
        c['decision']['witness'] = w
        assert not verify_certificate(m, c)

def test_input_source_and_rhs_edits_invalidate_binding():
    m = cases()['interval_conflict']
    c = farkas_witness(m)['certificate']
    for edit in ('rhs', 'source', 'scope'):
        mm = deepcopy(m)
        if edit == 'rhs':
            mm['claims'][0]['interval'][1] = 5
        elif edit == 'source':
            mm['claims'][0]['source']['id'] = 'other'
        else:
            mm['scope']['window'] = 'transient'
        assert not verify_certificate(mm, c)

def test_exact_json_preserves_decimal_and_python_float_is_unknown():
    m = model(['x'], [interval('A', '0.1', '0.2'), interval('B', '0.3', '0.4')])
    txt = json.dumps(m).replace('"0.1"', '0.1').replace('"0.2"', '0.2')
    exact = load_exact_json(txt)
    r = farkas_witness(exact)
    assert r['status'] == 'INCONSISTENT' and verify_certificate(exact, r['certificate'])
    assert farkas_witness(json.loads(txt))['status'] == 'UNKNOWN'
    m['claims'][0]['interval'][0] = True
    assert farkas_witness(m)['status'] == 'UNKNOWN'

def test_interval_unit_conversion_and_incompatible_units():
    m = model(['x'], [interval('percent', 0, 50, unit='%'), interval('unitless', '3/5', '4/5')])
    r = farkas_witness(m)
    assert r['status'] == 'INCONSISTENT' and verify_certificate(m, r['certificate'])
    m['claims'][0]['unit'] = 'N'
    assert farkas_witness(m)['status'] == 'UNKNOWN'

def test_grouped_equalities_are_claims_not_individual_rows():
    m = model(['x'], [claim('a', [row({'x': 1}, 0, '=')]), claim('b', [row({'x': 1}, 1, '=')])])
    r = farkas_witness(m)
    assert len(r['minimal_subset']) == 2 and verify_certificate(m, r['certificate'])

def test_malformed_scope_dimensions_and_duplicate_ids_abstain():
    for how in ('scope', 'duplicate', 'variable', 'unit', 'relation'):
        m = cases()['joint_triple']
        if how == 'scope':
            m['claims'][0]['scope'] = {}
        elif how == 'duplicate':
            m['claims'][1]['id'] = m['claims'][0]['id']
        elif how == 'variable':
            m['claims'][0]['constraints'][0]['coefficients'] = {'missing': 1}
        elif how == 'unit':
            m['variables'][0]['unit'] = 'N'
        else:
            m['claims'][0]['constraints'][0]['relation'] = 'strictly less'
        assert farkas_witness(m)['status'] == 'UNKNOWN'

def test_certified_oed_uses_noise_expanded_repair_predictions():
    m = cases()['interval_conflict']
    r = farkas_witness(m)
    o = r['suggested_discriminating_measurement']
    assert o['status'] == 'CERTIFIED_DISCRIMINATING_UNDER_MODEL'
    assert o['min_observable_gap'] == '9/5'
    assert len(o['hypotheses']) == 2 and o['where']['location'] == 'test domain'
    m['measurements'][0]['error_bound'] = 2
    assert farkas_witness(m)['suggested_discriminating_measurement']['status'] == 'UNKNOWN'

def test_oed_missing_location_overlap_and_unbounded_repair_abstain():
    m = cases()['interval_conflict']
    m['measurements'][0].pop('where')
    assert farkas_witness(m)['suggested_discriminating_measurement']['status'] == 'UNKNOWN'
    m = cases()['joint_triple']
    m['measurements'] = [measurement()]
    assert farkas_witness(m)['suggested_discriminating_measurement']['status'] == 'UNKNOWN'

def test_measurement_cost_and_noise_change_selected_sensor():
    m = cases()['interval_conflict']
    m['measurements'] = [measurement('slow', cost=10), measurement('fast', cost=1)]
    assert farkas_witness(m)['suggested_discriminating_measurement']['measurement_id'] == 'fast'

def test_minimum_budget_keeps_only_verified_inconsistency():
    m = cases()['joint_triple']
    r = farkas_witness(m, max_oracle_calls=1)
    assert r['status'] == 'INCONSISTENT' and r['minimal_subset'] is None
    assert r['minimality'] == 'UNKNOWN' and verify_certificate(m, r['certificate'])
    assert farkas_witness(m, max_claims=2)['status'] == 'UNKNOWN'

def test_scalar_self_contradiction_has_cardinality_one():
    m = model(['x'], [claim('impossible', [row({}, -1)])])
    r = farkas_witness(m)
    assert r['minimal_subset'] == ['impossible'] and verify_certificate(m, r['certificate'])

def test_cone_strong_and_weak_cases():
    import pytest
    pytest.importorskip('clarabel')
    m = cases()['strong_soc']
    r = farkas_witness(m)
    assert r['status'] == 'INCONSISTENT' and len(r['minimal_subset']) == 2
    assert verify_certificate(m, r['certificate'])
    assert farkas_witness(cases()['weak_soc'])['status'] == 'UNKNOWN'

def test_native_prose_and_sign_only_federation_abstain():
    assert farkas_witness(adapt_bodytwin({'id': 'native', 'claim': 'x is impossible'}))['status'] == 'UNKNOWN'
    g = {'graph_id': 'G', 'claims': [{'id': 'c', 'sign': -1, 'validity': {'x': [0, 1]}}]}
    assert farkas_witness(adapt_claim_federation(g))['status'] == 'UNKNOWN'

def test_federation_extension_namespaces_and_preserves_hard_contract():
    m = cases()['consistent']
    g = {'graph_id': 'G', 'farkas_context': {k: v for (k, v) in m.items() if k != 'claims'}, 'claims': [{'id': c['id'], 'sign': 1, 'farkas_claim': c} for c in m['claims']]}
    a = adapt_claim_federation(g)
    r = farkas_witness(a)
    assert r['status'] == 'CONSISTENT' and verify_certificate(a, r['certificate'])
    assert r['certificate']['decision']['claim_ids'] == ['G:one', 'G:two']

def test_g3_different_regimes_and_rejected_conventions_never_become_conflict():
    s = {'quantity': 'q', 'population': 'p', 'regime': {'time': 'acute'}, 'measurement_method': 'm', 'survival_status': 'SUPPORTED_UNDER_MODEL', 'value': 1, 'unit': '1', 'source': {'id': 's'}}
    t = deepcopy(s)
    t['value'] = 10
    t['regime'] = {'time': 'chronic'}
    e = {'contradiction_node_id': 'n', 'evidence': {'sides': [s, t]}}
    assert farkas_witness(adapt_g3_pair(e))['status'] == 'UNKNOWN'
    t['regime'] = s['regime']
    t['survival_status'] = 'REJECTED_CONVENTION'
    assert farkas_witness(adapt_g3_pair(e))['status'] == 'UNKNOWN'

def test_farkas_transitive_relation_inference_is_never_automatic():
    m = cases()['joint_triple']
    m['claims'][0]['constraints'][0]['rhs'] = '1/1000000000000'
    r = farkas_witness(m)
    assert r['status'] == 'INCONSISTENT' and verify_certificate(m, r['certificate'])

def test_malformed_adapters_abstain_without_exception():
    for adapted in (adapt_bodytwin({'id': 'n', 'farkas_model': 'not an object'}), adapt_claim_federation(['not a graph']), adapt_g3_pair({'evidence': {'sides': ['bad', 'bad']}})):
        assert farkas_witness(adapted)['status'] == 'UNKNOWN'

def test_invalid_measurement_keeps_verified_conflict_and_unknown_oed():
    for measurements in ('bad', [17]):
        m = cases()['interval_conflict']
        m['measurements'] = measurements
        r = farkas_witness(m)
        assert r['status'] == 'INCONSISTENT' and verify_certificate(m, r['certificate'])
        assert r['suggested_discriminating_measurement']['status'] == 'UNKNOWN'

@pytest.mark.parametrize('field', ['error_bound', 'resolution', 'cost'])
def test_review_invalid_sensor_rational_keeps_verified_decision(field):
    m = cases()['interval_conflict']
    m['measurements'][0][field] = '1/0'
    r = farkas_witness(m)
    assert r['status'] == 'INCONSISTENT' and r['minimality'] == 'CARDINALITY_MINIMUM'
    assert verify_certificate(m, r['certificate'])
    assert r['suggested_discriminating_measurement']['status'] == 'UNKNOWN'

def test_review_nonobject_sensor_coefficients_keep_verified_decision():
    m = cases()['interval_conflict']
    m['measurements'][0]['coefficients'] = ['x']
    r = farkas_witness(m)
    assert r['status'] == 'INCONSISTENT' and verify_certificate(m, r['certificate'])
    assert r['suggested_discriminating_measurement']['status'] == 'UNKNOWN'

def test_review_nonstring_sensor_ids_are_excluded_without_losing_conflict():
    m = cases()['interval_conflict']
    m['measurements'] = [measurement('sensor'), measurement(17), measurement(['sensor'])]
    r = farkas_witness(m)
    o = r['suggested_discriminating_measurement']
    assert r['status'] == 'INCONSISTENT' and verify_certificate(m, r['certificate'])
    assert o['status'] == 'CERTIFIED_DISCRIMINATING_UNDER_MODEL'
    assert o['measurement_id'] == 'sensor' and len(o['unusable_candidates']) == 2

def test_review_bad_sensor_does_not_hide_another_valid_sensor():
    m = cases()['interval_conflict']
    bad = measurement('bad')
    bad['cost'] = '1/0'
    m['measurements'] = [bad, measurement('good')]
    r = farkas_witness(m)
    o = r['suggested_discriminating_measurement']
    assert r['status'] == 'INCONSISTENT' and verify_certificate(m, r['certificate'])
    assert o['measurement_id'] == 'good' and len(o['unusable_candidates']) == 1

def test_review_sub_float_separation_is_exact():
    m = model(['x'], [interval('upper', 0, 1), interval('lower', '1.000000000000000000000000000001', 2)])
    r = farkas_witness(m)
    assert r['status'] == 'INCONSISTENT' and verify_certificate(m, r['certificate'])

def test_review_support_reuse_requires_nonzero_dual():
    from field_engine.experimental._farkas_check import support_reusable
    A = [[1]]
    b = [1]
    assert not support_reusable(A, b, A, b, [0])
    assert not support_reusable(A, b, [[0]], [99], [0])
    assert support_reusable(A, b, A, b, [1])
    assert not support_reusable(A, b, [[0]], [99], [1])
