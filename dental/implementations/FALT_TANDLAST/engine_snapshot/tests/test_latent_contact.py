from copy import deepcopy
from fractions import Fraction as Q
import pytest
from motion_engine.ncp.latent_contact import make_family, certify_family, verify_family, port_family

def problem(**kw):
    d = dict(G=[[2, 1], [1, 2]], q0=[-5, -6], q1=[7, -7], box=[-1, 1], step='1/100', source_id='published-Siconos-matrix:declared-extension', source_sha256='0' * 64, unit_system='EXTERNAL_UNSPECIFIED')
    d.update(kw)
    return make_family(**d)

def test_exact_external_matrix_and_interior_total_minimum():
    pb = problem()
    a = certify_family(pb)
    assert verify_family(pb, a)
    assert [(c['lo'], c['hi']) for c in a['cells']] == [(Q(-1), Q(-1, 3)), (Q(-1, 3), Q(4, 21)), (Q(4, 21), Q(1))]
    c = a['cells'][1]
    assert c['p0'] == [Q(4, 3), Q(7, 3)]
    assert a['bounds']['total_normal_force']['lo'] == Q(1100, 3)
    assert a['bounds']['total_normal_force']['hi'] == 650
    assert a['bounds']['total_normal_force']['lo'] < 600
    assert a['variation'] == 'VARIABLE' and a['uniqueness'] == 'UNIQUE_PER_LATENT_VALUE'

def test_identical_cause_preserves_difference_not_rectangular_marginals():
    pb = problem(q0=['-1/10', '-1/10'], q1=['1/10', '1/10'], unit_system='SI')
    a = certify_family(pb)
    b = a['bounds']['force_difference:0:1']
    assert b['lo'] == b['hi'] == 0
    assert a['bounds']['normal_force:0']['hi'] > a['bounds']['normal_force:0']['lo']
    assert a['unit'] == 'N'

@pytest.mark.parametrize('mutation', ['drop', 'duplicate', 'coefficient', 'source', 'bound', 'witness', 'unit', 'status', 'coverage', 'readout', 'gap'])
def test_solver_free_replay_rejects_one_sided_mutations(mutation):
    pb = problem()
    a = certify_family(pb)
    assert verify_family(pb, a)
    if mutation == 'drop':
        a['cells'].pop(1)
    elif mutation == 'duplicate':
        a['cells'].append(deepcopy(a['cells'][0]))
    elif mutation == 'coefficient':
        a['cells'][0]['p0'][0] += 1
    elif mutation == 'source':
        pb['source_sha256'] = '1' * 64
    elif mutation == 'bound':
        a['bounds']['total_normal_force']['hi'] += 1
    elif mutation == 'witness':
        a['bounds']['total_normal_force']['lo_witness']['z'] = Q(-1)
    elif mutation == 'unit':
        a['unit'] = 'N'
    elif mutation == 'status':
        a['uniqueness'] = 'UNIQUE_GLOBAL_CONSTANT'
    elif mutation == 'coverage':
        a['complete'] = False
    elif mutation == 'readout':
        a['readouts']['total_normal_force'] = (Q(2), Q(1))
    elif mutation == 'gap':
        a['cells'][1]['lo'] += Q(1, 100)
    assert not verify_family(pb, a)

@pytest.mark.parametrize('G', [[[1, 1], [1, 1]], [[1, 2], [2, 1]], [[2, 1], [0, 2]]])
def test_invalid_spd_license_is_rejected(G):
    with pytest.raises(ValueError):
        problem(G=G)

def test_input_budgets_nonfinite_boolean_and_empty_readouts():
    for value in (True, float('nan'), float('inf')):
        with pytest.raises(ValueError):
            problem(q0=[value, -1])
    with pytest.raises(ValueError):
        problem(box=[1, -1])
    with pytest.raises(ValueError):
        problem(step=0)
    with pytest.raises(ValueError):
        make_family([[1] * 9] * 9, [0] * 9, [0] * 9, [-1, 1], step=1, source_id='x', source_sha256='0' * 64)
    with pytest.raises(ValueError):
        certify_family(problem(), readouts={})

def test_known_coupled_primal_family_checks_actual_parameter_intervention():
    G = [[3, 1, 1], [1, 4, 1], [1, 1, 5]]
    p0 = [Q(1), Q(0), Q(2)]
    p1 = [Q(1, 4), Q(0), Q(-1, 4)]
    w0 = [Q(0), Q(1), Q(0)]
    w1 = [Q(0), Q(1, 5), Q(0)]
    q0 = [w0[i] - sum((x * y for (x, y) in zip(G[i], p0))) for i in range(3)]
    q1 = [w1[i] - sum((x * y for (x, y) in zip(G[i], p1))) for i in range(3)]
    pb = make_family(G, q0, q1, [-1, 1], step=1, source_id='known-original-KKT-primal', source_sha256='0' * 64)
    a = certify_family(pb)
    assert verify_family(pb, a)
    values = []
    for k in range(19):
        z = Q(k, 9) - 1
        c = next((c for c in a['cells'] if c['lo'] <= z <= c['hi']))
        p = [x + y * z for (x, y) in zip(c['p0'], c['p1'])]
        expected = [x + y * z for (x, y) in zip(p0, p1)]
        assert p == expected
        values.append(p)
    assert len(values) == 19 and values[0] != values[-1]

def test_singleton_domain_and_json_replay():
    import json
    from motion_engine.ncp.latent_contact import encode
    pb = problem(box=[0, 0])
    a = certify_family(pb)
    assert len(a['cells']) == 1 and verify_family(pb, a)
    assert verify_family(json.loads(json.dumps(encode(pb))), json.loads(json.dumps(encode(a))))

def port():
    return dict(schema='threefold-contact-port/v1', producer='scanner', latent_box=[['xi', '-1/1000', '1/1000']], contacts=[dict(id=str(i), gap=dict(unit='m', assurance='MODEL_BOUND', lower='-1/1000', upper='1/1000')) for i in range(2)], affine_gaps=[dict(candidate_id=str(i), constant_m='0', coefficients=[['xi', '1']]) for i in range(2)])

def test_port_exact_shared_cause_and_invalid_sigma_only():
    p = port()
    pb = port_family(p, [[2, 1], [1, 2]], ['-1/10'] * 2, step='1/100')
    a = certify_family(pb)
    assert verify_family(pb, a) and a['bounds']['force_difference:0:1']['lo'] == 0
    p['contacts'][0]['gap']['assurance'] = 'CONDITIONAL'
    with pytest.raises(ValueError):
        port_family(p, [[2, 1], [1, 2]], ['-1/10'] * 2, step='1/100')

@pytest.mark.parametrize('mutation', ['latent', 'band', 'duplicate', 'unit'])
def test_port_binding_rejects_invalid_inputs(mutation):
    p = port()
    if mutation == 'latent':
        p['affine_gaps'][0]['coefficients'][0][0] = 'independent'
    elif mutation == 'band':
        p['contacts'][0]['gap']['upper'] = '0'
    elif mutation == 'duplicate':
        p['contacts'][0]['id'] = '1'
    elif mutation == 'unit':
        p['contacts'][0]['gap']['unit'] = 'mm'
    with pytest.raises(ValueError):
        port_family(p, [[2, 1], [1, 2]], ['-1/10'] * 2, step='1/100')
