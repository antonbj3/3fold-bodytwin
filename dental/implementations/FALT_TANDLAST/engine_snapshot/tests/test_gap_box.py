from fractions import Fraction as Q
from dataclasses import replace
from copy import deepcopy
import pytest
from motion_engine.ncp import gap_box as G

@pytest.fixture
def scalar():
    return G.GapProblem.make([[1]], [1], [-1], [[0, '1/50']], '1/100')

@pytest.fixture
def scalar_answer(scalar):
    return G.enclose_gap_box(scalar)

def test_complete_scalar_witnessed_interval(scalar, scalar_answer):
    a = scalar_answer
    assert a['status'] == 'MULTIPLE' and a['complete'] and G.verify_gap_enclosure(scalar, a)
    assert (a['bounds']['normal_force:0']['lo'], a['bounds']['normal_force:0']['hi']) == (0, 100)
    assert (a['bounds']['velocity:0']['lo'], a['bounds']['velocity:0']['hi']) == (-1, 0)
    assert a['bounds']['net_force:0']['lo'] == -100

def test_closed_gap_nominal_not_box_answer(scalar):
    p = replace(scalar, gaps=((Q(1, 100), Q(1, 100)),))
    a = G.enclose_gap_box(p)
    assert a['status'] == 'UNIQUE' and G.verify_gap_enclosure(p, a)
    assert a['bounds']['velocity:0']['lo'] == -1
    assert a['bounds']['normal_force:0']['hi'] == 0

def test_exact_farkas_pruning(scalar):
    p = replace(scalar, gaps=((Q(1, 50), Q(3, 100)),))
    a = G.enclose_gap_box(p)
    assert a['status'] == 'UNIQUE' and G.verify_gap_enclosure(p, a)
    assert any((rec['kind'] == 'pruned' and rec['certificate']['kind'] == 'farkas' for rec in a['records']))
    assert a['bounds']['normal_force:0']['hi'] == 0

def test_prefix_budget_sticky(scalar):
    a = G.enclose_gap_box(scalar, max_nodes=2)
    assert a['status'] == 'UNCERTAIN' and (not a['complete']) and (a['bounds'] == {})
    assert not G.verify_gap_enclosure(scalar, a)

def test_zero_stub_rejected(scalar):
    a = G.enclose_gap_box(scalar)
    for b in a['bounds'].values():
        b['lo'] = b['hi'] = Q(0)
    assert not G.verify_gap_enclosure(scalar, a)

def test_missing_mode_rejected(scalar, scalar_answer):
    a = deepcopy(scalar_answer)
    a['records'].pop()
    assert not G.verify_gap_enclosure(scalar, a)

def test_forged_dual_rejected(scalar, scalar_answer):
    a = deepcopy(scalar_answer)
    leaf = next((r for r in a['records'] if r['prefix'] == ('closed',)))
    leaf['endpoints']['normal_force:0:lo']['y'] = [Q(0)] * len(leaf['endpoints']['normal_force:0:lo']['y'])
    assert not G.verify_gap_enclosure(scalar, a)

def test_input_hash_and_witness_binding(scalar, scalar_answer):
    assert not G.verify_gap_enclosure(replace(scalar, step=Q(1, 50)), scalar_answer)
    a = deepcopy(scalar_answer)
    a['bounds']['velocity:0']['lo_witness']['x'][0] = Q(-1)
    assert not G.verify_gap_enclosure(scalar, a)

def test_certificate_json_round_trip(scalar, scalar_answer):
    import json
    assert G.verify_gap_enclosure(scalar, json.loads(json.dumps(G.rational_json(scalar_answer))))

def test_planar_grip_entire_box_holds():
    p = G.GapProblem.make([[1, 0], [0, 1], [-1, 0], [0, 1]], [1, 1], [0, '-981/10000'], [['-1/500', '-3/2500']] * 2, '1/100', law='planar_coulomb', mu=['1/2'] * 2, compliance=[1, 1])
    a = G.enclose_gap_box(p, max_nodes=128)
    assert a['status'] == 'MULTIPLE' and G.verify_gap_enclosure(p, a)
    assert a['bounds']['velocity:1']['lo'] == a['bounds']['velocity:1']['hi'] == 0
    assert a['bounds']['contact_force:1']['lo'] == Q(981, 100)
    assert a['bounds']['net_force:1']['hi'] == 0

def test_planar_grip_can_slip_and_open():
    p = G.GapProblem.make([[1, 0], [0, 1], [-1, 0], [0, 1]], [1, 1], [0, '-981/10000'], [['-1/1000', '1/5000']] * 2, '1/100', law='planar_coulomb', mu=['1/2'] * 2, compliance=[1, 1])
    a = G.enclose_gap_box(p, max_nodes=128)
    assert a['status'] == 'MULTIPLE' and G.verify_gap_enclosure(p, a)
    assert a['bounds']['velocity:1']['lo'] == Q(-981, 10000) and a['bounds']['velocity:1']['hi'] == 0
    leaves = [rec for rec in a['records'] if rec['kind'] == 'leaf']
    assert leaves
    assert any(('open' in rec['prefix'] for rec in leaves))
    assert any(('slide-' in rec['prefix'] for rec in leaves))

def test_face_free_variable_and_force_fibre():
    p = G.GapProblem.make([[1], [1]], [1], [-1], [[0, 0], [0, 0]], 1)
    a = G.enclose_gap_box(p)
    assert a['status'] == 'MULTIPLE' and G.verify_gap_enclosure(p, a)
    assert a['bounds']['velocity:0']['lo'] == a['bounds']['velocity:0']['hi'] == 0
    assert a['bounds']['normal_force:0']['lo'] == 0 and a['bounds']['normal_force:0']['hi'] == 1
    assert a['bounds']['contact_force:0']['lo'] == a['bounds']['contact_force:0']['hi'] == 1

def test_box_corners_miss_interior_extremum():
    p = G.GapProblem.make([[1, 1], [1, 2]], [1, 1], [-1, 0], [[0, 1], [0, '1/10']], 1)
    a = G.enclose_gap_box(p)
    assert a['status'] == 'MULTIPLE' and G.verify_gap_enclosure(p, a)
    assert a['bounds']['velocity:1']['lo'] == Q(3, 10)
    witness = a['bounds']['velocity:1']['lo_witness']['x']
    assert witness[2] == Q(2, 5) and witness[3] == Q(1, 10)

def test_massratio_exact_control():
    p = G.GapProblem.make([[1, -1]], [1, '1/1000000'], [-1, 0], [[0, 0]], 1)
    a = G.enclose_gap_box(p)
    assert a['status'] == 'UNIQUE' and G.verify_gap_enclosure(p, a)
    assert a['bounds']['normal_force:0']['lo'] == Q(1000000, 1000001)

@pytest.mark.parametrize('value', [1.0, True, float('nan'), float('inf')])
def test_float_and_nonfinite_inputs_rejected(value):
    a = G.solve_gap_box([[1]], [1], [value], [[0, 1]], 1)
    assert a['status'] == 'UNCERTAIN' and (not a['complete'])

def test_unsupported_3d_coulomb_cannot_assert_completeness():
    a = G.solve_gap_box([[1]], [1], [-1], [[0, 1]], 1, law='coulomb3d', mu=[1])
    assert a['status'] == 'UNCERTAIN' and a['bounds'] == {}
    a = G.solve_gap_box([[1]], [1], [-1], [[0, 1]], 1, complete=True)
    assert a['status'] == 'UNCERTAIN'

def test_empty_inconsistent_gap_model_is_unknown():
    p = G.GapProblem.make([[0]], [1], [0], [[-1, -1]], 1)
    a = G.enclose_gap_box(p)
    assert a['status'] == 'UNCERTAIN' and a['bounds'] == {} and ('EMPTY_MODEL' in a['reason'])

def test_unbounded_force_fibre_cannot_give_finite_bounds():
    p = G.GapProblem.make([[0]], [1], [0], [[0, 0]], 1)
    a = G.enclose_gap_box(p)
    assert a['status'] == 'UNCERTAIN' and (not a['complete']) and (a['bounds'] == {})

def test_full_enumeration_same_exact_intervals(scalar, scalar_answer):
    control = G.enclose_gap_box(scalar, exhaustive=True)
    assert G.verify_gap_enclosure(scalar, control)
    assert {k: (v['lo'], v['hi']) for (k, v) in control['bounds'].items()} == {k: (v['lo'], v['hi']) for (k, v) in scalar_answer['bounds'].items()}

def test_all_eight_candidates_pruned_without_faking_coverage():
    p = G.GapProblem.make([[1]] * 8, [1], [-1], [['1/50', '3/100']] * 8, '1/100')
    a = G.enclose_gap_box(p, max_nodes=64)
    assert a['status'] == 'UNIQUE' and G.verify_gap_enclosure(p, a)
    assert a['stats']['visited_nodes'] == 17
    assert all((a['bounds'][f'normal_force:{i}']['hi'] == 0 for i in range(8)))

def test_candidate_load_and_conservative_measurement_contract(scalar, scalar_answer):
    c = scalar_answer['contact_candidates'][0]
    assert c['can_bear_load'] and (not c['must_bear_load']) and c['investigate_gap']
    a = deepcopy(scalar_answer)
    a['contact_candidates'][0]['must_bear_load'] = True
    assert not G.verify_gap_enclosure(scalar, a)

def test_arithmetic_dimension_budget():
    p = G.GapProblem.make([[1]] * 65, [1], [-1], [[0, 0]] * 65, 1)
    a = G.enclose_gap_box(p)
    assert a['status'] == 'UNCERTAIN' and a['reason'] == 'EXACT_ARITHMETIC_DIMENSION_BUDGET' and (a['bounds'] == {})

def test_forged_input_hash_rejected_even_if_equations_same(scalar, scalar_answer):
    a = deepcopy(scalar_answer)
    a['input_sha256'] = '0' * 64
    assert not G.verify_gap_enclosure(scalar, a)
