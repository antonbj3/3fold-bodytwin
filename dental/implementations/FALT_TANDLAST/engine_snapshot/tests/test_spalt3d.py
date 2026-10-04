from copy import deepcopy
from fractions import Fraction as Q
from dataclasses import replace
import pytest
from motion_engine.ncp.gap_box import GapProblem, enclose_gap_box, verify_gap_enclosure, rational_json
from motion_engine.ncp.polyhedral_gap import PolyhedralGapProblem, DIAMOND, SQUARE, solve_polyhedral_gap
from motion_engine.ncp.shared_gap import shared_gap_enclosure, verify_shared_gap
from motion_engine.ncp.coulomb_outer import SparseCoulombBox, energy_enclosure, verify_energy_enclosure, sqrt_upper

def single(polygon=DIAMOND):
    return PolyhedralGapProblem.make([[1, 0, 0], [0, 1, 0], [0, 0, 1]], [1, 1, 1], [-1, -1, -1], [[0, 0]], 1, mu=['1/2'], polygon=polygon)

def outer_pb():
    return SparseCoulombBox.make([[(0, 1)], [(1, 1)], [(2, 1)]], [1, 1, 1], [-1, -2, -1], [['-1/10', '1/10']], 1, ['1/2'])

def test_edge_sliding_interior_is_covered():
    pb = single()
    a = enclose_gap_box(pb, max_nodes=32)
    assert verify_gap_enclosure(pb, a)
    assert Q(a['bounds']['velocity:1']['lo']) == Q(-3, 4)
    assert Q(a['bounds']['velocity:1']['hi']) == Q(-3, 4)
    edges = [r for r in a['records'] if r['kind'] == 'leaf' and r['prefix'] == ('edge:0',)]
    assert len(edges) == 1
    b = deepcopy(a)
    b['records'] = [r for r in b['records'] if r['prefix'] != ('edge:0',)]
    assert not verify_gap_enclosure(pb, b)

@pytest.mark.parametrize('polygon', [DIAMOND, SQUARE])
def test_polygon_matches_full_and_rejects_corruption(polygon):
    pb = single(polygon)
    a = enclose_gap_box(pb)
    b = enclose_gap_box(pb, exhaustive=True)
    assert verify_gap_enclosure(pb, a) and verify_gap_enclosure(pb, b)
    assert {k: (v['lo'], v['hi']) for (k, v) in a['bounds'].items()} == {k: (v['lo'], v['hi']) for (k, v) in b['bounds'].items()}
    bad = deepcopy(a)
    bad['bounds']['velocity:1']['hi'] += 1
    assert not verify_gap_enclosure(pb, bad)
    bad = deepcopy(a)
    bad['input_sha256'] = '0' * 64
    assert not verify_gap_enclosure(pb, bad)

def test_polygon_validation_and_budget():
    for bad in [tuple(reversed(DIAMOND)), ((0, 0), (1, 0), (0, 1)), ((1, 1), (0, 1), (-1, 1), (-1, -1), (1, -1))]:
        with pytest.raises(ValueError):
            single(bad)

    class LargeSparse:
        shape = (2313, 1200)

        def toarray(self):
            raise AssertionError('no materialization before budget')
    out = solve_polyhedral_gap(LargeSparse(), [1] * 1200, [0] * 1200, [[0, 1]] * 771, 1, mu=[Q(1, 2)] * 771)
    assert out['status'] == 'UNCERTAIN' and out['bounds'] == {}

@pytest.mark.parametrize('mu', [Q(0), Q(1, 4), Q(1), Q(4)])
@pytest.mark.parametrize('polygon', [DIAMOND, SQUARE])
def test_polygon_friction_changes_the_solution(mu, polygon):
    pb = PolyhedralGapProblem.make([[1, 0, 0], [0, 1, 0], [0, 0, 1]], [1, 1, 1], [-1, -1, -1], [[0, 0]], 1, mu=[mu], polygon=polygon)
    a = enclose_gap_box(pb, max_nodes=32)
    assert verify_gap_enclosure(pb, a)
    expected = min(Q(0), -1 + (mu / 2 if polygon == DIAMOND else mu))
    assert (a['bounds']['velocity:1']['lo'], a['bounds']['velocity:1']['hi']) == (expected, expected)

@pytest.mark.parametrize('gaps', [[[0, 1], [0, Q(1, 10)]], [[-Q(1, 5), Q(1, 10)], [0, Q(1, 2)]], [[Q(1, 2), 1], [Q(1, 4), Q(1, 4)]]])
def test_shared_sharp_bounds_against_full(gaps):
    pb = GapProblem.make([[1], [1]], [Q(3, 2)], [-1], gaps, 1, compliance=[Q(1, 2), 2])
    a = shared_gap_enclosure(pb)
    b = enclose_gap_box(pb, exhaustive=True)
    assert verify_shared_gap(pb, rational_json(a)) and verify_gap_enclosure(pb, b)
    assert {k: (v['lo'], v['hi']) for (k, v) in a['bounds'].items()} == {k: (v['lo'], v['hi']) for (k, v) in b['bounds'].items()}
    bad = deepcopy(a)
    bad['bounds']['normal_force:0']['hi'] += Q(1, 100)
    assert not verify_shared_gap(pb, bad)
    assert not verify_shared_gap(replace(pb, J=((Q(1),), (Q(-1),))), a)

def test_shared_unsupported_license():
    pb = GapProblem.make([[1]], [1], [-1], [[0, 1]], 1)
    with pytest.raises(ValueError):
        shared_gap_enclosure(pb)

def test_outer_circle_point_and_rational_replay():
    pb = outer_pb()
    a = energy_enclosure(pb, [1, 0, 0])
    assert verify_energy_enclosure(pb, rational_json(a))
    stick = SparseCoulombBox.make([[(0, 1)], [(1, 1)], [(2, 1)]], [1, 1, 1], [-1, Q(-1, 4), Q(-1, 4)], [[0, 0]], 1, [Q(1, 2)])
    b = energy_enclosure(stick, [1, 0, 0])
    assert verify_energy_enclosure(stick, b)
    velocity_bounds = [v for (k, v) in b['bounds'].items() if k.startswith('velocity:')]
    assert len(velocity_bounds) == len(stick.inverse_mass)
    assert velocity_bounds
    assert all((v['lo'] <= 0 <= v['hi'] for v in velocity_bounds))
    assert a['existence'] == 'UNPROVED' and a['status'] == 'OUTER_ENCLOSURE'

def test_outer_friction_affects_virtual_separator():
    pb = outer_pb()
    a = energy_enclosure(pb, [2, 1, 0])
    assert a['alpha'] == Q(3, 2) and verify_energy_enclosure(pb, a)
    high = replace(pb, mu=(Q(1),))
    b = energy_enclosure(high, [2, 1, 0])
    assert b['alpha'] == 1 and verify_energy_enclosure(high, b)
    assert not verify_energy_enclosure(high, a)

@pytest.mark.parametrize('mutation', ['mu', 'escape', 'alpha', 'cap', 'velocity', 'hash'])
def test_outer_rejects_unsound_mutations(mutation):
    pb = outer_pb()
    a = energy_enclosure(pb, [1, 0, 0])
    bad = deepcopy(a)
    if mutation == 'mu':
        pb = replace(pb, mu=(Q(2),))
    elif mutation == 'escape':
        bad['escape'] = (Q(-1), 0, 0)
    elif mutation == 'alpha':
        bad['alpha'] *= 2
    elif mutation == 'cap':
        bad['total_normal_impulse_upper'] = 0
    elif mutation == 'velocity':
        bad['bounds']['velocity:0'] = {'lo': bad['velocity_center'][0], 'hi': bad['velocity_center'][0]}
    else:
        bad['input_sha256'] = '0' * 64
    assert not verify_energy_enclosure(pb, bad)

def test_outer_no_separator_and_empty_model():
    pb = outer_pb()
    a = energy_enclosure(pb, [-1, 0, 0])
    assert a['status'] == 'UNCERTAIN' and (not a['all_solutions_covered']) and (a['bounds'] == {})
    with pytest.raises(ValueError):
        SparseCoulombBox.make([], [], [], [], 1, [])
    with pytest.raises(TypeError):
        SparseCoulombBox.make([[(0.5, 1)], [(1, 1)], [(2, 1)]], [1, 1, 1], [-1, 0, 0], [[0, 0]], 1, [Q(1, 2)])
    assert not verify_energy_enclosure(replace(pb, compliance=(Q(-1),)), energy_enclosure(pb, [1, 0, 0]))

def test_sqrt_upper_exact_and_extreme():
    for q in [Q(0), Q(2), Q(10 ** 100, 3), Q(1, 2 ** 160)]:
        bound = sqrt_upper(q)
        assert bound >= 0 and bound ** 2 >= q

def test_reviewed_empty_certificate_rejected():
    pb = GapProblem.make([[1]], [1], [-1], [[0, 1]], 1)
    a = enclose_gap_box(pb)
    a['records'] = []
    a['bounds'] = {}
    a['status'] = 'ENTYDIG'
    assert not verify_gap_enclosure(pb, a)

def test_reviewed_enumeration_zero_budget_is_lazy():
    pb = GapProblem.make([[1]] * 12, [1], [-1], [[0, 1]] * 12, 1)
    a = enclose_gap_box(pb, exhaustive=True, max_nodes=0)
    assert a['status'] == 'UNCERTAIN' and a['stats']['visited_nodes'] == 0 and (not a['records'])
