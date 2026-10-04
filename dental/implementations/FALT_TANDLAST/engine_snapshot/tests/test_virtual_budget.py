from copy import deepcopy
from dataclasses import replace
from fractions import Fraction as Q
import pytest
from motion_engine.ncp.coulomb_outer import SparseCoulombBox, energy_enclosure
from motion_engine.ncp.virtual_budget import virtual_budget_enclosure, verify_virtual_budget
from motion_engine.ncp.gap_box import GapProblem, enclose_gap_box

def coupled():
    pb = SparseCoulombBox.make([[(0, 1)], [], [], [(0, -1), (1, 1)], [], []], [1, 1], [-1, -1], [['-1/10', '1/10']] * 2, 1, [0, 0], compliance=[1, 1])
    seed = energy_enclosure(pb, [1, 2])
    return (pb, seed)

def test_coupled_full_outcome_set_is_covered():
    (pb, seed) = coupled()
    answer = virtual_budget_enclosure(pb, seed, energy_fields=[['1/10', '1/5']], budget_fields=[[1, 1]], clusters=[('first', (0,)), ('both', (0, 1))])
    assert verify_virtual_budget(pb, answer)
    normal = GapProblem.make([[1, 0], [-1, 1]], [1, 1], [-1, -1], pb.gaps, 1, compliance=[1, 1])
    sharp = enclose_gap_box(normal, exhaustive=True, max_nodes=16)
    for (name, b) in sharp['bounds'].items():
        outer = answer['bounds'][name]
        assert outer['lo'] <= b['lo'] <= b['hi'] <= outer['hi']
    assert any((answer['bounds'][k]['hi'] < seed['bounds'][k]['hi'] for k in seed['bounds']))
    assert answer['existence'] == 'UNPROVED'
    assert answer['cluster_force_bounds']['both']['unit'] == 'N'

def test_exterior_work_cannot_be_dropped():
    (pb, seed) = coupled()
    with pytest.raises(ValueError, match='exterior'):
        virtual_budget_enclosure(pb, seed, budget_fields=[[1, 0]], clusters=[('first', (0,))])

def test_energy_field_must_pay_each_gap_budget():
    (pb, seed) = coupled()
    with pytest.raises(ValueError, match='gap budget'):
        virtual_budget_enclosure(pb, seed, energy_fields=[['1/10', '1/10']])

def test_zero_field_and_zero_outside_margin_are_valid():
    (pb, seed) = coupled()
    answer = virtual_budget_enclosure(pb, seed, budget_fields=[[0, 0], [1, 1]], clusters=[('first', (0,))])
    assert verify_virtual_budget(pb, answer)
    assert answer['budget_proofs'][-2]['upper'] == 0

@pytest.mark.parametrize('kind', ['force', 'velocity', 'margin', 'support', 'radius', 'cluster', 'status', 'existence', 'coverage', 'extra'])
def test_certificate_tampering_fails(kind):
    (pb, seed) = coupled()
    a = virtual_budget_enclosure(pb, seed, energy_fields=[['1/10', '1/5']], budget_fields=[[1, 1]], clusters=[('first', (0,))])
    b = deepcopy(a)
    if kind == 'force':
        b['bounds']['normal_force:0']['hi'] -= Q(1, 100)
    if kind == 'velocity':
        b['bounds']['velocity:0']['lo'] += Q(1, 100)
    if kind == 'margin':
        b['budget_proofs'][-1]['margins'] = list(b['budget_proofs'][-1]['margins'])
        b['budget_proofs'][-1]['margins'][1] = Q(1)
    if kind == 'support':
        b['budget_proofs'][-1]['upper'] -= Q(1, 100)
    if kind == 'radius':
        b['energy_proofs'][0]['radius'] -= Q(1, 100)
    if kind == 'cluster':
        b['cluster_force_bounds']['first']['hi'] -= Q(1, 100)
    if kind == 'status':
        b['status'] = 'ENTYDIG'
    if kind == 'existence':
        b['existence'] = 'PROVED'
    if kind == 'coverage':
        b['all_solutions_covered'] = False
    if kind == 'extra':
        b['unproved_claim'] = True
    assert not verify_virtual_budget(pb, b)

def test_model_and_reserve_changes_fail_closed():
    (pb, seed) = coupled()
    a = virtual_budget_enclosure(pb, seed)
    assert not verify_virtual_budget(replace(pb, mu=(Q(1), Q(0))), a)
    bad = deepcopy(seed)
    bad['total_normal_impulse_upper'] = Q(0)
    with pytest.raises(ValueError, match='reserve'):
        virtual_budget_enclosure(pb, bad)
    assert not verify_virtual_budget(pb, {})

@pytest.mark.parametrize('clusters', [[('x', (0, 0))], [('x', (2,))], [('x', ())], [('x', (0,)), ('x', (1,))]])
def test_invalid_cluster_queries_fail(clusters):
    (pb, seed) = coupled()
    with pytest.raises(ValueError, match='cluster'):
        virtual_budget_enclosure(pb, seed, clusters=clusters)

def test_public_function_produces_a_tighter_verified_answer():
    pb = SparseCoulombBox.make([[(0, 1)], [(1, 1)], [(2, 1)]], [1] * 3, [-1, -1, -1], [['-1/100', '1/100']], 1, ['1/2'])
    seed = energy_enclosure(pb, [3, 0, 0])
    a = virtual_budget_enclosure(pb, seed, energy_fields=[['1/100', 0, 0]], budget_fields=[[1, 0, 0]])
    assert a['all_solutions_covered'] and verify_virtual_budget(pb, a)
    assert a['bounds']['normal_force:0']['hi'] < seed['bounds']['normal_force:0']['hi']
    pb = replace(pb, free_velocity=(Q(-1), Q(-3, 5), Q(-4, 5)))
    seed = energy_enclosure(pb, [3, 0, 0])
    a = virtual_budget_enclosure(pb, seed, energy_fields=[['1/100', 0, 0]])
    for g in [Q(-1, 100), Q(0), Q(1, 100)]:
        pn = 1 - g
        v = [-g, Q(-3, 5) + pn * Q(3, 10), Q(-4, 5) + pn * Q(2, 5)]
        assert v[0] + g == 0 and pn > 0
        assert all((a['bounds'][f'velocity:{i}']['lo'] <= x <= a['bounds'][f'velocity:{i}']['hi'] for (i, x) in enumerate(v)))
        assert pn <= a['bounds']['normal_force:0']['hi']
