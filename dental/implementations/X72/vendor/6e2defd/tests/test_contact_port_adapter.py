from fractions import Fraction as Q
from dataclasses import replace
from types import SimpleNamespace
import pytest
import numpy as np
pytest.importorskip('field_engine', reason='optional field contact port package')
from field_engine.contact_port_v1 import *
from motion_engine.ncp.contact_port_adapter import *
S = Source('fixture', 'b' * 64, 'MODEL_ONLY')

@pytest.mark.parametrize('kind,g', [('SUPPORT', '0.00075'), ('VERTEX_TRIANGLE', '0.0005'), ('EDGE_EDGE', '0.0005')])
def test_cloth_thickness_convention(kind, g):
    gap = from_cloth_distance('0.001', thickness_m='0.0005', kind=kind, source=S)
    assert gap.lower == gap.upper == Q(g)

def test_mpm_requires_physical_cell_size():
    assert from_mpm_cells(-2, 3, cell_size_m='0.001', source=S).upper == Q(3, 1000)
    with pytest.raises(ValueError):
        from_mpm_cells(-2, 3, cell_size_m=0, source=S)
    with pytest.raises(TypeError):
        from_mpm_cells(-2, 3, source=S)

def test_motion_gap_is_not_scene_b():
    c = from_motion_contact({'gap': '1/3'}, source=S, participants=('a', 'b'), frame='world')
    assert c.gap.lower == Q(1, 3)
    with pytest.raises(ValueError):
        from_motion_contact({'b': 0.01}, source=S, participants=('a', 'b'), frame='world')

def test_normal_map_sigma_and_clip():
    g = Gap('-1/1000', '1/1000', Unit.M, S, sigma='19/100000', assurance=Assurance.CONDITIONAL, condition='z=1')
    r = solver_gap(g, dt_s='1/100', convention='RIGID_SPECULATIVE')
    assert r.normal_offset_m_s == (0, Q(1, 10)) and r.sigma_m_s == Q(19, 1000)
    assert solve_normal_gap(r, free_normal_m_s='-1/10', inverse_mass_per_kg=1) == (0, Q(1, 10))
    assert solver_gap(g, dt_s='1/100', convention='CLOTH_SIGNED').normal_offset_m_s == (-Q(1, 10), Q(1, 10))
    with pytest.raises(ValueError):
        replace(r, normal_offset_m_s=(0, 0))
    with pytest.raises(ValueError):
        replace(r, sigma_m_s=0)

@pytest.mark.parametrize('kwargs', [dict(mu='.4'), dict(inverse_mass_per_kg=-1)])
def test_unsupported_normal_law_rejected(kwargs):
    r = solver_gap(Gap(0, 1, Unit.M, S), dt_s=1, convention='RIGID_SPECULATIVE')
    args = dict(free_normal_m_s=-1, inverse_mass_per_kg=1)
    args.update(kwargs)
    with pytest.raises(ValueError):
        solve_normal_gap(r, **args)

def test_unknown_completion_rejected():
    r = solver_gap(Gap(0, 1, Unit.M, S, blocked=True), dt_s=1, convention='RIGID_SPECULATIVE')
    with pytest.raises(ValueError):
        solve_normal_gap(r, free_normal_m_s=-1, inverse_mass_per_kg=1)

def test_scene_realization_does_not_change_contract():
    from motion_engine.ncp.ncp_ref import Scene
    scene = Scene(np.eye(3), np.ones(3), np.array([-0.1, 0.0, 0.0]), np.array([0.0]), np.array([[0, -1]]), 0.01)
    g = Gap('-0.001', '0.001', Unit.M, S)
    req = solver_gap(g, dt_s=scene.dt, convention='RIGID_SPECULATIVE')
    p = point_scene(scene, req, contact_index=0, realized_gap_m='0.001')
    assert np.allclose(p.b(), np.zeros(3)) and scene.b_offset is None
    with pytest.raises(ValueError):
        point_scene(p, req, contact_index=0, realized_gap_m='0.001')
    with pytest.raises(ValueError):
        point_scene(scene, req, contact_index=0, realized_gap_m='0.002')

def test_reachability_preserves_merge_witnesses_and_unknown():
    from motion_engine.ncp.branch_set import propagate, Step, Child

    def transition(parent):
        return Step((Child(('merged',), ('from', parent)),), True)
    r = propagate([('a',), ('b',)], transition, 1, verify=lambda p, c: c.witness == ('from', p))
    bs = from_reachability(r, source=S)
    assert len(bs.branches) == 1
    data = __import__('json').loads(bs.coverage.payload_json)
    assert len(data['edges'][0]) == 2 and data['layers'][-1][0][1] == 2
    r = propagate([('a',)], lambda p: Step((), False, 'BUDGET'), 1, verify=lambda p, c: True)
    assert from_reachability(r, source=S).status == Status.UNKNOWN

def test_global_unknown_does_not_become_empty_complete():
    v = SimpleNamespace(status='UNKNOWN', reason='unsupported', velocity=(), impulse_witness=(), input_sha256='', scope='RIGID_ONLY')
    assert from_global_verdict(v, source=S).status == Status.UNKNOWN

def test_joint_query_refuses_marginal_or_unequal_relation():
    g = Gap('-1/1000', '1/1000', Unit.M, S)
    cs = tuple((ContactCandidate(k, ('body' + k, 'ground'), 'world', g, (0, 0, 1)) for k in ('a', 'b')))
    p = Port('field', 'motion', cs, BranchSet((), False, 'MODEL'))
    kwargs = dict(free_normal_m_s=-1, inverse_mass_per_kg=1, dt_s=1)
    with pytest.raises(ValueError):
        identical_normal_impulse_difference(p, 'a', 'b', **kwargs)
    p = replace(p, latent_box=(('xi', Q(-1, 1000), Q(1, 1000)),), affine_gaps=(AffineGap('a', 0, (('xi', 1),)), AffineGap('b', 0, (('xi', -1),))))
    with pytest.raises(ValueError):
        identical_normal_impulse_difference(p, 'a', 'b', **kwargs)

def test_global_velocity_unique_does_not_claim_unique_impulse():
    import json
    v = SimpleNamespace(status='UNIQUE', reason='fixed translation-model fixture', velocity=('1', '-2'), impulse_witness=('3', '4'), input_sha256='c' * 64, scope='RIGID_ONLY')
    bs = from_global_verdict(v, source=S)
    assert bs.status == Status.UNIQUE
    state = json.loads(bs.branches[0].state_json)
    assert set(state) == {'velocity'}
    assert state['velocity'] == ['1', '-2']
    payload = json.loads(bs.coverage.payload_json)
    assert payload['impulse_witness'] == ['3', '4']
    assert payload['velocity'] == ['1', '-2']
    assert payload['input_sha256'] == 'c' * 64

def test_identical_joint_affine_gap_has_zero_impulse_difference():
    g = Gap('-1/1000', '1/1000', Unit.M, S)
    cs = tuple((ContactCandidate(k, ('body' + k, 'ground'), 'world', g, (0, 0, 1)) for k in ('a', 'b')))
    p = Port('field', 'motion', cs, BranchSet((), False, 'MODEL'), latent_box=(('xi', Q(-1, 1000), Q(1, 1000)),), affine_gaps=(AffineGap('a', 0, (('xi', 1),)), AffineGap('b', 0, (('xi', 1),))))
    d = identical_normal_impulse_difference(p, 'a', 'b', free_normal_m_s=-1, inverse_mass_per_kg=2, dt_s=Q(1, 100))
    assert d == (Q(0), Q(0))

def test_fixed_width_fraction_cannot_shrink_the_scalar_impulse():
    r = solver_gap(Gap(0, 0, Unit.M, S), dt_s=1, convention='RIGID_SPECULATIVE')
    d = solve_normal_gap(r, free_normal_m_s=Q(np.int64(-2 ** 62)), inverse_mass_per_kg=Q(1, 4))
    assert d == (Q(2 ** 64), Q(2 ** 64))
    assert type(d[0].numerator) is int

def test_gap_is_normalized_before_unit_and_time_conversion():
    v = Q(np.int64(2 ** 62))
    r = solver_gap(Gap(v, v, Unit.M, S), dt_s=Q(1, 4), convention='RIGID_SPECULATIVE')
    assert r.normal_offset_m_s == (Q(2 ** 64), Q(2 ** 64))
    assert type(r.gap.lower.numerator) is int
