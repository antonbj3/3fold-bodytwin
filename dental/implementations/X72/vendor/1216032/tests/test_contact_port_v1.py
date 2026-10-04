import json
from dataclasses import replace
from types import SimpleNamespace
from fractions import Fraction as Q
import pytest
from field_engine.contact_port_v1 import *
from field_engine.contact_port_adapter import from_contact_band, from_univariate, from_aligner_port, from_shell_equilibria
S = Source('fixture', 'a' * 64, 'MODEL_ONLY')

def sample():
    w = Witness('test', canonical({'input': 1}), digest({'input': 1}), S)
    b = Branch('a', canonical({'q': ['1/3', '1/2']}), (w,))
    g = Gap('-1/1000', '1/1000', Unit.M, S, sigma='1/10000')
    c = ContactCandidate('a', ('shell', 'tooth'), 'world', g, (0, 0, 1))
    return Port('field', 'motion', (c,), BranchSet((b,), True, 'SCALAR_MODEL', w), (('scan', Q(-1, 1000), Q(1, 1000)),), (AffineGap('a', 0, (('scan', 1),)),), evidence=(w,))

@pytest.mark.parametrize('value,expected', [('OSAKER', Status.UNKNOWN), ('UNCERTAIN', Status.UNKNOWN), ('UNKNOWN', Status.UNKNOWN), ('UNCERTAIN', Status.UNKNOWN), ('AMBIGUOUS', Status.MULTIPLE), ('ENTYDIG', Status.UNIQUE)])
def test_status_preserves_ambiguity(value, expected):
    assert status(value) == expected

def test_exact_units_and_sigma():
    g = Gap('-0.19', '0.19', Unit.MM, S, sigma='0.19', assurance=Assurance.CONDITIONAL, condition='z=1')
    assert g.to(Unit.M).sigma == Q(19, 100000)
    assert g.to(Unit.M).to(Unit.MM) == g
    assert Gap(0.1, 0.2, Unit.MM, S).to(Unit.M).lower == Q(0.1) / 1000

@pytest.mark.parametrize('low,high,expected', [(1, 2, Status.OPEN), (-2, -1, Status.CLOSED), (-1, 1, Status.UNKNOWN), (0, 0, Status.UNKNOWN)])
def test_gap_classes(low, high, expected):
    assert Gap(low, high, Unit.M, S).classify() == expected

@pytest.mark.parametrize('args', [dict(lower=2, upper=1), dict(lower=float('nan'), upper=1), dict(lower=0, upper=1, sigma=-1), dict(lower=0, upper=1, unit='cell'), dict(lower=0, upper=1, assurance=Assurance.CONDITIONAL)])
def test_bad_gap_rejected(args):
    d = dict(lower=0, upper=1, unit=Unit.M, source=S)
    d.update(args)
    with pytest.raises((ValueError, TypeError)):
        Gap(**d)

def test_nonzero_threshold_and_guard():
    g = Gap(1, 2, Unit.MM, S)
    assert g.classify(3) == Status.CLOSED
    assert replace(g, blocked=True).classify(3) == Status.UNKNOWN

def test_roundtrip_exact_witnesses_sources_relations():
    p = sample()
    assert Port.from_json(p.to_json()) == p
    assert p.branch_set.status == Status.UNIQUE
    assert replay_witnesses(p.branch_set, {'test': lambda x: x == {'input': 1}})
    assert not replay_witnesses(p.branch_set, {})

@pytest.mark.parametrize('mutate', [lambda d: d.update(schema='v2'), lambda d: d.update(physical_status='UNIQUE'), lambda d: d['contacts'][0]['gap'].update(sigma='-1'), lambda d: d['affine_gaps'][0].update(constant_m='1'), lambda d: d['contacts'][0].update(frame=''), lambda d: d.update(unsupported=1)])
def test_bad_wire_rejected(mutate):
    d = json.loads(sample().to_json())
    mutate(d)
    with pytest.raises((ValueError, TypeError)):
        Port.from_json(json.dumps(d))

def test_tampered_witness_replay_rejected():
    p = sample()
    w = replace(p.branch_set.coverage, payload_json=canonical({'input': 2}))
    bs = replace(p.branch_set, coverage=w)
    assert not replay_witnesses(bs, {'test': lambda d: True})

def test_unknown_sticky_and_empty_not_unique():
    p = sample()
    assert replace(p.branch_set, complete=False).status == Status.UNKNOWN
    assert BranchSet((), True, 'MODEL', p.branch_set.coverage).status == Status.EMPTY
    with pytest.raises(ValueError):
        BranchSet(p.branch_set.branches, True, 'MODEL')

def test_band_keeps_sigma_and_completion():
    band = SimpleNamespace(lower_mm=[1.0, -2.0], upper_mm=[2.0, -1.0], sigma_tot_mm=[0.19, 0.19], classes=[1, 0], completion_unknown=[False, True], sign_conflict=[False, False])
    c = from_contact_band(band, source=S, condition='joint gap sigma, z=3', participants=[('a', 'b'), ('c', 'b')], frame='world')
    assert c[0].gap.sigma == Q(0.19) / 1000
    assert c[1].gap.classify() == Status.UNKNOWN
    with pytest.raises(ValueError):
        from_contact_band(band, source=S, condition='', participants=[('a', 'b'), ('c', 'b')], frame='world')

def test_polynomial_replay_and_missing_root():
    from field_engine.univariate_equilibrium import certify
    cert = certify(['0', '-1', '0', '1'])
    bs = from_univariate(cert, source=S)
    assert bs.complete and len(bs.branches) == 3
    cert['roots'].pop()
    with pytest.raises(ValueError):
        from_univariate(cert, source=S)

def test_legacy_clinical_unknown_is_preserved():
    d = {'schema': 'dental-aligner-branch/v1', 'status': 'OSAKER', 'model_branches': {'natural': {'q_mm': ['1', '2']}, 'everted': {'q_mm': ['-2', '-1']}}}
    bs = from_aligner_port(d, source=S)
    assert bs.status == Status.UNKNOWN and len(bs.branches) == 2
    assert json.loads(bs.branches[0].witnesses[0].payload_json) == d

def test_wire_rejects_implicit_numeric_units():
    d = json.loads(sample().to_json())
    d['contacts'][0]['gap']['lower'] = -0.001
    with pytest.raises(ValueError):
        Port.from_json(json.dumps(d))

def test_covariance_payload_is_transported_without_independence_claim():
    p = sample()
    payload = {'unit': 'mm^2', 'covariance': [['1/25', '1/25'], ['1/25', '1/25']], 'projection': ['1', '-1'], 'joint_source': 'one scan'}
    w = Witness('joint_covariance', canonical(payload), digest(payload), S)
    p = replace(p, evidence=(w,))
    assert Port.from_json(p.to_json()).evidence == p.evidence
    assert json.loads(p.evidence[0].payload_json) == payload

def test_reduced_shell_adapter_preserves_all_states_and_full_witness():
    cert = {'model_scope': 'REDUCED_MODEL_ONLY', 'complete': True, 'status': 'MULTIPLE', 'states': [{'q': '-1', 'stability': 'stable'}, {'q': '0', 'stability': 'unstable'}, {'q': '1', 'stability': 'stable'}], 'supports': {'open': 'verified', 'closed': 'verified'}, 'force_unit': 'N'}

    def verify(d):
        return d['model_scope'] == 'REDUCED_MODEL_ONLY' and d['complete'] is True and ([s['q'] for s in d['states']] == ['-1', '0', '1']) and (d['supports'] == {'open': 'verified', 'closed': 'verified'})
    bs = from_shell_equilibria(cert, source=S, replay=verify)
    assert bs.status == Status.MULTIPLE
    assert len(bs.branches) == 3
    assert [json.loads(b.state_json)['stability'] for b in bs.branches] == ['stable', 'unstable', 'stable']
    for b in bs.branches:
        assert b.witnesses[0].kind == 'reduced_shell_contact_set'
        assert json.loads(b.witnesses[0].payload_json) == cert
        assert b.witnesses[0].source == S
    assert replay_witnesses(bs, {'reduced_shell_contact_set': verify})
