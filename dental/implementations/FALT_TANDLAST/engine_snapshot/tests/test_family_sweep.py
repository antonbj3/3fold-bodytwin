from copy import deepcopy
from fractions import Fraction as Q
import pytest
from motion_engine.ncp.latent_contact import make_family, certify_family
from motion_engine.ncp.family_sweep import certify_sweep, verify_sweep

def fixture(clearance=Q(1, 5), coupled=True):
    pb = make_family([[2, 1], [1, 2]], [-5, -6], [7, -7], [-1, 1], step='1/100', source_id='declared-LCP-motion-fixture', source_sha256='0' * 64, unit_system='EXTERNAL_UNSPECIFIED')
    a = certify_family(pb)
    cells = []
    for c in a['cells']:
        init = [[Q(1, 4), Q(1, 4), clearance], [0, 0, 0], [1, 0, 0], [0, 1, 0]]
        base = [[list(map(Q, p)) for p in init] for _ in range(2)]
        slope = [[[Q(0), Q(0), Q(1, 10)] for _ in range(4)] for _ in range(2)]
        if coupled:
            base[1][0][2] += Q(1, 100) * c['p0'][0]
            slope[1][0][2] += Q(1, 100) * c['p1'][0]
            for j in (1, 2, 3):
                base[1][j][2] += Q(1, 100) * c['p0'][1]
                slope[1][j][2] += Q(1, 100) * c['p1'][1]
        cells.append(dict(base=base, slope=slope))
    return (pb, a, [dict(id='PT', kind='PT', skin=0, cells=cells)])

def test_all_latent_values_with_exact_differential_response():
    (pb, a, g) = fixture()
    s = certify_sweep(pb, a, g, complete=True)
    assert verify_sweep(pb, a, g, s, complete=True)
    assert s['status'] == 'SAFE_PRESCRIBED_LINEAR_PATH_FAMILY'
    assert min((c['margin'] for c in s['certificates'])) == Q(27, 200)

def test_response_intervention_changes_decision_and_is_real_crossing():
    (pb, a, g) = fixture(Q(1, 100))
    s = certify_sweep(pb, a, g, complete=True)
    assert s['status'] == 'UNKNOWN'
    c = g[0]['cells'][-1]
    z = Q(1)
    initial = c['base'][0][0][2] - c['base'][0][1][2]
    final = c['base'][1][0][2] + z * c['slope'][1][0][2] - c['base'][1][1][2] - z * c['slope'][1][1][2]
    assert initial > 0 and final < 0
    (_, _, without) = fixture(Q(1, 100), False)
    assert certify_sweep(pb, a, without, complete=True)['status'] == 'SAFE_PRESCRIBED_LINEAR_PATH_FAMILY'

@pytest.mark.parametrize('mutation', ['drop', 'duplicate', 'normal', 'margin', 'family', 'geometry', 'hash', 'coverage'])
def test_sweep_replay_rejects_one_sided_mutations(mutation):
    (pb, a, g) = fixture()
    s = certify_sweep(pb, a, g, complete=True)
    assert verify_sweep(pb, a, g, s, complete=True)
    if mutation == 'drop':
        s['certificates'].pop()
    elif mutation == 'duplicate':
        s['certificates'][1] = deepcopy(s['certificates'][0])
    elif mutation == 'normal':
        s['certificates'][0]['direction'] = (0, 0, 0)
    elif mutation == 'margin':
        s['certificates'][0]['margin'] += 1
    elif mutation == 'family':
        a['cells'][0]['p0'][0] += 1
    elif mutation == 'geometry':
        g[0]['cells'][0]['base'][1][0][2] -= 1
    elif mutation == 'hash':
        s['input_sha256'] = '1' * 64
    elif mutation == 'coverage':
        s['coverage'] = 'ALL_UNDISCOVERED_PAIRS'
    assert not verify_sweep(pb, a, g, s, complete=True)

def test_no_vacuous_or_incomplete_acceptance():
    (pb, a, g) = fixture()
    assert certify_sweep(pb, a, [], complete=True)['status'] == 'UNKNOWN'
    assert certify_sweep(pb, a, g, complete=False)['status'] == 'UNKNOWN'
    s = certify_sweep(pb, a, g, complete=True)
    s['certificates'] = []
    assert not verify_sweep(pb, a, g, s, complete=True)

def test_time_interior_collision_with_uncertain_endpoint_family_is_not_corner_distance():
    (pb, a, g) = fixture(coupled=False)
    for c in g[0]['cells']:
        c['base'][1][0][2] = -Q(1, 5)
    assert certify_sweep(pb, a, g, complete=True)['status'] == 'UNKNOWN'

def test_strict_skin_bound_and_degenerate_triangle_conservatism():
    (pb, a, g) = fixture(coupled=False)
    g[0]['skin'] = Q(1, 5)
    assert certify_sweep(pb, a, g, complete=True)['status'] == 'UNKNOWN'
    g[0]['skin'] = Q(1, 10)
    for c in g[0]['cells']:
        for t in (0, 1):
            c['base'][t][2] = c['base'][t][1][:]
            c['base'][t][3] = c['base'][t][1][:]
    s = certify_sweep(pb, a, g, complete=True)
    assert verify_sweep(pb, a, g, s, complete=True)

def test_oblique_fixed_direction_has_norm_squared_skin_license():
    (pb, a, g) = fixture(coupled=False)
    for c in g[0]['cells']:
        for t in (0, 1):
            c['base'][t][0][0] = 2
    g[0]['skin'] = Q(1, 10)
    s = certify_sweep(pb, a, g, complete=True, directions=[(1, 0, 1)])
    assert verify_sweep(pb, a, g, s, complete=True)

def test_certificate_generator_is_not_called_by_verifier(monkeypatch):
    import motion_engine.ncp.latent_contact as lc
    import motion_engine.ncp.family_sweep as fs
    (pb, a, g) = fixture()
    s = fs.certify_sweep(pb, a, g, complete=True)

    def fail(*args, **kwargs):
        raise RuntimeError('discovery forbidden during replay')
    monkeypatch.setattr(lc, 'certify_family', fail)
    monkeypatch.setattr(fs, 'certify_sweep', fail)
    assert lc.verify_family(pb, a) and fs.verify_sweep(pb, a, g, s, complete=True)
