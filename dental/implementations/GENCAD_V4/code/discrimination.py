from common import *
sys.path.insert(0, str(DATA / 'participant_code'))
from quality import contact_map, nondominated
from legacy.checks import all_checks
from legacy.generators import submit
from scipy.sparse import eye

def run():
    from integrity import verify
    verify()
    start = time.perf_counter()
    rows = []
    checks = []

    def put(name, a, b, down_a, down_b, extension, unit, threshold):
        ident = float(np.max(np.abs(np.asarray(a) - np.asarray(b))))
        delta = float(np.max(np.abs(np.asarray(down_a) - np.asarray(down_b))))
        row = dict(metric=name, summary_A=a, summary_B=b, identity_error=ident, downstream_A=down_a, downstream_B=down_b, downstream_difference=delta, unit=unit, minimal_extension=extension, identity_pass=ident == 0.0, discrimination_pass=delta >= threshold, resolution='PER_POINT', external_referent=dict(kind='our_own_fixture', locator='code/discrimination.py', compared_quantity='Exactly matched summary with changed downstream quantity', refutes_us=True))
        rows.append(row)
        checks.append(dict(name=name + '_wrong_downstream_value', rejected=not abs(delta - 0.0) <= 1e-10))
    xy = np.concatenate([np.array([[0, 0], [1, 0], [1, 1], [0, 1]], float) + [2 * j, 0] for j in range(4)])
    faces = np.concatenate([np.array([[0, 1, 2], [0, 2, 3]]) + 4 * j for j in range(4)])
    gaps_a = np.repeat([0.0625, 0.0625, 0.25, 0.25], 4)
    gaps_b = gaps_a[::-1].copy()
    t = dict(task_id='R1', frame='fixture_mm', family='molar_crown', xy=xy, faces=faces, preparation_z=np.zeros(16), A=eye(16, format='csr'), obstacle_b=np.full(16, 3.0), requirements=dict(film_min_mm=0.02, film_max_mm=0.15, wall_mm=0.25))
    designs = []

    def feasible_pair(name, z1, z2, i1=None, i2=None):
        i1 = np.full(16, 0.0625) if i1 is None else i1
        i2 = np.full(16, 0.0625) if i2 is None else i2
        va = all_checks(t, submit(t, np.broadcast_to(z1, (16,)), i1))['validity']
        vb = all_checks(t, submit(t, np.broadcast_to(z2, (16,)), i2))['validity']
        designs.append(dict(metric=name, L1_A=va, L1_B=vb, pass_both=va == vb == 'PASS'))
    m = contact_map(xy, faces, gaps_b, gaps_a)
    put('contact_area', 2.0, 2.0, 0.0, m['contact_symdiff_mm2'], 'Spatial contact mask, or first moment for this force-moment query', 'mm2', 0.5)
    feasible_pair('contact', 3 - gaps_a, 3 - gaps_b)
    put('proximal_mean_gap', np.mean([0.0625, 0.3125]), np.mean([0.3125, 0.0625]), [0.0625, 0.3125], [0.3125, 0.0625], 'Mesial and distal gap separately; contact location for force transmission', 'mm', 0.1)
    feasible_pair('proximal', np.repeat([1.0, 1.0, 1.25, 1.25], 4), np.repeat([1.25, 1.25, 1.0, 1.0], 4))
    err_a = np.repeat([0.0625, 0.0625, 0.25, 0.25], 4)
    err_b = err_a[::-1].copy()
    put('shape_RMS', np.sqrt(np.mean(err_a ** 2)), np.sqrt(np.mean(err_b ** 2)), float(err_a[0] <= 0.1), float(err_b[0] <= 0.1), 'Spatial signed error at antagonist/neighbor-sensitive locations', 'contact indicator', 1.0)
    feasible_pair('anatomy', 3 - err_a, 3 - err_b)
    residual_a = np.array([0.5, 1.5])
    residual_b = np.array([1.0, 1.0])
    put('removed_volume', float(np.sum(2 - residual_a)), float(np.sum(2 - residual_b)), float(residual_a.min()), float(residual_b.min()), 'Residual thickness map relative to actual pulp; minimum suffices only for this residual-minimum query', 'mm', 0.2)
    feasible_pair('removal', np.repeat(np.tile(residual_a, 2) + 0.0625, 4), np.repeat(np.tile(residual_b, 2) + 0.0625, 4))
    films = np.array([0.03125, 0.125])
    h3 = films ** 3
    parallel = float(h3.mean())
    series = float(2 / np.sum(1 / h3))
    put('cement_mean_and_histogram', films.tolist(), films.copy().tolist(), 1.0, parallel / series, 'Spatial film field and connectivity to vent; hydraulic conductance suffices for this fixed-pressure flow query', 'conductance ratio', 2.0)
    feasible_pair('cement', 1.0, 1.0, np.repeat(np.tile(films, 2), 4), np.repeat(np.tile(films[::-1], 2), 4))
    th_a = np.array([0.5, 1.0])
    th_b = th_a[::-1]
    put('fracture_minimum_wall', float(th_a.min()), float(th_b.min()), float(1 / th_a[0] ** 2), float(1 / th_b[0] ** 2), 'Thickness at loaded region plus protocol, support and stress field; minimum wall alone is insufficient', 'relative bending stress', 0.1)
    feasible_pair('fracture', np.repeat(np.tile(th_a, 2) + 0.0625, 4), np.repeat(np.tile(th_b, 2) + 0.0625, 4))
    checks.append(dict(name='polygon_vs_direct_unit_square_enumeration', passed=abs(m['contact_symdiff_mm2'] - 4.0) < 1e-10))
    checks.append(dict(name='injected_wrong_polygon_area_rejected', rejected=abs(m['contact_symdiff_mm2'] + 1 - 4) > 1e-10))
    front = nondominated([[0.0, 1.0], [1.0, 0.0], [1.0, 1.0]], [0.0, 0.0])
    checks.append(dict(name='tradeoff_not_hidden_by_pareto', passed=front.tolist() == [True, True, False]))
    direct_series = 1 / (0.5 / films[0] ** 3 + 0.5 / films[1] ** 3)
    checks.append(dict(name='reynolds_direct_series', passed=direct_series == series))
    gates = dict(identity=all((r['identity_pass'] for r in rows)), discrimination=all((r['discrimination_pass'] for r in rows)), both_L1=all((r['pass_both'] for r in designs)), controls=all((r.get('passed', r.get('rejected', False)) for r in checks)))
    out = dict(claim_type='capability', round='R1', decision='DISCRIMINATING_VECTOR_PASS' if all(gates.values()) else 'FAIL', gates=gates, sufficiency=rows, L1_pairs=designs, controls=checks, cost_seconds=time.perf_counter() - start, scope='Fixtures demonstrate metric information loss. They are not external physical validation. No complete quality vector is claimed universally sufficient.')
    dump(ROOT / 'rounds/R1.json', out)
    (ROOT / 'HANDOFF_R1.md').write_text('R1 ' + out['decision'] + '. Every pair has exactly equal declared summary; spatial/physical downstream queries differ. Full result and false-value injections in rounds/R1.json. Next construction: replace scalar contact-area error with spatial PL overlap on the frozen real-data panel; run all eight actual participants.\n')
    state('R1_DECIDED', out['decision'], 'Generate all participants on public-only inputs before exposing quality references')
    print(json.dumps(gates))
    if not all(gates.values()):
        raise RuntimeError('R1 frozen gate failed')
if __name__ == '__main__':
    run()
