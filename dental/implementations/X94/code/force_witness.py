from dental_release.paths import expand as _release_expand
import sys, json, hashlib, time
from fractions import Fraction as Q
from pathlib import Path
from geometry import R, PARENT, write, sha
sys.dont_write_bytecode = True
sys.path.insert(0, str(PARENT / 'code'))
from model import edges_from_x21, exact_point, replay_point

def stringify(x):
    if isinstance(x, Q):
        return str(x)
    if isinstance(x, dict):
        return {str(k): stringify(v) for (k, v) in x.items()}
    if isinstance(x, (list, tuple)):
        return list(map(stringify, x))
    return x

def main():
    start = time.perf_counter()
    (e, p) = edges_from_x21(_release_expand('@DENTAL_CASE_ID@'))
    g = [(v[2] + v[3]) / 2 for v in e]
    teeth = sorted({t for edge in e for t in edge[:2]})
    T = Q(205)
    states = []
    failed = []
    for (name, c) in [('equal750', {t: Q(1, 750) for t in teeth}), ('stiff_anterior_soft_posterior', {t: Q(1, 750) if t % 10 <= 3 else Q(1) for t in teeth})]:
        try:
            a = exact_point(e, g, T, c)
            assert replay_point(e, stringify(a))
            states.append((name, a))
        except Exception as err:
            failed.append(dict(name=name, error=repr(err)))
    if len(states) != 2:
        write('raw/FORCE_WITNESS_FAILED.json', dict(failures=failed))
        raise ArithmeticError('Witness discovery failed; preserve failure')
    (A, B) = (states[0][1], states[1][1])
    summary = lambda a: dict(source_sha256=sha(p), gaps=list(map(str, a['gap_mm'])), pose='unchanged_X21_source', total_N=str(a['total_N']))
    (s1, s2) = (summary(A), summary(B))
    assert s1 == s2
    ant = lambda a, j: sum((a['tooth_force_N'].get(t, Q(0)) for t in teeth if (t < 30) == (j == 'upper') and t % 10 <= 3)) / T
    diff = {j: float(abs(ant(A, j) - ant(B, j))) for j in ['upper', 'lower']}
    bad = stringify(A)
    bad['edge_force_N'][0] = str(Q(bad['edge_force_N'][0]) + 1)
    rejected = not replay_point(e, bad)
    out = dict(case=_release_expand('@DENTAL_CASE_ID@'), geometry_summary_identical=True, identity_error_mm=0.0, total_load_identity_error_N=0.0, summary_sha256=hashlib.sha256(json.dumps(s1, sort_keys=True).encode()).hexdigest(), anterior_share_difference=diff, states=[dict(name=n, anterior_share={j: float(ant(a, j)) for j in ['upper', 'lower']}, exact_KKT_replay=True, answer=stringify(a)) for (n, a) in states], smallest_extension='For this regional output: measured anterior force share at the same pose/load/time. Individual tooth output requires additional resolved force constraints; geometry or MUST classes do not supply them', status='GEOMETRY_ONLY_IDENTIFIABILITY_REFUTED_IN_DECLARED_MODEL', support_status='PHENOMENOLOGICAL counterexamples, not physiological stiffness measurements', external_referent=dict(kind='closed_form', locator=str(PARENT / 'QUANTITY_LICENSE.md'), compared_quantity='KKT balance/complementarity of a declared convex contact model', refutes_us=True), injected_plus1N_rejected=rejected, wall_s=time.perf_counter() - start)
    write('raw/FORCE_WITNESS.json', out)
    print(json.dumps({k: v for (k, v) in out.items() if k != 'states'}))
if __name__ == '__main__':
    main()
