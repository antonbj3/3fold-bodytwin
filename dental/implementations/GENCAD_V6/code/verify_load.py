from common import *
from loadsets import M, F, classes, edges
import copy

def run():
    native = {r['case_key']: r for r in read(ROOT / 'raw/R3_NATIVE.json')}
    rows = read(ROOT / 'raw/R3_ROWS.json')
    checked = 0
    mutants = []

    def verify(r):
        if 'classes' not in r:
            return True
        e = [v for v in native[r['case_key']]['edges'] if v['lower'] != r['source_fdi']] + r['incident_edges']
        for (delta, key) in [(0.0, 'classes'), (0.05, 'robust_classes')]:
            expected = classes(e, delta)
            if clean(expected) != r[key]:
                return False
        return True
    for r in rows:
        if 'classes' not in r:
            continue
        if not verify(r):
            raise AssertionError('Changed output ' + r['uid'])
        checked += 1
        if len(mutants) < 12:
            bad = copy.deepcopy(r)
            t = next(iter(bad['classes']))
            bad['classes'][t]['share_hull_pp'][1] = 123.0
            mutants.append(dict(uid=r['uid'], mutant='changed reported share upper to123', rejected=not verify(bad)))
    A = [dict(upper=16, lower=36, gap_bounds_mm=[0.0, 0.0]), dict(upper=26, lower=46, gap_bounds_mm=[1.0, 1.0])]
    B = [dict(upper=16, lower=36, gap_bounds_mm=[1.0, 1.0]), dict(upper=26, lower=46, gap_bounds_mm=[0.0, 0.0])]
    (ca, cb) = (classes(A), classes(B))
    identity = max((abs(x - y) for (x, y) in zip(sorted((r['gap_bounds_mm'][0] for r in A)), sorted((r['gap_bounds_mm'][0] for r in B)))))
    changes = [t for t in ca if ca[t]['classification'] != cb[t]['classification']]
    suff = dict(summary='sorted minimum gaps and number of contact edges', identity_error=identity, changed_tooth_classes=changes, A={t: ca[t] for t in changes}, B={t: cb[t] for t in changes}, minimum_extension='tooth-labelled incidence of minimum-gap edges; complete labelled gap intervals for future changes', resolution='PER_TOOTH', external_referent=dict(kind='our_own_fixture', locator='code/verify_load.py', compared_quantity='mathematical summary insufficiency; not measured tooth load', refutes_us=True))
    direct = []
    for r in [A, B]:
        for t in M.FDI:
            inc = [x for x in r if t in [x['upper'], x['lower']]]
            other = [x for x in r if t not in [x['upper'], x['lower']]]
            cls = 'NEVER' if not inc else 'MUST' if not other or min((x['gap_bounds_mm'][1] for x in inc)) <= min((x['gap_bounds_mm'][0] for x in other)) else 'CAN'
            direct.append(classes(r)[str(t)]['classification'] == cls)
    controls = read(ROOT / 'raw/R3_CONTROLS.json')
    kernel = all((c.get('unpruned_error_mm', 0) <= 1e-07 and c.get('lower_bound_error_mm', 0) <= 1e-12 and c.get('injected_plus1_rejected', True) and c.get('nominal_pass', True) and c.get('injected_rejected', True) for c in controls))
    out = dict(checked_actual_crown_outputs=checked, changed_output_mutants=mutants, direct_formula_checks=len(direct), direct_all=all(direct), kernel_and_KKT_checks=len(controls), kernel_all=kernel, sufficiency_exact=identity == 0 and len(changes) == 4, all_pass=all(direct) and kernel and all((x['rejected'] for x in mutants)) and (identity == 0) and (len(changes) == 4))
    dump(ROOT / 'raw/SUFFICIENCY_LOAD.json', suff)
    dump(ROOT / 'raw/LOAD_CONTROLS.json', out)
    if not out['all_pass']:
        raise AssertionError(out)
    return out
if __name__ == '__main__':
    print(clean(run()))
