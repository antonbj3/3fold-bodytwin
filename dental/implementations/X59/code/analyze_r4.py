"""Exact observational identifiability; physical closures remain explicit."""
from fractions import Fraction as F
from extract_sources import ROOT, save

def rank(matrix):
    a = [[F(x) for x in row] for row in matrix]
    i = 0
    for j in range(len(a[0])):
        p = next((k for k in range(i, len(a)) if a[k][j]), None)
        if p is None:
            continue
        (a[i], a[p]) = (a[p], a[i])
        pivot = a[i][j]
        a[i] = [x / pivot for x in a[i]]
        for k in range(len(a)):
            if k != i:
                v = a[k][j]
                a[k] = [x - v * y for (x, y) in zip(a[k], a[i])]
        i += 1
        if i == len(a):
            break
    return i

def run():
    rows = {'D': [1, 0, 0, 1, 0], 'W': [1, 1, 0, 1, 0], 'E': [1, 1, 1, 1, 0], 'O': [1, 1, 1, 0, 1]}
    original = [rows['D'], rows['O']]
    extended = list(rows.values())
    contrasts = {'seating': [0, 1, 0, 0, 0], 'handling': [0, 0, 1, 0, 0], 'detector_minus_CT': [0, 0, 0, -1, 1], 'absolute_true_gap': [1, 0, 0, 0, 0]}
    tests = {k: dict(original=rank(original + [v]) == rank(original), four_readings=rank(extended + [v]) == rank(extended)) for (k, v) in contrasts.items()}
    assert all((tests[k]['four_readings'] for k in ['seating', 'handling', 'detector_minus_CT']))
    assert not tests['absolute_true_gap']['four_readings']
    worlds = [[80, 20, 0, 0, 0], [80, 0, 0, 0, 20]]
    readings = [{k: sum((F(a) * F(b) for (a, b) in zip(row, w))) for (k, row) in rows.items()} for w in worlds]
    identity = max((abs(readings[0][k] - readings[1][k]) for k in ['D', 'O']))
    assert identity == 0
    expanded = [[1, 0, 0, 1, 0, 0, 0], [1, 1, 0, 0, 1, 0, 0], [1, 1, 1, 0, 0, 1, 0], [1, 1, 1, 0, 0, 0, 1]]
    variable_CT = {k: rank(expanded + [v]) == rank(expanded) for (k, v) in {'seating': [0, 1, 0, 0, 0, 0, 0], 'handling': [0, 0, 1, 0, 0, 0, 0], 'same_film_detector_minus_CT': [0, 0, 0, 0, 0, -1, 1]}.items()}
    result = dict(round='R4', claim_type='capability', outcome='FOUR_OBSERVATION_CONTRASTS_IDENTIFIABLE_CONDITIONAL_ON_CT_STATE_CALIBRATION', external_referent=dict(kind='closed_form', locator='https://ocw.mit.edu/courses/18-06sc-linear-algebra-fall-2011/62a9db9eeab190694d40afe4734068ca_MIT18_06SCF11_Ses1.10sum.pdf ; DECOMPOSITION_R4.json declared exact observation matrix', compared_quantity='contrast identifiability under declared additive observation model', refutes_us=True), resolution='PER_POINT', columns=['g', 's', 'f', 'c', 'o'], original_rank=rank(original), four_reading_rank=rank(extended), contrast_identifiability=tests, common_gauge=[1, 0, 0, -1, -1], state_specific_CT_identifiability=variable_CT, sufficiency_test=dict(kind='our_own_fixture', latent_worlds=worlds, readings=[{k: str(v) for (k, v) in r.items()} for r in readings], summary=['D', 'O'], identity_max_abs_error_um=str(identity), downstream_in_situ_CT_difference_um=str(abs(readings[0]['W'] - readings[1]['W'])), minimum_extension='Add in-situ and sectioned-film CT readings, plus independently bounded common CT state bias or state-specific phantom calibration'), physical_measurements_performed=0, physical_validation='UNKNOWN', absolute_true_gap='UNKNOWN_TRACEABLE_THICKNESS_REFERENCE_REQUIRED', executable_successor='code/lab_bridge.py freeze / score on actual paid CSV ; no synthetic data fitted as external reference', full_cost_unknown='physical lab time, radiation protocol validation, contrast-specific phantom, registered sectioning repeatability and actual calibration sample acquisition')
    save('RESULTS_R4.json', result)
    print('R4 exact rank4; absolute gap refused; state-specific CT offsets refute seating identification')
    return result
if __name__ == '__main__':
    run()
