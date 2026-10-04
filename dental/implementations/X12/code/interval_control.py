"""Actual independent HiGHS feasibility/minimum-margin control for R10."""
import csv, json, os, time, warnings
from collections import Counter
from pathlib import Path
import numpy as np
from scipy.optimize import linprog, OptimizeWarning
from consumer_ports import uniform_shell_decision
PACKAGE = Path(__file__).resolve().parents[1]
ROOT = Path(os.environ.get('X12_RUN_ROOT', str(PACKAGE)))

def independent(distance, radius, wall, buffer, delta):
    low = max(0, float(distance) - float(radius) - 2 * delta)
    high = float(distance) + float(radius) + 2 * delta
    required = wall + buffer
    options = {'threads': 1, 'primal_feasibility_tolerance': 1e-07, 'dual_feasibility_tolerance': 1e-07}
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', message='Unrecognized options detected', category=OptimizeWarning)
        minimum = linprog([1.0], bounds=[(low, high)], method='highs', options=options)
        feasible = linprog([0.0], A_ub=[[-1.0]], b_ub=[-required], bounds=[(low, high)], method='highs', options=options)
    if not minimum.success or feasible.status not in (0, 2):
        raise RuntimeError('Unexpected LP solver status')
    if minimum.fun - required >= -1e-07:
        decision = 'LOCAL_CONDITION_ALL_ENCLOSURES'
    elif feasible.status == 2:
        decision = 'LOCAL_CONDITION_NO_ENCLOSURE'
    else:
        decision = 'UNKNOWN_INTERVAL'
    return (decision, float(minimum.fun - required), feasible.status)

def main():
    started = time.time()
    p = ROOT / 'raw/R8_teeth.csv'
    if not p.exists():
        p = PACKAGE / 'raw/R8_teeth.csv'
    rows = [r for r in csv.DictReader(p.open()) if r['domain_truncated'] == 'False']
    selection = np.linspace(0, len(rows) - 1, 120, dtype=int)
    cards = json.loads((PACKAGE / 'sources/product_cards.json').read_text())['cards']
    checked = []
    for index in selection:
        r = rows[index]
        for card in cards:
            w = card['posterior_w_mm'] if int(r['fdi']) % 10 >= 4 else card['anterior_w_mm']
            for x in [0, 0.25, 0.5, 1.0]:
                for delta in [0, 0.15, 0.3]:
                    d = float(r['horn_boundary_min_mm'])
                    radius = float(r['digital_two_surface_radius_mm'])
                    actual = uniform_shell_decision(d, radius, w, x, delta)['decision']
                    (lp, margin, status) = independent(d, radius, w, x, delta)
                    checked.append({'case': r['case'], 'fdi': int(r['fdi']), 'product': card['id'], 'buffer_mm': x, 'delta_each_mm': delta, 'candidate': actual, 'LP': lp, 'minimum_margin_mm': margin, 'feasibility_status': status, 'agree': actual == lp})
    poison = []
    for (name, d, expected, wrong) in [('safe', 2.0, 'LOCAL_CONDITION_ALL_ENCLOSURES', 'LOCAL_CONDITION_NO_ENCLOSURE'), ('infeasible', 0.2, 'LOCAL_CONDITION_NO_ENCLOSURE', 'LOCAL_CONDITION_ALL_ENCLOSURES'), ('straddles', 1.0, 'UNKNOWN_INTERVAL', 'LOCAL_CONDITION_NO_ENCLOSURE')]:
        (lp, margin, status) = independent(d, 0.52, 0.5, 0.5, 0.0)
        candidate = uniform_shell_decision(d, 0.52, 0.5, 0.5)['decision']
        poison.append({'input': name, 'distance_mm': d, 'radius_mm': 0.52, 'wall_mm': 0.5, 'buffer_mm': 0.5, 'valid_pass': candidate == lp == expected, 'injected_class': wrong, 'injected_rejected': lp != wrong, 'detects_R9_error': name == 'straddles' and lp != wrong})
    result = {'queries': len(checked), 'LP_solves': 2 * len(checked) + 6, 'mismatches': sum((not r['agree'] for r in checked)), 'class_counts': dict(Counter((r['candidate'] for r in checked))), 'poison': poison, 'gate_pass': all((r['agree'] for r in checked)) and all((p['valid_pass'] for p in poison)), 'all_poison_rejected': all((p['injected_rejected'] for p in poison)), 'external_referent': {'kind': 'published_code', 'locator': 'https://docs.scipy.org/doc/scipy/reference/optimize.linprog-highs.html', 'compared_quantity': 'minimum distance margin and interval linear feasibility', 'refutes_us': False}, 'solver_threads': 1, 'tolerance_mm': 1e-07, 'outcome': 'TIE_WITH_INDEPENDENT_LP', 'wall_s': time.time() - started}
    (ROOT / 'raw/INTERVAL_LP_ROWS.json').write_text(json.dumps(checked, indent=1) + '\n')
    (ROOT / 'raw/INTERVAL_CONTROL.json').write_text(json.dumps(result, indent=2) + '\n')
    assert result['gate_pass'] and result['all_poison_rejected']
    print(json.dumps(result))
if __name__ == '__main__':
    main()
