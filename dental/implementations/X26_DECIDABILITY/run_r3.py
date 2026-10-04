import copy
import csv
import hashlib
import importlib.util
import json
import math
import sys
import time
from pathlib import Path
import numpy as np
from force_budget import peak, upper_stress
from decidability import ScalarPort, interval_decision
ROOT = Path(__file__).resolve().parent

def dump(name, obj):
    (ROOT / name).write_text(json.dumps(obj, indent=2, ensure_ascii=False) + '\n')

def checks(rows):
    return {'source_mixed_FE_match': all((abs(r['nominal_peak_MPa'] - r['source_mixed_peak_MPa']) / max(1, r['source_mixed_peak_MPa']) < 1e-08 for r in rows)), 'simplex_upper_match': all((abs(r['curve'][0]['upper_MPa'] - r['source_simplex_upper_MPa']) / max(1, r['source_simplex_upper_MPa']) < 1e-08 for r in rows)), 'exact_convex_upper_LP_match': all((abs(q['rayleigh_LP_MPa'] - q['upper_MPa']) / max(1, q['upper_MPa']) < 1e-08 for r in rows for q in r['curve'])), 'all_witnesses_feasible': all((abs(sum(q['witness_fraction']) - 1) < 1e-10 and all((l - 1e-10 <= a <= u + 1e-10 for (a, l, u) in zip(q['witness_fraction'], q['lower'], q['upper']))) for r in rows for q in r['curve'])), 'monotone_information': all((r['curve'][i + 1]['upper_MPa'] <= r['curve'][i]['upper_MPa'] + 1e-08 * max(1, r['curve'][i]['upper_MPa']) for r in rows for i in range(len(r['curve']) - 1))), 'spectral_control_is_upper': all((q['upper_MPa'] <= q['conservative_spectral_upper_MPa'] + 1e-08 * max(1, q['upper_MPa']) for r in rows for q in r['curve'])), 'no_physical_certificate': all((q['physical_status'] == 'UNKNOWN' for r in rows for q in r['curve']))}

def load_module(cell):
    path = ROOT / 'inputs/bodytwin' / cell / (cell + '_model.py')
    spec = importlib.util.spec_from_file_location(cell, str(path))
    mod = importlib.util.module_from_spec(spec)
    sys.modules[cell] = mod
    spec.loader.exec_module(mod)
    return mod

def hand_cases():
    hemo = load_module('SURG_HEMOSTASIS')
    coll = load_module('SURG_COLLAGEN')
    h = []
    for layer in hemo.LAYERS:
        if layer.vessel_density_per_mm2 == 0:
            continue
        f = layer.upstream_fraction_f
        port = ScalarPort(math.log(2), 'log ratio', math.log(1 / f), (), coefficient=1)
        klass = port.evaluate(0)['class']
        hand = 1 if f > 0.5 else 2
        h.append({'layer': layer.name, 'f': f, 'factor_1_over_f': 1 / f, 'hand_factor2_class': hand, 'typed_class': klass, 'agreement': klass == hand, 'parameter_level': 'PHENOMENOLOGICAL', 'physical_status': 'UNKNOWN', 'binding_quantity': 'upstream resistance fraction', 'analogy': 'cement film resolves local gap; upstream escape resistance remains unobserved'})
    p = coll.Parameters()
    band = [coll.r_recruit(p, k) for k in p.kappa_slide_over_rupture_band]
    (lo, hi) = (min(band), max(band))
    classification = interval_decision(lo, hi, 2)
    c = {'hand_band': [lo, hi], 'engineering_ratio_threshold': 2, 'typed_class': 2 if classification == 'ABSTAIN' else 1, 'hand_class': 2, 'binding_quantity': 'sliding-to-rupture work ratio kappa', 'agreement': classification == 'ABSTAIN', 'physical_status': 'UNKNOWN', 'analogy': 'cement interface dissipative work versus ceramic rupture work', 'parameter_level': 'PHENOMENOLOGICAL', 'external_scope': 'source hand model, not measured physiological facit'}
    equality = {'f': 0.5, 'hand_predictable': True, 'strict_finite_class': ScalarPort(math.log(2), 'log ratio', math.log(2)).evaluate(0.01)['class'], 'difference': 'same equality convention issue as atlas v2; no finite strict error slack'}
    return {'bleeding': h, 'collagen': c, 'equality_difference': equality, 'all_nonboundary_classes_agree': all((x['agreement'] for x in h)) and c['agreement'], 'physical_hand_parameters_unmeasured': True}

def main():
    start = time.perf_counter()
    pr = json.loads((ROOT / 'PREREG_R3.json').read_text())
    assert hashlib.sha256((ROOT / 'PREREG_R3.json').read_bytes()).hexdigest() == (ROOT / 'PREREG_R3.sha256').read_text().strip()
    manifest = json.loads((ROOT / 'SOURCE_MANIFEST_R3.json').read_text())
    assert all((hashlib.sha256((ROOT / x['local']).read_bytes()).hexdigest() == x['sha256'] for x in manifest))
    parents = json.loads((ROOT / 'inputs/X18_results.json').read_text())['force_uncertainty_capability']['rows']
    rows = []
    for source in parents:
        f = next((x for x in manifest if x.get('fdi') == source['fdi']))
        basis = np.load(ROOT / f['local'])['stress_tensors_MPa']
        center = source['independent_mixed_solve']['alpha']
        mixed = peak(np.einsum('j,jea->ea', center, basis))
        curve = [upper_stress(basis, center, eps) for eps in pr['epsilon_fraction_grid']]
        passes = [q for q in curve if q['upper_MPa'] < pr['threshold_MPa']]
        record = {'case': source['case'], 'fdi': source['fdi'], 'patches': len(center), 'channels_sufficient': len(center) - 1, 'basis_sha256': f['sha256'], 'center': center, 'center_origin': 'synthetic prospective pressure; existing mixed FE test vector', 'source_mixed_peak_MPa': source['independent_mixed_solve']['peak_MPa'], 'source_simplex_upper_MPa': source['sharp_upper_peak_MPa'], 'mixed_reference_relative_error': abs(mixed - source['independent_mixed_solve']['peak_MPa']) / max(1, mixed), 'simplex_relative_error': abs(curve[0]['upper_MPa'] - source['sharp_upper_peak_MPa']) / max(1, source['sharp_upper_peak_MPa']), 'largest_tested_sufficient_fraction_error': max((q['epsilon_fraction'] for q in passes), default=None), 'nominal_peak_MPa': mixed, 'curve': curve, 'physical_status': 'UNKNOWN'}
        rows.append(record)
        print(source['fdi'], record['largest_tested_sufficient_fraction_error'], flush=True)
    gates = checks(rows)
    mutations = {}
    witnesses = {}
    bad = copy.deepcopy(rows)
    bad[0]['nominal_peak_MPa'] *= 2
    mutations['source_mixed_FE_match'] = not checks(bad)['source_mixed_FE_match']
    witnesses['source_mixed_FE_match'] = {'nominal_peak_MPa_multiplier': 2}
    bad = copy.deepcopy(rows)
    bad[0]['curve'][0]['upper_MPa'] *= 2
    mutations['simplex_upper_match'] = not checks(bad)['simplex_upper_match']
    witnesses['simplex_upper_match'] = {'unmeasured_upper_multiplier': 2}
    for (gate, changes) in {'exact_convex_upper_LP_match': {'rayleigh_LP_MPa': -1}, 'all_witnesses_feasible': {'witness_fraction': [1.1] + [0] * (rows[0]['patches'] - 1)}, 'monotone_information': {'upper_MPa': rows[0]['curve'][0]['upper_MPa'] + 1000}, 'spectral_control_is_upper': {'conservative_spectral_upper_MPa': -1}, 'no_physical_certificate': {'physical_status': 'CERTIFIED'}}.items():
        bad = copy.deepcopy(rows)
        idx = 1 if gate == 'monotone_information' else 0
        bad[0]['curve'][idx].update(changes)
        mutations[gate] = not checks(bad)[gate]
        witnesses[gate] = changes
    hc = hand_cases()
    dump('BODYTWIN_HAND_CASE_COMPARISON.json', hc)
    result = {'round': 'R3', 'claim_type': 'capability', 'outcome': 'PROSPECTIVE_FORCE_MEASUREMENT_ACCURACY_BUDGETS', 'gates': gates, 'mutation_rejections': mutations, 'mutation_injected_values': witnesses, 'rows': rows, 'bodytwin_hand_comparison': hc, 'dropout': {'source_teeth': len(parents), 'kept': len(rows), 'rejected': 0, 'fraction': 0}, 'external_referent': pr['external_referent'], 'numerical_control': 'independent Rayleigh support LP matches convex vertex upper; no algorithm superiority claimed', 'physical_sensor_measurements': 0, 'physical_stress_accuracy': 'UNKNOWN', 'mesh_convergence': 'UNKNOWN', 'cost': {'wall_seconds': time.perf_counter() - start, 'fit_seconds': 0, 'basis_preparation_bytes': sum((x['bytes'] for x in manifest)), 'inherited_FE_and_lab_cost': 'UNKNOWN', 'discovery_and_review': 'UNKNOWN'}, 'next_operation': 'real registered pressure intervals plus independent strain/load-displacement on same crown; mesh refinement with finite contact/material length'}
    dump('RESULTS_R3.json', result)
    with (ROOT / 'FORCE_MEASUREMENT_BUDGET.csv').open('w') as f:
        fields = ['case', 'fdi', 'patches', 'channels_sufficient', 'largest_tested_sufficient_fraction_error', 'nominal_peak_MPa', 'physical_status', 'basis_sha256']
        w = csv.DictWriter(f, fieldnames=fields, extrasaction='ignore')
        w.writeheader()
        w.writerows(rows)
    assert all(gates.values()) and all(mutations.values()) and hc['all_nonboundary_classes_agree']
    print('R3 all declared gates and same-validator injected-error probes pass', flush=True)
if __name__ == '__main__':
    main()
