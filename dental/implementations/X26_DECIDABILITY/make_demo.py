"""Render the reader artifact; persistent results use frozen scientific quantities."""
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT = Path(__file__).resolve().parent

def get(name):
    return json.loads((ROOT / name).read_text())

def dump(name, x):
    (ROOT / name).write_text(json.dumps(x, indent=2, ensure_ascii=False) + '\n')

def main():
    r1 = get('RESULTS_R1.json')
    r2 = get('RESULTS_R2.json')
    r3 = get('RESULTS_R3.json')
    r4 = get('RESULTS_R4.json')
    raw = get('RAW_R2_CANAL.json')
    facit = get('EXTERNAL_FACIT_CHECK.json')
    plt.rcParams.update({'font.size': 9, 'axes.titlesize': 11})
    (fig, axs) = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
    ax = axs[0, 0]
    m = np.sort(np.abs([x['margin'] for x in raw]))
    cdf = np.arange(1, len(m) + 1) / len(m)
    ax.plot(m, cdf, color='#215a86')
    floor = r2['canal_floor']['conditional_floor_mm']
    ax.axvline(floor, color='#b95730', ls='--', label=f'Conditional floor {floor:.3f} mm')
    ax.axvline(floor + np.sqrt(3) * 0.3, color='#527a50', ls=':', label='Envelope at h = 0.3 mm')
    ax.set(xlim=(0, 5), ylim=(0, 1), xlabel='Absolute archived canal margin [mm]', ylabel='Fraction of question instances', title='A. Fixed-image conditional budgets | PER_SURFACE_REGION')
    ax.legend(loc='lower right', fontsize=8)
    ax.text(0.02, 0.92, 'All physical surface certificates: UNKNOWN', transform=ax.transAxes, fontsize=8)
    ax = axs[0, 1]
    ids = {35, 37, 44, 47}
    for r in r3['rows']:
        if r['fdi'] in ids:
            ax.plot(range(len(r['curve'])), [q['upper_MPa'] for q in r['curve']], marker='o', label=f"FDI {r['fdi']}")
    ax.axhline(100, color='black', ls='--', label='Engineering screen 100 MPa')
    ax.set(xticks=range(7), xticklabels=['1', '.2', '.1', '.05', '.02', '.01', '0'], xlabel='Assumed absolute fraction error per patch', ylabel='Conditional upper tensile stress [MPa]', title='B. Different information | fixed FE model, PER_POINT')
    ax.legend(fontsize=8)
    ax.text(0.38, 0.65, 'Synthetic sensor centres', transform=ax.transAxes, fontsize=8)
    ax = axs[1, 0]
    cells = facit['rows']
    settings = [int(x['arm'].split()[0][1:]) for x in cells]
    ax.errorbar(settings, [x['source_mean_um'] for x in cells], yerr=[x['source_sd_um'] for x in cells], fmt='o-', capsize=4, label='Measured group mean ± specimen SD')
    ax.plot(settings, settings, '--', color='#999999', label='CAD value = reported gap')
    ax.set(xlabel='CAD spacer setting [µm]', ylabel='Measured dry marginal gap [µm]', title='C. Independent Laboratory reference | POPULATION')
    ax.legend(fontsize=8)
    ax.text(0.38, 0.4, 'Cureus 2023, Table 2\nDOI 10.7759/cureus.38688', transform=ax.transAxes, fontsize=8)
    ax = axs[1, 1]
    h = np.asarray(r4['h_over_ell'])
    ax.loglog(h, r4['peak_normalized'], 'o-', label='Ideal sharp-crack point stress: class 3')
    ax.axhline(2, color='#527a50', ls='--', label='Fixed-support average / σ(ell) = 2')
    ax.invert_xaxis()
    ax.set(xlabel='h / ell (refinement →)', ylabel='Stress normalized by σ(ell)', title='D. External closed form; actual crown divergence UNKNOWN')
    ax.legend(fontsize=8)
    ax.text(0.03, 0.75, 'Irwin 1957\nell must be physically measured', transform=ax.transAxes, fontsize=8)
    fig.suptitle('Decision before resolution: what changes the answer?', fontsize=15)
    fig.savefig(ROOT / 'DECIDABILITY_DEMO.png', dpi=160)
    fig.savefig(ROOT / 'DECIDABILITY_DEMO.pdf', metadata={'CreationDate': None, 'ModDate': None})
    plt.close(fig)
    file = ROOT / 'DECIDABILITY_TABLE.csv'
    with file.open() as f:
        reader = csv.DictReader(f)
        fields = reader.fieldnames
        table = list(reader)
    table = [x for x in table if not x['decision'].startswith('REFERENCE:')]
    table.append(dict(decision='REFERENCE:ideal_sharp_crack_peak', quantity='point stress as h/ell -> 0', resolution_level='PER_POINT', dominant_sigma_or_bound='sigma_peak proportional to h^(-1/2); not random sigma', **{'class': 3}, binding_quantity='physical support/process-zone/contact length and different observable', certifying_resolution_or_budget='no finite convergent point-peak budget', physical_status='UNKNOWN', claim_type='capability', timescale='SIMULTANEOUS', certifying_scope='ideal LEFM reference only; actual dental class3 UNKNOWN', external_locator='doi:10.1115/1.4011547', confidence='closed form, not empirical confidence'))
    with file.open('w') as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(table)
    x5 = get('inputs/X5_results.json')
    snapshots = get('SOURCE_MANIFEST.json') + get('SOURCE_MANIFEST_R3.json')
    first_cost_file = ROOT / 'FIRST_RUN_COST.json'
    if not first_cost_file.exists():
        dump('FIRST_RUN_COST.json', {'legacy_seconds': sum((x['seconds'] for x in get('LEGACY_REPRODUCTION.json'))), 'R2_seconds': r2['cost']['wall_seconds'], 'R3_seconds': r3['cost']['wall_seconds'], 'preparation_fit_discovery_validation_questions_fallback': {'preparation': 'selected frozen files, hashes and copies; human source selection time UNKNOWN', 'fit': 0, 'discovery': 'reading, coding and inherited FE/acquisition time UNKNOWN', 'validation': 'six legacy entrypoints, endpoint checks, mixed-FE and LP, primary table re-extraction, injected errors', 'questions': '4132 canal queries, 56 pressure-box queries and six normalized singular queries', 'fallback': 'UNKNOWN physical quantities -> named measurement; empty measurement port rejects'}, 'measured_lab_cost': 'NOT_PERFORMED', 'independent_review_cost': 'UNKNOWN'})
    counts = r2['canal_counts']
    result = {'lane': 'X26-decidability', 'claim_type': 'capability', 'review_state': 'PENDING_INDEPENDENT_REVIEW', 'capability': 'Return a strict finite conditional numerical budget or the different measurement needed for each scoped dental decision; prospective fixed-basis force accuracy is computable.', 'outcome': 'CONDITIONAL_BUDGETS_DELIVERED_PHYSICAL_CALIBRATION_UNKNOWN', 'external_referent': {'kind': 'independent_measurement', 'locator': 'doi:10.7759/cureus.38688 Table2; inputs/external_PMC10246932.xml', 'compared_quantity': 'reported marginal gap group means and SD at four CAD spacer settings; these reject equality of CAD setting and observed gap, not validate a local error bound', 'refutes_us': True}, 'external_referents': [facit['external_referent'], r3['external_referent'], r4['external_referent'], {'kind': 'published_dataset', 'locator': 'https://ditto.ing.unimore.it/toothfairy2/; SOURCE_MANIFEST.json', 'compared_quantity': 'annotated per-surface-region canal clearance margins; physical boundary error not available', 'refutes_us': False}, {'kind': 'independent_measurement', 'locator': 'doi:10.1016/j.ajodo.2011.06.021; PMID22051495', 'compared_quantity': 'buccal CBCT/dissection thickness MAE0.13mm and LoA -0.32/+0.38mm at0.3mm acquisition; does not calibrate our canal or extrema', 'refutes_us': True}, {'kind': 'published_dataset', 'locator': get('inputs/X24_results_R1.json')['external_referent']['locator'], 'compared_quantity': 'scanner-matched polymer phantom ROI responses; max heldout calibration error146.508HU_ref, not bone modulus', 'refutes_us': True}], 'legacy_reproduction': {'entrypoints': 6, 'complete_gates_match': True, 'naive_refine_misfires': 3, 'atlas_instances': 5, 'scope': 'curated per-instance branch classifier; not a generic empirical classifier'}, 'canal_conditional': {'class1': int(counts['1']), 'class2': int(counts['2']), 'query_instances': len(raw), 'represented_scans': x5['scan_count_archived'], 'two_candidate_policies': 'dependent alternatives, not independent observations', 'fixed_acquisition_mm': 0.3, 'numerical_operator': 'two-query 1-Lipschitz resampling, error <= sqrt(3)*h; not the archived solver mesh', 'floor_mm': floor, 'sigma_w_mm': r2['canal_floor']['sigma_w_mm'], 'tau_mm': r2['canal_floor']['tau_mm'], 'parameter_resolution_level': 'PHENOMENOLOGICAL', 'query_resolution_level': 'PER_SURFACE_REGION', 'conditional_scalar_coverage': 0.9973002039367398, 'physical_surface_confidence': 'UNKNOWN', 'physical_certificates': 0, 'replacement_measurement': 'same-site anatomical signed edge errors and simultaneous full-surface maximum coverage', 'dropout': r2['dropout']}, 'parent_decisions': {'covered': 8, 'structural_binding_prerequisite_class2': 8, 'physical_certificates': 0}, 'force_measurement': {'fixed_model_threshold_MPa': 100, 'threshold_type': 'PHENOMENOLOGICAL engineering screen; not measured strength', 'resolution_level': 'PER_POINT', 'measurement_resolution_level': 'PER_SURFACE_REGION', 'sensitive_teeth': [{k: r[k] for k in ['case', 'fdi', 'patches', 'channels_sufficient', 'largest_tested_sufficient_fraction_error']} for r in r3['rows'] if r['curve'][0]['upper_MPa'] >= 100], 'already_below_unmeasured_upper': sum((r['curve'][0]['upper_MPa'] < 100 for r in r3['rows'])), 'real_pressure_measurements': 0, 'centres': 'synthetic prospective source mixed-load vectors', 'physical_stress_accuracy': 'UNKNOWN', 'meshing_error': 'UNKNOWN', 'not_a_proven_minimum_channel_count': True}, 'divergence': {'actual_dental_class3': 'UNKNOWN', 'external_ideal_crack_class3': True, 'normalized_fixed_support_average': 2, 'physical_support_length': 'UNKNOWN', 'peak_and_regularized_mean_are_different_observables': True}, 'hand_cases': get('BODYTWIN_HAND_CASE_COMPARISON.json'), 'negative_results': ['No physical surface/stress/failure/biological certificate', 'Population image-label sigma does not certify a regional minimum', 'Original atlas-v2 reachable() says True at equality while n_req() returns None; strict finite classifier returns class2', 'Graph experiment dispatch rejected stale historical result hash; new governor port registered as evidence definition', 'No measured dental class3 or lab calibration gain'], 'same_information_control': 'All frozen numerical tolerances pass: endpoint arithmetic and Rayleigh LP match; no algorithm superiority claim.', 'full_cost': get('FIRST_RUN_COST.json'), 'full_cost_status': 'UNKNOWN: FIRST_RUN_COST is the first successfully collected component timing, not total project/lab/failed-attempt cost; see COST_LEDGER.json', 'source_manifests': {'R2': hashlib.sha256((ROOT / 'SOURCE_MANIFEST.json').read_bytes()).hexdigest(), 'R3': hashlib.sha256((ROOT / 'SOURCE_MANIFEST_R3.json').read_bytes()).hexdigest()}, 'frozen_predictions': {x.name: hashlib.sha256(x.read_bytes()).hexdigest() for x in sorted(ROOT.glob('FROZEN_PREDICTIONS*.json'))}, 'clinical_recommendation': False, 'next_construction': 'same-crown registered pressure intervals and independent mechanical response; same-site anatomical maximum-error reference; refine with physically measured contact/process-zone length', 'run_command': './run_all.sh'}
    dump('results.json', result)
    print('Demo figure and stable scientific results written')
if __name__ == '__main__':
    main()
