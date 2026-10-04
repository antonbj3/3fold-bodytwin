from dental_release.paths import expand as _release_expand
from pathlib import Path
import os, json, hashlib, csv, datetime
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT = Path(os.environ.get('X33_RUN_ROOT', str(Path(__file__).resolve().parent)))
ORIGINAL = Path(__file__).resolve().parent

def read(p):
    return json.loads((ROOT / p).read_text())

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def save(p, x):
    (ROOT / p).write_text(json.dumps(x, indent=2, ensure_ascii=False, allow_nan=False) + '\n')

def main():
    r1 = read('raw/R1_RESULTS.json')
    r2 = read('raw/R2_RESULTS.json')
    r3 = read('raw/R3_RESULTS.json')
    (fig, ax) = plt.subplots(2, 3, figsize=(15, 9), layout='constrained')
    g = np.load(ROOT / 'raw/GEOMETRY_SMALL.npz')
    z = 44
    img = g['prepared'][z].astype(float)
    img[g['pulp'][z]] = 2
    ax[0, 0].imshow(img, origin='lower', extent=[-0.15, (img.shape[1] - 0.5) * 0.3, -0.15, (img.shape[0] - 0.5) * 0.3], cmap='cividis', vmin=0, vmax=2)
    path = read('raw/SCAN_PATH.json')
    p = np.asarray(path['positions_zyx_mm'])
    ax[0, 0].scatter(p[:, 2], p[:, 1], c='orangered', s=13, label='scan on exposed face')
    ax[0, 0].set(xlabel='local x (mm)', ylabel='local y (mm)', title='Published tooth/pulp masks; virtual preparation')
    ax[0, 0].legend(fontsize=8)
    ax[0, 0].text(0.02, 0.98, 'Slice z=13.2mm\nEnamel/dentin split UNKNOWN', transform=ax[0, 0].transAxes, va='top', fontsize=9, color='white')
    pitches = [r['pitch_mm'] for r in r1['refinement']]
    uniform_fine = [r['maximum_C'] for r in r3['cross_resolution_uniform_history']]
    ax[0, 1].plot(pitches, uniform_fine, 'o-', color='#254f77')
    ax[0, 1].invert_xaxis()
    ax[0, 1].axhline(5.5, color='#9c3131', ls='--')
    ax[0, 1].set(xlabel='numerical pitch (mm)', ylabel='max pulp rise (C)', title='Same uniform retained heat: 0.5J')
    ax[0, 1].text(0.04, 0.06, 'PER_POINT, fixed voxel union\nFinite differences are not error bounds', transform=ax[0, 1].transAxes, fontsize=9)
    colors = {'uniform': '#555555', 'low_peak': '#3b8053', 'high_peak': '#bb4939'}
    for name in colors:
        h = np.loadtxt(ROOT / f'raw/R2_HISTORY_{name}.csv', delimiter=',', skiprows=1)
        ax[0, 2].plot(h[:, 0], h[:, 1], label=name.replace('_', ' '), color=colors[name])
    ax[0, 2].axhline(5.5, color='#9c3131', ls='--')
    ax[0, 2].legend(fontsize=8)
    ax[0, 2].set(xlabel='time (s)', ylabel='max pulp rise (C)', title='Equal total heat; opposite decisions')
    ax[0, 2].text(0.04, 0.06, 'Simulated heat allocations\nNo calorimeter measurement occurred', transform=ax[0, 2].transAxes, fontsize=9)
    for name in colors:
        eta = np.asarray(r2['histories_eta'][name])
        ax[1, 0].plot(np.arange(42) / 6, eta * 0.25 * 1000, drawstyle='steps-mid', color=colors[name], label=name.replace('_', ' '))
    ax[1, 0].set(xlabel='pulse time (s)', ylabel='retained heat per pulse (mJ)', title='Source history is information')
    ax[1, 0].legend(fontsize=8)
    data = r1['external']['rows']
    obs = [r['reported_rise_C'] for r in data]
    pred = [r['proxy_predicted_peak_C'] for r in data]
    for (split, marker, color) in [('TRAIN', 'o', '#254f77'), ('HELD_ENERGY', 's', '#bb4939')]:
        take = [j for (j, r) in enumerate(data) if r['split'] == split]
        ax[1, 1].scatter(np.array(obs)[take], np.array(pred)[take], marker=marker, c=color, label=split.lower())
    ax[1, 1].plot([0, 3], [0, 3], '--', color='#777777')
    ax[1, 1].legend(fontsize=8, loc='lower right')
    ax[1, 1].set(xlabel='2005 reported signed rise (C)', ylabel='half-space proxy peak (C)', title='Published referent rejects source transfer')
    ax[1, 1].text(0.04, 0.97, f"Held energy MAE={r1['external']['held_proxy_MAE_C']:.3f}C; gate0.5C FAIL\nStatistic/sensor matching UNKNOWN", transform=ax[1, 1].transAxes, va='top', fontsize=9)
    for (rel, style) in [(0.0, '-'), (0.01, '--'), (0.02, ':')]:
        rows = [r for r in r3['rows'] if r['relative_measurement_error'] == rel]
        ax[1, 2].plot([r['bins'] for r in rows], [r['worst_peak_pulp_C'] for r in rows], style, marker='o', label=f'bin error {100 * rel:.0f}%')
    ax[1, 2].axhline(5.5, color='#9c3131', ls='--')
    ax[1, 2].set_xscale('log', base=2)
    ax[1, 2].legend(fontsize=8)
    ax[1, 2].set(xlabel='number of retained-heat time bins', ylabel='adverse pulp rise (C)', title='Acquisition budget on 0.3mm grid only')
    ax[1, 2].text(0.04, 0.06, 'All 3 coarse-grid passes fail finer-grid check', transform=ax[1, 2].transAxes, fontsize=9)
    fig.suptitle('X33: pulsed laser heating needs physical source information and resolved anatomy\nConstitutive-model probes; physical certification UNKNOWN', fontsize=14)
    for fmt in ['png', 'pdf']:
        fig.savefig(ROOT / f'figures/laser_pulp_demo.{fmt}', dpi=180)
    plt.close(fig)
    costs = {name: r['cost'] for (name, r) in [('R1', r1), ('R2', r2), ('R3', r3)]}
    measured_wall = sum((r['wall_s'] for r in costs.values()))
    dropout = {'discovery_set': {'records': 27, 'selected_for_matched_direct_ablation_validation': 0, 'rejected': 27, 'rejected_fraction': 1.0, 'reason': 'wrong laser/material operation, review, nonthermal outcome, or missing matched pulse-history/geometry/sensor contract; this is selection scope, not proof of literature absence'}, 'additional_primary_protocols': {'records': 2, 'retained_as_temperature_referents': 2, 'fully_matched_for_physical_prediction': 0, 'rejected_for_physical_prediction': 2, 'rejected_fraction': 1.0, 'reason': '2005 signed statistic/point sensor and 2mm stationary bovine preparation; 2008 unknown duration/path and nominal-power conflict'}, 'numerical_partitions': {'records': 24, 'rejected_on_native_grid': 21, 'native_rejection_fraction': 21 / 24, 'native_passes': 3, 'passes_rejected_by_finer_grid': 3, 'physical_certificates': 0}, 'metadata_correction': 'raw/R1_RESULTS.json has an incorrect25/27 retained-primary denominator. Both retained protocols lie outside those27; correct denominators are here. Original hashed record preserved.'}
    physical_debts = [{'quantity': 'eta_p: heat remaining after each pulse, reflection/transmission/ejection/vaporization', 'resolution_level': 'PER_SURFACE_REGION', 'status': 'UNKNOWN', 'replacement_measurement': 'same-specimen time/position resolved net retained heat; optical absorptance alone is insufficient'}, {'quantity': 'wet coverage, h(x,t), Tw(t)', 'resolution_level': 'PER_SURFACE_REGION', 'status': 'PHENOMENOLOGICAL', 'replacement_measurement': 'laser-off spray transient plus actual inlet/outlet water temperature/flow, calibrated sensor and source footprint'}, {'quantity': 'remaining dentin/enamel and ablation front', 'resolution_level': 'PER_POINT', 'status': 'PHENOMENOLOGICAL', 'replacement_measurement': 'same specimen micro-CT and pre/post scan with enamel/dentin boundary; pulse-by-pulse depth or validated ablation law'}, {'quantity': 'k(x,T), rho c(x,T)', 'resolution_level': 'PER_POINT', 'status': 'PHENOMENOLOGICAL', 'replacement_measurement': 'same specimen thermal transient with independently known net reference energy; no population constant treated as local bound'}, {'quantity': 'point temperature -> paste/sensor temperature', 'resolution_level': 'PER_SURFACE_REGION', 'status': 'UNKNOWN', 'replacement_measurement': 'sensor position, paste contact, time-response calibration and heldout transient'}, {'quantity': 'CEM43 -> pulp injury', 'resolution_level': 'PHENOMENOLOGICAL', 'status': 'UNKNOWN', 'replacement_measurement': 'independent duration-temperature injury data; neither collagen equilibrium Tm nor5.5C implies injury probability'}]
    artifacts = []
    for p in sorted((ROOT / 'raw').glob('*')):
        if p.is_file():
            artifacts.append({'path': str(p), 'bytes': p.stat().st_size, 'sha256': sha(p)})
    artifacts.append({'path': r2['matrix_path'], 'bytes': r2['matrix_bytes'], 'sha256': r2['matrix_sha256']})
    outcome = {'lane': 'X33-laser-pulp', 'claim_type': 'capability', 'review_state': 'PENDING_INDEPENDENT_REVIEW', 'scientific_admission': False, 'capability_delivered': 'Runnable pulse/time/position energy routing on published real tooth/pulp annotation geometry; identifies why resolution and total calorimetry alone fail the pulp maximum query.', 'physical_decision': {'question': 'DeltaT_pulp<5.5C', 'class': 'KLASS_2_INFORMATION_LIMITED', 'binding_quantity': 'pulse-specific net retained heat and spray/ablation/observation closure', 'dominant_sigma': None, 'sigma_status': 'UNKNOWN_NO_DISTRIBUTION_FITTED', 'resolution_level': 'PER_POINT', 'X26_conditional_class': r1['decision']['X26']['class'], 'physical_certification': 'UNKNOWN', 'certifying_resolution_mm': None, 'certifying_budget': None, 'reason': 'Numerical refinement does not bound physical source/geometry/sensor errors; R1 also fails its numerical convergence gate.'}, 'external_referent': r1['external_referent'], 'external_validation': {'status': 'UNKNOWN_UNMATCHED', 'proxy_transfer_MAE_C': r1['external']['held_proxy_MAE_C'], 'proxy_gate': 'FAIL', 'independent_primary_XML_peak_C': 0.84, 'sample_SD_C': 0.55, 'resolution_level': 'POPULATION', 'no_SD_to_local_sigma_conversion': True, 'negative_water_changes_refute_zero_coolant_drive': True}, 'R1': {'refinement': r1['refinement'], 'numerical_tests': r1['numerical_tests'], 'refinement_gate': 'FAIL', 'physical_validity': 'UNKNOWN'}, 'R2': {'same_total_heat_J': 0.5, 'low_peak_C': r2['minimum_possible_maximum_C'], 'high_peak_C': r2['maximum_possible_maximum_C'], 'uniform_peak_C': r2['uniform_peak_C'], 'resolution_level': 'PER_POINT', 'grid_mm': 0.3, 'controls': r2['gates'], 'synthetic_heat_input': True, 'external_referent': r2['external_referent']}, 'R3': {'native_min_bins_only': r3['minimum_bins_for_native_grid_only'], 'native_acceptances': r3['coarse_acceptances'], 'native_acceptances_refuted_by_finer_history': r3['coarse_acceptances_refuted_by_finer_feasible_history'], 'uniform_cross_resolution': r3['cross_resolution_uniform_history'], 'resolution_consistency_gate': 'FAIL', 'external_referent': r3['external_referent']}, 'uncertainty_debts': physical_debts, 'dropout': dropout, 'cost': {'measured_successful_round_wall_s_sum': measured_wall, 'rounds': costs, 'maxRSS_KiB': max((c['maxRSS_KiB'] for c in costs.values())), 'threads_max': 4, 'GPU': False, 'preparation_fit_discovery_validation_queries_fallback': 'PREREG_R1/R2/R3 and COMMANDS.md separate each operation. Timings include numerical fit/validation/LP/replay, not agent reasoning.', 'failed_R2_v1_cost': 'UNKNOWN_UNINSTRUMENTED; serialization failure and exact code preserved', 'total_research_cost': 'UNKNOWN'}, 'control_outcomes': {'slab_and_rectangular_pulse': 'PASS; deliberate flux/depth error rejects', 'energy_balance': 'PASS; extra energy rejects', 'primary_XML_value': 'PASS; +1C rejects', 'source_interval_X26': 'ABSTAIN, conditional class2; erroneous below label rejects', 'pulse_matrix_and_LP': 'PASS; no algorithm novelty claimed', 'resolution': 'FAIL; coarse acceptance is rejected by fine feasible history'}, 'graph_status': 'DISPATCH_FAILED_STALE_INPUT; feedback pending independent review and receipt currently unavailable', 'artifacts': artifacts, 'limitations': ['No actual heat/calorimeter measurement performed', 'No same-specimen thermal validation', 'Fixed virtual preparation; no predictive pulse-by-pulse ablation', 'Constant-property model invalid for extreme surface temperature', 'Finite numerical refinements are not certified continuum error bounds', 'Anatomical segmentation error unknown', '5.5C reporting threshold and CEM43 are not calibrated injury endpoints'], 'next_construction': 'Same-specimen pulsed source/pause/spray protocol with known retained heat history, ablation geometry and independent local sensor transient; derive a bounded continuum/model error before any physical certificate.'}
    save('results.json', outcome)
    for p in ['FROZEN_PREDICTIONS.json', 'FROZEN_PREDICTIONS_R3.json']:
        if (ROOT / p).exists():
            (ROOT / (p + '.sha256')).write_text(sha(ROOT / p) + '  ' + p + '\n')
    table = f"""| Construction | Frozen test | Outcome |
|---|---|---|
| R1 pulsed heat in real voxel geometry | Energy1e-7; refinement0.15 degrees C; external proxy MAE0.5 degrees C | Energy/slab PASS; refinement FAIL; proxy0.602 degrees C FAIL; physical reference UNKNOWN |
| R2 one total heat measurement | Same0.5 J, different pulse histories | Maximum4.076-6.938 degrees C: decision can reverse |
| R3 more time-resolved heat measurement | 1-42 bins; 0-2% measurement error; finer-grid countertest | Only 42 bins among eight tested partitions pass 0.3 mm model; all 3 approvals fail refinement |
"""
    explanation = _release_expand(f"""Pulp temperature is now computable as a pulse-, location- and time-bound response in a real published tooth mask. The experiment establishes that **contactless laser cannot obtain a physical CLASS1 certificate by buying numerical resolution**. CLASS2 binds through UNKNOWN residual heat per pulse, cooling, ablation, local anatomy and sensor mapping. Physical temperature prediction and certifying resolution/budget UNKNOWN.

{table}
On same0.3 mm domain and same0.5 J, three admissible simulated heat histories give maxima {r2['uniform_peak_C']:.3f}, {r2['minimum_possible_maximum_C']:.3f} and {r2['maximum_possible_maximum_C']:.3f} degrees C, PER_POINT. Total calorimetry alone loses decision-relevant information. Frozen uniform history gives {uniform_fine[0]:.3f}/{uniform_fine[1]:.3f}/{uniform_fine[2]:.3f} degrees C at0.3/0.15/0.075 mm. Finer grids share same anatomical voxel union; no new anatomical measurements. 0.5 J and pulse distributions are scenarios, not laboratory measurements.

External reference: [Geraldo-Martins2005](https://doi.org/10.1089/pho.2005.23.182), p184 Table2 reports positive dry and four negative wet temperature changes at2 mm residual dentin. A source with only positive laser heating and cooling water at initial temperature cannot produce negative values. Dry half-space proxy with shared heat fraction fails held-out350 mJ; signed measured quantity and model peak also unmatched. Original PDF download returned403; primary author-uploaded full-text table manually transcribed, requiring original review.

Separately stored local [2008 study](https://doi.org/10.1590/S1678-77572008000300009), Table1 Er:YAG states explicit maximum sensor temperature increase0.84+/-0.55 degrees C, POPULATION. XML control rereads value and rejects+1 degree C. Exact laser trajectory/duration absent; stated3.5 W differs from250 mJ times4 Hz=1 W. Sensor in thermal paste is another observation than maximum over all pulp voxels. No matched reference for our human tooth/local maximum. SD is not used as local sigma.

Controls: analytical slab solution and rectangular pulse quadrature pass; energy conserved; independent scalar replays and HiGHS match pulse matrix responses. Flux, depth, energy, table value, pulse order and refinement decision faults rejected by respective controls. No algorithmic superiority.

**What fails:** pulse sequence does not yet change geometry; virtual0.5 mm preparation fixed. Separate dentin/enamel, physical segmentation uncertainty, temperature-dependent material properties, perfusion and MEASURED wet surface absent. Large eta=1 responses are linear coefficients, not physical forecasts: extreme surface temperatures violate material validity. Finite refinement differences do not prove a numerical error band. CEM43 includes time but has no validated pulp injury reference; equilibrium Tm is not an injury criterion. No clinical recommendation.

Run `./run_all.sh`. Recomputes all three constructions in a new replay directory and compares frozen first outcomes. [Figure](figures/laser_pulp_demo.png), [raw data](raw/R1_RESULTS.json) and [results](results.json) state levels, gates and uncertainties. Successful first numeric runs total {measured_wall:.2f} s; maxRSS {max((c['maxRSS_KiB'] for c in costs.values())) / 1024:.1f} MiB, at most4 threads, no GPU. Full research cost/failed R2v1 time UNKNOWN. Response matrix {r2['matrix_bytes'] / 1000000.0:.1f} MB under @DENTAL_EXTERNAL_ROOT@/storage, hash-bound in results.json. No large archive extraction.

Next construction: same-specimen experiment with known pulse-specific net energy flux, separate spray history, MEASURED ablation and independent local temperature trajectory. [MEASUREMENT_SPEC.md](MEASUREMENT_SPEC.md) separates proposed measurements from our insufficient model. All results PENDING_INDEPENDENT_REVIEW. Graph rank/packet/define dispatch refused by STALE_INPUT; no generated graph rebuilt.
""")
    (ROOT / 'RESULTS.md').write_text('# Pulse-bound heat and limits of resolution certification\n\n' + explanation)
    demo = f"""# Pulsed laser and pulp temperature — executable research demo

Follow residual energy from each Er:YAG pulse through actual tooth/pulp geometry and query largest local temperature increase. Contact force disappears, but ablation residual energy must be determined. Query uses5.5 degrees C as reporting threshold.

```bash
cd {ORIGINAL}
./run_all.sh
```

Command simulates42 pulses on virtually prepared X12 P1/FDI36 tooth, refines heat solution, runs external observation/model gates, solves histories with identical total energy and tests1,2,3,6,7,14,21,42 measurement bins. New table, figure and replay check written to replay/<UTC>/. First outcomes/preregistrations retained. Python3 with local NumPy/SciPy/Numba required; figure uses system Matplotlib. No installation/GPU/dataset download.

{table}
Useful capability rejects temperature certificates losing pulse history or accepting coarse grids. Physical certification CLASS2/UNKNOWN. Same0.5 J gives4.08-6.94 degrees C in frozen0.3 mm model, so total energy alone has no unique answer. Finest tested uniform history gives6.36 degrees C while coarse grid gives5.15 degrees C. Simulations with explicit closures, never an actual laser safety statement.

References: published Er:YAG measurements DOI10.1089/pho.2005.23.182 p184 Table2 and DOI10.1590/S1678-77572008000300009 Table1. They reject simplified cooling and reveal remaining observation gap. Temperature references are bovine lab protocols at POPULATION level; anatomical mask and maximum simulated pulp field are different inputs. Reference does not calibrate our tooth. [RESULTS.md](RESULTS.md) gives full comparison and failed transfer.

What fails: known net energy per pulse, changing ablation geometry, dentin/enamel labels, local material fields and mapped actual sensor absent. 42 bins/1% precision prospective model contract, not MEASURED instrument capability. Three finite grids give no rigorous continuum error band. No physical measurement performed; injected faults test software. [Measurement contract](MEASUREMENT_SPEC.md) and [frozen predictions](FROZEN_PREDICTIONS.json) make next test reviewable.

Data/license: same local published ToothFairy2/Pulpy3D CT pair as X12. ToothFairy2 CC BY-SA4.0; original ToothFairy data CC BY-SA. Explicit separate license for Pulpy3D pulp extension UNKNOWN in reviewed primary sources; no dataset distributed. X12 FRAME_CONTRACT/SHA256 bind0.3 mm scale and coordinates. Local 2008 XML states CC BY-NC3.0; 2005 author upload has no verified reuse terms and only derived table values used. No Bits2Bites, Teeth3DS or mandibular defect data used. No personal data.

Budget: about45 s successful first runs; under0.6 GiB MEASURED RSS, at most4 threads. 93 MB pulse matrix needed for all local maximum queries stored on data disk with hash. Replay produces a new matrix; lane refuses if3 GB limit at risk. Full research preparation time/cost UNKNOWN. Status PENDING_INDEPENDENT_REVIEW.
"""
    (ROOT / 'README_DEMO.md').write_text(demo)
    with (ROOT / 'DECIDABILITY_TABLE.csv').open('w') as f:
        w = csv.DictWriter(f, fieldnames=['quantity', 'resolution_level', 'class', 'dominant_sigma', 'binding_quantity', 'certifying_resolution_mm', 'physical_status'])
        w.writeheader()
        w.writerow({'quantity': 'max_t,x_in_pulp DeltaT<5.5C', 'resolution_level': 'PER_POINT', 'class': 2, 'dominant_sigma': 'UNKNOWN; retained heat is an interval, not SD', 'binding_quantity': 'pulse-specific net retained heat; spray and moving ablation; local anatomy/sensor', 'certifying_resolution_mm': 'UNKNOWN', 'physical_status': 'UNKNOWN'})
    print(json.dumps({'report_root': str(ROOT), 'outcome': 'KLASS_2_PHYSICAL_CERTIFICATE_UNKNOWN', 'figure': 'figures/laser_pulp_demo.png'}), flush=True)
if __name__ == '__main__':
    main()
