import numpy as np, time
from common import *
import legacy_drill_thermal as dt
from r1 import sindel, SD, OBS, metrics

def build_basis():
    datasets = []
    curves = []

    def add(study, group, y, e, w, debt, locator):
        datasets.append(dict(study=study, group=group, observed_rise_C=float(y), locator=locator, resolution='POPULATION', support='PER_POINT', physical_input_debt=debt))
        curves.append((np.asarray(e), np.asarray(w)))
    for irrig in [False, True]:
        edge = sindel(irrig, eta=1.0, chi=0.0)
        wall = sindel(irrig, eta=1.0, chi=1.0)
        for (i, s) in enumerate(edge['steps']):
            t = edge['t']
            i0 = np.searchsorted(t, s['t_start_s'])
            i1 = np.searchsorted(t, s['t_end_drill_s'] + 30.0)
            e = edge['probes']['tc'][i0:i1] - edge['probes']['tc'][max(0, i0 - 1)]
            w = wall['probes']['tc'][i0:i1] - wall['probes']['tc'][max(0, i0 - 1)]
            for key in ['12', '30'] if irrig else ['none']:
                add('Sindel_legacy', f'{key}_D{SD[i]}', OBS[key][i], e, w, ['feed, drilling endpoint, cortical fraction and sensor depth nominal', 'primary fulltext missing; old table extraction only'], f'PMC5549525 Table1 drill{i + 1} {key}; inherited anchors_thermal.py')
    records = load('raw/measurements.json')
    for row in records:
        if not row['calibration_use']:
            continue
        sid = row['study']
        if sid == 'PMC10173627':
            feed = 2.0
            depth = 9.0
            tc = 1.5
            irrig = True
            wall_irrig = row['guide'] == 'Non-limiting'
            T0 = 21.0
            debt = ['feed=2 mm/s and cortical thickness=1.5 mm nominal; heat source transferred from unrelated torque study']
        else:
            feed = row['feed_mm_s']
            depth = 12.0
            tc = 1.5
            irrig = False
            wall_irrig = True
            T0 = 22.0
            debt = ['cooling assumed absent, drill endpoint=12 mm, layered ovine material=K3 nominal; measured force does not identify torque or heat partition']
        outputs = []
        for chi in [0.0, 1.0]:
            r = dt.simulate([dict(D=row['diameter_mm'], d_prev=0.0, depth=depth, rpm=row['rpm'], feed=feed, irrig=irrig, Tcool=T0)], t_cort_mm=tc, T0=T0, perfusion=False, h_mm=0.1, probes=[('tc', row['diameter_mm'] / 2 + row['sensor_wall_distance_mm'], row['sensor_depth_mm'])], cool_after_s=30.0, eta=1.0, chi=chi, wall_irrig_during=wall_irrig, R_extra_mm=8.0, Z_extra_mm=6.0)
            outputs.append(r['probes']['tc'] - T0)
        add(sid, row['group'], row['value_C'], *outputs, debt, row['table_locator'])
    nmax = max((len(e) for (e, w) in curves))
    E = np.zeros((len(curves), nmax))
    W = E.copy()
    lengths = []
    for (i, (e, w)) in enumerate(curves):
        assert len(e) == len(w)
        lengths.append(len(e))
        E[i, :len(e)] = e
        W[i, :len(w)] = w
    np.savez_compressed(ROOT / 'raw/source_basis.npz', edge=E, wall=W, lengths=lengths)
    dump('raw/R2_dataset.json', datasets)
    return (datasets, E, W)

def candidate_grid(bounded):
    lim = 20 if bounded else 40
    return [(i * 0.05, j * 0.05) for i in range(lim + 1) for j in range(lim + 1) if (not bounded or i + j <= 20) and i + j > 0]

def predict(E, W, a, b):
    return np.max(a * E + b * W, axis=1)

def main():
    start = time.perf_counter()
    verify('PREREG_R2.json')
    contract = dict(legacy_code_sha256=sha(ROOT / 'code/legacy_drill_thermal.py'), builder_sha256=sha(ROOT / 'code/r2.py'), prereg_sha256=sha(ROOT / 'PREREG_R2.json'), measurements_sha256=sha(ROOT / 'raw/measurements.json'))
    valid_cache = (ROOT / 'raw/source_basis.npz').exists() and (not (ROOT / 'raw/basis_contract.json').exists() or load('raw/basis_contract.json') == contract)
    if valid_cache:
        z = np.load(ROOT / 'raw/source_basis.npz')
        E = z['edge']
        W = z['wall']
        rows = load('raw/R2_dataset.json')
    else:
        (rows, E, W) = build_basis()
    dump('raw/basis_contract.json', contract)
    studies = sorted(set((r['study'] for r in rows)))
    s = np.array([r['study'] for r in rows])
    y = np.array([r['observed_rise_C'] for r in rows])
    base = predict(E, W, 0.25, 0.25)
    fits = []
    held_predictions = np.zeros(len(y))
    ctrl_predictions = held_predictions.copy()
    grids = {b: candidate_grid(b) for b in [True, False]}
    predictions = {b: np.array([predict(E, W, *ab) for ab in grids[b]]) for b in [True, False]}
    for held in studies:
        train = s != held
        test = ~train
        weights = np.array([1.0 / np.sum(s[train] == study) for study in s[train]])
        fold = dict(held_study=held, training_studies=[v for v in studies if v != held], training_rows=int(train.sum()), held_rows=int(test.sum()))
        for bounded in [True, False]:
            pred = predictions[bounded]
            scores = np.sum(abs(np.log(np.maximum(pred[:, train], 1e-12) / y[train])) * weights, axis=1)
            k = int(np.argmin(scores))
            (a, b) = grids[bounded][k]
            (held_predictions if bounded else ctrl_predictions)[test] = pred[k, test]
            label = 'bounded_source' if bounded else 'same_information_unbounded_control'
            fold[label] = dict(a_edge=a, b_wall=b, eta=a + b, chi_wall=b / (a + b), physically_admissible=a + b <= 1 + 1e-15, training_weighted_loss=float(scores[k]))
            yf = y.copy()
            yf[test] += 10
            fault_scores = np.sum(abs(np.log(np.maximum(pred[:, train], 1e-12) / yf[train])) * weights, axis=1)
            assert np.array_equal(scores, fault_scores)
        fits.append(fold)
    payload = dict(round='R2', claim_type='information_link', prereg_sha256=sha(ROOT / 'PREREG_R2.json'), basis_sha256=sha(ROOT / 'raw/source_basis.npz'), dataset_sha256=sha(ROOT / 'raw/R2_dataset.json'), fit_description='retrospective whole-study held-out nominal-operator fit; no parameter or metric changes after scoring', folds=fits, predictions=[dict(study=r['study'], group=r['group'], baseline_C=float(p), bounded_predicted_rise_C=float(q), unbounded_control_rise_C=float(c)) for (r, p, q, c) in zip(rows, base, held_predictions, ctrl_predictions)])
    if not (ROOT / 'FROZEN_PREDICTIONS_R2.json').exists():
        freeze('FROZEN_PREDICTIONS_R2.json', payload)
    else:
        verify('FROZEN_PREDICTIONS_R2.json')
        f = load('FROZEN_PREDICTIONS_R2.json')
        f.pop('created_at_utc')
        assert f == payload, 'Frozen predictions drifted'
    measures = []
    for fold in fits:
        m = s == fold['held_study']
        b = metrics(base[m], y[m])
        q = metrics(held_predictions[m], y[m])
        c = metrics(ctrl_predictions[m], y[m])
        wrong = metrics(held_predictions[m], y[m] + 10.0)
        wrong_control = metrics(ctrl_predictions[m], y[m] + 10.0)
        measures.append(dict(**fold, baseline=b, bounded=q, unbounded_control=c, injected_plus10_scoring=wrong, injected_plus10_unbounded_control=wrong_control, frozen_transfer_gate=q['within_factor2_fraction'] >= 0.6 and q['mae_C'] <= 2.0, fault_rejected=wrong['within_factor2_fraction'] < 0.6 or wrong['mae_C'] > 2.0, unbounded_control_fault_rejected=wrong_control['within_factor2_fraction'] < 0.6 or wrong_control['mae_C'] > 2.0))
    allb = metrics(base, y)
    allq = metrics(held_predictions, y)
    primary_complete = load('raw/extraction_checks.json')['strict_complete_rows']
    outcome = dict(claim_type='information_link', folds=measures, overall_baseline=allb, overall_bounded=allq, overall_unbounded=metrics(ctrl_predictions, y), aggregate_MAE_reduction_fraction=1 - allq['mae_C'] / allb['mae_C'], transfer_gate_pass=all((d['frozen_transfer_gate'] for d in measures)) and allq['mae_C'] <= 0.7 * allb['mae_C'], physical_primary_evidence_gate_pass=False, complete_primary_rows=primary_complete, primary_studies_in_fit=2, legacy_unverified_studies_in_fit=1, protocol_by_human_bone_class='UNKNOWN', CEM43_validation='UNKNOWN: no full published same-region trajectories or local bone dose-to-injury closure', numerical_guarantee='Finite candidate grid only; rigorous continuum/float enclosure absent', faults=dict(no_heldout_leakage=True, all_injected_plus10_rejected=all((d['fault_rejected'] for d in measures))), runtime_s=time.perf_counter() - start)
    dump('attempts/R2_results.json', outcome)
    (ROOT / 'attempts/HANDOFF_R2.md').write_text('R2: bounded spatial-source repair and same-information control executed with each whole study held out. Predictions frozen before scoring. See R2_results.json. Inputs remain nominal, one legacy original absent, no complete primary physical calibration row; human bone-class protocol and CEM43 remain UNKNOWN. Next construction adds a known-power thermal calibration and retained per-specimen histories instead of another peak fit.\n')
    state('R2_DECIDED', 'LOSO scored; physical evidence gate FAIL independently of numerical fit', 'R3: known-power heater breaks thermal scale ambiguity; quantify dose sufficiency and export lab contract')
    print('R2', allq, 'gate', outcome['transfer_gate_pass'])
if __name__ == '__main__':
    main()
