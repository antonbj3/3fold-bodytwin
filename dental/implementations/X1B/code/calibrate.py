"""Registered REML study random intercept and equally informed conventional control."""
import numpy as np
from scipy.linalg import cho_factor, cho_solve
from scipy.optimize import minimize
from scipy.stats import t
from common import R, read, dump, sha, now, state

def feature(r):
    return np.array([1.0, r['material'] == '4Y', r['material'] == '5Y', np.log(r['thickness_mm']), r['angle_deg'] / 30.0, r['cement'] == 'resin'], float)

def basis(X):
    (_, s, vh) = np.linalg.svd(X, full_matrices=False)
    rank = int((s > s[0] * 1e-10).sum())
    return (vh[:rank].T, rank)

def sampling(r):
    return np.log1p((r['sd_N'] / r['mean_N']) ** 2) / r['n']

def fit(rows, kind='mixed'):
    X = np.stack([feature(r) for r in rows])
    y = np.log([r['mean_N'] for r in rows])
    v = np.array([sampling(r) for r in rows])
    (B, rank) = basis(X)
    A = X @ B
    ids = sorted({r['study'] for r in rows})
    Z = np.array([[r['study'] == s for s in ids] for r in rows], float)

    def solve(z):
        (tau, sig) = np.exp(z)
        V = np.diag(v + sig) + tau * (Z @ Z.T)
        C = cho_factor(V, lower=True)
        ViA = cho_solve(C, A)
        Cbeta = np.linalg.inv(A.T @ ViA)
        beta = Cbeta @ (ViA.T @ y)
        e = y - A @ beta
        obj = 0.5 * (2 * np.log(np.diag(C[0])).sum() + np.linalg.slogdet(A.T @ ViA)[1] + e @ cho_solve(C, e) + (len(y) - rank) * np.log(2 * np.pi))
        return (obj, beta, Cbeta, tau, sig)
    if kind == 'mixed':
        op = minimize(lambda z: solve(z)[0], np.log([0.15, 0.005]), method='L-BFGS-B', bounds=[(np.log(1e-08), np.log(4))] * 2)
        (_, beta, Cbeta, tau, sig) = solve(op.x)
        success = bool(op.success)
    else:
        counts = {s: sum((r['study'] == s for r in rows)) for s in ids}
        w = np.array([1 / counts[r['study']] for r in rows])
        H = A.T @ (w[:, None] * A)
        beta = np.linalg.solve(H, A.T @ (w * y))
        e = y - A @ beta
        means = np.array([np.mean(e[[r['study'] == s for r in rows]]) for s in ids])
        tau = float(np.var(means, ddof=1)) if len(ids) > 1 else 1.0
        sig = float(np.sum(w * e ** 2) / w.sum())
        Cbeta = np.linalg.inv(H) * (tau + sig)
        success = True
    cvlog = np.mean([np.log1p((r['sd_N'] / r['mean_N']) ** 2) for r in rows])
    return dict(beta=B @ beta, cov=B @ Cbeta @ B.T, basis=B, rank=rank, tau2=float(tau), sigma2=float(sig), specimen_log_variance=float(cvlog), studies=ids, optimizer_success=success, training_row_ids=[r['row_id'] for r in rows])

def predict(model, r, n=6):
    x = feature(r)
    estimable = float(np.linalg.norm(x - model['basis'] @ (model['basis'].T @ x))) < 1e-07
    mu = float(x @ model['beta'])
    var = float(x @ model['cov'] @ x + model['tau2'] + model['sigma2'] + model['specimen_log_variance'] / n)
    q = float(t.ppf(0.95, max(1, len(model['studies']) - 1)))
    half = q * np.sqrt(max(0, var))
    (lo, hi) = np.exp(np.clip([mu - half, mu + half], -30, 30))
    return dict(point_N=float(np.exp(mu)), interval_N=[float(lo), float(hi)], interval_level=0.9, interval_kind='provisional new-study group arithmetic mean, log-normal approximation', target_n=n, estimable=estimable, upper_lower_ratio=float(hi / lo), log_standard_error=float(np.sqrt(var)))

def serial(model):
    return {k: v.tolist() if isinstance(v, np.ndarray) else v for (k, v) in model.items()}

def calibrate():
    rows = [r for r in read('raw/LITERATURE_GROUPS.json') if r['eligible_R1']]
    studies = sorted({r['study'] for r in rows})
    folds = []
    for study in studies:
        train = [r for r in rows if r['study'] != study]
        test = [r for r in rows if r['study'] == study]
        models = {kind: fit(train, kind) for kind in ['mixed', 'ordinary']}
        preds = [dict(row_id=r['row_id'], study=study, material=r['material'], thickness_mm=r['thickness_mm'], angle_deg=r['angle_deg'], locator=r['locator'], doi=r['doi'], predictions={k: predict(m, r, r['n']) for (k, m) in models.items()}) for r in test]
        fp = R / f'raw/FROZEN_LOSO_{study}.json'
        content = dict(protocol_sha256=sha(R / 'PREREG_LITERATURE_R1.json'), code_sha256=sha(__file__), train_row_ids=[r['row_id'] for r in train], held_out_study=study, predictions=preds, retrospective=True, prospective_lab_measurement=False)
        if fp.exists():
            old = read(str(fp.relative_to(R)))
            timestamp = old['frozen_utc']
            if {k: v for (k, v) in old.items() if k != 'frozen_utc'} != content:
                raise ValueError('Frozen LOSO drift; create a new preregistered version')
        else:
            timestamp = now()
            dump(str(fp.relative_to(R)), content | dict(frozen_utc=timestamp))
        fp.with_suffix('.sha256').write_text(sha(fp) + '\n')
        observed = {r['row_id']: r['mean_N'] for r in test}
        for p in preds:
            p['observed_mean_N'] = observed[p['row_id']]
        folds.append(dict(study=study, frozen_prediction_sha256=sha(fp), models={k: serial(v) for (k, v) in models.items()}, rows=preds))
    scores = {}
    for kind in ['mixed', 'ordinary']:
        per = []
        coverage = []
        valid = []
        for fold in folds:
            rr = fold['rows']
            per.append(float(np.mean([np.log(r['predictions'][kind]['point_N'] / r['observed_mean_N']) ** 2 for r in rr])))
            coverage.append(float(np.mean([r['predictions'][kind]['interval_N'][0] <= r['observed_mean_N'] <= r['predictions'][kind]['interval_N'][1] for r in rr])))
            valid.extend((r['predictions'][kind]['estimable'] for r in rr))
        scores[kind] = dict(study_balanced_log_RMSE=float(np.sqrt(np.mean(per))), study_balanced_interval_coverage=float(np.mean(coverage)), all_fold_targets_estimable=bool(all(valid)), fold_MSE=per, fold_coverage=coverage)
    full = fit(rows)
    base = read('history/X1/FROZEN_PREDICTIONS.json')
    designs = []
    for d in base['designs']:
        preds = {}
        for angle in [0, 30]:
            r = dict(material=d['material'], thickness_mm=d['geometry_parameters']['t_occ'], angle_deg=angle, cement='resin')
            p = predict(full, r)
            support = [x for x in rows if x['material'] == d['material'] and x['angle_deg'] == angle and (x['cement'] == 'resin')]
            span = [min((x['thickness_mm'] for x in support)), max((x['thickness_mm'] for x in support))] if support else None
            in_span = bool(span and span[0] <= r['thickness_mm'] <= span[1])
            p.update(validity='UNKNOWN', thickness_support_span_mm=span, in_setup_span=in_span, scope='Cross-study empirical diagnostic; exact X1 molar geometry and 18GPa support are not matched by these premolar/stylized crown measurements')
            preds[str(angle)] = p
        designs.append(dict(design=d['design'], material=d['material'], geometry_parameters=d['geometry_parameters'], gap=d['gap'], load_predictions=preds, original_crown_stl_sha256=d['crown_stl_sha256']))
    gates = dict(min_five_studies=len(studies) >= 5, LOSO_RMSE=scores['mixed']['study_balanced_log_RMSE'] <= 0.3, LOSO_coverage=scores['mixed']['study_balanced_interval_coverage'] >= 0.8, all_targets_estimable=scores['mixed']['all_fold_targets_estimable'], interval_width=all((p['upper_lower_ratio'] <= 3 for d in designs for p in d['load_predictions'].values())))
    gain = 1 - scores['mixed']['study_balanced_log_RMSE'] / scores['ordinary']['study_balanced_log_RMSE']
    gates['beats_same_information_control_10percent'] = gain >= 0.1
    frozen = dict(schema='X1b-literature-R1-diagnostic', frozen_utc=now(), measurement_status='NO_LAB_MEASUREMENT', prereg_sha256=sha(R / 'PREREG_LITERATURE_R1.json'), code_sha256=sha(__file__), table_sha256=sha(R / 'raw/LITERATURE_GROUPS.json'), designs=designs, all_gates=gates)
    fp = R / 'FROZEN_PREDICTIONS_R1_DIAGNOSTIC.json'
    if fp.exists():
        old = read(fp.name)
        frozen['frozen_utc'] = old['frozen_utc']
        if old != frozen:
            raise ValueError('R1 prediction drift')
    else:
        dump(fp.name, frozen)
    fp.with_suffix('.sha256').write_text(sha(fp) + '\n')
    result = dict(round='R1', status='PASS' if all(gates.values()) else 'FAIL_OR_UNKNOWN', n_eligible_groups=len(rows), n_studies=len(studies), studies=studies, models={'mixed': serial(full), 'ordinary': serial(fit(rows, 'ordinary'))}, scores=scores, gates=gates, relative_RMSE_improvement=gain, folds=folds, new_design_predictions=designs, external_referent=dict(kind='independent_measurement', locator=[r['doi'] + '; ' + r['locator'] for r in rows], compared_quantity='held-out primary arithmetic mean crown fracture force N', refutes_us=True), no_novelty_claim=True)
    dump('raw/CALIBRATION_R1.json', result)
    state('R1_LITERATURE_DECIDED', gates, 'R2: matched specimen family and conditional thickness contrast, remove unknowable common lab scale')
    (R / 'HANDOFF_R1.md').write_text('R1 executes primary extraction and real REML/ordinary LOSO. Only three primary studies satisfy explicit class/thickness/mean/SD/n. Other values are retained as UNKNOWN, not borrowed from preparation reductions. ' + json_string(gates) + '\nAbsolute X1 transfer is UNKNOWN; numeric diagnostic does not justify physical narrow bands. Next change: matched same-tooth thickness contrasts and one batch-scale measurement, with frozen conditional predictions.\n')
    return result

def json_string(x):
    import json
    return json.dumps(x)
if __name__ == '__main__':
    z = calibrate()
    print(z['scores'])
    print(z['gates'])
