"""Study-random-intercept calibration instrument, never a new physics model.

The first convenience subset has missing thickness and nonexchangeable setups.
Accordingly the fitted intercept is an explicitly unqualified transport stress
test; a separate compatibility gate denies design-specific calibration.
"""
import json, math, time
import numpy as np
from scipy.optimize import minimize_scalar
from .io import ROOT, dump, digest, freeze, sha
from .scoring import interval_score

def fit(train):
    y = np.log([r['mean'] for r in train])
    var = np.array([(r['sd'] / r['mean']) ** 2 / r['n'] for r in train])
    groups = np.array([r['study'] for r in train])
    one = np.ones(len(y))
    Z = (groups[:, None] == groups[None, :]).astype(float)

    def profile(tau2, full=False):
        V = np.diag(var) + tau2 * Z
        iv1 = np.linalg.solve(V, one)
        ivy = np.linalg.solve(V, y)
        info = one @ iv1
        mu = one @ ivy / info
        e = y - mu
        obj = np.linalg.slogdet(V)[1] + np.log(info) + e @ np.linalg.solve(V, e)
        return (mu, 1 / info) if full else obj
    opt = minimize_scalar(profile, bounds=(0, 25), method='bounded', options={'xatol': 1e-12})
    tau = float(opt.x) if profile(opt.x) < profile(0) else 0.0
    (mu, muvar) = profile(tau, True)
    v = tau + muvar + float(np.median(var))
    (lo, hi) = np.exp([mu - 1.95996398454 * np.sqrt(v), mu + 1.95996398454 * np.sqrt(v)])
    studym = [float(np.mean(y[groups == s])) for s in sorted(set(groups))]
    bm = float(np.mean(studym))
    bv = float(np.var(studym, ddof=1))
    bpredvar = bv * (1 + 1 / len(studym)) + float(np.median(var))
    (bl, bh) = np.exp([bm - 1.95996398454 * np.sqrt(bpredvar), bm + 1.95996398454 * np.sqrt(bpredvar)])
    return dict(model=dict(prediction_N=float(np.exp(mu)), lower_N=float(lo), upper_N=float(hi), study_variance_log=tau, mean_variance_log=float(muvar), fixed_effects=['intercept'], missing_covariates='thickness/cement unknown in some studies; no imputation'), control=dict(prediction_N=float(np.exp(bm)), lower_N=float(bl), upper_N=float(bh), study_variance_log=bv), n_training_studies=len(studym), n_training_groups=len(train))

def predict_loso():
    tic = time.perf_counter()
    source = ROOT / 'data/published_measurements.json'
    rows = json.loads(source.read_text())['rows']
    folds = []
    for heldout in sorted({r['study'] for r in rows}):
        train = [r for r in rows if r['study'] != heldout]
        targets = [{k: r[k] for k in ['row_id', 'study', 'material_product', 'specimen', 'thickness_mm', 'support', 'cement', 'load_angle_deg', 'ageing']} for r in rows if r['study'] == heldout]
        fitrow = fit(train)
        folds.append(dict(heldout=heldout, training_studies=sorted({r['study'] for r in train}), training_data_sha256=digest(train), targets=targets, **fitrow, design_specific_calibration='UNKNOWN', reason='No matching full setup across independent studies'))
    payload = dict(measurements_sha256=sha(source), policy='retrospective LOSO; entire study held out; intercept-only stress test', folds=folds, quantity='crown_fracture_force_N', interval_level=0.95)
    frozen = freeze(ROOT / 'FROZEN_PREDICTIONS.json', payload)
    dump(ROOT / 'raw/calibration_fit_cost.json', dict(wall_s=time.perf_counter() - tic))
    return frozen

def score_loso():
    frozen = json.loads((ROOT / 'FROZEN_PREDICTIONS.json').read_text())
    payload = frozen['payload']
    if digest(payload) != frozen['payload_sha256']:
        raise ValueError('Prediction freeze drift')
    source = ROOT / 'data/published_measurements.json'
    if sha(source) != payload['measurements_sha256']:
        raise ValueError('Measurement drift')
    data = json.loads(source.read_text())
    lookup = {r['row_id']: r for r in data['rows']}
    scores = []
    for f in payload['folds']:
        assert f['heldout'] not in f['training_studies']
        for target in f['targets']:
            r = lookup[target['row_id']]
            for name in ['model', 'control']:
                p = f[name]
                y = r['mean']
                scores.append(dict(row_id=r['row_id'], study=r['study'], method=name, measured_N=y, sd_N=r['sd'], prediction_N=p['prediction_N'], lower_N=p['lower_N'], upper_N=p['upper_N'], interval_score_N=interval_score(p['lower_N'], p['upper_N'], y), standardized_error=(p['prediction_N'] - y) / r['sd'], covered=p['lower_N'] <= y <= p['upper_N'], locator=r['locator']))
    summary = {}
    for name in ['model', 'control']:
        groups = []
        for s in sorted({r['study'] for r in scores}):
            rr = [r for r in scores if r['study'] == s and r['method'] == name]
            groups.append(dict(study=s, interval_score_N=float(np.mean([r['interval_score_N'] for r in rr])), coverage=float(np.mean([r['covered'] for r in rr]))))
        summary[name] = dict(by_study=groups, study_balanced_interval_score_N=float(np.mean([r['interval_score_N'] for r in groups])), study_balanced_coverage=float(np.mean([r['coverage'] for r in groups])))
    out = dict(rows=scores, summary=summary, calibration_admission='UNKNOWN_SETUP_MISMATCH', n_studies=3, n_groups=len(data['rows']), frozen_predictions_sha256=sha(ROOT / 'FROZEN_PREDICTIONS.json'), inference_limits=['Convenience extraction, not systematic selection', 'Only two training studies per fold', 'Support, veneer, cement, angle and ageing are confounded with study', 'Normal log-response closure not empirically established', 'No prediction for generated crowns; no FE validation'], external_referent=dict(kind='independent_measurement', locator='; '.join(('https://doi.org/' + s['doi'] for s in data['sources'])), compared_quantity='Held-out published group mean crown fracture force, N', refutes_us=True))
    dump(ROOT / 'raw/calibration_results.json', out)
    return out
