from common import *
from data import reports, treatment_rows, load_geometry
import numpy as np, pickle, collections, sys
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.pipeline import make_pipeline
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression

def matrix(rows, g, cols):
    return np.array([[g[r['case_id']][r['fdi']].get(k, np.nan) for k in cols] for r in rows])

def metrics(rows, pred):
    y = np.array([r['y'] for r in rows], int)
    p = np.array(pred, int)
    tp = int(np.sum((y == 1) & (p == 1)))
    fn = int(np.sum((y == 1) & (p == 0)))
    tn = int(np.sum((y == 0) & (p == 0)))
    fp = int(np.sum((y == 0) & (p == 1)))
    se = tp / (tp + fn) if tp + fn else None
    sp = tn / (tn + fp) if tn + fp else None
    return dict(n=len(rows), positive_teeth=tp + fn, negative_teeth=tn + fp, positive_patients=len({r['case_id'] for r in rows if r['y'] == 1}), negative_patients=len({r['case_id'] for r in rows if r['y'] == 0}), confusion=dict(TN=tn, FP=fp, FN=fn, TP=tp), sensitivity=se, specificity=sp, balanced_accuracy=(se + sp) / 2 if se is not None and sp is not None else None, resolution='POPULATION over PER_TOOTH observations')

def threshold(rows, scores):
    candidates = sorted(set([0.0, 1.0000001] + list(scores)))
    return max(candidates, key=lambda t: (metrics(rows, np.array(scores) >= t)['balanced_accuracy'], metrics(rows, np.array(scores) >= t)['specificity'], t))

def fit():
    st = time.perf_counter()
    cpu = time.process_time()
    g = load_geometry()
    m = read(X7 / 'raw/DATA_MANIFEST.json')
    train = read('raw/R1_train_ROWS.json')
    cal = read('raw/R1_calibration_ROWS.json')
    cols = sorted({k for c in g.values() for row in c.values() for k in row})
    xx = matrix(train, g, cols)
    y = np.array([r['y'] for r in train])
    hgb = HistGradientBoostingClassifier(max_iter=150, max_depth=3, learning_rate=0.05, min_samples_leaf=20, l2_regularization=2, random_state=79)
    logit = make_pipeline(SimpleImputer(strategy='median', add_indicator=True), StandardScaler(), LogisticRegression(C=0.2, max_iter=2000, random_state=79))
    hgb.fit(xx, y)
    logit.fit(xx, y)
    priors = {pos: float(np.mean([r['y'] for r in train if r['fdi'] % 10 == pos])) for pos in range(1, 9) if any((r['fdi'] % 10 == pos for r in train))}
    models = {'shape': hgb, 'logistic': logit}
    thresholds = {}
    calx = matrix(cal, g, cols)
    for (name, model) in models.items():
        thresholds[name] = float(threshold(cal, model.predict_proba(calx)[:, 1]))
    thresholds['tooth_position_prior'] = float(threshold(cal, [priors.get(r['fdi'] % 10, float(y.mean())) for r in cal]))
    candidate = [dict(case_id=c, fdi=f) for c in m['test'] if c in g for f in sorted(g[c])]
    testx = matrix(candidate, g, cols)
    scores = {name: model.predict_proba(testx)[:, 1] for (name, model) in models.items()}
    scores['tooth_position_prior'] = np.array([priors.get(r['fdi'] % 10, float(y.mean())) for r in candidate])
    preds = [dict(**r, scores={n: float(s[i]) for (n, s) in scores.items()}, predicted={n: int(s[i] >= thresholds[n]) for (n, s) in scores.items()}) for (i, r) in enumerate(candidate)]
    write('raw/R1_PREDICTIONS.json', preds)
    write('raw/R1_CALIBRATION.json', dict(thresholds=thresholds, tooth_position_priors=priors, training=metrics(train, [0] * len(train)), calibration={name: metrics(cal, models[name].predict_proba(calx)[:, 1] >= thresholds[name] if name in models else np.array([priors.get(r['fdi'] % 10, float(y.mean())) for r in cal]) >= thresholds[name]) for name in thresholds}))
    DATA.mkdir(parents=True, exist_ok=True)
    modelp = DATA / 'restoration_models.pkl'
    modelp.write_bytes(pickle.dumps(dict(columns=cols, models=models, thresholds=thresholds, prior=priors)))
    freeze('FROZEN_PREDICTIONS_R1.json', dict(frozen_utc=now(), phase='Predictions on every test tooth before numerical test-label parsing and comparison; test text keyword feasibility was exposed, older targets already analyzed', predictions_path=str(ROOT / 'raw/R1_PREDICTIONS.json'), predictions_sha256=sha(ROOT / 'raw/R1_PREDICTIONS.json'), model_path=str(modelp), model_sha256=sha(modelp), parser_sha256=sha(ROOT / 'code/report_parser.py'), geometry_sha256=sha(ROOT / 'raw/TOOTH_FEATURES.json'), prereg_sha256=sha(ROOT / 'PREREG_R1_RESTORATION_GEOMETRY.json'), split_sha256=sha(X7 / 'raw/DATA_MANIFEST.json')))
    write('raw/R1_FIT_COST.json', cost(st, cpu))
    state('R1_PREDICTIONS_FROZEN', 'NOT_EVALUATED', 'Read test photo reports and source-gold after immutable prediction freeze')
    print('R1 frozen', flush=True)

def bootstrap(rows, model, draws=500):
    rng = np.random.default_rng(79)
    groups = collections.defaultdict(list)
    for r in rows:
        groups[r['case_id']].append(r)
    ids = list(groups)
    vals = []
    for _ in range(draws):
        sample = [r for c in rng.choice(ids, len(ids), replace=True) for r in groups[c]]
        mm = metrics(sample, [r['predicted'][model] for r in sample])
        if mm['balanced_accuracy'] is not None:
            vals.append([mm['sensitivity'], mm['specificity'], mm['balanced_accuracy']])
    return dict(sensitivity95=np.quantile(np.array(vals)[:, 0], [0.025, 0.975]).tolist(), specificity95=np.quantile(np.array(vals)[:, 1], [0.025, 0.975]).tolist(), BA95=np.quantile(np.array(vals)[:, 2], [0.025, 0.975]).tolist(), draws=draws, unit='patient-cluster bootstrap; descriptive, annotation/site bias excluded') if vals else None

def evaluate():
    st = time.perf_counter()
    cpu = time.process_time()
    fr = read('FROZEN_PREDICTIONS_R1.json')
    for (p, key) in [(fr['predictions_path'], 'predictions_sha256'), (ROOT / 'code/report_parser.py', 'parser_sha256'), (ROOT / 'raw/TOOTH_FEATURES.json', 'geometry_sha256')]:
        assert sha(p) == fr[key]
    m = read(X7 / 'raw/DATA_MANIFEST.json')
    g = load_geometry()
    rr = reports(m['test'])
    (rows, ex) = treatment_rows(rr, g)
    pred = {(r['case_id'], r['fdi']): r for r in read(fr['predictions_path'])}
    for r in rows:
        r.update(pred[r['case_id'], r['fdi']])
    write('raw/R1_TEST_OBSERVATIONS.json', rr)
    write('raw/R1_TEST_ROWS.json', rows)
    write('raw/R1_TEST_EXCLUSIONS.json', ex)
    results = {}
    for name in ['shape', 'logistic', 'tooth_position_prior']:
        results[name] = metrics(rows, [r['predicted'][name] for r in rows])
        results[name]['bootstrap'] = bootstrap(rows, name)
    qa = read('raw/PARSER_GOLD_VALIDATION.json') if (ROOT / 'raw/PARSER_GOLD_VALIDATION.json').exists() else None
    p = read('PREREG_R1_RESTORATION_GEOMETRY.json')
    mm = results['shape']
    lift = mm['balanced_accuracy'] - results['tooth_position_prior']['balanced_accuracy']
    support = mm['positive_patients'] >= 20 and mm['negative_patients'] >= 20
    numeric = mm['balanced_accuracy'] >= 0.65 and lift >= 0.1
    outcome = 'UNKNOWN_INSUFFICIENT_PATIENT_SUPPORT' if not support else 'PASS' if numeric and qa and qa['pass_gate'] else 'FAIL'
    result = dict(claim_type=p['claim_type'], external_referent=p['external_referent'], outcome=outcome, metrics=results, lift_over_no_individual_geometry_prior=lift, support_gate=support, numeric_gate=numeric, parser_gate=qa['pass_gate'] if qa else None, excluded_records=len(ex), exclusion_reasons=dict(collections.Counter((r['reason'] for r in ex))), evaluable_test_patients=len({r['case_id'] for r in rows}), total_test_patients=len(m['test']), patient_dropout_fraction=(len(m['test']) - len({r['case_id'] for r in rows})) / len(m['test']), spectrum_bias='Positive named teeth in treated patients vs all represented teeth in explicitly untreated patients; cannot estimate general population accuracy.', cost=cost(st, cpu))
    write('raw/R1_RESULTS.json', result)
    write('raw/R1_ERROR_CASES.json', [r for r in rows if r['y'] != r['predicted']['shape']])
    (ROOT / 'HANDOFF_R1.md').write_text(f"R1 {outcome}: shape BA={mm['balanced_accuracy']:.6f}; sensitivity={mm['sensitivity']:.6f}; specificity={mm['specificity']:.6f}; lift={lift:.6f}. Explicit photo-report truth only. No intrinsic material inference. Change operation next: source-linked material applicability port, not more same-shape tuning. Failed criterion remains frozen.\n")
    state('R1_DECIDED', outcome, 'R2 add observed restoration state to exact tooth consumer before physics')
    print(json.dumps(result), flush=True)
if __name__ == '__main__':
    globals()[sys.argv[1]]()
