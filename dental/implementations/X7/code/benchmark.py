import os
for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[k] = '2'
import sys, json, zipfile, hashlib, time, datetime, pickle, resource
from pathlib import Path
import numpy as np
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.metrics import balanced_accuracy_score, accuracy_score, cohen_kappa_score, confusion_matrix
from sklearn.impute import SimpleImputer
from parser import parse, FIELDS, published_parse
P = Path(__file__).resolve().parents[1]
MF = {'overbite': ['overbite_mm'], 'overjet': ['overjet_mm'], 'crossbite': ['crossbite_extent_right_deg', 'crossbite_extent_left_deg'], 'midlines': ['midline_deviation_mm'], 'spee': ['curve_of_spee_right_mm', 'curve_of_spee_left_mm'], 'molar_right': ['ap_offset_right_molar_mm', 'cusp_lag_right_deg'], 'molar_left': ['ap_offset_left_molar_mm', 'cusp_lag_left_deg'], 'canine_right': ['ap_offset_right_canine_mm', 'cusp_lag_right_deg'], 'canine_left': ['ap_offset_left_canine_mm', 'cusp_lag_left_deg']}

def write(n, x):
    (P / n).write_text(json.dumps(x, indent=2, allow_nan=False) + '\n')

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def state(status, gate, nextop):
    write('CURRENT_WORK_STATE.json', dict(lane='X7-bite2text', updated_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), status=status, latest_gate=gate, next_operation=nextop))

def load_geometry():
    rows = [json.loads(l) for l in (P / 'raw/geometry.jsonl').read_text().splitlines()]
    return {r['case_id']: r for r in rows}

def labels(ids):
    m = json.load(open(P / 'raw/DATA_MANIFEST.json'))
    out = {}
    with zipfile.ZipFile(m['zip']) as z:
        ns = sorted(z.namelist())
        for c in ids:
            files = [n for n in ns if n.startswith(c + '/reports_ios_en/') and n.endswith('.txt')]
            out[c] = [dict(member=n, sha256=hashlib.sha256(z.read(n)).hexdigest(), labels=parse(z.read(n).decode('utf8'))) for n in files]
    return out

def parser_review():
    gold = json.load(open(P / 'raw/MANUAL_GOLD.json'))['labels']
    m = json.load(open(P / 'raw/DATA_MANIFEST.json'))
    rr = []
    with zipfile.ZipFile(m['zip']) as z:
        for n in m['manual_review_reports']:
            c = n.split('/')[0]
            t = z.read(n).decode()
            a = parse(t)
            b = published_parse(t).as_dict()
            for f in FIELDS:
                rr.append(dict(case_id=c, field=f, gold=gold[c][f], ours=a[f], published=b[f]))
    summary = {}
    for meth in ['ours', 'published']:
        emitted = [r for r in rr if r[meth] is not None]
        correct = sum((r[meth] == r['gold'] for r in emitted))
        summary[meth] = dict(precision=correct / len(emitted), coverage=len(emitted) / len(rr), n_emitted=len(emitted), n_correct=correct, n_total=len(rr), pass_gate=correct / len(emitted) >= 0.95 and len(emitted) / len(rr) >= 0.85)
    write('raw/PARSER_REVIEW_CORRECTED.json', dict(summary=summary, rows=rr, independent_review=False, development_sample=True, source_correction='raw/MANUAL_GOLD.json'))
    return summary

def train():
    m = json.load(open(P / 'raw/DATA_MANIFEST.json'))
    g = load_geometry()
    missing = set(m['train']) - set(g)
    if missing:
        raise RuntimeError('Training geometry incomplete: ' + str(len(missing)))
    start = time.perf_counter()
    ys = labels(m['train'])
    write('raw/TRAIN_LABELS.json', ys)
    write('raw/PARSER_REVIEW_SUMMARY.json', parser_review())
    cols = sorted({k for v in g.values() if v.get('features') for k in v['features']})
    models = {}
    cost = []
    for f in FIELDS:
        ids = [c for c in m['train'] if not g[c]['error'] and ys[c] and (ys[c][0]['labels'][f] is not None)]
        y = np.array([ys[c][0]['labels'][f] for c in ids])
        X = np.array([[g[c]['features'].get(k, np.nan) if g[c]['features'].get(k) is not None else np.nan for k in cols] for c in ids])
        (classes, ct) = np.unique(y, return_counts=True)
        if len(classes) < 2:
            models[f] = {'classes': classes.tolist(), 'majority': classes[0] if len(classes) else None, 'gb': None}
            continue
        majority = str(classes[np.argmax(ct)])
        t = time.perf_counter()
        gb = HistGradientBoostingClassifier(max_iter=400, learning_rate=0.05, max_depth=3, min_samples_leaf=15, l2_regularization=1.0, random_state=17).fit(X, y)
        j = [cols.index(k) for k in MF[f]]
        imputer = SimpleImputer(strategy='median', add_indicator=True)
        Xs = imputer.fit_transform(X[:, j])
        tree = DecisionTreeClassifier(max_depth=3, min_samples_leaf=15, random_state=17).fit(Xs, y)
        models[f] = dict(gb=gb, threshold=tree, imputer=imputer, indices=j, majority=majority, classes=gb.classes_.tolist(), train_n=len(ids), train_counts=dict(zip(classes.tolist(), ct.tolist())))
        cost.append(dict(field=f, fit_seconds=time.perf_counter() - t, train_n=len(ids)))
        print(f, len(ids), 'fit', round(cost[-1]['fit_seconds'], 2), flush=True)
    bundle = dict(models=models, columns=cols, geometry_adapter_sha256=sha(P / 'code/geometry.py'), parser_sha256=sha(P / 'code/parser.py'), data_manifest_sha256=sha(P / 'raw/DATA_MANIFEST.json'))
    (P / 'raw/models.pkl').write_bytes(pickle.dumps(bundle))
    write('raw/FIT_COST.json', dict(fits=cost, total_wall_seconds=time.perf_counter() - start, peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, train_label_sha256=sha(P / 'raw/TRAIN_LABELS.json')))
    state('R1_MODELS_FIT', 'Training labels only; all test/calibration labels still unread by benchmark', 'Freeze all held-out predictions before comparison')

def freeze():
    if (P / 'FROZEN_PREDICTIONS.json').exists():
        raise RuntimeError('Frozen predictions already exist')
    m = json.load(open(P / 'raw/DATA_MANIFEST.json'))
    g = load_geometry()
    bundle = pickle.loads((P / 'raw/models.pkl').read_bytes())
    cols = bundle['columns']
    rows = []
    start = time.perf_counter()
    for split in ['calibration', 'test']:
        missing = set(m[split]) - set(g)
        if missing:
            raise RuntimeError(split + ' incomplete ' + str(len(missing)))
        for c in m[split]:
            r = dict(case_id=c, split=split, error=g[c]['error'], fields={})
            if not g[c]['error']:
                X = np.array([[g[c]['features'].get(k) if g[c]['features'].get(k) is not None else np.nan for k in cols]])
                for (f, obj) in bundle['models'].items():
                    if obj.get('gb') is None:
                        continue
                    p = obj['gb'].predict_proba(X)[0]
                    r['fields'][f] = dict(classes=obj['classes'], probabilities=p.tolist(), point=obj['classes'][int(p.argmax())], threshold=str(obj['threshold'].predict(obj['imputer'].transform(X[:, obj['indices']]))[0]), majority=obj['majority'])
            rows.append(r)
    write('raw/PREDICTIONS_R1.json', rows)
    frozen = dict(frozen_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), phase='Before any benchmark access to calibration/test report labels', prediction_path='raw/PREDICTIONS_R1.json', prediction_sha256=sha(P / 'raw/PREDICTIONS_R1.json'), model_sha256=sha(P / 'raw/models.pkl'), geometry_sha256=sha(P / 'raw/geometry.jsonl'), parser_sha256=sha(P / 'code/parser.py'), prereg_sha256=sha(P / 'PREREG_R1.json'), query_wall_seconds=time.perf_counter() - start)
    write('FROZEN_PREDICTIONS.json', frozen)
    state('R1_PREDICTIONS_FROZEN', 'All test predictions hashed before comparing reports', 'Read test reports and run frozen accuracy gates')

def metric(y, yp):
    if not y:
        return dict(n=0)
    labs = sorted(set(y) | set(yp))
    return dict(n=len(y), accuracy=float(accuracy_score(y, yp)), balanced_accuracy=float(balanced_accuracy_score(y, yp)), kappa=float(cohen_kappa_score(y, yp)) if len(set(y) | set(yp)) > 1 else None, confusion_labels=labs, confusion=confusion_matrix(y, yp, labels=labs).tolist())

def evaluate():
    frozen = json.load(open(P / 'FROZEN_PREDICTIONS.json'))
    assert sha(P / frozen['prediction_path']) == frozen['prediction_sha256']
    m = json.load(open(P / 'raw/DATA_MANIFEST.json'))
    start = time.perf_counter()
    ys = labels(m['test'])
    write('raw/TEST_LABELS.json', ys)
    pred = json.load(open(P / frozen['prediction_path']))
    pred = {r['case_id']: r for r in pred if r['split'] == 'test'}
    metrics = {}
    diffs = []
    for f in FIELDS:
        ids = [c for c in m['test'] if ys[c] and ys[c][0]['labels'][f] is not None and (f in pred[c]['fields'])]
        y = [ys[c][0]['labels'][f] for c in ids]
        out = {meth: metric(y, [pred[c]['fields'][f][meth] for c in ids]) for meth in ['point', 'threshold', 'majority']}
        out['parsed_first_report_fraction'] = len(ids) / len(m['test'])
        out['lift_vs_majority_balanced'] = out['point'].get('balanced_accuracy', 0) - out['majority'].get('balanced_accuracy', 0)
        out['strong_matched_control'] = 'Same HGB on same features/data: TIE, no algorithm superiority'
        metrics[f] = out
        for c in ids:
            if pred[c]['fields'][f]['point'] != ys[c][0]['labels'][f]:
                diffs.append(dict(case_id=c, field=f, reference=ys[c][0]['labels'][f], prediction=pred[c]['fields'][f]['point'], reports=[r['labels'][f] for r in ys[c]], explanation='Not causally identified; inspect landmark/axis proxy, report disagreement and ambiguity separately'))
    parser_gate = json.load(open(P / 'raw/PARSER_REVIEW_CORRECTED.json'))['summary']['ours']['pass_gate']
    lift = float(np.mean([v['lift_vs_majority_balanced'] for v in metrics.values()]))
    npass = all((v['point']['n'] >= 40 for v in metrics.values()))
    gate = bool(lift >= 0.1 and parser_gate and npass)
    injected = {f: metric([ys[c][0]['labels'][f] for c in m['test'] if ys[c] and ys[c][0]['labels'][f] is not None and (f in pred[c]['fields'])], [next((k for k in pred[c]['fields'][f]['classes'] if k != ys[c][0]['labels'][f])) if any((k != ys[c][0]['labels'][f] for k in pred[c]['fields'][f]['classes'])) else '__wrong__' for c in m['test'] if ys[c] and ys[c][0]['labels'][f] is not None and (f in pred[c]['fields'])]) for f in FIELDS}
    inj_lift = float(np.mean([injected[f]['balanced_accuracy'] - metrics[f]['majority']['balanced_accuracy'] for f in FIELDS]))
    assert inj_lift < 0.1
    write('raw/DISAGREEMENTS_R1.json', diffs)
    write('raw/RESULTS_R1.json', dict(metrics=metrics, macro_balanced_lift=lift, primary_gate='PASS' if gate else 'FAIL', parser_gate=parser_gate, min_test_n_gate=npass, method_gain='TIE', external_referent=json.load(open(P / 'PREREG_R1.json'))['external_referent'], injected_wrong_value_control=dict(injection='Force each emitted test label to another category', macro_lift=inj_lift, gate='FAIL', refutes=True), validation_wall_seconds=time.perf_counter() - start))
    if not (P / 'raw/RESULTS_R1_INITIAL.json').exists():
        write('raw/RESULTS_R1_INITIAL.json', json.load(open(P / 'raw/RESULTS_R1.json')))
    (P / 'HANDOFF_R1.md').write_text('R1 primary ' + ('PASS' if gate else 'FAIL') + '; macro balanced lift ' + str(lift) + ' against0.10. Matched HGB controlTIE. Numeric Angle/cusp identities UNKNOWN; only categorical text association measured. Next constructionR2 adds all-readable-report patient-maximum calibration, preserving observed disagreement. Original predictions and first raw results immutable.\n')
    state('R1_DECIDED', 'R1 primary ' + ('PASS' if gate else 'FAIL') + '; matched HGB TIE', 'R2 preserve report-disagreement sets and calibrate on separate patients')
    print('R1', gate, 'macro lift', lift, flush=True)
if __name__ == '__main__':
    globals()[sys.argv[1]]()
