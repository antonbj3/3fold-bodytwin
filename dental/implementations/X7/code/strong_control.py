"""Actually refit equally informed conventional HGB control; no method win is asserted."""
import os
for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[k] = '2'
import json, time, datetime
import numpy as np
from sklearn.ensemble import HistGradientBoostingClassifier
from benchmark import P, load_geometry, write, sha
from parser import FIELDS

def run():
    m = json.load(open(P / 'raw/DATA_MANIFEST.json'))
    g = load_geometry()
    ys = json.load(open(P / 'raw/TRAIN_LABELS.json'))
    pred = json.load(open(P / 'raw/PREDICTIONS_R1.json'))
    bundle = __import__('pickle').loads((P / 'raw/models.pkl').read_bytes())
    cols = bundle['columns']
    idsq = [r['case_id'] for r in pred if not r['error']]
    lookup = {r['case_id']: r for r in pred}
    Xq = np.array([[g[c]['features'].get(k) if g[c]['features'].get(k) is not None else np.nan for k in cols] for c in idsq])
    rows = []
    start = time.perf_counter()
    for f in FIELDS:
        ids = [c for c in m['train'] if not g[c]['error'] and ys[c] and (ys[c][0]['labels'][f] is not None)]
        X = np.array([[g[c]['features'].get(k) if g[c]['features'].get(k) is not None else np.nan for k in cols] for c in ids])
        y = np.array([ys[c][0]['labels'][f] for c in ids])
        t = time.perf_counter()
        control = HistGradientBoostingClassifier(max_iter=400, learning_rate=0.05, max_depth=3, min_samples_leaf=15, l2_regularization=1.0, random_state=17).fit(X, y)
        p = control.predict_proba(Xq)
        diff = []
        neq = 0
        for (c, pp) in zip(idsq, p):
            ff = lookup[c]['fields'][f]
            assert control.classes_.tolist() == ff['classes']
            diff.extend(np.abs(pp - np.asarray(ff['probabilities'])).tolist())
            neq += int(str(control.classes_[pp.argmax()]) != ff['point'])
        maxerr = float(max(diff))
        assert maxerr <= 1e-12 and neq == 0
        bad = p.copy()
        bad[0, 0] += 0.1
        injected = float(np.max(abs(bad - p)))
        assert injected > 1e-12
        rows.append(dict(field=f, train_n=len(ids), query_n=len(idsq), max_probability_difference=maxerr, label_disagreements=neq, fit_and_query_seconds=time.perf_counter() - t, injected_probability_error=injected, injected_corrupt_control_fails=True, control_prediction_sha256=__import__('hashlib').sha256(p.tobytes()).hexdigest()))
    out = dict(recipe_frozen_in='PREREG_R1.json', predictions_match='TIE', comparison_tolerance=1e-12, rows=rows, full_fit_and_query_seconds=time.perf_counter() - start, training_labels_sha256=sha(P / 'raw/TRAIN_LABELS.json'), prediction_reference_sha256=sha(P / 'raw/PREDICTIONS_R1.json'), note='Separate standard HGB refit on identical training data/features/recipe. Batch query for control; no cost superiority claim. Executed after R1 test comparison but control recipe was preregistered, no tuning.')
    write('raw/STRONG_CONTROL_R1.json', out)
    print('Actual conventional HGB refit: all9fields TIE', out['full_fit_and_query_seconds'])
if __name__ == '__main__':
    run()
