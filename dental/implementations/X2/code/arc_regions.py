"""Monotone arch-coordinate classifier; fit on independent local FDI labels."""
import json, time
import numpy as np
from occlusion_operator import H, features, load_obj, teeth_files, predict_regions, TYPE_TO_REGION, sha, dump, state

def angle(v, ras=False):
    (X, frame) = features(v, ras)
    return (np.arctan2(X[:, 0], 0.25 + (X[:, 1] + 1) / 2), frame)

def predict(v, model, ras=False):
    (t, f) = angle(v, ras)
    r = np.searchsorted(model['thresholds_rad'], t)
    return (r + 4 * (v[:, 0] < f['x_center_mm']), f)

def fit():
    models = {}
    records = []
    start = time.perf_counter()
    edges = np.linspace(0, np.pi, 513)
    for jaw in ['upper', 'lower']:
        hist = np.zeros((4, 512))
        inputs = []
        for (pid, p) in teeth_files(jaw, 'training')[:50]:
            v = load_obj(p)
            lab = np.asarray(json.loads(p.with_suffix('.json').read_text())['labels'])
            (theta, fr) = angle(v)
            for r in range(4):
                hist[r] += np.histogram(theta[(lab > 0) & (TYPE_TO_REGION[lab % 10] == r)], edges)[0]
            inputs.append(dict(patient=pid, mesh=str(p), mesh_sha256=sha(p), labels_sha256=sha(p.with_suffix('.json'))))
        hist /= np.maximum(hist.sum(1)[:, None], 1)
        cost = hist.sum(0)[None, :] - hist
        prefix = np.c_[np.zeros(4), cost.cumsum(1)]
        dp = np.full((5, 513), np.inf)
        prev = np.zeros((5, 513), int)
        dp[0, 0] = 0
        for n in range(1, 5):
            for end in range(n, 513):
                starts = np.arange(n - 1, end)
                values = dp[n - 1, starts] + prefix[n - 1, end] - prefix[n - 1, starts]
                best = values.argmin()
                dp[n, end] = values[best]
                prev[n, end] = starts[best]
        cuts = []
        end = 512
        for n in range(4, 0, -1):
            end = int(prev[n, end])
            cuts.append(end)
        cuts = sorted(cuts)[1:]
        model = dict(thresholds_rad=edges[cuts].tolist(), training_arches=50, training_balanced_error=float(dp[4, 512] / 4), inputs=inputs)
        conf = np.zeros((4, 4), int)
        control = np.zeros((4, 4), int)
        r1 = json.loads((H / 'raw/segmentation_model.json').read_text())[jaw]
        for (pid, p) in teeth_files(jaw, 'testing')[30:60]:
            v = load_obj(p)
            (prediction, fr) = predict(v, model)
            (r0, _, _) = predict_regions(v, r1)
            ph = __import__('hashlib').sha256(prediction.astype('int8').tobytes()).hexdigest()
            lab = np.asarray(json.loads(p.with_suffix('.json').read_text())['labels'])
            ok = lab > 0
            true = TYPE_TO_REGION[lab[ok] % 10]
            c = np.bincount(true * 4 + prediction[ok] % 4, minlength=16).reshape(4, 4)
            c0 = np.bincount(true * 4 + r0[ok] % 4, minlength=16).reshape(4, 4)
            conf += c
            control += c0
            records.append(dict(patient=pid, jaw=jaw, prediction_sha256=ph, mesh=str(p), mesh_sha256=sha(p), labels_sha256=sha(p.with_suffix('.json')), confusion=c.tolist(), control_confusion=c0.tolist()))

        def metrics(c):
            recall = np.diag(c) / np.maximum(c.sum(1), 1)
            ap = (c[:2, 2:].sum() + c[2:, :2].sum()) / c.sum()
            return dict(confusion=c.tolist(), per_region_recall=recall.tolist(), macro_recall=float(recall.mean()), anterior_posterior_error=float(ap), gate=bool(recall.mean() >= 0.85 and ap <= 0.05))
        model['validation'] = metrics(conf)
        model['R1_matched_control'] = metrics(control)
        models[jaw] = model
        print('ARC', jaw, model['validation'], flush=True)
    dump(H / 'raw/segmentation_r2_model.json', models)
    dump(H / 'raw/segmentation_r2_validation.json', dict(models=models, records=records, wall_s=time.perf_counter() - start))
    state('R2_REGION_FIT_FINISHED', str({j: m['validation']['gate'] for (j, m) in models.items()}), 'Reuse fixed-pose maps, compare fresh-region predictions and fine geometry')
    return models
if __name__ == '__main__':
    fit()
