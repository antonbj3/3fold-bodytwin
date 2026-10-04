"""Fit only explicit development donors; serialize numeric arrays, no pickle."""
from common import *
from generator import features, predict
from scipy.spatial import cKDTree

def dataset(keys):
    ds = {}
    rejected = []
    access = []
    for key in keys:
        refs = dev_reference(key)
        access.append(key)
        for t in tasks(key)[::4]:
            if t['status'] != 'READY':
                rejected.append({'key': key, 'family': t['family'], 'reason': t['site_error']})
                continue
            t = scene(t)
            ref = np.array(refs[t['family']])
            ok = np.isfinite(ref)
            if ok.mean() < 0.5:
                rejected.append({'key': key, 'family': t['family'], 'reason': 'less than 50% native height coverage'})
                continue
            if not ok.all():
                (_, j) = cKDTree(t['uv'][ok]).query(t['uv'][~ok])
                ref[~ok] = ref[ok][j]
            gap = t['ceiling'] - ref
            contact = (np.isfinite(gap) & (gap >= 0) & (gap <= 0.1)).astype(float)
            ds.setdefault(t['family'], []).append(dict(key=key, x=features(t), y=ref - t['prior'], contact=contact, ref_valid_mask=ok))
    return (ds, rejected, access)

def fit_arrays(records, kind, param, output_modes=12):
    X = np.array([r['x'] for r in records])
    Y = np.array([r['y'] for r in records])
    xm = X.mean(0)
    xs = np.maximum(X.std(0), 0.1)
    X = (X - xm) / xs
    ym = Y.mean(0)
    Yr = Y - ym
    model = {'kind': np.array(kind), 'xmean': xm, 'xscale': xs, 'ymean': ym, 'training_keys': np.array([r['key'] for r in records])}
    if kind == 'kernel':
        (gamma, ridge) = param
        (_, _, V) = np.linalg.svd(Yr, full_matrices=False)
        basis = V[:min(output_modes, len(V))] if output_modes else np.eye(Y.shape[1])
        out = Yr @ basis.T
        D = np.mean((X[:, None, :] - X[None, :, :]) ** 2, axis=2)
        K = np.exp(-gamma * D)
        model.update(xtrain=X, gamma=np.array(gamma), basis=basis, coef=np.linalg.solve(K + ridge * np.eye(len(X)), out))
    elif kind == 'ridge':
        model['coef'] = X.T @ np.linalg.solve(X @ X.T + float(param) * np.eye(len(X)), Yr)
    elif kind == 'forest':
        from sklearn.ensemble import RandomForestRegressor
        rf = RandomForestRegressor(n_estimators=96, min_samples_leaf=int(param), max_features=1.0, n_jobs=2, random_state=6142)
        rf.fit(X, Y)
        model['ntrees'] = np.array(len(rf.estimators_))
        for (i, tr) in enumerate(rf.estimators_):
            a = tr.tree_
            prefix = 'tree' + str(i) + '_'
            model.update({prefix + 'left': a.children_left, prefix + 'right': a.children_right, prefix + 'feature': a.feature, prefix + 'threshold': a.threshold, prefix + 'values': a.value[:, :, 0]})
        for j in [0, len(records) - 1]:
            x = X[j]
            v = []
            for tr in rf.estimators_:
                a = tr.tree_
                node = 0
                xx = x.astype(np.float32)
                while a.children_left[node] != -1:
                    node = int(a.children_left[node] if xx[a.feature[node]] <= a.threshold[node] else a.children_right[node])
                v.append(a.value[node, :, 0])
            if np.max(np.abs(np.mean(v, 0) - rf.predict(x[None, :])[0])) > 1e-12:
                raise ValueError('forest export differs')
    return model

def loadmodel(path):
    with np.load(path, allow_pickle=False) as a:
        return dict(a)
