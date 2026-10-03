"""geomgr.evaluate: LOSO / external evaluation of conditional instantiation (numpy/scipy only; cloud shards).

  python3 -m geomgr.evaluate --data frozen --task femur_imp|femur_vsd|femur_pool|tibia --shard i --nshards n --out DIR

Writes DIR/<task>_<i>.json (rows) and DIR/<task>_<i>.npz (coverage curves on KAPPA_GRID, landmark z2).
Scenarios and noise follow PREREG.md / PREREG_ADD1.md.
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np

from . import core as C
from . import certify as CE

SEED = 20260924
N_DRAWS = 5
LM5 = ['HC6', 'SGT', 'LT', 'MEC', 'LEC']
FEMUR_SCEN = {'S0': ([], []), 'S1': (['L_mech'], []), 'S2': (C.X5, []), 'S3': (C.T5, []),
              'S4': ([], LM5), 'S5': (['L_mech'], ['MEC', 'LEC'])}
# S6 (added after PREREG, secondary, not graded): X5 + LM5 = the demo request
FEMUR_SCEN_EXTRA = {'S6': (C.X5, LM5)}
SCEN_FILTER = None


def scen_items():
    allsc = dict(FEMUR_SCEN, **FEMUR_SCEN_EXTRA)
    if SCEN_FILTER:
        return [(k, allsc[k]) for k in SCEN_FILTER]
    return list(FEMUR_SCEN.items())
TIBIA_SCEN = {'T0': [], 'TH': ['height'], 'THW': ['height', 'weight'], 'TL': ['L_tib'], 'THL': ['height', 'L_tib']}
TIBIA_NOISE_SD = {'height': 1.0, 'weight': 1.0, 'age': 1.0, 'L_tib': 2.0, 'W_plat': 2.0}
TIBIA_FEATS = ['height', 'weight', 'age', 'L_tib', 'W_plat']


# ------------------------------------------------------------------ femur setup
class FemurData:
    def __init__(self, d):
        P = np.load(Path(d) / 'femur_pop.npz')
        O = np.load(Path(d) / 'vsd_obs.npz')
        self.S, self.F = P['S'], P['F']
        self.names = [str(x) for x in P['names']]
        self.groups = np.array([str(x) for x in P['groups']])
        self.dataset = np.array([str(x) for x in P['dataset']])
        self.lm_tri, self.lm_bc = P['lm_tri'], P['lm_bc']
        self.loo_tri, self.loo_bc = P['loo_tri'], P['loo_bc']
        self.subj = [str(x) for x in O['subj']]
        self.obs_lm, self.obs_feat, self.feat_mean20 = O['obs_lm'], O['obs_feat'], O['feat_mean20']
        self.dev_anat, self.ndef = O['dev_anat'], O['ndef']
        self.valid = np.isfinite(self.obs_lm[:, :, 0, 0])

    def anchors(self, tri, bc):
        bcn = np.clip(bc, 0, None)
        bcn = bcn / bcn.sum(1, keepdims=True)
        return C.Anchors(C.LM_NAMES, [self.F[t] for t in tri], list(bcn), self.S.shape[1])

    def features(self, anchors, idx):
        return np.array([C.femur_features_from_landmarks(C.femur_points(anchors, self.S[i]), self.S[i]) for i in idx])

    def feat_pool(self, exclude=None):
        rows = []
        for s in range(len(self.subj)):
            if s == exclude:
                continue
            for r in np.flatnonzero(self.valid[s]):
                rows.append(self.obs_feat[s, r] - self.feat_mean20[s])
        return np.array(rows)

    def lm_cov(self, exclude=None, ndef_for=None):
        """Rater covariance (anatomical frame) per point + definition second moment for VSD tests."""
        out = {}
        for j, n in enumerate(C.PT_NAMES):
            D = [self.dev_anat[s, r, j] for s in range(len(self.subj)) if s != exclude
                 for r in np.flatnonzero(self.valid[s])]
            D = np.array(D)
            Cr = np.cov(D.T)
            if ndef_for is not None:
                Nd = np.array([self.ndef[ndef_for, q, j] for q in range(len(self.subj)) if q != ndef_for])
                Cr = Cr + Nd.T @ Nd / len(Nd)
            out[n] = Cr
        return out

    def noise_draw(self, rng, exclude=None):
        pairs = [(s, r) for s in range(len(self.subj)) if s != exclude for r in np.flatnonzero(self.valid[s])]
        s, r = pairs[rng.integers(len(pairs))]
        return self.obs_feat[s, r] - self.feat_mean20[s], self.dev_anat[s, r]


def frame_fn(P):
    return C.knee_frame(P)


def run_case(model, anchors, F, truth, feat_obs, R_feat, lm_obs, C_lm, placed_truth=True):
    """Condition, predict, evaluate against truth (N,3) in the observation frame. Returns row, curves."""
    t0 = time.process_time()
    post = C.condition(model, anchors, C.femur_points, feat_obs, R_feat,
                       lm_obs, C_lm, frame_fn=frame_fn)
    X, Cv = C.predictive(model, post)
    t_cond = time.process_time() - t0
    s, R, t = C.umeyama(X, truth, scale=False)
    Xa = C.apply(1.0, R, t, X)
    Ca = np.einsum('ij,njk,lk->nil', R, Cv, R)
    e = truth - Xa
    nrm = C.vertex_normals(Xa, F)
    cur = C.coverage_curves(e, Ca, nrm)
    row = dict(E_shape=float(np.sqrt((e ** 2).sum(1).mean())), z2_median=cur['z2_median'],
               sig_n_median=cur['sig_n_median'], t_condition_cpu_s=t_cond)
    lp = C.landmark_predictive(model, anchors, C.femur_points, post, C.PT_NAMES)
    Pt = C.femur_points(anchors, truth)
    z2 = []
    for n in C.PT_NAMES:
        p, Cl = lp[n]
        pa = R @ p + t
        Cla = R @ Cl @ R.T
        d = Pt[n] - pa
        z2.append(float(d @ np.linalg.solve(Cla, d)))
        if n == 'HC6':
            row['E_HJC'] = float(np.linalg.norm(d))
    curves = dict(cov3=cur['cov3'], cov1=cur['cov1'], lmz2=np.array(z2))
    cert = CE.numeric_certificate(model, post, X, F, feat_obs)
    row['cert_accepted'] = cert['accepted']
    row['cert_reasons'] = cert['reasons']
    row['extrapolation_inside'] = cert['checks']['extrapolation']['inside']
    row['stage_a_p'] = post.stage_a['p']
    if post.stage_b is not None:
        row['stage_b_p'] = post.stage_b['p']
        row['pose_identified'] = bool(post.pose_identified)
        row['gn_iters'] = post.stage_b['iters']
        if post.pose_identified and placed_truth:
            Xp, Cp = C.predictive(model, post, placed=True)
            ep = truth - Xp
            curp = C.coverage_curves(ep, Cp, C.vertex_normals(Xp, F))
            row['E_placed'] = float(np.sqrt((ep ** 2).sum(1).mean()))
            lpp = C.landmark_predictive(model, anchors, C.femur_points, post, ['HC6'], placed=True)
            row['E_HJC_placed'] = float(np.linalg.norm(Pt['HC6'] - lpp['HC6'][0]))
            dd = Pt['HC6'] - lpp['HC6'][0]
            row['HJC_placed_z2'] = float(dd @ np.linalg.solve(lpp['HC6'][1], dd))
            curves['cov3_placed'] = curp['cov3']
            curves['cov1_placed'] = curp['cov1']
    return row, curves


def femur_model(D, tr_idx, anchors):
    M = D.features(anchors, tr_idx)
    m = C.ShapeModel(D.S[tr_idx], M, C.FEMUR_FEATURES, D.groups[tr_idx])
    m.set_oos(C.oos_variance(D.S[tr_idx], D.groups[tr_idx]))
    return m


def synth_obs(D, anchors, truth, fnames, lnames, rng, exclude=None):
    fn, dev = D.noise_draw(rng, exclude)
    ftrue = C.femur_features_from_landmarks(C.femur_points(anchors, truth), truth)
    fo = {k: float(ftrue[C.FEMUR_FEATURES.index(k)] + fn[C.FEMUR_FEATURES.index(k)]) for k in fnames}
    Pt = C.femur_points(anchors, truth)
    _, Ax = C.knee_frame(Pt)
    lo = {k: Pt[k] + dev[C.PT_NAMES.index(k)] @ Ax for k in lnames}
    return fo, lo


def real_obs(D, s, r, fnames, lnames):
    fo = {k: float(D.obs_feat[s, r, C.FEMUR_FEATURES.index(k)]) for k in fnames}
    lo = {k: D.obs_lm[s, r, C.PT_NAMES.index(k)] for k in lnames}
    return fo, lo


def _append(rows, curves_acc, meta, row, cur):
    row.update(meta)
    rows.append(row)
    for k, v in cur.items():
        curves_acc.setdefault(k, []).append(v)
        curves_acc.setdefault(k + '_row', []).append(len(rows) - 1)


def task_femur_imp(D, persons, rows, acc):
    anchors = D.anchors(D.lm_tri, D.lm_bc)
    imp = np.flatnonzero(D.dataset == 'imperial')
    pool = D.feat_pool()
    Clm = D.lm_cov()
    all_persons = sorted(set(D.groups[imp]))
    for p in persons:
        rng = np.random.default_rng([SEED, int(p[1:])])
        tr = np.array([i for i in imp if D.groups[i] != p])
        te = [i for i in imp if D.groups[i] == p]
        t0 = time.time()
        model = femur_model(D, tr, anchors)
        t_model = time.time() - t0
        for i in te:
            truth = D.S[i]
            for dr in range(N_DRAWS):
                for sc, (fn, ln) in scen_items():
                    fo, lo = synth_obs(D, anchors, truth, fn, ln, rng)
                    R = np.cov(pool[:, [C.FEMUR_FEATURES.index(k) for k in fn]].T).reshape(len(fn), len(fn)) if fn else None
                    row, cur = run_case(model, anchors, D.F, truth, fo, R, lo, Clm)
                    if sc == 'S1':
                        mm = C.femur_features_from_landmarks(C.femur_points(anchors, model.mu), model.mu)
                        c = model.mu.mean(0)
                        B0 = c + fo['L_mech'] / mm[C.FEMUR_FEATURES.index('L_mech')] * (model.mu - c)
                        s_, Rb, tb = C.umeyama(B0, truth, scale=False)
                        row['E_shape_B0'] = float(np.sqrt(((truth - C.apply(1, Rb, tb, B0)) ** 2).sum(1).mean()))
                    _append(rows, acc, dict(task='femur_imp', person=p, member=D.names[i], draw=dr, scenario=sc,
                                            t_model_s=t_model), row, cur)
            if SCEN_FILTER:
                continue
            # placebo: S2 with the features of the next person's femur (same side index)
            q = all_persons[(all_persons.index(p) + 1) % len(all_persons)]
            j = [k for k in imp if D.groups[k] == q][te.index(i)]
            fo, _ = synth_obs(D, anchors, D.S[j], C.X5, [], rng)
            R = np.cov(pool[:, [C.FEMUR_FEATURES.index(k) for k in C.X5]].T)
            row, cur = run_case(model, anchors, D.F, truth, fo, R, {}, Clm)
            _append(rows, acc, dict(task='femur_imp', person=p, member=D.names[i], draw=0, scenario='S2_placebo',
                                    placebo_features_from=D.names[j]), row, cur)


def task_femur_vsd(D, subjects, rows, acc):
    imp = np.flatnonzero(D.dataset == 'imperial')
    sig2 = C.oos_variance(D.S[imp], D.groups[imp])
    for sid in subjects:
        s = D.subj.index(sid)
        rng = np.random.default_rng([SEED, 1000 + s])
        anchors = D.anchors(D.loo_tri[s], D.loo_bc[s])
        M = D.features(anchors, imp)
        model = C.ShapeModel(D.S[imp], M, C.FEMUR_FEATURES, D.groups[imp])
        model.set_oos(sig2)
        pool = D.feat_pool(exclude=s)
        Clm = D.lm_cov(exclude=s, ndef_for=s)
        truth = D.S[D.names.index(f'vsd:{sid}')]
        ratings = rng.choice(np.flatnonzero(D.valid[s]), size=N_DRAWS, replace=False)
        for dr, r in enumerate(ratings):
            for sc, (fn, ln) in scen_items():
                fo, lo = real_obs(D, s, r, fn, ln)
                R = np.cov(pool[:, [C.FEMUR_FEATURES.index(k) for k in fn]].T).reshape(len(fn), len(fn)) if fn else None
                row, cur = run_case(model, anchors, D.F, truth, fo, R, lo, Clm)
                _append(rows, acc, dict(task='femur_vsd', person=sid, member=f'vsd:{sid}', draw=dr, rating=int(r),
                                        scenario=sc), row, cur)


def task_femur_pool(D, persons, rows, acc):
    allidx = np.arange(len(D.S))
    for p in persons:
        rng = np.random.default_rng([SEED, 5000 + sorted(set(D.groups)).index(p)])
        is_vsd = p in D.subj
        s = D.subj.index(p) if is_vsd else None
        anchors = D.anchors(D.loo_tri[s], D.loo_bc[s]) if is_vsd else D.anchors(D.lm_tri, D.lm_bc)
        tr = allidx[D.groups != p]
        te = allidx[D.groups == p]
        model = femur_model(D, tr, anchors)
        pool = D.feat_pool(exclude=s)
        Clm = D.lm_cov(exclude=s, ndef_for=s) if is_vsd else D.lm_cov()
        for i in te:
            truth = D.S[i]
            if is_vsd:
                draws = [(dr, r) for dr, r in enumerate(rng.choice(np.flatnonzero(D.valid[s]), N_DRAWS, replace=False))]
            else:
                draws = [(dr, None) for dr in range(N_DRAWS)]
            for dr, r in draws:
                for sc, (fn, ln) in scen_items():
                    if is_vsd:
                        fo, lo = real_obs(D, s, r, fn, ln)
                    else:
                        fo, lo = synth_obs(D, anchors, truth, fn, ln, rng)
                    R = np.cov(pool[:, [C.FEMUR_FEATURES.index(k) for k in fn]].T).reshape(len(fn), len(fn)) if fn else None
                    row, cur = run_case(model, anchors, D.F, truth, fo, R, lo, Clm)
                    _append(rows, acc, dict(task='femur_pool', person=p, member=D.names[i], dataset=str(D.dataset[i]),
                                            draw=dr, scenario=sc, real_obs=is_vsd), row, cur)


# ------------------------------------------------------------------ tibia
class TibiaData:
    def __init__(self, d):
        import json as _j
        P = np.load(Path(d) / 'tibia_pop.npz')
        self.S, self.F = P['S'], P['F']
        self.case = [str(x) for x in P['case']]
        reg = _j.loads((Path(d) / 'registry_tibia_r.json').read_text())
        lms = [e for e in reg['entities'] if e['kind'] == 'landmark']
        self.anchors = C.Anchors([e['name'] for e in lms], [e['indices'] for e in lms], [e['weights'] for e in lms],
                                 self.S.shape[1])
        age = P['age_yr'].copy()
        age[np.isnan(age)] = np.nanmean(age)          # declared: 1 missing age imputed by mean
        self.cov = dict(height=P['height_cm'], weight=P['weight_kg'], age=age)

    def feats(self, i, V=None):
        V = self.S[i] if V is None else V
        g = C.tibia_geom_features(self.anchors.dict(V))
        return np.array([self.cov['height'][i], self.cov['weight'][i], self.cov['age'][i], g['L_tib'], g['W_plat']])


def task_tibia(T, persons, rows, acc):
    n = len(T.S)
    for p in persons:
        i = T.case.index(p)
        rng = np.random.default_rng([SEED, 9000 + i])
        tr = np.array([k for k in range(n) if k != i])
        M = np.array([T.feats(k) for k in tr])
        model = C.ShapeModel(T.S[tr], M, TIBIA_FEATS, np.array(T.case)[tr])
        model.set_oos(C.oos_variance(T.S[tr], np.array(T.case)[tr]))
        truth = T.S[i]
        ftrue = T.feats(i)
        for dr in range(N_DRAWS):
            for sc, fn in TIBIA_SCEN.items():
                fo = {k: float(ftrue[TIBIA_FEATS.index(k)] + rng.normal(0, TIBIA_NOISE_SD[k])) for k in fn}
                R = np.diag([TIBIA_NOISE_SD[k] ** 2 for k in fn]) if fn else None
                post = C.condition(model, T.anchors, C.tibia_points, fo, R, None, None)
                X, Cv = C.predictive(model, post)
                s_, Rr, tt = C.umeyama(X, truth, scale=False)
                Xa = C.apply(1, Rr, tt, X)
                Ca = np.einsum('ij,njk,lk->nil', Rr, Cv, Rr)
                e = truth - Xa
                cur = C.coverage_curves(e, Ca, C.vertex_normals(Xa, T.F))
                lp = C.landmark_predictive(model, T.anchors, C.tibia_points, post, T.anchors.names)
                Pt = T.anchors.dict(truth)
                z2 = []
                for nme in T.anchors.names:
                    pp, Cl = lp[nme]
                    d = Pt[nme] - (Rr @ pp + tt)
                    z2.append(float(d @ np.linalg.solve(Rr @ Cl @ Rr.T, d)))
                cert = CE.numeric_certificate(model, post, X, T.F, fo)
                row = dict(E_shape=float(np.sqrt((e ** 2).sum(1).mean())), z2_median=cur['z2_median'],
                           sig_n_median=cur['sig_n_median'], cert_accepted=cert['accepted'], cert_reasons=cert['reasons'],
                           stage_a_p=post.stage_a['p'], extrapolation_inside=cert['checks']['extrapolation']['inside'])
                _append(rows, acc, dict(task='tibia', person=p, member=p, draw=dr, scenario=sc), row,
                        dict(cov3=cur['cov3'], cov1=cur['cov1'], lmz2=np.array(z2)))


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--data', required=True)
    ap.add_argument('--task', required=True)
    ap.add_argument('--shard', type=int, default=0)
    ap.add_argument('--nshards', type=int, default=1)
    ap.add_argument('--out', required=True)
    ap.add_argument('--limit', type=int, default=0)
    ap.add_argument('--scen', default='')
    a = ap.parse_args(argv)
    global SCEN_FILTER
    SCEN_FILTER = [x for x in a.scen.split(',') if x] or None
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    rows, acc = [], {}
    if a.task.startswith('femur'):
        D = FemurData(a.data)
        if a.task == 'femur_imp':
            units = sorted(set(D.groups[D.dataset == 'imperial']))
        elif a.task == 'femur_vsd':
            units = list(D.subj)
        else:
            units = sorted(set(D.groups))
    else:
        T = TibiaData(a.data)
        units = list(T.case)
    mine = units[a.shard::a.nshards]
    if a.limit:
        mine = mine[:a.limit]
    fn = {'femur_imp': task_femur_imp, 'femur_vsd': task_femur_vsd, 'femur_pool': task_femur_pool}.get(a.task)
    if fn is not None:
        fn(D, mine, rows, acc)
    else:
        task_tibia(T, mine, rows, acc)
    tag = f'{a.task}_{a.shard:02d}' + (f"_{'-'.join(SCEN_FILTER)}" if SCEN_FILTER else '')
    (out / f'{tag}.json').write_text(json.dumps(dict(task=a.task, shard=a.shard, nshards=a.nshards, units=mine,
                                                     wall_s=time.time() - t0, rows=rows), default=float))
    np.savez_compressed(out / f'{tag}.npz', **{k: np.asarray(v) for k, v in acc.items()})
    print(tag, len(rows), round(time.time() - t0, 1), flush=True)


if __name__ == '__main__':
    main()
