"N39 post hoc (EJ preregistered, logged as outlier): is A2's gain the muscle term or just group-k?\nN1g = k_g·|GRF| per activity group (LOSO), A2_placebo = A2 with MB from another trial (same rule as L1 placebo).\nReusing n39_score.py (same directory, own code). Typing ../posthoc.json."
import json
import numpy as np
import n39_score as S

T, excl = S.load()
keys = sorted(T)
picks = json.load(open(S.N12B / 'scores.json'))['LOSO_picks']['cprime']
for k in keys:
    z = T[k]['z']
    X = dict(np.load(S.XOUT / f'{k}.npz'))
    m = np.isin(np.round(X['t'], 9), np.round(z['t'], 9))
    ae = np.maximum(X['Msag'][m] / X['ext_M'][m], 0); af = np.maximum(X['Msag'][m] / X['flx_M'][m], 0)
    z['MB'] = ae * X['ext_C'][m] + af * X['flx_C'][m]
    z['G'] = np.linalg.norm(z['grf'], axis=1)
for P in S.PERSONS:
    ks = [k for k in keys if T[k]['person'] == P]
    for i, k in enumerate(ks):
        src = T[ks[(i + 1) % len(ks)]]['z']['MB']; n = len(T[k]['z']['t'])
        T[k]['z']['MB_pl'] = np.interp(np.linspace(0, 1, n), np.linspace(0, 1, len(src)), src)
grp = {k: ('altered' if T[k]['act'] in S.ALTERED else 'normal') for k in keys}
cat = lambda ks, f: np.concatenate([f(T[k]['z']) for k in ks])
pred, par = {}, {}
for P in S.PERSONS:
    tr = [k for k in keys if T[k]['person'] != P]
    par[P] = {}
    for g in ('normal', 'altered'):
        ks = [k for k in tr if grp[k] == g]
        y = cat(ks, lambda z: z['meas_tot']); G = cat(ks, lambda z: z['G'])
        kg = float((G * y).sum() / (G * G).sum())
        pl = [float(x) for x in S.lsq([G, cat(ks, lambda z: z['MB_pl'])], y)]
        par[P][g] = dict(k_g=kg, A2pl=pl)
    for k in keys:
        if T[k]['person'] != P:
            continue
        z, q = T[k]['z'], par[P][grp[k]]
        pred[k] = dict(N1g=q['k_g'] * z['G'], A2_placebo=q['A2pl'][0] * z['G'] + q['A2pl'][1] * z['MB_pl'])
base = json.load(open(S.ROOT / 'scores.json'))
out = dict(params=par, summary={}, COMAK_common={}, B24_peak_rmse_pctBW={})
for mn in ('N1g', 'A2_placebo'):
    per = {k: S.metrics(pred[k][mn], T[k]['z']['meas_tot'], float(T[k]['z']['bwN'])) for k in keys}
    pp = {P: float(np.median([per[k]['rmse_BW'] for k in keys if T[k]['person'] == P])) for P in S.PERSONS}
    png = {P: float(np.median([per[k]['rmse_BW'] for k in keys if T[k]['person'] == P and T[k]['act'] == 'ngait'])) for P in S.PERSONS}
    pa = {a: float(np.median([per[k]['rmse_BW'] for k in keys if T[k]['act'] == a])) for a in sorted({T[k]['act'] for k in keys})}
    out['summary'][mn] = dict(per_person=pp, median_over_persons=float(np.median(list(pp.values()))), per_person_ngait=png,
                              median_over_persons_ngait=float(np.median(list(png.values()))), per_activity_median=pa)
    cs = dict(np.load(S.N12 / 'comak_series.npz')); cm = {}
    for sess, trials in (('DM6', ['ngait_og1', 'ngait_og3', 'ngait_og4', 'ngait_og5', 'ngait_og6']), ('JW1', ['ngait_2', 'ngait_3', 'ngait_4', 'ngait_5', 'ngait_6'])):
        v = []
        for tr in trials:
            key = f'{sess}__{tr}'
            if key not in T or f'{key}__V1' not in cs: continue
            z = T[key]['z']; tc, _ = cs[f'{key}__V1']; m = (z['t'] >= tc[0]) & (z['t'] <= tc[-1])
            if m.sum() < 10: continue
            v.append(S.metrics(pred[key][mn][m], z['meas_tot'][m], float(z['bwN']))['rmse_BW'])
        cm[sess] = float(np.median(v))
    out['COMAK_common'][mn] = cm
    pk = {}
    for P in S.PERSONS:
        ks = [k for k in keys if T[k]['person'] == P and T[k]['act'] == 'ngait']
        pk[P] = (np.median([per[k]['peak_pred_pctBW'] for k in ks]), np.median([per[k]['peak_meas_pctBW'] for k in ks]))
    out['B24_peak_rmse_pctBW'][mn] = float(np.sqrt(np.mean([(a - b) ** 2 for a, b in pk.values()])))
json.dump(out, open(S.ROOT / 'posthoc.json', 'w'), indent=1)
for mn, s in out['summary'].items():
    print(mn, round(s['median_over_persons'], 3), {P: round(v, 3) for P, v in s['per_person'].items()}, 'ngait', round(s['median_over_persons_ngait'], 3),
          {P: round(v, 3) for P, v in s['per_person_ngait'].items()}, out['COMAK_common'][mn], round(out['B24_peak_rmse_pctBW'][mn], 1))
    print('   act', {a: round(v, 3) for a, v in s['per_activity_median'].items()})
print(json.dumps(par))
