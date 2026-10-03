"N39 scoring according to PREREG.md (sha256 in ../PREREG.sha256).\nCopied with source attribution from results/N12b/code/n12b_score.py: PERSON, metrics(), plate_mixed(), load() (unchanged),\nN1 adaptation and the COMAK frame rule. c′ = N12b's LOSO selection (N12b/scores.json -> LOSO_picks.cprime).\nUsage: python3 n39_score.py. Writes ../scores.json."
import glob
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
N12 = ROOT.parent / 'N12'
N12B = ROOT.parent / 'N12b'
OUT = N12B / 'cloud_out' / 'out'
XOUT = ROOT / 'cloud_out' / 'out'
PERSON = {'JW1': 'JW', 'JW4': 'JW', 'DM2': 'DM', 'DM6': 'DM', 'SC3': 'SC', 'PS5': 'PS'}
PERSONS = ['JW', 'DM', 'SC', 'PS']
MIX_MIN = 60
ALTERED = {'crouch', 'mildcrouch', 'bouncy', 'medthrust', 'mtgait', 'mtpgait'}
NORMAL = {'ngait', 'smooth', 'tsgait', 'trunksway', 'wpgait', 'rightturn'}
THETAS = [20, 25, 30, 35, 40]
SEED = 20260924


# ---- copied from N12b n12b_score.py
def metrics(p, m, bw):
    e = p - m
    r = float(np.corrcoef(p, m)[0, 1]) if np.std(p) > 0 and np.std(m) > 0 else float('nan')
    return dict(rmse_N=float(np.sqrt(np.mean(e ** 2))), rmse_BW=float(np.sqrt(np.mean(e ** 2)) / bw), r=r,
                peak_err_pct=float((p.max() - m.max()) / m.max() * 100), peak_pred_pctBW=float(p.max() / bw * 100),
                peak_meas_pctBW=float(m.max() / bw * 100))


def plate_mixed(pu):
    for s in pu:
        if s == 'none':
            continue
        r, l = s[1:].split('/L')
        if int(r) > MIX_MIN and int(l) > MIX_MIN:
            return True
    return False


def load():
    T, excl = {}, {}
    for f in sorted(glob.glob(str(OUT / '*_summary.json'))):
        S = json.load(open(f))
        sess = S['sess']
        for t, r in S['trials'].items():
            key = f'{sess}__{t}'
            if r.get('status') != 'ok':
                excl[key] = r.get('status') + (': ' + r.get('error', '') if r.get('error') else '')
                continue
            if r['nan_frac'] > 0.10:
                excl[key] = f"nan_frac {r['nan_frac']:.2f}"
                continue
            if plate_mixed(r['plate_use']):
                excl[key] = f"platta delad {r['plate_use']}"
                continue
            z = dict(np.load(OUT / f'{key}.npz'))
            inr = z['inr'].astype(bool)
            n = len(inr)
            if inr.mean() < 0.9:
                excl[key] = f'measuring signal covers {inr.mean():.2f} by window'
                continue
            z['inr_mask'] = inr.copy()
            z = {k: (v[..., inr] if (hasattr(v, 'shape') and v.ndim >= 1 and v.shape[-1] == n) else
                     (v[inr] if (hasattr(v, 'shape') and v.ndim >= 1 and v.shape[0] == n) else v)) for k, v in z.items()}
            T[key] = dict(sess=sess, person=PERSON[sess], trial=t, act=r['activity'], z=z, info=r)
    return T, excl
# ---- slut kopia


def lsq(cols, y):
    return np.linalg.lstsq(np.stack(cols, 1), y, rcond=None)[0]


def main():
    T, excl = load()
    keys = sorted(T)
    picks = json.load(open(N12B / 'scores.json'))['LOSO_picks']['cprime']
    res = dict(n_trials=len(T), n_excluded=len(excl), cprime_picks=picks)
    # ---- extraherade ID-storheter, samma rutor
    chk = []
    for k in keys:
        z = T[k]['z']
        X = dict(np.load(XOUT / f'{k}.npz'))
        tfull = X['t']
        m = np.isin(np.round(tfull, 9), np.round(z['t'], 9))
        assert m.sum() == len(z['t']), k
        for f in ('Msag', 'Mfront', 'N0_tot', 'ext_M', 'ext_C', 'flx_M', 'flx_C'):
            z['X_' + f] = X[f][m]
        chk.append(float(np.abs(z['X_N0_tot'] - z['N0_tot']).max()))
        ae = np.maximum(z['X_Msag'] / z['X_ext_M'], 0)
        af = np.maximum(z['X_Msag'] / z['X_flx_M'], 0)
        z['MB'] = ae * z['X_ext_C'] + af * z['X_flx_C']
        z['both_pos'] = int(((ae > 0) & (af > 0)).sum())
        z['Mf_d'] = np.abs(z['X_Mfront']) / float(X['d'])
        z['G'] = np.linalg.norm(z['grf'], axis=1)
        z['cp'] = z[f"{picks[T[k]['person']]}_tot"]
        z['flex_med'] = float(np.median(z['knee_flex']))
    res['extract_check'] = dict(N0_max_abs_diff_N=max(chk), frames_both_groups_active=int(sum(T[k]['z']['both_pos'] for k in keys)))
    grp = {k: ('altered' if T[k]['act'] in ALTERED else 'normal') for k in keys}
    assert all(T[k]['act'] in ALTERED | NORMAL for k in keys), {T[k]['act'] for k in keys}

    def cat(ks, f):
        return np.concatenate([f(T[k]['z']) for k in ks])

    Y = lambda z: z['meas_tot']  # noqa: E731
    params, pred = {}, {k: {} for k in keys}
    # placebo-MB: next trial (cyclic) by the same person, time-normalized
    for P in PERSONS:
        ks = [k for k in keys if T[k]['person'] == P]
        for i, k in enumerate(ks):
            src = T[ks[(i + 1) % len(ks)]]['z']['MB']
            n = len(T[k]['z']['t'])
            T[k]['z']['MB_pl'] = np.interp(np.linspace(0, 1, n), np.linspace(0, 1, len(src)), src)
    rng = np.random.default_rng(SEED)
    grp_pl = {}
    for P in PERSONS:
        ks = [k for k in keys if T[k]['person'] == P]
        lab = [grp[k] for k in ks]
        perm = rng.permutation(len(ks))
        for i, k in enumerate(ks):
            grp_pl[k] = lab[perm[i]]

    def fitA1(tr, k1, G_):
        w = {}
        for g in ('normal', 'altered'):
            ks = [k for k in tr if G_[k] == g]
            if not ks:
                w[g] = 0.0
                continue
            n1 = cat(ks, lambda z: k1 * z['G'])
            dv = cat(ks, lambda z: z['cp']) - n1
            w[g] = float(np.clip((dv * (cat(ks, Y) - n1)).sum() / (dv * dv).sum(), 0, 1))
        return w

    for P in PERSONS:
        tr = [k for k in keys if T[k]['person'] != P]
        te = [k for k in keys if T[k]['person'] == P]
        y = cat(tr, Y)
        G = cat(tr, lambda z: z['G'])
        k1 = float((G * y).sum() / (G * G).sum())
        MB = cat(tr, lambda z: z['MB'])
        pp = dict(k_N1=k1)
        pp['A1_w'] = fitA1(tr, k1, grp)
        pp['A1pl_w'] = fitA1(tr, k1, grp_pl)
        pp['A2'] = {}
        for g in ('normal', 'altered'):
            ks = [k for k in tr if grp[k] == g]
            pp['A2'][g] = [float(x) for x in lsq([cat(ks, lambda z: z['G']), cat(ks, lambda z: z['MB'])], cat(ks, Y))]
        pp['L1'] = [float(x) for x in lsq([G, MB], y)]
        pp['L1pl'] = [float(x) for x in lsq([G, cat(tr, lambda z: z['MB_pl'])], y)]
        pp['L2'] = [float(x) for x in lsq([G, MB, cat(tr, lambda z: z['Mf_d'])], y)]
        # L3: θ* is selected on the training persons
        best = None
        for th in THETAS:
            hi = [k for k in tr if T[k]['z']['flex_med'] > th]
            if hi:
                n1 = cat(hi, lambda z: k1 * z['G'])
                dv = cat(hi, lambda z: z['cp']) - n1
                w = float(np.clip((dv * (cat(hi, Y) - n1)).sum() / (dv * dv).sum(), 0, 1))
            else:
                w = 0.0
            per_p = []
            for Q_ in PERSONS:
                if Q_ == P:
                    continue
                r_ = []
                for k in tr:
                    if T[k]['person'] != Q_:
                        continue
                    z = T[k]['z']
                    ww = w if z['flex_med'] > th else 0.0
                    pr = (1 - ww) * k1 * z['G'] + ww * z['cp']
                    r_.append(np.sqrt(np.mean((pr - z['meas_tot']) ** 2)) / float(z['bwN']))
                per_p.append(np.median(r_))
            sc = float(np.median(per_p))
            if best is None or sc < best[0]:
                best = (sc, th, w)
        pp['L3'] = dict(theta=best[1], w=best[2], train_score=best[0])
        params[P] = pp
        for k in te:
            z = T[k]['z']
            n1 = k1 * z['G']
            w = pp['A1_w'][grp[k]]
            wpl = pp['A1pl_w'][grp_pl[k]]
            a2 = pp['A2'][grp[k]]
            ww = pp['L3']['w'] if z['flex_med'] > pp['L3']['theta'] else 0.0
            pred[k] = {'N1': n1, 'N0': z['N0_tot'], 'cprime': z['cp'], 'Hg_N12b': None,
                       'A1': (1 - w) * n1 + w * z['cp'], 'A1_placebo': (1 - wpl) * n1 + wpl * z['cp'],
                       'A2': a2[0] * z['G'] + a2[1] * z['MB'],
                       'L1': pp['L1'][0] * z['G'] + pp['L1'][1] * z['MB'],
                       'L1_placebo': pp['L1pl'][0] * z['G'] + pp['L1pl'][1] * z['MB_pl'],
                       'L2': pp['L2'][0] * z['G'] + pp['L2'][1] * z['MB'] + pp['L2'][2] * z['Mf_d'],
                       'L3': (1 - ww) * n1 + ww * z['cp'],
                       'P0': z['N0_tot'] + z['MB']}
            del pred[k]['Hg_N12b']
    res['params_LOSO'] = params
    # kontroll mot N12b:s N1-k
    n12b = json.load(open(N12B / 'scores.json'))
    res['N1_k_vs_N12b_max_abs'] = max(abs(params[P]['k_N1'] - n12b['N1_params'][P]['k']) for P in PERSONS)
    MODELS = list(pred[keys[0]].keys())
    per = {k: {mn: metrics(pred[k][mn], T[k]['z']['meas_tot'], float(T[k]['z']['bwN'])) for mn in MODELS} for k in keys}

    def pmed(mn, P, act=None, field='rmse_BW'):
        v = [per[k][mn][field] for k in keys if T[k]['person'] == P and (act is None or T[k]['act'] in act)]
        return float(np.median(v)) if v else None

    summ = {}
    for mn in MODELS:
        d = dict(per_person={P: pmed(mn, P) for P in PERSONS},
                 per_person_ngait={P: pmed(mn, P, {'ngait'}) for P in PERSONS},
                 per_person_peak_err={P: pmed(mn, P, None, 'peak_err_pct') for P in PERSONS},
                 per_person_r={P: float(np.nanmedian([per[k][mn]['r'] for k in keys if T[k]['person'] == P])) for P in PERSONS})
        d['median_over_persons'] = float(np.median(list(d['per_person'].values())))
        d['median_over_persons_ngait'] = float(np.median(list(d['per_person_ngait'].values())))
        d['per_activity_median'] = {a: float(np.median([per[k][mn]['rmse_BW'] for k in keys if T[k]['act'] == a]))
                                    for a in sorted({T[k]['act'] for k in keys})}
        summ[mn] = d
    res['summary'] = summ
    res['activity_n'] = {a: sum(T[k]['act'] == a for k in keys) for a in sorted({T[k]['act'] for k in keys})}
    # ---- COMAK, same frames (N12's rule, copied from N12b)
    cs = dict(np.load(N12 / 'comak_series.npz'))
    cm = {}
    for sess, trials in (('DM6', ['ngait_og1', 'ngait_og3', 'ngait_og4', 'ngait_og5', 'ngait_og6']),
                         ('JW1', ['ngait_2', 'ngait_3', 'ngait_4', 'ngait_5', 'ngait_6'])):
        rows = {}
        for tr in trials:
            key = f'{sess}__{tr}'
            if key not in T or f'{key}__V1' not in cs:
                rows[tr] = 'saknas: ' + str(excl.get(key, 'ingen COMAK'))
                continue
            z = T[key]['z']
            tc, fc_ = cs[f'{key}__V1']
            m = (z['t'] >= tc[0]) & (z['t'] <= tc[-1])
            if m.sum() < 10:
                rows[tr] = f'too few common boxes ({int(m.sum())})'
                continue
            bw = float(z['bwN'])
            meas = z['meas_tot'][m]
            r = dict(n=int(m.sum()), COMAK=metrics(np.interp(z['t'][m], tc, fc_), meas, bw))
            for mn in MODELS:
                r[mn] = metrics(pred[key][mn][m], meas, bw)
            rows[tr] = r
        ok = [v for v in rows.values() if isinstance(v, dict)]
        cm[sess] = dict(trials=rows, median_rmse_BW={mn: float(np.median([v[mn]['rmse_BW'] for v in ok])) for mn in ['COMAK'] + MODELS},
                        median_peak_err={mn: float(np.median([v[mn]['peak_err_pct'] for v in ok])) for mn in ['COMAK'] + MODELS})
    res['COMAK_common'] = cm
    # ---- B24-like peak null (walking, person level), copied rule from N12b
    pk = {}
    for P in PERSONS:
        ks = [k for k in keys if T[k]['person'] == P and T[k]['act'] == 'ngait']
        pk[P] = dict(meas=float(np.median([per[k]['N1']['peak_meas_pctBW'] for k in ks])), n=len(ks),
                     **{mn: float(np.median([per[k][mn]['peak_pred_pctBW'] for k in ks])) for mn in MODELS})
    for P in pk:
        pk[P]['B24null'] = float(np.median([pk[Q_]['meas'] for Q_ in pk if Q_ != P]))
    res['B24_peak_ngait'] = dict(per_person=pk, rmse_pctBW={mn: float(np.sqrt(np.mean([(pk[P][mn] - pk[P]['meas']) ** 2 for P in pk])))
                                                            for mn in MODELS + ['B24null']})
    # ---- kriterier
    crit = {}
    N1s = summ['N1']
    for mn in MODELS:
        if mn in ('N1',):
            continue
        s = summ[mn]
        k1 = s['median_over_persons'] < N1s['median_over_persons']
        k2 = {sess: cm[sess]['median_rmse_BW'][mn] < cm[sess]['median_rmse_BW']['COMAK'] for sess in cm}
        dng = {P: s['per_person_ngait'][P] - N1s['per_person_ngait'][P] for P in PERSONS}
        k3 = s['median_over_persons_ngait'] <= N1s['median_over_persons_ngait'] and max(dng.values()) <= 0.02
        k4 = res['B24_peak_ngait']['rmse_pctBW'][mn] < res['B24_peak_ngait']['rmse_pctBW']['B24null']
        crit[mn] = dict(K1=bool(k1), K1_delta_BW=s['median_over_persons'] - N1s['median_over_persons'],
                        persons_beating_N1=[P for P in PERSONS if s['per_person'][P] < N1s['per_person'][P]],
                        K2=k2, K2_all=bool(all(k2.values())), K3=bool(k3), K3_ngait_delta_BW=dng,
                        K4=bool(k4), beats_N1_and_COMAK=bool(k1 and all(k2.values()) and k3),
                        approved_claim=bool(k1 and all(k2.values()) and k3 and k4))
    res['criteria'] = crit
    res['per_trial'] = {k: dict(person=T[k]['person'], act=T[k]['act'], group=grp[k], flex_med=T[k]['z']['flex_med'],
                                rmse_BW={mn: per[k][mn]['rmse_BW'] for mn in MODELS},
                                peak_err_pct={mn: per[k][mn]['peak_err_pct'] for mn in MODELS}) for k in keys}
    json.dump(res, open(ROOT / 'scores.json', 'w'), indent=1)
    print('check', res['extract_check'], 'N1 k diff', res['N1_k_vs_N12b_max_abs'])
    for mn in MODELS:
        s = summ[mn]
        print(f"{mn:12s} med {s['median_over_persons']:.3f} ngait {s['median_over_persons_ngait']:.3f}",
              {P: round(v, 3) for P, v in s['per_person'].items()},
              'COMAK', {x: round(cm[x]['median_rmse_BW'][mn], 3) for x in cm}, 'B24', round(res['B24_peak_ngait']['rmse_pctBW'][mn], 1))
    print('COMAK', {x: round(cm[x]['median_rmse_BW']['COMAK'], 3) for x in cm}, 'B24null', round(res['B24_peak_ngait']['rmse_pctBW']['B24null'], 1))
    print(json.dumps(params, indent=0)[:2500])
    for mn, c in crit.items():
        print(mn, {k: v for k, v in c.items() if k != 'K3_ngait_delta_BW'})


if __name__ == '__main__':
    main()
