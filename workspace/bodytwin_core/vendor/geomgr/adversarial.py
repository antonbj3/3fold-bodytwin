"""geomgr.adversarial: the 10 pre-registered adversarial cases (PREREG K4) + genuine controls.

  python3 -m geomgr adversarial [--subject z001] [--out results/N7a/adversarial.json]
Base = a real single rating of the subject (external, LOO manager). Each case must be REJECTED by the certificate.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from . import core as C
from .api import GeometryManager

L5 = ['HC6', 'SGT', 'LT', 'MEC', 'LEC']


def run(subject='z001', out=None, data=None, gm=None):
    root = Path(__file__).resolve().parents[1]
    data = Path(data or '/media/anton/sdc1-tmp/bodytwin/N7a/frozen')
    gm = gm or GeometryManager(data, population=('imperial',), exclude_vsd=subject,
                               cache_dir='/media/anton/sdc1-tmp/bodytwin/N7a')
    O = np.load(data / 'vsd_obs.npz')
    s = [str(x) for x in O['subj']].index(subject)
    r0 = int(np.flatnonzero(np.isfinite(O['obs_lm'][s, :, 0, 0]))[0])
    f = {k: float(O['obs_feat'][s, r0, C.FEMUR_FEATURES.index(k)]) for k in C.FEMUR_FEATURES}
    lm = {k: O['obs_lm'][s, r0, C.PT_NAMES.index(k)].copy() for k in C.PT_NAMES}
    U = lambda ks: {k: C.feature_units(k) for k in ks}  # noqa: E731
    X5 = {k: f[k] for k in C.X5}
    L = {k: lm[k] for k in L5}
    _, Ax = C.knee_frame(lm)
    cases = {}

    def req(name, **kw):
        r = gm.instantiate(subject=f'{subject}:{name}', with_points=False, **kw)
        cases[name] = dict(accepted=bool(r.accepted), reasons=r.certificate['reasons'],
                           flags=r.certificate.get('flags', []))
        return r

    # genuine controls (must be accepted)
    req('control_X5', features=X5, feature_units=U(X5))
    req('control_L5', landmarks=L)
    req('control_X5+L5', features=X5, feature_units=U(X5), landmarks=L)
    # 1 unit error: L_mech in metres labelled mm
    x = dict(X5, L_mech=X5['L_mech'] / 1000)
    req('A1_unit_error_m_as_mm', features=x, feature_units=U(x))
    # 2 absurd value
    x = dict(X5, CCD=175.0)
    req('A2_absurd_CCD_175deg', features=x, feature_units=U(x))
    # 3 contradictory measures: mechanical length -40 mm, trochanter-epicondyle length +40 mm
    x = {'L_mech': f['L_mech'] - 40, 'L_tl': f['L_tl'] + 40}
    req('A3_contradictory_lengths', features=x, feature_units=U(x))
    # 4 swapped MEC/LEC
    y = dict(L, MEC=L['LEC'], LEC=L['MEC'])
    req('A4_swapped_MEC_LEC', landmarks=y)
    # 5 left given as right (mirror x)
    y = {k: v * np.array([-1.0, 1, 1]) for k, v in L.items()}
    req('A5_left_as_right_mirrored', landmarks=y)
    # 6 one landmark outlier 30 mm (anterior)
    y = dict(L, SGT=L['SGT'] + 30 * Ax[1])
    req('A6_outlier_SGT_30mm', landmarks=y)
    # 7 extreme direct coefficient: first fold on modes 1..5, scanned in SD units; control at +-3 SD accepted
    scan = []
    hit = None
    for k in range(5):
        for c in (3, -3, 5, -5, 8, -8, 12, -12, 16, -16, 20, -20, 30, -30):
            b = np.zeros(gm.model.r)
            b[k] = c * np.sqrt(gm.model.lam[k])
            r = gm.instantiate_coefficients(b, subject=f'mode{k + 1}_{c}sd')
            nf = r.certificate['checks']['fold']['n_folded_faces']
            scan.append(dict(mode=k + 1, sd=c, folded=nf, accepted=bool(r.accepted)))
            if nf > 0 and hit is None:
                hit = (k, c, r)
    if hit is not None:
        k, c, r = hit
        cases['A7_extreme_coefficient_folds'] = dict(accepted=bool(r.accepted), reasons=r.certificate['reasons'],
                                                    mode=k + 1, sd=c)
    else:
        cases['A7_extreme_coefficient_folds'] = dict(accepted=None, reasons=['no fold found up to 30 SD'])
    ctrl3 = [x for x in scan if abs(x['sd']) == 3]
    # 8 NaN
    x = dict(X5, AV=float('nan'))
    req('A8_nan_value', features=x, feature_units=U(x))
    # 9 unknown measure
    req('A9_unknown_measure', features={'femur_length_cm': 43.0}, feature_units={'femur_length_cm': 'cm'})
    # 10 wrong units string
    req('A10_wrong_unit_string', features={'L_mech': X5['L_mech'] / 10}, feature_units={'L_mech': 'cm'})
    adv = {k: v for k, v in cases.items() if k.startswith('A')}
    ctl = {k: v for k, v in cases.items() if k.startswith('control')}
    res = dict(subject=subject, rating_index=r0, n_adversarial=len(adv),
               n_rejected=int(sum(v['accepted'] is False for v in adv.values())),
               controls_accepted=int(sum(v['accepted'] for v in ctl.values())), n_controls=len(ctl),
               fold_scan_control_3sd_accepted=int(sum(x['accepted'] for x in ctrl3)), n_control_3sd=len(ctrl3),
               cases=cases, fold_scan=scan)
    out = Path(out or root / 'adversarial.json')
    out.write_text(json.dumps(res, indent=1, default=float))
    return res
