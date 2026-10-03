"""geomgr.summarize: pool cloud shards -> eval_summary.json + calibration.json (kappa per scenario).

  python3 -m geomgr.summarize --in DIR [DIR ...] --out results/N7a
kappa: Imperial LOSO cross-fitted per person (kappa_p from the other persons' curves); VSD uses kappa from all
Imperial LOSO curves of the same scenario (no parameter chosen on VSD). Same scheme within femur_pool / tibia.
"""
from __future__ import annotations

import argparse
import glob
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy import stats

from . import core as C


def load(dirs):
    cases = []
    for d in dirs:
        for jf in sorted(glob.glob(str(Path(d) / '**' / '*.json'), recursive=True)):
            try:
                J = json.loads(Path(jf).read_text())
            except Exception:
                continue
            if 'rows' not in J or 'task' not in J:
                continue
            Z = np.load(jf[:-5] + '.npz')
            rows = J['rows']
            cur = [dict() for _ in rows]
            for k in Z.files:
                if k.endswith('_row'):
                    continue
                for v, ri in zip(Z[k], Z[k + '_row']):
                    cur[int(ri)][k] = v
            for r, c in zip(rows, cur):
                r['_curves'] = c
                r['_file'] = Path(jf).name
                cases.append(r)
    # de-duplicate (a shard can appear in several collected dirs)
    seen, out = set(), []
    for r in cases:
        key = (r['task'], r['person'], r['member'], r['draw'], r['scenario'])
        if key not in seen:
            seen.add(key)
            out.append(r)
    return out


def cov_at(cases, key, kappa):
    v = [C.at_kappa(r['_curves'][key], kappa) for r in cases if key in r['_curves']]
    return float(np.mean(v)) if v else None


def lm_cov(cases, kappa, idx=None):
    z = np.concatenate([np.asarray(r['_curves']['lmz2'])[idx if idx is not None else slice(None)]
                        for r in cases if 'lmz2' in r['_curves']])
    return float(np.mean(z <= C.CHI2_3_90 * kappa ** 2))


def kappa_of(cases, key):
    cs = [r['_curves'][key] for r in cases if key in r['_curves']]
    return C.kappa_for(cs) if cs else None


def crossfit(cases, key):
    """Per person: kappa from the other persons; returns mean coverage and kappa range."""
    by = defaultdict(list)
    for r in cases:
        if key in r['_curves']:
            by[r['person']].append(r)
    covs, ks, lmc = [], [], []
    for p, rs in by.items():
        others = [r for q, x in by.items() if q != p for r in x]
        k = kappa_of(others, key)
        ks.append(k)
        covs.extend(C.at_kappa(r['_curves'][key], k) for r in rs)
        if key == 'cov3':
            for r in rs:
                z = np.asarray(r['_curves']['lmz2'])
                lmc.append(np.mean(z <= C.CHI2_3_90 * k ** 2))
    if not ks:
        return None
    return dict(coverage=float(np.mean(covs)), kappa_min=float(min(ks)), kappa_max=float(max(ks)),
                landmark_coverage=float(np.mean(lmc)) if lmc else None)


def person_mean(cases, field):
    by = defaultdict(list)
    for r in cases:
        if r.get(field) is not None:
            by[r['person']].append(r[field])
    return {p: float(np.mean(v)) for p, v in by.items()}


def summarize_task(cases, task, ref_cases=None):
    out = {}
    scen = sorted({r['scenario'] for r in cases})
    for sc in scen:
        cs = [r for r in cases if r['scenario'] == sc]
        pm = person_mean(cs, 'E_shape')
        d = dict(n_cases=len(cs), n_persons=len(pm), E_shape_mean_mm=float(np.mean(list(pm.values()))),
                 E_shape_p95_person_mm=float(np.percentile(list(pm.values()), 95)))
        for f in ('E_placed', 'E_HJC', 'E_HJC_placed', 'E_shape_B0', 'sig_n_median'):
            pmf = person_mean(cs, f)
            if pmf:
                d[f + '_mean'] = float(np.mean(list(pmf.values())))
        d['coverage90_kappa1'] = cov_at(cs, 'cov3', 1.0)
        d['coverage90_normal_kappa1'] = cov_at(cs, 'cov1', 1.0)
        d['landmark_coverage90_kappa1'] = lm_cov(cs, 1.0)
        d['cert_false_reject_rate'] = float(np.mean([not r['cert_accepted'] for r in cs]))
        d['extrapolation_flag_rate'] = float(np.mean([not r.get('extrapolation_inside', True) for r in cs]))
        rej = defaultdict(int)
        for r in cs:
            for x in r['cert_reasons']:
                rej[x.split('(')[0].strip()] += 1
        d['cert_reject_reasons'] = dict(rej)
        if 'HJC_placed_z2' in cs[0] or any('HJC_placed_z2' in r for r in cs):
            z = np.array([r['HJC_placed_z2'] for r in cs if 'HJC_placed_z2' in r])
            d['HJC_placed_coverage90_kappa1'] = float(np.mean(z <= C.CHI2_3_90))
        if ref_cases is None:      # calibrate within the task (cross-fit)
            for key in ('cov3', 'cov1', 'cov3_placed'):
                cf = crossfit(cs, key)
                if cf:
                    d[f'crossfit_{key}'] = cf
            d['kappa_all'] = {k: kappa_of(cs, k) for k in ('cov3', 'cov1', 'cov3_placed') if kappa_of(cs, k)}
        else:
            ref = [r for r in ref_cases if r['scenario'] == sc]
            ks = {k: kappa_of(ref, k) for k in ('cov3', 'cov1', 'cov3_placed') if kappa_of(ref, k)}
            d['kappa_from_reference'] = ks
            if 'cov3' in ks:
                d['coverage90_kappa_ref'] = cov_at(cs, 'cov3', ks['cov3'])
                d['landmark_coverage90_kappa_ref'] = lm_cov(cs, ks['cov3'])
            if 'cov1' in ks:
                d['coverage90_normal_kappa_ref'] = cov_at(cs, 'cov1', ks['cov1'])
            if 'cov3_placed' in ks:
                d['coverage90_placed_kappa_ref'] = cov_at(cs, 'cov3_placed', ks['cov3_placed'])
        out[sc] = d
    return out


def wilcoxon_order(cases, a, b):
    pa, pb = person_mean([r for r in cases if r['scenario'] == a], 'E_shape'), \
        person_mean([r for r in cases if r['scenario'] == b], 'E_shape')
    ps = sorted(set(pa) & set(pb))
    x = np.array([pa[p] for p in ps]) - np.array([pb[p] for p in ps])
    return dict(n=len(ps), mean_diff_mm=float(x.mean()), better=int((x < 0).sum()),
                p_one_sided=float(stats.wilcoxon(x, alternative='less').pvalue))


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--in', dest='inp', nargs='+', required=True)
    ap.add_argument('--out', required=True)
    a = ap.parse_args(argv)
    cases = load(a.inp)
    by = defaultdict(list)
    for r in cases:
        by[r['task']].append(r)
    imp = [r for r in by['femur_imp'] if r['scenario'] != 'S2_placebo']
    S = {}
    S['femur_imp'] = summarize_task(imp, 'femur_imp')
    if by['femur_vsd']:
        S['femur_vsd'] = summarize_task(by['femur_vsd'], 'femur_vsd', ref_cases=imp)
    if by['femur_pool']:
        S['femur_pool'] = summarize_task(by['femur_pool'], 'femur_pool')
        for ds in ('imperial', 'vsd', 'tlem'):
            sub = [r for r in by['femur_pool'] if r.get('dataset') == ds]
            if sub:
                S[f'femur_pool[{ds}]'] = {sc: dict(E_shape_mean_mm=v['E_shape_mean_mm'], n_persons=v['n_persons'],
                                                   coverage90_kappa1=v['coverage90_kappa1'])
                                          for sc, v in summarize_task(sub, 'x', ref_cases=sub).items()}
    if by['tibia']:
        S['tibia'] = summarize_task(by['tibia'], 'tibia')
        S['tibia_order'] = {'THW<T0': wilcoxon_order(by['tibia'], 'THW', 'T0'),
                            'TH<T0': wilcoxon_order(by['tibia'], 'TH', 'T0'),
                            'TL<T0': wilcoxon_order(by['tibia'], 'TL', 'T0'),
                            'THL<TL': wilcoxon_order(by['tibia'], 'THL', 'TL')}
    S['K2_order_imperial'] = {'S2<S1': wilcoxon_order(imp, 'S2', 'S1'), 'S1<S0': wilcoxon_order(imp, 'S1', 'S0'),
                              'S3<S1': wilcoxon_order(imp, 'S3', 'S1'), 'S4<S1': wilcoxon_order(imp, 'S4', 'S1')}
    pl = [r for r in by['femur_imp'] if r['scenario'] == 'S2_placebo']
    if pl:
        k = kappa_of([r for r in imp if r['scenario'] == 'S2'], 'cov3')
        S['placebo_S2'] = dict(n=len(pl), kappa_S2=k, coverage90_at_kappa=cov_at(pl, 'cov3', k),
                               coverage90_kappa1=cov_at(pl, 'cov3', 1.0),
                               E_shape_mean_mm=float(np.mean([r['E_shape'] for r in pl])),
                               stage_a_reject_rate=float(np.mean([not r['cert_accepted'] for r in pl])))
    gen = [r for r in imp + by['femur_vsd'] if r['scenario'] in ('S1', 'S2', 'S3', 'S4', 'S5')]
    S['K4_false_reject_genuine'] = dict(n=len(gen), rate=float(np.mean([not r['cert_accepted'] for r in gen])),
                                        imperial=float(np.mean([not r['cert_accepted'] for r in gen if r['task'] == 'femur_imp'])),
                                        vsd=float(np.mean([not r['cert_accepted'] for r in gen if r['task'] == 'femur_vsd'])) if by['femur_vsd'] else None)
    S['n_cases'] = {k: len(v) for k, v in by.items()}
    S['sources'] = sorted({r['_file'] for r in cases})
    calib = {'femur_r': {sc: v['kappa_all'] for sc, v in S['femur_imp'].items()}}
    if 'tibia' in S:
        calib['tibia_r'] = {sc: v['kappa_all'] for sc, v in S['tibia'].items()}
    out = Path(a.out)
    (out / 'eval_summary.json').write_text(json.dumps(S, indent=1))
    (out / 'calibration.json').write_text(json.dumps(dict(
        source='Imperial LOSO (femur_imp) / Keast LOSO (tibia), all persons; kappa scales the covariance so pooled '
               '90 % coverage = 0.90', kappa=calib), indent=1))
    print(json.dumps({k: v for k, v in S.items() if k.startswith(('K', 'placebo', 'n_cases', 'tibia_order'))}, indent=1))


if __name__ == '__main__':
    main()
