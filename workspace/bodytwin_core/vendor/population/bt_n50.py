"""BT-N50: LOSO-kappa-kalibrering av N40b:s hoftkraftsband pa kraftniva.
Ren efterbearbetning av frysta N40b K=10-dragningar. Se PREREG.md.
Kallor: n40b_collect.py (band90_coverage), BT-N46/geomgr/bandcal.py (B86-regeln).
"""
import numpy as np, glob, json, os

S = '/mnt/shared_data/bodytwin_work/N40b/cloud_k10/solve_out/'
G = '/mnt/shared_data/bodytwin_work/N40b/cloud_k10/geom_out/'
Z90 = 1.6448536269514722
JOINTS = ['hip_r', 'hip_l', 'knee_r', 'knee_l', 'l5s1', 'gh_r', 'gh_l']

Fs = []; ids = []; var = []; geo = []
for f in sorted(glob.glob(S + 'F_shard*.npz')):
    z = np.load(f); Fs.append(z['F'])
    ids += list(map(str, z['ids'])); var += list(map(str, z['var'])); geo += list(map(str, z['geo']))
F = np.concatenate(Fs); ids = np.array(ids); var = np.array(var); geo = np.array(geo)
peak = F.max(1)
assert F.shape[1] == 141 and len(JOINTS) == F.shape[2]

SYN = ['syn%03d' % i for i in range(100)]
VSD = ['vsd_z%03d' % i for i in [1, 9, 13, 19, 23, 27, 35, 42, 46, 49]]
VSDS6 = ['vsdS6_z%03d' % i for i in [1, 9, 13, 19, 23, 27, 35, 42, 46, 49]]
COHORTS = {'syn': SYN, 'vsd': VSD, 'vsdS6': VSDS6}


def val(i, v, g, j):
    m = (ids == i) & (var == v) & (geo == g)
    return float(peak[m, j][0]) if m.any() else None


def q(x, p):
    return float(np.quantile(np.asarray(x, float), p))


def summarise(v):
    v = np.asarray(v, float)
    return dict(min=float(v.min()), median=float(np.median(v)), mean=float(v.mean()), max=float(v.max()))


def fit_scope(rows, target=0.90):
    """B86/BT-N46 geomgr.bandcal.fit_scope (copied semantics)."""
    r = [x['r'] for x in rows]
    c_train = float(np.mean([x['cov'] for x in rows]))
    b = q(r, min(c_train, 1.0))
    q_target = q(r, target)
    t = max(b, q_target)
    return dict(n=len(rows), C_train=c_train, b=b, q_target=q_target, t=t, lam=t / b, kappa=t / Z90)


def build_rows(j, cohort=None):
    rows = []
    for coh, lst in COHORTS.items():
        if cohort and coh != cohort:
            continue
        for i in lst:
            pt = val(i, 'A', 'point', j); tr = val(i, 'A', 'truth', j)
            dr = np.array([val(i, 'A', 'draw%02d' % k, j) for k in range(10)])
            if pt is None or tr is None or not np.all(np.isfinite(dr)):
                continue
            sd = dr.std(ddof=1); lo, hi = np.percentile(dr, [5, 95]); e = pt - tr
            rows.append(dict(id=i, cohort=coh, point=pt, truth=tr, e=e, sd=sd, lo=lo, hi=hi,
                             r=abs(e) / sd, cov=bool(lo <= tr <= hi),
                             cov_normal=bool(abs(e) <= Z90 * sd),
                             draw_z=list(np.abs(dr - dr.mean()) / sd)))
    return rows


def loso(rows):
    """BT-N46 bandcal.leave_one_subject_out semantics, individuals as subjects."""
    subjects = sorted({x['id'] for x in rows})
    lam, kap, covered = {}, {}, 0
    for s in subjects:
        fs = fit_scope([x for x in rows if x['id'] != s])
        lam[s], kap[s] = fs['lam'], fs['kappa']
        covered += sum(1 for x in rows if x['id'] == s and x['r'] <= fs['t'])
    return dict(n=len(rows), n_subjects=len(subjects),
                coverage_before=float(np.mean([x['cov'] for x in rows])),
                coverage_before_normal=float(np.mean([x['cov_normal'] for x in rows])),
                covered_after=int(covered), coverage_after=covered / len(rows),
                lam=summarise(list(lam.values())), kappa=summarise(list(kap.values())),
                kappa_by_subject=kap)


def band_stats(rows, kappa_by_subject=None, kappa=None):
    emp, nor, after = [], [], []
    for x in rows:
        a = abs(x['point'])
        emp.append((x['hi'] - x['lo']) / a)
        nor.append(2 * Z90 * x['sd'] / a)
        k = kappa if kappa is not None else kappa_by_subject[x['id']]
        after.append(2 * k * Z90 * x['sd'] / a)
    return dict(empirical_before_pct_median=float(100 * np.median(emp)),
                normal_before_pct_median=float(100 * np.median(nor)),
                calibrated_after_pct_median=float(100 * np.median(after)),
                width_change_vs_empirical_pct=float(100 * (np.median(after) / np.median(emp) - 1)))


def geom_ratio(ids_list):
    pr, sdg = [], []
    for i in ids_list:
        p = G + 'geom_%s.json' % i
        if not os.path.exists(p):
            continue
        d = json.load(open(p))['delta_attach_mm']
        pr.append(d['point_vs_truth_median']); sdg.append(d['draw_sd_median'])
    pr, sdg = np.array(pr), np.array(sdg)
    return dict(n=len(pr), point_vs_truth_median_mm=float(np.median(pr)),
                draw_sd_median_mm=float(np.median(sdg)), ratio_median=float(np.median(pr / sdg)))


results = dict(job='BT-N50', z90=Z90, n_rows=int(len(F)), n_individuals=120,
               joints=JOINTS, sources=dict(F=S, geom=G))

for j, jn in enumerate(JOINTS):
    rows = build_rows(j)
    if not rows:
        continue
    common = loso(rows)
    full = fit_scope(rows)
    per_coh = {}
    for coh in COHORTS:
        cr = [x for x in rows if x['cohort'] == coh]
        r = loso(cr)
        r.update(band_stats(cr, kappa_by_subject=r['kappa_by_subject']))
        per_coh[coh] = r
    # tail diagnostic: actual r vs the draws' own 90th percentile
    zd = np.concatenate([np.array(x['draw_z']) for x in rows])
    rr = np.array([x['r'] for x in rows])
    tail = dict(actual_r_q50=float(np.percentile(rr, 50)), actual_r_q90=float(np.percentile(rr, 90)),
                actual_r_q95=float(np.percentile(rr, 95)),
                draw_z_q90=float(np.percentile(zd, 90)), draw_z_q95=float(np.percentile(zd, 95)),
                frac_r_gt_draw_q90=float(np.mean(rr > np.percentile(zd, 90))))
    # geometry vs measurement, force level (relative)
    e_rel = np.array([x['e'] / abs(x['point']) for x in rows])
    sd_rel = np.array([x['sd'] / abs(x['point']) for x in rows])
    var_e = float(np.var(e_rel, ddof=1)); mean_sd2 = float(np.mean(sd_rel ** 2))
    e2 = float(np.mean(e_rel ** 2))
    # placebo (syn): wrong individual's measures; no placebo draws exist, compare to nominal band
    pl = []
    for i in SYN:
        pt = val(i, 'A', 'placebo', j); tr = val(i, 'A', 'truth', j)
        dr = np.array([val(i, 'A', 'draw%02d' % k, j) for k in range(10)])
        if pt is None or not np.all(np.isfinite(dr)):
            continue
        pl.append(dict(err=abs(pt - tr) / abs(pt), cov_normal=bool(abs(pt - tr) <= Z90 * dr.std(ddof=1))))
    # variant controls (point-truth error, no draws)
    ctrl = {}
    for v in ('A', 'B', 'C'):
        er = []
        for lst in COHORTS.values():
            for i in lst:
                pt = val(i, v, 'point', j); tr = val(i, v, 'truth', j)
                if pt is not None and tr is not None:
                    er.append(abs(pt - tr) / abs(pt))
        ctrl[v] = dict(point_err_median_pct=float(100 * np.median(er)), n=len(er))
    results[jn] = dict(
        common=common,
        full_data=dict(kappa=full['kappa'], lam=full['lam'], t=full['t'], b=full['b'],
                       q_target=full['q_target'], C_train=full['C_train']),
        deploy_band_stats=band_stats(rows, kappa=full['kappa']),
        loso_band_stats=band_stats(rows, kappa_by_subject=common['kappa_by_subject']),
        per_cohort=per_coh,
        tail=tail,
        placebo=dict(n=len(pl), coverage_normal=float(np.mean([x['cov_normal'] for x in pl])),
                     err_median_pct=float(100 * np.median([x['err'] for x in pl]))),
        geometry_vs_measurement=dict(
            var_e_rel=var_e, mean_sd2_rel=mean_sd2, E_e2_rel=e2,
            geometry_variance_share=max(0.0, 1.0 - mean_sd2 / e2) if e2 > 0 else None,
            e_median_abs_pct=float(100 * np.median(np.abs(e_rel))),
            sd_median_pct=float(100 * np.median(sd_rel)),
            geom=geom_ratio(SYN + VSD + VSDS6)),
        variant_point_err=ctrl)

# derived: decomposition of the undercoverage gap for the primary joint
for jn in JOINTS:
    if jn not in results:
        continue
    c = results[jn]['common']; t = results[jn]['tail']
    ideal_meas = 1.0 - t['frac_r_gt_draw_q90']
    gap = 0.90 - c['coverage_before']
    meas_part = ideal_meas - c['coverage_before']
    geo_part = 0.90 - ideal_meas
    results[jn]['derived'] = dict(
        ideal_measurement_band_coverage=ideal_meas,
        undercoverage_gap_to_nominal=gap,
        gap_from_measurement_band_construction=meas_part,
        gap_from_geometry_tail=geo_part,
        measurement_share_of_gap=(meas_part / gap) if gap > 0 else None,
        geometry_share_of_gap=(geo_part / gap) if gap > 0 else None)

json.dump(results, open('results.json', 'w'), indent=1, default=float)

print('=== BT-N50 hip_r (variant A) ===')
for jn in ('hip_r', 'hip_l', 'knee_r'):
    r = results[jn]; c = r['common']; bs = r['loso_band_stats']
    print(f"{jn}: before_emp={c['coverage_before']:.3f} before_normal={c['coverage_before_normal']:.3f} "
          f"after={c['coverage_after']:.3f} kappa_med={c['kappa']['median']:.3f} lam_med={c['lam']['median']:.3f}")
    print(f"   width emp {bs['empirical_before_pct_median']:.2f}% -> cal {bs['calibrated_after_pct_median']:.2f}% "
          f"({bs['width_change_vs_empirical_pct']:+.1f}%)  normal {bs['normal_before_pct_median']:.2f}%")
    for coh in ('syn', 'vsd', 'vsdS6'):
        pc = r['per_cohort'][coh]
        print(f"   {coh:6s} before={pc['coverage_before']:.3f} after={pc['coverage_after']:.3f} "
              f"kappa_med={pc['kappa']['median']:.3f}")
    print('   tail:', {k: round(v, 3) for k, v in r['tail'].items()})
    gm = r['geometry_vs_measurement']
    print(f"   geom ratio {gm['geom']['ratio_median']:.3f} (p_vs_t {gm['geom']['point_vs_truth_median_mm']:.2f}mm / "
          f"draw_sd {gm['geom']['draw_sd_median_mm']:.2f}mm) geo_var_share={gm['geometry_variance_share']}")
    print('   variant point err:', {k: round(v['point_err_median_pct'], 3) for k, v in r['variant_point_err'].items()})
