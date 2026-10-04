"Geometry pipeline with numeric evaluation and optional local OpenSim and signed-distance-field exports."
from __future__ import annotations

import json
import os
import time
from pathlib import Path

import numpy as np

from . import core as C
from . import identity as I
from . import leg as LG
from . import measures as MS
from . import ops as OP
from . import certify2 as C2
from .api import GeometryManager

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DATA = Path(os.environ.get('GEOMGR_DATA', 'external_media'))
EXTRA = Path(os.environ.get('GEOMGR_EXTRA', ROOT / 'frozen_extra'))
SCEN = {'S2': (C.X5, []), 'S4': ([], ['HC6', 'SGT', 'LT', 'MEC', 'LEC']), 'S6': (C.X5, ['HC6', 'SGT', 'LT', 'MEC', 'LEC'])}
DEMO10 = ['z009', 'z013', 'z019', 'z023', 'z027', 'z035', 'z036', 'z042', 'z046', 'z049']


class Timer:
    def __init__(self):
        self.t = {}

    def __call__(self, name):
        tm = self

        class _C:
            def __enter__(s):
                s.w, s.c = time.time(), time.process_time()

            def __exit__(s, *a):
                tm.t[name] = dict(wall_s=round(time.time() - s.w, 4), cpu_s=round(time.process_time() - s.c, 4))
        return _C()


def _js(o):
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, (np.floating, np.integer, np.bool_)):
        return o.item()
    if isinstance(o, set):
        return sorted(o)
    return str(o)


def kappa_table():
    f = EXTRA / 'calibration_N7a.json'
    if not f.exists():
        f = ROOT.parent / 'N7a' / 'calibration.json'
    k = json.loads(f.read_text())['kappa']['femur_r']
    return {sc: v.get('cov3_placed', v.get('cov3')) for sc, v in k.items()}


def observations(O, s, obs, rating=None):
    r0 = int(np.flatnonzero(np.isfinite(O['obs_lm'][s, :, 0, 0]))[0]) if rating is None else rating
    fn, ln = SCEN[obs]
    feats = {k: float(O['obs_feat'][s, r0, C.FEMUR_FEATURES.index(k)]) for k in fn}
    lms = {k: O['obs_lm'][s, r0, C.PT_NAMES.index(k)] for k in ln}
    return r0, feats, lms


def tlem_scaled_points(leg, gm, L_mech, W_epi):
    "landmark-based baseline: TLEM femur scaled about the TLEM knee centre, diag(s_w, s_l, s_w) in the TLEM femur"
    P = np.load(DATA / 'femur_pop.npz')
    k = leg.d['tlem_member']['index']
    R, t = np.array(leg.d['tlem_member']['R_native_to_tlem']), np.array(leg.d['tlem_member']['t_native_to_tlem'])
    V = P['S'][k] @ R.T + t
    vals = I.evaluate(gm.registry, V)
    L = MS.landmarks_of(vals)
    m = C.measures(L, V)
    sl, sw = L_mech / m['L_mech'], W_epi / m['W_epi']
    kc = leg.J['Knee|Femur']['center']
    S = np.diag([sw, sl, sw])
    sc = lambda x: kc + (np.asarray(x, float) - kc) @ S.T  # noqa: E731
    fp = {'_hjc': sc(leg.J['Hip|Femur']['center'])}
    for n, x in LG.tlem_construct_points(leg).items():
        fp[n] = sc(x)
    for el in sum(LG.GROUPS.values(), []):
        for kk, (kind, seg, xyz) in enumerate(leg.muscles[el]['points']):
            if seg == 'femur':
                key = f'{el.replace(" ", "_")}|{kind}' if kind in ('Origin', 'Insertion') else f'TLEM.via.{el.replace(" ", "_")}.{kk}'
                fp[key] = sc(xyz)
    return fp, dict(s_l=sl, s_w=sw, tlem_L_mech=m['L_mech'], tlem_W_epi=m['W_epi'])


def identity_chain(gm, res):
    inst = res.instance
    out = {}
    d = gm.sample(res, 1, seed=1)[0]
    out['posterior_draw'] = I.check_identity(inst, d, exact=False)
    m = I.op_mirror(inst, inst.reg.part.replace('_r', '_l'))
    out['mirror'] = I.check_identity(inst, m, transform=lambda x: x * np.array([-1.0, 1, 1]), tol=1e-9)
    mm = I.op_mirror(m, inst.reg.part)
    out['mirror_twice_bitwise'] = dict(ok=bool(np.array_equal(mm.V, inst.V) and np.array_equal(mm.F, inst.F)
                                               and mm.reg.manifest_hash() == inst.reg.manifest_hash()))
    fr = inst.values()[f'{inst.reg.part}/frame/ISB']
    Rm, t = fr['axes'], -fr['axes'] @ fr['origin']
    r = I.op_rigid(inst, Rm, t, 'ISB')
    out['rigid_to_ISB'] = I.check_identity(inst, r, transform=lambda x: x @ Rm.T + t, tol=1e-9)
    u = I.op_units(inst, 'm')
    out['units_mm_to_m'] = I.check_identity(inst, u, transform=lambda x: x * 1e-3, tol=1e-12)
    s = I.op_subdivide(inst)
    out['subdivide'] = I.check_identity(inst, s, tol=1e-9)
    V2 = I.laplacian_smooth(s.V, s.F, iters=2, lam=0.3)
    rt, disp = I.op_retopologize(inst, V2, s.F)
    ck = I.check_identity(inst, rt, exact=False)
    ck['max_anchor_displacement_mm'] = float(max(disp.values()))
    ck['version_bumped'] = all(rt.reg.entities[b].version == inst.reg.entities[b].version + 1
                               for b in inst.reg.bases()
                               if inst.reg.entities[b].kind in ('landmark', 'attachment', 'region', 'construct', 'wrap'))
    out['retopologize_smoothed'] = ck
    n_ids = len(inst.reg.entities)
    lost = sum(len(v.get('lost', [])) for v in out.values())
    exact = [k for k in ('mirror', 'rigid_to_ISB', 'units_mm_to_m', 'subdivide')]
    return dict(chain=out, n_ids=n_ids, total_lost=lost, exact_max_dev=max(out[k]['max_point_dev'] for k in exact),
                exact_ok=all(out[k]['ok'] for k in exact))


def numeric(subject, obs='S6', out=None, draws=50, register=False, all23=False, seed=0):
    tm = Timer()
    out = Path(out or ROOT / 'pipeline' / subject)
    out.mkdir(parents=True, exist_ok=True)
    rep = dict(subject=subject, obs=obs, draws=draws, measure_definition=MS.MEASURE_VERSION)
    kap = kappa_table()
    with tm('1_load_population'):
        cache = os.environ.get('GEOMGR_CACHE', 'external_media')
        gm = GeometryManager(DATA, population=('imperial',), exclude_vsd=subject, kappa=kap,
                             cache_dir=cache if Path(cache).is_dir() else EXTRA)
    leg = LG.LegFrozen.load(EXTRA / 'tlem_leg.json')
    O = {k: v for k, v in np.load(DATA / 'vsd_obs.npz').items()}
    subj = [str(x) for x in O['subj']]
    s = subj.index(subject)
    r0, feats, lms = observations(O, s, obs)
    rep['observations'] = dict(rating_index=r0, features=feats, landmarks={k: v.tolist() for k, v in lms.items()})
    with tm('2_instantiate'):
        res = gm.instantiate(features=feats, feature_units={k: C.feature_units(k) for k in feats}, landmarks=lms,
                             subject=subject, scenario=obs)
    rep['instantiate'] = res.summary()
    with tm('3_outlier_B35'):
        o = {}
        if lms and len(lms) >= 4:
            o['observed'] = OP.op_outlier(gm, lms, list(lms))
        if all23:
            full = {n: O['obs_lm'][s, r0, j] for j, n in enumerate(C.PT_NAMES)}
            o['ALL23'] = OP.op_outlier(gm, full, list(C.PT_NAMES))
            shifted = dict(full)
            _, ax = C.knee_frame(full)
            shifted['SGT'] = full['SGT'] + 30.0 * ax[1]
            o['ALL23_SGT_plus30mm_anterior'] = OP.op_outlier(gm, shifted, list(C.PT_NAMES))
            if lms:
                l5 = {k: (shifted[k] if k == 'SGT' else v) for k, v in lms.items()}
                o['L5_SGT_plus30mm_anterior'] = OP.op_outlier(gm, l5, list(lms))
        rep['outlier'] = o
    mu_inst = I.Instance(gm.registry, gm.model.mu, gm.F)
    with tm('4_certificate_v2_numeric'):
        shape_inst = I.Instance(res.instance.reg, res.X_shape, gm.F, subject=subject)
        rep['certificate_v2'] = C2.instance_certificate(mu_inst, shape_inst, 'instantiate', local=False)
        rep['certificate_v1_numeric_N7a'] = res.certificate
    with tm('5_morphometrics'):
        L = MS.landmarks_of(res.instance.values())
        rep['morphometrics'] = dict(chosen=MS.ccd_av(L, res.instance.V, gm.F),
                                    b42_signed=MS.ccd_av(L, res.instance.V, gm.F, 'b42'))
    with tm('6_edit_AV+10'):
        ed, erep = OP.op_edit_osteotomy(res.instance, 'AV', 10.0)
        rep['edit'] = erep
    with tm('7_tissue'):
        th, src = None, 'UNKNOWN (no CT thickness for this individual)'
        f = EXTRA / f'e_layers_{subject}.npz'
        if f.exists():
            th, src = np.load(f)['thickness_mm'], f'CT half-max thickness of {subject} (results/N7b e_regions, reused)'
        tis, trep = OP.op_tissue(res.instance, th, src)
        rep['tissue'] = trep
    with tm('8_wraps'):
        rep['wraps'] = OP.wraps_summary(res.instance)
    with tm('9a_arms_point_and_band'):
        band, A = LG.arm_band(gm, res, leg, n=draws, seed=seed, kappa=res.kappa or 1.0)
    with tm('9b_arms_truth_baseline_edit_placebo'):
        names = [str(x) for x in np.load(DATA / 'femur_pop.npz')['names']]
        Vt = np.load(DATA / 'femur_pop.npz')['S'][names.index(f'vsd:{subject}')]
        vt = I.evaluate(gm.registry, Vt)
        truth = LG.moment_arms(vt, leg, V=Vt)
        truth_nowrap = LG.moment_arms(vt, leg, V=Vt, wraps=False)
        edited = LG.moment_arms(ed.values(), leg, V=ed.V)
        Lm = feats.get('L_mech', C.measures(L, res.instance.V)['L_mech'])
        We = feats.get('W_epi', C.measures(L, res.instance.V)['W_epi'])
        fp, sinfo = tlem_scaled_points(leg, gm, Lm, We)
        base = LG.moment_arms(None, leg, femur_points=fp)
        # S2 (measures only) and placebo (S2 measures of the next subject)
        r2, f2, _ = observations(O, s, 'S2')
        res2 = gm.instantiate(features=f2, feature_units={k: C.feature_units(k) for k in f2}, subject=subject, scenario='S2')
        band2, _ = LG.arm_band(gm, res2, leg, n=draws, seed=seed + 1, kappa=res2.kappa or 1.0)
        sp = (s + 1) % len(subj)
        _, fpl, _ = observations(O, sp, 'S2')
        resp = gm.instantiate(features=fpl, feature_units={k: C.feature_units(k) for k in fpl}, subject=subject, scenario='S2')
        plac = LG.moment_arms(I.evaluate(gm.registry, resp.X_shape), leg, V=resp.X_shape)
        fpp, _ = tlem_scaled_points(leg, gm, fpl['L_mech'], fpl['W_epi'])
        base_plac = LG.moment_arms(None, leg, femur_points=fpp)
        mean_shape = LG.moment_arms(I.evaluate(gm.registry, gm.model.mu), leg, V=gm.model.mu)
    rep['arms'] = dict(quantities=[LG.qname(*q) for q in LG.QUANTITIES], band_S6=band, band_S2=band2, truth=truth,
                       truth_nowrap=truth_nowrap, edited_AV10=edited, baseline_tlem_scaled=base, baseline_scaling=sinfo,
                       placebo_S2_next_subject=dict(subject=subj[sp], geomgr=plac, baseline=base_plac),
                       population_mean_shape=mean_shape, kappa=res.kappa, kappa_S2=res2.kappa)
    with tm('10_identity_chain'):
        rep['identity'] = identity_chain(gm, res)
    if register:
        with tm('11_register_surface'):
            dz = np.load(EXTRA / f'vsd_{subject}.npz')
            ri, rrep = OP.op_register_surface(gm, dz['dense'], dz['normals'], subject)
            anc = gm.anchors.dict(ri.V)
            tr = gm.anchors.dict(Vt)
            rrep['landmark_error_vs_P1_registration_mm'] = dict(
                median=float(np.median([np.linalg.norm(anc[n] - tr[n]) for n in C.LM_NAMES])),
                max=float(np.max([np.linalg.norm(anc[n] - tr[n]) for n in C.LM_NAMES])))
            m20 = O['obs_lm'][s].mean(0) if np.all(np.isfinite(O['obs_lm'][s])) else None
            if m20 is not None:
                e = [np.linalg.norm(anc[n] - m20[j]) for j, n in enumerate(C.LM_NAMES)]
                rrep['landmark_error_vs_20rating_mean_mm'] = dict(median=float(np.median(e)), p90=float(np.percentile(e, 90)))
            rep['register'] = rrep
            np.save(out / 'registered_V.npy', ri.V)
    rep['timings'] = tm.t
    np.savez(out / 'instance.npz', V=res.instance.V, X_shape=res.X_shape, F=gm.F, V_edit=ed.V,
             sigma=res.sigma_mm, draws_S6=A)
    (out / 'registry.json').write_text(json.dumps(res.instance.reg.to_json()))
    (out / 'numeric.json').write_text(json.dumps(rep, indent=1, default=_js))
    return rep


def local(subject, out=None, sdf=None, validate_full=None):
    """Local stage on a saved numeric result (certificate v1 + LV1, exports + id checks). Timed."""
    from types import SimpleNamespace
    from . import export as X
    from . import certify as CE
    tm = Timer()
    out = Path(out or ROOT / 'pipeline' / subject)
    num = json.loads((out / 'numeric.json').read_text())
    d = np.load(out / 'instance.npz')
    F = d['F']
    reg = I.Registry.from_json(json.loads((out / 'registry.json').read_text()), faces=F)
    inst = I.Instance(reg, d['V'], F, frame=num['instantiate']['frame'], subject=subject)
    th = EXTRA / f'e_layers_{subject}.npz'
    inst, _ = OP.op_tissue(inst, np.load(th)['thickness_mm'] if th.exists() else None, 'see numeric.json')
    rep = dict(subject=subject)
    with tm('L1_certificate_v1_bodytwin'):
        rep['bodytwin_v1'] = CE.geometry_certificate_local(inst.V, F)
    with tm('L2_certificate_LV1_interval'):
        rep['lv1_interval'] = CE.interval_volume(inst.V, F)
    with tm('L3_certificate_v2_local_edit'):
        ed = inst.copy()
        ed.V = d['V_edit']
        rep['certificate_v2_edit_local'] = C2.instance_certificate(inst, ed, 'edit_osteotomy:AV+10', local=True)
    res = SimpleNamespace(instance=inst, subject=subject, frame=inst.frame, accepted=num['instantiate']['accepted'],
                          certificate=num['certificate_v2'], kappa=num['arms']['kappa'], points=None,
                          sigma_mm=d['sigma'])
    ids_all = {e.id for e in reg.entities.values()}
    ex = {}
    with tm('L4_export_json'):
        X.to_json(res, out / f'{subject}_instance.json')
        got = X.ids_in_json(out / f'{subject}_instance.json')
    ex['json'] = dict(expected=len(ids_all), found=len(got & ids_all), lost=sorted(ids_all - got))
    pts = set(X.point_ids(res))
    with tm('L8_export_opensim_and_load'):
        eo = X.to_opensim(res, out / 'opensim')
        co = X.check_opensim([eo])
        lost = co['rows'][0].get('lost', []) if co.get('rows') else ['<load failed>']
        ex['opensim'] = dict(n_markers=eo['n_markers'], load_ok=co.get('ok'), ids_expected=len(pts), ids_lost=lost)
    if sdf if sdf is not None else subject == 'z001':
        with tm('L9_export_sdf_field_engine'):
            ex['sdf'] = X.to_sdf([(f'{subject}_femur_r', inst.V, F)], out / 'sdf')
    rep['exports'] = ex
    rep['timings'] = tm.t
    (out / 'local.json').write_text(json.dumps(rep, indent=1, default=_js))
    return rep


def main(argv):
    import argparse
    ap = argparse.ArgumentParser(prog='geomgr pipeline')
    ap.add_argument('--subject', default='z001')
    ap.add_argument('--subjects', default=None)
    ap.add_argument('--obs', default='S6')
    ap.add_argument('--stage', choices=('numeric', 'local', 'all'), default='all')
    ap.add_argument('--draws', type=int, default=50)
    ap.add_argument('--register', action='store_true')
    ap.add_argument('--all23', action='store_true')
    ap.add_argument('--out', default=None)
    a = ap.parse_args(argv)
    subs = (DEMO10 if a.subjects == 'demo10' else a.subjects.split(',')) if a.subjects else [a.subject]
    base = Path(a.out) if a.out else ROOT / 'pipeline'
    summ = {}
    for sid in subs:
        o = base / sid
        if a.stage in ('numeric', 'all'):
            r = numeric(sid, a.obs, o, a.draws, register=a.register, all23=a.all23 or sid == 'z001')
            summ[sid] = dict(accepted=r['instantiate']['accepted'], lost_ids=r['identity']['total_lost'],
                             timings={k: v['wall_s'] for k, v in r['timings'].items()})
        if a.stage in ('local', 'all'):
            r = local(sid, o)
            summ.setdefault(sid, {})['local'] = {k: v['wall_s'] for k, v in r['timings'].items()}
        print(json.dumps({sid: summ[sid]}, default=_js), flush=True)
    return 0
