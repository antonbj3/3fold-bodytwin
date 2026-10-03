"Batch geometry processing with optional local OpenSim and signed-distance-field exports."
from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np

from . import core as C
from . import certify as CE
from . import export as X
from . import identity as I
from .api import GeometryManager, Result
from .demo import identity_chain

L5 = ['HC6', 'SGT', 'LT', 'MEC', 'LEC']
SEED = 424242


def _rot(rng, deg):
    w = rng.normal(size=3)
    w = w / np.linalg.norm(w) * np.radians(rng.uniform(-deg, deg))
    return C._expm_so3(w)


def individuals(gm, data, n):
    O = np.load(Path(data) / 'vsd_obs.npz')
    valid = np.isfinite(O['obs_lm'][:, :, 0, 0])
    pairs = [(s, r) for s in range(valid.shape[0]) for r in np.flatnonzero(valid[s])]
    rng = np.random.default_rng(SEED)
    out = []
    for i in range(n):
        b = rng.standard_normal(gm.model.r) * np.sqrt(gm.model.lam)
        R, t = _rot(rng, 20), rng.uniform(-50, 50, 3)
        V = gm.model.shape(b) @ R.T + t
        s, r = pairs[rng.integers(len(pairs))]
        fn = O['obs_feat'][s, r] - O['feat_mean20'][s]
        dev = O['dev_anat'][s, r]
        P = C.femur_points(gm.anchors, V)
        f = C.femur_features_from_landmarks(P, V)
        feats = {k: float(f[C.FEMUR_FEATURES.index(k)] + fn[C.FEMUR_FEATURES.index(k)]) for k in C.X5}
        _, Ax = C.knee_frame(P)
        lms = {k: P[k] + dev[C.PT_NAMES.index(k)] @ Ax for k in L5}
        out.append(dict(subject=f'syn{i:03d}', truth=V, features=feats, landmarks=lms))
    return out


def stage_numeric(data, out, n, kappa=None):
    out = Path(out)
    (out / 'ind').mkdir(parents=True, exist_ok=True)
    T = {}
    w0, c0 = time.time(), time.process_time()
    gm = GeometryManager(data, population=('imperial',), kappa=kappa, cache_dir=out)
    T['load_population_and_fit'] = dict(wall_s=time.time() - w0, cpu_s=time.process_time() - c0, per='once')
    inds = individuals(gm, data, n)
    steps = {}
    rows = []
    for ind in inds:
        st = {}

        def tick(name, fn):
            w, c = time.time(), time.process_time()
            v = fn()
            st[name] = dict(wall_s=time.time() - w, cpu_s=time.process_time() - c)
            return v
        res = tick('instantiate', lambda: gm.instantiate(features=ind['features'],
                                                         feature_units={k: C.feature_units(k) for k in ind['features']},
                                                         landmarks=ind['landmarks'], subject=ind['subject'],
                                                         scenario='S6'))
        lv1 = tick('certificate_LV1_interval', lambda: CE.interval_volume(res.instance.V, res.instance.F))
        idc = tick('identity_chain', lambda: identity_chain(gm, res, out))
        jp = out / 'ind' / f"{ind['subject']}.json"
        tick('export_json', lambda: X.to_json(res, jp))
        got = X.ids_in_json(jp)
        ids_all = {e.id for e in res.instance.reg.entities.values()}
        eo = tick('export_opensim_text', lambda: X.to_opensim(res, out / 'opensim'))
        e = ind['truth'] - res.instance.V
        z2 = np.einsum('ni,nij,nj->n', e, np.linalg.inv(res.cov), e)
        np.savez_compressed(out / 'ind' / f"{ind['subject']}.npz", V=res.instance.V.astype(np.float64),
                            sigma=res.sigma_mm.astype(np.float32),
                            pts_names=np.array(list(res.points)), pts=np.array([p for p, _ in res.points.values()]),
                            pts_cov=np.array([c for _, c in res.points.values()]))
        lost = {k: v.get('n_lost', 0) for k, v in idc.items() if isinstance(v, dict) and 'n_lost' in v}
        rows.append(dict(subject=ind['subject'], accepted=res.accepted, reasons=res.certificate['reasons'],
                         E_placed_mm=float(np.sqrt((e ** 2).sum(1).mean())),
                         coverage90_placed=float(np.mean(z2 <= C.CHI2_3_90)), pose_identified=res.pose_identified,
                         lv1_float_inside=lv1['float_inside'], lv1_width_mm3=lv1['width_mm3'],
                         identity_ok={k: bool(v.get('ok', True)) for k, v in idc.items()}, identity_lost=lost,
                         identity_max_dev=max(v.get('max_point_dev', 0.0) for v in idc.values()),
                         retopo_max_disp_mm=idc['retopologize_smoothed']['max_anchor_displacement_mm'],
                         json_lost=len(ids_all - got), osim=eo['osim'], osim_expected=eo['expected'],
                         timings=st))
        for k, v in st.items():
            steps.setdefault(k, []).append(v)
    T['per_step'] = {k: dict(sum_wall_s=float(sum(x['wall_s'] for x in v)), sum_cpu_s=float(sum(x['cpu_s'] for x in v)),
                             first_wall_s=float(v[0]['wall_s']), first_cpu_s=float(v[0]['cpu_s']), n=len(v))
                     for k, v in steps.items()}
    rep = dict(stage='numeric', n=n, timings=T, rows=rows,
               totals=dict(accepted=int(sum(r['accepted'] for r in rows)),
                           identity_ids_lost_total=int(sum(sum(r['identity_lost'].values()) + r['json_lost'] for r in rows)),
                           identity_all_ok=bool(all(all(r['identity_ok'].values()) for r in rows)),
                           max_exact_dev=float(max(r['identity_max_dev'] for r in rows)),
                           lv1_all_inside=bool(all(r['lv1_float_inside'] for r in rows)),
                           E_placed_mean_mm=float(np.mean([r['E_placed_mm'] for r in rows])),
                           coverage90_placed_mean=float(np.mean([r['coverage90_placed'] for r in rows]))),
               registry_manifest=gm.registry.manifest_hash(), kappa=gm.kappa)
    (out / 'batch_numeric.json').write_text(json.dumps(rep, indent=1, default=float))
    return rep


def _result_from_saved(gm, npz, subject):
    d = np.load(npz)
    inst = I.Instance(gm.registry, d['V'], gm.F, frame=f'observation@{subject}', subject=subject)
    pts = {str(n): (p, c) for n, p, c in zip(d['pts_names'], d['pts'], d['pts_cov'])}
    return Result(subject=subject, accepted=True, certificate={}, instance=inst, sigma_mm=d['sigma'], cov=None,
                  points=pts, frame=inst.frame, kappa=None, pose_identified=True, post=None, timings={})


def stage_local(data, out, n, full_validate_n=10):
    out = Path(out)
    num = json.loads((out / 'batch_numeric.json').read_text())
    gm = GeometryManager(data, population=('imperial',), cache_dir=out)
    T = {}
    res = [_result_from_saved(gm, out / 'ind' / f"{r['subject']}.npz", r['subject']) for r in num['rows'][:n]]
    import resource

    def timed(name, fn, per=None):
        w, c = (time.time(), time.process_time())
        v = fn()
        T[name] = dict(wall_s=time.time() - w, cpu_s=time.process_time() - c, n=per, maxrss_mb_self=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, maxrss_mb_children=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / 1024)
        (out / 'batch_local_progress.json').write_text(json.dumps(T, indent=1))
        return v
    geo = timed('certificate_bodytwin_geometry', lambda: [CE.geometry_certificate_local(r.instance.V, r.instance.F) for r in res], len(res))
    exps = [dict(osim=str(out / 'opensim' / Path(r['osim']).name), expected=r['osim_expected']) for r in num['rows'][:n]]
    osim = timed('load_opensim_4.6', lambda: X.check_opensim(exps), len(exps))
    sdf = timed('export_sdf_field_engine', lambda: X.to_sdf([(r.subject, r.instance.V, r.instance.F) for r in res], out / 'sdf'), len(res))
    rep = dict(stage='local', n=len(res), timings=T, bodytwin_geometry_accepted=int(sum((g['accepted'] is True for g in geo))), opensim_ok=int(sum((x['ok'] for x in osim.get('rows', [])))), opensim_ids_lost_total=int(sum((len(x.get('lost', [])) for x in osim.get('rows', [])))), opensim_max_marker_dev_m=float(max((x.get('max_marker_dev_m', np.inf) for x in osim.get('rows', [])), default=np.inf)), sdf_ok=sdf.get('ok'), sdf_n=sdf.get('n'), sdf_volume_rel_err_max=float(max((x['volume_relative_error'] for x in sdf.get('rows', [{'volume_relative_error': np.nan}])))), sdf_rows_sample=sdf.get('rows', [])[:2])
    (out / 'batch_local.json').write_text(json.dumps(rep, indent=1, default=float))
    return rep
