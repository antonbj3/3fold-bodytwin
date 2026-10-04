"Geometry demonstration with numeric checks and optional local exports."
from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np

from . import core as C
from . import certify as CE
from . import export as X
from . import identity as I
from .api import GeometryManager

SCEN = {'S2': (C.X5, []), 'S4': ([], ['HC6', 'SGT', 'LT', 'MEC', 'LEC']), 'S6': (C.X5, ['HC6', 'SGT', 'LT', 'MEC', 'LEC']),
        'S1': (['L_mech'], []), 'S3': (C.T5, [])}


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


def load_calibration(root):
    f = Path(root) / 'calibration.json'
    if not f.exists():
        return {}
    k = json.loads(f.read_text())['kappa'].get('femur_r', {})
    return {sc: v for sc, v in k.items()}


def identity_chain(gm, res, tmpdir):
    """All operations on one instance; returns per-op identity checks (lost ids, max deviation)."""
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
                               for b in inst.reg.bases() if inst.reg.entities[b].kind in ('landmark', 'attachment', 'region'))
    out['retopologize_smoothed'] = ck
    return out


def export_roundtrips(res, outdir, weight_kg, validate_full=True, tm=None):
    tm = tm or Timer()
    outdir = Path(outdir)
    ids_all = {e.id for e in res.instance.reg.entities.values()}
    pts = set(X.point_ids(res))
    r = {}
    with tm('export_json'):
        X.to_json(res, outdir / f'{res.subject}_instance.json')
        got = X.ids_in_json(outdir / f'{res.subject}_instance.json')
    r['json'] = dict(expected=len(ids_all), found=len(got & ids_all), lost=sorted(ids_all - got))
    with tm('export_opensim'):
        eo = X.to_opensim(res, outdir / 'opensim')
    with tm('load_opensim_4.6'):
        co = X.check_opensim([eo])
    idm = json.loads(Path(eo['osim'].replace('.osim', '.idmap.json')).read_text())
    got = {idm[n] for n in co['rows'][0].get('markers_found', [])} if False else None
    lost = co['rows'][0].get('lost', []) if co.get('rows') else ['<load failed>']
    r['opensim'] = dict(osim=eo['osim'], n_markers=eo['n_markers'], mass_kg=eo['mass_kg'], load=co, ids_expected=len(pts), ids_lost=[idm.get(x, x) for x in lost])
    return (r, tm)


def run(subject='z001', obs='S6', out=None, data=None, validate_full=True):
    root = Path(__file__).resolve().parents[1]
    data = Path(data or 'external_media')
    out = Path(out or root / f'demo_{subject}')
    out.mkdir(parents=True, exist_ok=True)
    tm = Timer()
    kappa = {sc: v.get('cov3_placed', v.get('cov3')) for sc, v in load_calibration(root).items()}
    with tm('load_population_and_fit'):
        gm = GeometryManager(data, population=('imperial',), exclude_vsd=subject, kappa=kappa,
                             cache_dir='external_media')
    O = np.load(data / 'vsd_obs.npz')
    subj = [str(x) for x in O['subj']]
    s = subj.index(subject)
    r0 = int(np.flatnonzero(np.isfinite(O['obs_lm'][s, :, 0, 0]))[0])
    fn, ln = SCEN[obs]
    feats = {k: float(O['obs_feat'][s, r0, C.FEMUR_FEATURES.index(k)]) for k in fn}
    lms = {k: O['obs_lm'][s, r0, C.PT_NAMES.index(k)] for k in ln}
    with tm('instantiate'):
        res = gm.instantiate(features=feats, feature_units={k: C.feature_units(k) for k in feats}, landmarks=lms,
                             subject=subject, scenario=obs)
    cert = dict(res.certificate)
    with tm('certificate_bodytwin_geometry'):
        cert['bodytwin_geometry_certificate'] = CE.geometry_certificate_local(res.instance.V, res.instance.F)
    with tm('certificate_LV1_interval_volume'):
        cert['lv1_interval'] = CE.interval_volume(res.instance.V, res.instance.F)
    cert['accepted'] = bool(res.accepted and cert['bodytwin_geometry_certificate']['accepted'] is True
                            and cert['lv1_interval']['float_inside'])
    res.certificate = cert
    with tm('identity_chain'):
        idc = identity_chain(gm, res, out)
    weight = float(O['meta_h_w_age'][s, 1])
    exp, _ = export_roundtrips(res, out, weight, validate_full=validate_full, tm=tm)
    with tm('export_sdf_field_engine'):
        sdf = X.to_sdf([(f'{subject}_femur_r', res.instance.V, res.instance.F)], out / 'sdf')
    # C1 contract quantity for the hip centre
    with tm('c1_contract'):
        q = c1_quantities(res, subject)
    # error vs truth (subject's CT registration, LOO anchors)
    P = np.load(data / 'femur_pop.npz')
    names = [str(x) for x in P['names']]
    truth = P['S'][names.index(f'vsd:{subject}')]
    ev = {}
    if res.pose_identified:
        e = truth - res.instance.V
        Ci = np.linalg.inv(res.cov)
        z2 = np.einsum('ni,nij,nj->n', e, Ci, e)
        ev['E_placed_rms_mm'] = float(np.sqrt((e ** 2).sum(1).mean()))
        ev['coverage90_placed'] = float(np.mean(z2 <= C.CHI2_3_90))
    sR = C.umeyama(res.X_shape, truth, scale=False)
    Xa = C.apply(1, sR[1], sR[2], res.X_shape)
    ev['E_shape_rms_mm'] = float(np.sqrt(((truth - Xa) ** 2).sum(1).mean()))
    hc_t = C.femur_points(gm.anchors, truth)['HC6']
    hc_p, hc_c = res.points['HC6']
    ev['HJC_error_mm'] = float(np.linalg.norm(hc_t - hc_p)) if res.pose_identified else None
    ev['HJC_sigma_max_mm'] = float(np.sqrt(np.linalg.eigvalsh(hc_c)[-1]))
    ev['sigma_per_vertex_mm'] = dict(median=float(np.median(res.sigma_mm)), p95=float(np.percentile(res.sigma_mm, 95)))
    summary = dict(subject=subject, observations=dict(scenario=obs, rating_index=r0, features=feats,
                                                      landmarks={k: v.tolist() for k, v in lms.items()}),
                   population=dict(n=gm.model.n, members='Imperial 70 (35 persons); subject excluded'),
                   kappa_applied=res.kappa, accepted=cert['accepted'], certificate=cert, identity=idc, exports=exp,
                   sdf=sdf, c1=q, error_vs_truth=ev, timings=tm.t,
                   registry=dict(manifest_sha256=res.instance.reg.manifest_hash(), n_entities=len(res.instance.reg.entities)))
    (out / f'demo_{subject}.json').write_text(json.dumps(summary, indent=1, default=_js))
    return summary


def _js(o):
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, (np.floating, np.integer, np.bool_)):
        return o.item()
    return str(o)


def c1_quantities(res, subject):
    K = __import__('geomgr._local', fromlist=['contract']).contract()
    fid = f'{res.frame}'
    if fid not in K.FRAMES:
        K.register_frame(K.Frame(id=fid, body=subject, up=None, note='geomgr instance frame (CT of the subject when placed)'))
    p, Cc = res.points['HC6']
    post = res.post
    unc = dict(shape_posterior_and_incompleteness=float(np.sqrt(np.trace(Cc) / 3)),
               registration_floor_P1=0.46, pose=('in covariance' if res.pose_identified else K.UNKNOWN),
               landmark_definition='in covariance (VSD LOO second moment)', calibration_kappa=res.kappa)
    unc = {k: v for k, v in unc.items() if not isinstance(v, str) or v == K.UNKNOWN}
    q = K.Q(id=f'Point:{subject}/HJC@geomgr-0.1', value=np.asarray(p, float), units='mm', frame=fid, body=subject,
            status='calibrated_public_model', kind='point', uncertainty=unc,
            assumptions=('Gaussian conditioning in PCA space of 70 Imperial femurs', 'rigid pose, flat prior',
                         'VSD rater noise model'),
            derivation=('geomgr.core.condition', 'geomgr.core.landmark_predictive'),
            covariance={'matrix': np.asarray(Cc).tolist(), 'units': 'mm^2'},
            provenance_key=K.provenance_key({'subject': subject, 'obs': res.instance.meta}, 'geomgr-0.1/condition'))
    return q.to_json()
