"""geomgr.api: the geometry manager as a Python API.

    from geomgr.api import GeometryManager
    gm = GeometryManager(data_dir)                          # femur_r, population = Imperial 70
    res = gm.instantiate(features={'L_mech': 431.0, ...}, feature_units={'L_mech': 'mm', ...},
                         landmarks={'MEC': [...], ...}, subject='z001')
    res.accepted, res.certificate, res.instance (identity.Instance), res.sigma_mm, res.points
    gm.sample(res, n) ; gm.batch(list_of_requests)

numpy/scipy only (cloud-safe). Exports live in geomgr.export (local-only validators).
"""
from __future__ import annotations

import hashlib
import json
import time
from pathlib import Path

import numpy as np

from . import core as C
from . import certify as CE
from . import identity as I


def _sha(a):
    return hashlib.sha256(np.ascontiguousarray(a).tobytes()).hexdigest()


class Result:
    def __init__(self, **kw):
        self.__dict__.update(kw)

    def summary(self):
        return dict(subject=self.subject, accepted=self.accepted, reasons=self.certificate.get('reasons'),
                    flags=self.certificate.get('flags'), frame=self.frame, kappa=self.kappa,
                    sigma_mm_median=None if self.sigma_mm is None else float(np.median(self.sigma_mm)),
                    pose_identified=self.pose_identified)


class GeometryManager:
    """A managed population in correspondence (one bone) with a versioned identity registry."""

    def __init__(self, data_dir, bone='femur_r', population=('imperial',), exclude_vsd=None, kappa=None,
                 cache_dir=None, extra_dir='auto'):
        """extra_dir: N7c frozen_extra (extended registry: tissue regions, wraps, TLEM constructs; B30 noise).
        'auto' = results/N7c/frozen_extra next to the package (or $GEOMGR_EXTRA); None = N7a registry only."""
        import os
        self.data_dir = Path(data_dir)
        if extra_dir == 'auto':
            extra_dir = os.environ.get('GEOMGR_EXTRA') or (Path(__file__).resolve().parents[1] / 'frozen_extra')
        self.extra_dir = Path(extra_dir) if extra_dir and (Path(extra_dir) / 'registry_femur_r_n7c.json').exists() else None
        self.bone = bone
        self.kappa = dict(kappa or {})          # scenario -> kappa (from LOSO calibration); default 1.0
        t0 = time.time()
        if bone == 'femur_r':
            P = np.load(self.data_dir / 'femur_pop.npz')
            O = np.load(self.data_dir / 'vsd_obs.npz')
            self.F = P['F']
            ds = np.array([str(x) for x in P['dataset']])
            groups = np.array([str(x) for x in P['groups']])
            sel = np.isin(ds, population)
            if exclude_vsd is not None:
                sel &= groups != exclude_vsd
            self.members = [str(x) for x, s in zip(P['names'], sel) if s]
            S = P['S'][sel]
            subj = [str(x) for x in O['subj']]
            s_ex = subj.index(exclude_vsd) if exclude_vsd in subj else None
            if s_ex is not None:        # landmark definitions without the excluded subject (LOO consensus)
                tri, bc = P['loo_tri'][s_ex], P['loo_bc'][s_ex]
            else:
                tri, bc = P['lm_tri'], P['lm_bc']
            bc = np.clip(bc, 0, None)
            bc = bc / bc.sum(1, keepdims=True)
            self.anchors = C.Anchors(C.LM_NAMES, [self.F[t] for t in tri], list(bc), S.shape[1])
            reg = json.loads(((self.extra_dir / 'registry_femur_r_n7c.json') if self.extra_dir is not None
                              else (self.data_dir / 'registry_femur_r.json')).read_text())
            self.registry = I.Registry.from_json(reg, faces=self.F)
            if s_ex is not None:
                for n, t, w in zip(C.LM_NAMES, tri, bc):
                    e = self.registry.entities[f'femur_r/landmark/{n}']
                    e.indices, e.weights = [int(i) for i in self.F[t]], [float(x) for x in w]
                    e.version += 1
                    e.note = f'LOO consensus without {exclude_vsd} (external demo)'
            self.point_fn = C.femur_points
            self.frame_fn = C.knee_frame
            self.feat_names = list(C.FEMUR_FEATURES)
            M = np.array([C.femur_features_from_landmarks(C.femur_points(self.anchors, V), V) for V in S])
            valid = np.isfinite(O['obs_lm'][:, :, 0, 0])
            keep = [s for s in range(len(subj)) if s != s_ex]
            pool = np.array([O['obs_feat'][s, r] - O['feat_mean20'][s] for s in keep for r in np.flatnonzero(valid[s])])
            self.R_pool = pool
            self.C_lm = {}
            for j, n in enumerate(C.PT_NAMES):
                D = np.array([O['dev_anat'][s, r, j] for s in keep for r in np.flatnonzero(valid[s])])
                Nd = O['ndef_all'][keep, j]
                self.C_lm[n] = np.cov(D.T) + Nd.T @ Nd / len(Nd)
            self.known_points = list(C.PT_NAMES)
        elif bone == 'tibia_r':
            from .evaluate import TibiaData, TIBIA_FEATS, TIBIA_NOISE_SD
            T = TibiaData(self.data_dir)
            self.F = T.F
            self.members = list(T.case)
            S = T.S
            self.anchors = T.anchors
            self.registry = I.Registry.from_json(json.loads((self.data_dir / 'registry_tibia_r.json').read_text()),
                                                 faces=self.F)
            self.point_fn = C.tibia_points
            self.frame_fn = None
            self.feat_names = list(TIBIA_FEATS)
            M = np.array([T.feats(i) for i in range(len(S))])
            self.R_pool = None
            self.noise_sd = TIBIA_NOISE_SD
            self.C_lm = None
            self.known_points = list(T.anchors.names)
        else:
            raise ValueError(bone)
        self.model = C.ShapeModel(S, M, self.feat_names, None)
        key = _sha(np.asarray(S, np.float64))[:16]
        cache = Path(cache_dir or self.data_dir) / f'oos_{bone}_{key}.npy'
        if cache.exists():
            self.model.set_oos(np.load(cache))
        else:
            grp = [m.split(':')[-1][:3] if m.startswith('imp:') else m for m in self.members]
            self.model.set_oos(C.oos_variance(S, np.array(grp)))
            try:
                np.save(cache, self.model.sig2_oos)
            except OSError:
                pass
        self.t_load_s = time.time() - t0
        # attachment anchors (for per-point covariance of attachments)
        att = [e for e in self.registry.entities.values() if e.kind == 'attachment']
        self.att_anchors = C.Anchors([e.name for e in att], [e.indices for e in att], [e.weights for e in att],
                                     self.model.N) if att else None

    # ------------------------------------------------------------------ requests
    def feature_R(self, names):
        if not names:
            return None
        if self.R_pool is not None:
            idx = [self.feat_names.index(k) for k in names]
            return np.cov(self.R_pool[:, idx].T).reshape(len(idx), len(idx))
        return np.diag([self.noise_sd[k] ** 2 for k in names])

    def instantiate(self, features=None, feature_units=None, landmarks=None, landmark_units='mm', subject='?',
                    scenario=None, with_points=True):
        t = {}
        t0 = time.process_time()
        features, landmarks = dict(features or {}), {k: np.asarray(v, float) for k, v in (landmarks or {}).items()}
        try:
            CE.check_inputs(self.model, features, feature_units, landmarks, landmark_units, self.known_points)
        except CE.InputError as ex:
            cert = dict(op='instantiate', accepted=False, reasons=[f'input contract: {ex}'], checks={}, flags=[])
            return Result(subject=subject, accepted=False, certificate=cert, instance=None, sigma_mm=None,
                          cov=None, points=None, frame=None, kappa=None, pose_identified=False, post=None, timings=t)
        post = C.condition(self.model, self.anchors, self.point_fn, features, self.feature_R(list(features)),
                           landmarks, self.C_lm, frame_fn=self.frame_fn)
        t['condition_cpu_s'] = time.process_time() - t0
        t1 = time.process_time()
        placed = post.stage_b is not None and post.pose_identified
        X_shape, _ = C.predictive(self.model, post)
        X, Cv = C.predictive(self.model, post, placed=placed)
        kappa = float(self.kappa.get(scenario, 1.0))
        Cv = Cv * kappa ** 2
        t['predictive_cpu_s'] = time.process_time() - t1
        t2 = time.process_time()
        cert = CE.numeric_certificate(self.model, post, X_shape, self.F, features)
        t['certificate_numeric_cpu_s'] = time.process_time() - t2
        frame = f'observation@{subject}' if placed else f'model@{self.bone}'
        inst = I.Instance(self.registry, X, self.F, units='mm', side='right', frame=frame, subject=subject,
                          meta=dict(features=features, landmarks={k: v.tolist() for k, v in landmarks.items()}))
        inst.log.append(dict(op='instantiate', accepted=cert['accepted'], scenario=scenario,
                             n_features=len(features), n_landmarks=len(landmarks)))
        # per-point sigma: sqrt of the largest eigenvalue of the 3x3 covariance (conservative radius)
        sig = np.sqrt(np.linalg.eigvalsh(Cv)[:, -1])
        pts = None
        if with_points:
            t3 = time.process_time()
            pts = {}
            lp = C.landmark_predictive(self.model, self.anchors, self.point_fn, post,
                                       self.known_points, placed=placed)
            for n, (p, Cl) in lp.items():
                pts[n] = (p, Cl * kappa ** 2)
            if self.att_anchors is not None:
                lpa = C.landmark_predictive(self.model, self.att_anchors, C.tibia_points, post,
                                            self.att_anchors.names, placed=placed)
                for n, (p, Cl) in lpa.items():
                    pts['att:' + n] = (p, Cl * kappa ** 2)
            t['points_cpu_s'] = time.process_time() - t3
        return Result(subject=subject, accepted=cert['accepted'], certificate=cert, instance=inst, sigma_mm=sig,
                      cov=Cv, points=pts, frame=frame, kappa=kappa, pose_identified=bool(placed), post=post,
                      X_shape=X_shape, timings=t)

    def instantiate_coefficients(self, b, subject='direct'):
        """Direct shape edit in PCA space (b in units of mm, model basis). Certified like any instance."""
        post = C.Posterior()
        post.b = np.asarray(b, float)
        post.Sbb = np.zeros((self.model.r, self.model.r))
        post.stage_a = dict(d2=0.0, dof=0, p=1.0, z={})
        post.stage_b = None
        post.pose_identified = False
        X = self.model.shape(post.b)
        cert = CE.numeric_certificate(self.model, post, X, self.F, {})
        inst = I.Instance(self.registry, X, self.F, frame=f'model@{self.bone}', subject=subject)
        return Result(subject=subject, accepted=cert['accepted'], certificate=cert, instance=inst, sigma_mm=None,
                      cov=None, points=None, frame=inst.frame, kappa=None, pose_identified=False, post=post,
                      X_shape=X, timings={})

    def sample(self, res, n=1, seed=0):
        """Posterior draws as Instances (same registry, same ids)."""
        rng = np.random.default_rng(seed)
        L = C._psd_sqrt(res.post.Sbb)
        out = []
        for _ in range(n):
            b = res.post.b + L @ rng.standard_normal(L.shape[1])
            X = self.model.shape(b)
            if res.pose_identified:
                R, t = res.post.pose
                X = X @ R.T + t
            inst = I.Instance(self.registry, X, self.F, frame=res.frame, subject=res.subject + ':draw')
            inst.log.append(dict(op='posterior_draw', seed=seed))
            out.append(inst)
        return out

    def batch(self, requests):
        return [self.instantiate(**rq) for rq in requests]

    def landmark_predictions(self, landmarks, names=None):
        from .outlier import landmark_predictions
        return landmark_predictions(self, landmarks, names)

    def landmark_outlier_test(self, landmarks, names=None):
        from .outlier import landmark_outlier_test
        return landmark_outlier_test(self, landmarks, names)

    def leave_one_landmark_out_test(self, landmarks, names=None):
        from .outlier import leave_one_landmark_out_test
        return leave_one_landmark_out_test(self, landmarks, names)
