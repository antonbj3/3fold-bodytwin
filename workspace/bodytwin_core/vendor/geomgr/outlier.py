"""Leave-one-landmark-out Mahalanobis testing for a GeometryManager."""
from __future__ import annotations

import argparse
import copy
import json
import time
from pathlib import Path

import numpy as np
from scipy import stats

from . import core as C
from .api import GeometryManager

D2_THRESHOLD = float(stats.chi2.ppf(0.999, 3))
Z_THRESHOLD = float(np.sqrt(D2_THRESHOLD))
L5 = ['HC6', 'SGT', 'LT', 'MEC', 'LEC']


def _as_point(value, name):
    value = np.asarray(value, float)
    if value.shape != (3,) or not np.all(np.isfinite(value)):
        raise ValueError(f'landmark {name}: expected finite 3D point')
    return value


def _names(gm, names):
    selected = list(gm.known_points if names is None else names)
    if len(selected) < 4:
        raise ValueError('at least four landmarks are required')
    if len(set(selected)) != len(selected):
        raise ValueError('landmark names must be unique')
    unknown = [n for n in selected if n not in gm.known_points]
    if unknown:
        raise ValueError(f'unknown landmarks: {unknown}')
    return selected


def _measurement_covariance(gm, post, target, rotation, translation):
    if target not in gm.C_lm:
        return np.zeros((3, 3))
    shape = gm.model.shape(post.b) @ np.asarray(rotation).T + np.asarray(translation)
    points = gm.point_fn(gm.anchors, shape)
    _, axes = gm.frame_fn(points)
    anatomical = np.asarray(gm.C_lm[target], float)
    return axes.T @ anatomical @ axes


def _fast_landmark_prediction(gm, post, target, placed):
    shape = gm.model.shape(post.b)
    points = gm.point_fn(gm.anchors, shape)
    mean = np.asarray(points[target], float)
    if target in gm.anchors.index:
        row = gm.anchors.M[gm.anchors.index[target]]
        jacobian = np.einsum('v,vkr->kr', row.toarray()[0], gm.model.G())
        s2 = gm.model.sig2_oos
        s2i = float(row.data @ s2[row.indices])
    else:
        h = 1e-3 * np.sqrt(np.maximum(gm.model.lam, 1e-12))
        jacobian = np.zeros((3, gm.model.r))
        for k in range(gm.model.r):
            delta = np.zeros(gm.model.r)
            delta[k] = h[k]
            plus = gm.point_fn(gm.anchors, gm.model.shape(post.b + delta))[target]
            minus = gm.point_fn(gm.anchors, gm.model.shape(post.b - delta))[target]
            jacobian[:, k] = (plus - minus) / (2.0 * h[k])
        s2 = gm.model.sig2_oos
        s2i = float(np.mean([s2[gm.anchors.M[gm.anchors.index[h_name]].indices].mean() for h_name in C.HEAD6]))
    inflation = 1.0 + 1.0 / gm.model.n
    if not placed:
        covariance = jacobian @ post.Sbb @ jacobian.T * inflation + s2i * np.eye(3)
        return mean, covariance
    rotation, translation = post.pose
    mean = rotation @ mean + translation
    jacobian = np.c_[rotation @ jacobian, -C._skew(rotation @ points[target]), np.eye(3)]
    covariance = jacobian @ post.Sigma @ jacobian.T * inflation + s2i * np.eye(3)
    return mean, covariance


def _observation_prediction(gm, post, target, conditioning, observations):
    if post.pose_identified:
        mean, covariance = _fast_landmark_prediction(gm, post, target, True)
        rotation, translation = post.pose
        return mean, covariance, True, rotation, translation
    shape = gm.model.shape(post.b)
    points = gm.point_fn(gm.anchors, shape)
    model_points = np.asarray([points[n] for n in conditioning], float)
    observed_points = np.asarray([observations[n] for n in conditioning], float)
    _, rotation, translation = C.umeyama(model_points, observed_points, scale=False)
    mean, covariance = _fast_landmark_prediction(gm, post, target, False)
    return rotation @ mean + translation, rotation @ covariance @ rotation.T, False, rotation, translation


def landmark_predictions(gm, landmarks, names=None):
    """Fit each leave-one-out posterior and return its held-out prediction."""
    if gm.bone != 'femur_r':
        raise ValueError('landmark leave-one-out testing currently requires femur_r')
    if gm.frame_fn is None or gm.C_lm is None:
        raise ValueError('landmark leave-one-out testing requires anatomical landmark noise')
    selected = _names(gm, names)
    observations = {n: _as_point(landmarks[n], n) for n in selected}
    predictions = []
    for target in selected:
        conditioning = [n for n in selected if n != target]
        post = C.condition(gm.model, gm.anchors, gm.point_fn,
                           lm_obs={n: observations[n] for n in conditioning},
                           C_lm={n: gm.C_lm[n] for n in conditioning}, frame_fn=gm.frame_fn,
                           max_iter=8, tol=1e-4)
        mean, covariance, aligned_pose, rotation, translation = _observation_prediction(
            gm, post, target, conditioning, observations)
        covariance = covariance + _measurement_covariance(gm, post, target, rotation, translation)
        covariance = 0.5 * (np.asarray(covariance, float) + np.asarray(covariance, float).T)
        eigenvalues = np.linalg.eigvalsh(covariance)
        if not np.all(np.isfinite(covariance)) or not np.all(np.isfinite(eigenvalues)) or eigenvalues[0] <= 0:
            raise ValueError(f'{target}: predictive covariance is not positive definite')
        stage_b = post.stage_b or {}
        predictions.append(dict(landmark=target, n_conditioning=len(conditioning), pose_identified=bool(post.pose_identified),
                                aligned_pose=bool(aligned_pose), condition_iters=int(stage_b.get('iters', 0)),
                                pose_min_eig=float(stage_b.get('pose_min_eig', np.nan)),
                                covariance_min_eig=float(eigenvalues[0]), mean=np.asarray(mean, float),
                                covariance=covariance))
    return predictions


def score_prediction(prediction, observation):
    """Return the held-out Mahalanobis score for one observation."""
    y = _as_point(observation, prediction['landmark'])
    residual = y - prediction['mean']
    covariance = prediction['covariance']
    try:
        d2 = float(residual @ np.linalg.solve(covariance, residual))
    except np.linalg.LinAlgError as ex:
        raise ValueError(f'{prediction["landmark"]}: covariance solve failed') from ex
    if not np.isfinite(d2):
        raise ValueError(f'{prediction["landmark"]}: non-finite Mahalanobis residual')
    if d2 < -1e-8:
        raise ValueError(f'{prediction["landmark"]}: negative Mahalanobis residual {d2}')
    d2 = max(0.0, d2)
    return dict(landmark=prediction['landmark'], n_conditioning=prediction['n_conditioning'],
                pose_identified=prediction['pose_identified'], aligned_pose=prediction['aligned_pose'],
                condition_iters=prediction['condition_iters'], pose_min_eig=prediction['pose_min_eig'],
                covariance_min_eig=prediction['covariance_min_eig'], residual_mm=float(np.linalg.norm(residual)),
                mahalanobis2=d2, standardized_residual=float(np.sqrt(d2)),
                threshold_d2=D2_THRESHOLD, flag=bool(d2 > D2_THRESHOLD))


def landmark_outlier_test(gm, landmarks, names=None):
    """Run the frozen leave-one-landmark-out test and return one row per target."""
    return [score_prediction(prediction, landmarks[prediction['landmark']])
            for prediction in landmark_predictions(gm, landmarks, names)]


def leave_one_landmark_out_test(gm, landmarks, names=None):
    return landmark_outlier_test(gm, landmarks, names)


def _jsonable(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.floating, np.integer, np.bool_)):
        return value.item()
    if isinstance(value, dict):
        return {str(k): _jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(v) for v in value]
    return value


def _summary(rows, protocol, displacement, landmark=None):
    chosen = [r for r in rows if r['protocol'] == protocol and r['displacement_mm'] == displacement
              and (landmark is None or r['landmark'] == landmark)]
    valid = [r for r in chosen if not r.get('error')]
    z = np.asarray([r['standardized_residual'] for r in valid], float)
    return dict(protocol=protocol, displacement_mm=displacement, n_cases=len(chosen), n_valid=len(valid),
                n_errors=len(chosen) - len(valid), n_subjects=len({r['subject'] for r in chosen}),
                n_targets=len({r['landmark'] for r in chosen}), n_flagged=int(sum(r['flag'] for r in valid)),
                flag_rate=float(np.mean([r['flag'] for r in valid])) if valid else None,
                standardized_residual_median=float(np.median(z)) if len(z) else None,
                standardized_residual_p95=float(np.percentile(z, 95)) if len(z) else None,
                standardized_residual_min=float(np.min(z)) if len(z) else None,
                standardized_residual_max=float(np.max(z)) if len(z) else None,
                pose_identified_rate=float(np.mean([r['pose_identified'] for r in valid])) if valid else None,
                source='rows with the stated protocol and displacement',
                derivation='counts, flag_rate, and standardized-residual quantiles over valid rows')


def _per_landmark(rows, protocol, displacement):
    chosen = [r for r in rows if r['protocol'] == protocol and r['displacement_mm'] == displacement]
    names = sorted({r['landmark'] for r in chosen})
    return {n: _summary(chosen, protocol, displacement, n) for n in names}


def _write_result(path, result):
    path = Path(path)
    path.write_text(json.dumps(_jsonable(result), indent=1, sort_keys=True))


def _subject_manager(base, population, observations, subject, subject_index):
    tri = population['loo_tri'][subject_index]
    bc = np.clip(population['loo_bc'][subject_index], 0, None)
    bc = bc / bc.sum(1, keepdims=True)
    anchors = C.Anchors(C.LM_NAMES, [base.F[t] for t in tri], list(bc), base.model.N)
    valid = np.isfinite(observations['obs_lm'][:, :, 0, 0])
    keep = [i for i in range(len(observations['subj'])) if i != subject_index]
    covariances = {}
    for j, name in enumerate(C.PT_NAMES):
        deviations = np.asarray([observations['dev_anat'][i, r, j] for i in keep
                                  for r in np.flatnonzero(valid[i])])
        definition = observations['ndef_all'][keep, j]
        covariances[name] = np.cov(deviations.T) + definition.T @ definition / len(definition)
    manager = copy.copy(base)
    manager.anchors = anchors
    manager.C_lm = covariances
    return manager


def run(data, out, subjects=None, displacements=(0, 5, 10, 20, 30), cache_dir=None):
    """Run the preregistered genuine and shifted-landmark comparison."""
    data = Path(data)
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    cache_dir = Path(cache_dir or out.parent)
    cache_dir.mkdir(parents=True, exist_ok=True)
    with np.load(data / 'vsd_obs.npz') as archive:
        p = {key: archive[key] for key in archive.files}
    with np.load(data / 'femur_pop.npz') as archive:
        population = {key: archive[key] for key in archive.files}
    subject_names = [str(x) for x in p['subj']]
    selected_subjects = subject_names if subjects is None else [str(x) for x in subjects]
    rows = []
    t0 = time.time()
    cpu0 = time.process_time()
    base_manager = GeometryManager(data, population=('imperial',), cache_dir=cache_dir)
    for sid in selected_subjects:
        s = subject_names.index(sid)
        valid = np.flatnonzero(np.isfinite(p['obs_lm'][s, :, 0, 0]))
        if not len(valid):
            continue
        rating = int(valid[0])
        base = {n: np.asarray(p['obs_lm'][s, rating, C.PT_NAMES.index(n)], float) for n in C.PT_NAMES}
        gm = _subject_manager(base_manager, population, p, sid, s)
        _, axes = C.knee_frame(base)
        for protocol, names in [('L5', L5), ('ALL23', C.PT_NAMES)]:
            try:
                predictions = landmark_predictions(gm, base, names)
            except Exception as ex:
                for name in names:
                    for displacement in displacements:
                        rows.append(dict(subject=sid, rating_index=rating, protocol=protocol, landmark=name,
                                         displacement_mm=int(displacement), target_excluded=True, error=str(ex), flag=False,
                                         source='frozen VSD first-valid rating; no usable prediction'))
                continue
            by_name = {p['landmark']: p for p in predictions}
            for name in names:
                prediction = by_name[name]
                for displacement in displacements:
                    y = base[name] + float(displacement) * axes[1]
                    try:
                        row = score_prediction(prediction, y)
                        row.update(subject=sid, rating_index=rating, protocol=protocol,
                                   displacement_mm=int(displacement), target_excluded=True,
                                   source='frozen VSD first-valid rating; shifted target only',
                                   derivation='D2=(y-mu)^T Sigma^-1(y-mu), flag D2>chi2_3(0.999)')
                        rows.append(row)
                    except Exception as ex:
                        rows.append(dict(subject=sid, rating_index=rating, protocol=protocol, landmark=name,
                                         displacement_mm=int(displacement), target_excluded=True, error=str(ex), flag=False,
                                         source='frozen VSD first-valid rating; shifted target only'))
    displacements = [int(x) for x in displacements]
    protocols = ['L5', 'ALL23']
    summary = {protocol: {str(d): _summary(rows, protocol, d) for d in displacements} for protocol in protocols}
    per_landmark = {protocol: {str(d): _per_landmark(rows, protocol, d) for d in displacements} for protocol in protocols}
    comparison = {}
    for d in displacements:
        a = summary['L5'][str(d)]
        b = summary['ALL23'][str(d)]
        comparison[str(d)] = dict(L5=a, ALL23=b,
                                  detection_rate_delta_ALL23_minus_L5=(None if a['flag_rate'] is None or b['flag_rate'] is None
                                                                      else b['flag_rate'] - a['flag_rate']))
    prereg_path = Path(__file__).resolve().parents[1] / 'PREREG.sha256'
    prereg_sha256 = prereg_path.read_text().split()[0] if prereg_path.exists() else None
    result = dict(id='BT-B35', threshold_d2=D2_THRESHOLD, threshold_standardized_residual=Z_THRESHOLD,
                  prereg_sha256=prereg_sha256,
                  protocols={p: (L5 if p == 'L5' else list(C.PT_NAMES)) for p in protocols},
                  subjects=selected_subjects, rating_rule='first valid rating per subject', displacements_mm=displacements,
                  summary=summary, per_landmark=per_landmark, comparison=comparison, rows=rows,
                  errors=[r for r in rows if r.get('error')], timings=dict(wall_s=time.time() - t0,
                                                                             cpu_s=time.process_time() - cpu0),
                  source_files=[str(data / 'femur_pop.npz'), str(data / 'vsd_obs.npz')],
                  model=dict(members=int(base_manager.model.n), shape_vertices=int(base_manager.model.N),
                             shape_modes=int(base_manager.model.r)),
                  conditioning='target excluded from a cold start; at most eight Gauss-Newton updates',
                  predictive_covariance='core landmark posterior + target anatomical measurement covariance',
                  implementation='geomgr/outlier.py; frozen D2 threshold from PREREG.md')
    _write_result(out, result)
    return result


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('--data', required=True)
    parser.add_argument('--out', required=True)
    parser.add_argument('--subjects', default='')
    parser.add_argument('--displacements', default='0,5,10,20,30')
    parser.add_argument('--cache-dir', default=None)
    args = parser.parse_args(argv)
    subjects = [x for x in args.subjects.split(',') if x] or None
    displacements = tuple(int(x) for x in args.displacements.split(',') if x)
    result = run(args.data, args.out, subjects, displacements, args.cache_dir)
    print(json.dumps({'n_rows': len(result['rows']), 'n_errors': len(result['errors']),
                      'summary': result['summary'], 'wall_s': result['timings']['wall_s']}, indent=1))
    return 0 if not result['errors'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
