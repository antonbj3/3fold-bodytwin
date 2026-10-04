"""Check both endpoints of nearest-surface witnesses; reject coordinate forgery."""
import argparse
import copy
import json
import sys
from pathlib import Path
import numpy as np
import trimesh

def distance(mesh, p):
    return float(trimesh.proximity.closest_point_naive(mesh, np.asarray(p)[None, :])[1][0])

def original_valid(w, mesh):
    (p, q, d) = (np.array(w['source_point_mm']), np.array(w['target_point_mm']), w['distance_mm'])
    return bool(np.isfinite(p).all() and np.isfinite(q).all() and np.isfinite(d) and (abs(np.linalg.norm(p - q) - abs(d)) <= 1e-07) and (abs(distance(mesh, p) - abs(d)) <= 1e-07))

def corrected_valid(w, target, source):
    try:
        return bool(original_valid(w, target) and distance(target, w['target_point_mm']) <= 1e-07 and (distance(source, w['source_point_mm']) <= 1e-07))
    except (ValueError, KeyError, TypeError):
        return False

def subset(mesh, ids):
    t = mesh.triangles[ids]
    return trimesh.Trimesh(t.reshape(-1, 3), np.arange(3 * len(t)).reshape(-1, 3), process=False)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--demo', type=Path, required=True)
    ap.add_argument('--data', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    a = ap.parse_args()
    sys.path.insert(0, str(a.demo / 'code'))
    from designgate.geometry import load
    report = json.loads((a.demo / 'reports/real_designgate.json').read_text())
    contract = json.loads((a.data / 'real_X18_110_34/contract_R2.json').read_text())
    crown = load(a.data / 'real_X18_110_34/crown.stl', 1)['mesh']
    target = subset(crown, contract['regions']['exterior'])
    source = subset(crown, contract['regions']['intaglio'])
    w = report['rules']['material_wall']['witness']
    (p, q) = (np.array(w['source_point_mm']), np.array(w['target_point_mm']))
    forged = copy.deepcopy(w)
    forged['target_point_mm'] = (2 * p - q).tolist()
    bad_distance = copy.deepcopy(w)
    bad_distance['distance_mm'] += 0.2
    bad_nan = copy.deepcopy(w)
    bad_nan['distance_mm'] = float('nan')
    measured = distance(target, p)
    result = {'wall_distance_mm_independent': measured, 'reported_distance_mm': w['distance_mm'], 'original_valid_witness_passes': original_valid(w, target), 'original_accepts_forged_target_endpoint': original_valid(forged, target), 'forged_endpoint_distance_to_target_mm': distance(target, forged['target_point_mm']), 'corrected_valid_witness_passes': corrected_valid(w, target, source), 'corrected_forged_endpoint_rejected': not corrected_valid(forged, target, source), 'corrected_bad_distance_rejected': not corrected_valid(bad_distance, target, source), 'corrected_nan_rejected': not corrected_valid(bad_nan, target, source), 'conclusion': 'Regional wall FAIL unchanged; complete commercial crown validation UNKNOWN.', 'resolution': 'PER_POINT'}
    ball = json.loads((a.demo / 'reports/stepped_cavity_closed.json').read_text())['rules']['ball_milling']
    ball_mesh = load(a.data / 'stepped_cavity_closed/crown.stl', 1)['mesh']
    bw = ball['witness']
    center = np.array(bw['ball_center_mm'])
    obstacle = np.array(bw['obstacle_point_mm'])

    def old_ball(x):
        return bool(np.isfinite(x['ball_center_mm']).all() and abs(np.linalg.norm(np.array(x['ball_center_mm']) - np.array(x['obstacle_point_mm'])) - x['obstacle_distance_to_center_mm']) <= 1e-07 and (abs(distance(ball_mesh, x['ball_center_mm']) - x['obstacle_distance_to_center_mm']) <= 1e-07) and (distance(ball_mesh, x['ball_center_mm']) < ball['radius_final_design_mm'] - 1e-05))

    def fixed_ball(x):
        return old_ball(x) and distance(ball_mesh, x['obstacle_point_mm']) <= 1e-07
    forged_ball = copy.deepcopy(bw)
    candidates = [2 * center - obstacle, center + np.roll(obstacle - center, 1), center + np.roll(obstacle - center, 2)]
    forged_ball['obstacle_point_mm'] = max(candidates, key=lambda x: distance(ball_mesh, x)).tolist()
    result.update(ball_original_accepts_forged_obstacle=old_ball(forged_ball), ball_corrected_positive_passes=fixed_ball(bw), ball_corrected_forged_obstacle_rejected=not fixed_ball(forged_ball), forged_obstacle_distance_to_mesh_mm=distance(ball_mesh, forged_ball['obstacle_point_mm']))
    a.output.write_text(json.dumps(result, indent=2) + '\n')
    assert abs(measured - w['distance_mm']) <= 1e-07
    assert result['corrected_valid_witness_passes']
    assert all((result[k] for k in ['corrected_forged_endpoint_rejected', 'corrected_bad_distance_rejected', 'corrected_nan_rejected']))
    assert result['ball_corrected_positive_passes'] and result['ball_corrected_forged_obstacle_rejected']
    print(json.dumps(result))
if __name__ == '__main__':
    main()
