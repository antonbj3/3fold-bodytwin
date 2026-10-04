"""Freeze successor expectations and geometry before distance/collision queries."""
import datetime as dt
import json
import sys
from pathlib import Path
import numpy as np
import trimesh
sys.path.insert(0, str(Path(__file__).resolve().parent))
from initialize import H, dump, sha
from prepare_demo import D
from designgate.vendor.geometry import export_stl, parse_stl

def main():
    if (H / 'FROZEN_PREDICTIONS_R2.json').exists():
        raise RuntimeError('Refuse refreeze')
    for case in ['valid', 'thin', 'film_low', 'film_high', 'penetration', 'real_X18_110_34']:
        m = json.loads((D / case / 'contract.json').read_text())
        m['preparation_outward_axis'] = [0, 0, 1]
        dump(D / case / 'contract_R2.json', m)
    for name in ['rotation', 'micrometers', 'wrong_scale']:
        m = json.loads((D / (name + '.json')).read_text())
        m['preparation_outward_axis'] = m['contact_axis']
        dump(D / (name + '_R2.json'), m)
    out = D / 'stepped_cavity'
    out.mkdir(exist_ok=True)
    crown = trimesh.creation.revolve(np.array([[2, 0], [2, 3], [1.5, 3], [1.5, 1], [0.8, 1], [0.8, 0]], float), sections=32)
    prep = trimesh.creation.revolve(np.array([[0, 0], [0.5, 0], [0.5, 1.3], [1.2, 1.3], [1.2, 2.1], [0, 2.1]], float), sections=32)
    paths = {k: out / (k + '.stl') for k in ('prep', 'crown', 'antagonist')}
    for (k, m) in [('prep', prep), ('crown', crown)]:
        export_stl(paths[k], m.vertices, m.faces)
    export_stl(paths['antagonist'], np.array([[-3, -3, 5], [3, -3, 5], [0, 3, 5]]), np.array([[0, 1, 2]]))
    (v, f) = parse_stl(paths['crown'].read_bytes())
    m = trimesh.Trimesh(v, f, process=False)
    ct = m.triangles_center
    rad = np.linalg.norm(ct[:, :2], axis=1)
    intag = np.flatnonzero((rad < 1.6) & (np.sum(m.face_normals[:, :2] * ct[:, :2], axis=1) < -0.1))
    meta = dict(units='mm', common_frame_confirmed=True, contact_axis=[0, 0, 1], input_sha256={k: sha(p) for (k, p) in paths.items()}, regions={'intaglio': intag.tolist()}, geometry_scope='closed_3D_material', crown_solid_integrity_asserted=True, prep_solid_integrity_asserted=True, extraction_axis=[0, 0, 1], extraction_travel_mm=2.0, ball_radius_final_mm=1.0, physical_geometry_status='Analytic revolved stepped cavity and mushroom prep; synthetic witness fixture', occlusal_gap_limits_mm=[-0.03, 0.1])
    dump(out / 'contract_R2.json', meta)
    for (name, gap) in [('proximal_valid', 0.02), ('proximal_penetration', -0.2)]:
        out = D / name
        out.mkdir(exist_ok=True)
        xy = np.array([[-1, -1], [1, -1], [-1, 1], [1, 1]], float)
        ff = np.array([[0, 1, 2], [1, 3, 2]])
        paths = {k: out / (k + '.stl') for k in ('prep', 'crown', 'antagonist', 'neighbor')}
        export_stl(paths['crown'], np.c_[np.zeros(4), xy], ff)
        export_stl(paths['neighbor'], np.c_[np.full(4, gap), xy], ff)
        export_stl(paths['prep'], np.c_[np.full(4, -1), xy], ff)
        export_stl(paths['antagonist'], np.c_[np.full(4, 1), xy], ff)
        dump(out / 'contract_R2.json', dict(units='mm', common_frame_confirmed=True, contact_axis=[1, 0, 0], neighbor_contact_axis=[1, 0, 0], input_sha256={k: sha(p) for (k, p) in paths.items()}, regions={'approximal': [0, 1], 'neighbor': [0, 1]}, approximal_gap_limits_mm=[-0.03, 0.1], physical_geometry_status='closed-form proximal patch fixture'))
    import importlib.metadata as im
    import trimesh.proximity as proximity
    files = []
    for p in sorted(D.rglob('*')):
        if p.is_file() and (p.suffix == '.stl' or p.name.endswith('_R2.json')):
            files.append(dict(path=str(p), sha256=sha(p), bytes=p.stat().st_size))
    freeze = dict(round='R2', claim_type='capability', frozen_utc=dt.datetime.now(dt.timezone.utc).isoformat(), prereg_sha256=sha(H / 'PREREG_R2.json'), expected=dict(valid={'material_wall': 'PASS', 'film_min': 'PASS', 'film_max': 'PASS'}, thin={'material_wall': 'FAIL'}, film_low={'film_min': 'FAIL'}, film_high={'film_max': 'FAIL'}, rotation={'material_wall': 'PASS', 'film_min': 'PASS', 'film_max': 'PASS'}, micrometers={'material_wall': 'PASS', 'film_min': 'PASS', 'film_max': 'PASS'}, wrong_scale={'film_max': 'FAIL'}, stepped_cavity={'ball_milling': 'FAIL', 'insertion': 'FAIL'}, proximal_valid={'approximal_contact': 'PASS'}, proximal_penetration={'approximal_contact': 'FAIL'}), value_tolerance_mm=0.0001, independent_naive_parity_mm=1e-07, real_prediction={'contact_gap_mm': 0.025, 'tolerance_mm': 1e-05, 'wall_decision': 'No foreknown sign claim; answer interval or UNKNOWN; actual closest witness must match independent full-face control', 'preparation': 'virtual plane, not external measured prep'}, files=files, code={str(p.relative_to(H)): sha(p) for p in (H / 'code').rglob('*.py')}, runtime=dict(python=sys.version, packages={x: im.version(x) for x in ['numpy', 'scipy', 'trimesh', 'rtree', 'shapely']}, proximity_source_sha256=sha(proximity.__file__)), physical_measurement_performed=False)
    dump(H / 'FROZEN_PREDICTIONS_R2.json', freeze)
    (H / 'FROZEN_PREDICTIONS_R2.sha256').write_text(sha(H / 'FROZEN_PREDICTIONS_R2.json') + '\n')
    dump(H / 'CURRENT_WORK_STATE.json', dict(lane='X34-stl-design-gate', status='R2_PREDICTIONS_FROZEN', latest_gate='R1 failed gate preserved; 3D successor expectations frozen', next_operation='Evaluate 3D distance envelopes and obstruction witnesses', updated_utc=freeze['frozen_utc']))
    print('Frozen R2:', len(files), 'inputs', sum((x['bytes'] for x in files)), 'bytes')
if __name__ == '__main__':
    main()
