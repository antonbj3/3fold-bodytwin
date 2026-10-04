"""Paired crown/die research geometry and a pre-measurement manifest; no physical predictions invented."""
import numpy as np
import trimesh
from .common import *
from .geometry import export_stl

def run():
    rows = read(ROOT / 'raw/SCORED_ROWS_R3.json')
    eligible = [r for r in rows if r.get('arm') == 'axial_flare' and r.get('wall', {}).get('status') == 'PASS']
    exports = []
    for r in eligible:
        a = np.load(r['file'], allow_pickle=False)
        v = a['preparation_vertices']
        f = a['preparation_faces']
        n = 24
        center = np.array([[0.0, 0.0, v[0, 2]]])
        V = np.r_[v, center]
        F = np.r_[f, np.array([[i, len(v), (i + 1) % n] for i in range(n)])]
        m = trimesh.Trimesh(V, F, process=False)
        if not m.is_watertight or not m.is_winding_consistent or m.volume <= 0:
            raise ValueError('closed die export failed')
        path = DATA / 'lab_pairs' / r['task_id']
        path.mkdir(parents=True, exist_ok=True)
        die = path / 'preparation_die.stl'
        export_stl(die, V, F)
        m2 = trimesh.load_mesh(die, process=True)
        if not m2.is_watertight:
            raise ValueError('die readback failed')
        crown = Path(r['file']).with_suffix('.stl')
        exports.append(dict(task_id=r['task_id'], units='mm', frame='same local crown/die frame', crown_path=str(crown), crown_sha256=sha(crown), die_path=str(die), die_sha256=sha(die), die_volume_relative_export_error=float(abs(m2.volume - m.volume) / m.volume), digital_prediction=dict(cap_wall_status='PASS', minimum_cap_wall_mm=0.5, nominal_film_min_mm=0.04, nominal_film_max_mm=0.12, antagonist_min_projected_gap_mm=r['occlusion']['measurement']['minimum_gap_mm']), physical_fracture_force_N=None, physical_seated_gap_um=None, physical_prediction_status='UNKNOWN_NO_SETUP_MATCHED_CALIBRATION'))
    freeze(ROOT / 'FROZEN_PREDICTIONS_FOR_LAB.json', dict(claim_type='capability', prospective_physical_measurement_performed=False, geometry_pairs=exports, scope='Frozen paired research geometry and digital predicates before any future manufacture/measurement. Physical force and actual fit predictions deliberately null.', minimum_measurement='Manufacture one declared material/process pair, record full batch/process/cementation metadata, independently scan both parts and measure seated regional and spatial gaps; report failure mode, not just an agreeing force.'))
    return exports
if __name__ == '__main__':
    print('paired research exports', len(run()))
