"""Preserve R2; new closed fixture and robust harness, unchanged decision thresholds."""
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
    if (H / 'PREREG_R2B.json').exists():
        raise RuntimeError('Refuse overwrite')
    p = json.loads((H / 'PREREG_R2.json').read_text())
    p['round'] = 'R2B'
    p['frozen_utc'] = dt.datetime.now(dt.timezone.utc).isoformat()
    p['parent'] = 'PREREG_R2.json; rounds/R2.json; raw/R2_validation.log'
    p['changed_operation'] = 'Retain R2 distance/obstruction operators; close the intended synthetic annular material with the missing bottom ring in a NEW input. Harness handles absent witnesses as failed tests. No numerical gate changed.'
    dump(H / 'PREREG_R2B.json', p)
    (H / 'PREREG_R2B.sha256').write_text(sha(H / 'PREREG_R2B.json') + '\n')
    out = D / 'stepped_cavity_closed'
    out.mkdir(exist_ok=True)
    crown = trimesh.creation.revolve(np.array([[2, 0], [2, 3], [1.5, 3], [1.5, 1], [0.8, 1], [0.8, 0], [2, 0]], float), sections=32)
    export_stl(out / 'crown.stl', crown.vertices, crown.faces)
    for k in ['prep', 'antagonist']:
        (out / (k + '.stl')).write_bytes((D / 'stepped_cavity' / (k + '.stl')).read_bytes())
    (v, f) = parse_stl((out / 'crown.stl').read_bytes())
    m = trimesh.Trimesh(v, f, process=False)
    ct = m.triangles_center
    rad = np.linalg.norm(ct[:, :2], axis=1)
    ids = np.flatnonzero((rad < 1.6) & (np.sum(m.face_normals[:, :2] * ct[:, :2], axis=1) < -0.1))
    meta = json.loads((D / 'stepped_cavity/contract_R2.json').read_text())
    meta['regions']['intaglio'] = ids.tolist()
    meta['input_sha256'] = {k: sha(out / (k + '.stl')) for k in ['prep', 'crown', 'antagonist']}
    meta['physical_geometry_status'] += '; R2B bottom annular boundary explicitly closed'
    dump(out / 'contract_R2.json', meta)
    f = json.loads((H / 'FROZEN_PREDICTIONS_R2.json').read_text())
    f['round'] = 'R2B'
    f['frozen_utc'] = dt.datetime.now(dt.timezone.utc).isoformat()
    f['prereg_sha256'] = sha(H / 'PREREG_R2B.json')
    f['expected']['stepped_cavity_closed'] = f['expected'].pop('stepped_cavity')
    f['files'] += [dict(path=str(q), sha256=sha(q), bytes=q.stat().st_size) for q in sorted(out.glob('*'))]
    f['code'] = {str(q.relative_to(H)): sha(q) for q in (H / 'code').rglob('*.py')}
    f['revision_reason'] = 'Original open fixture and aborted validation preserved. New fixture completes same explicitly intended analytic closed ring volume. Thresholds unchanged.'
    dump(H / 'FROZEN_PREDICTIONS_R2B.json', f)
    (H / 'FROZEN_PREDICTIONS_R2B.sha256').write_text(sha(H / 'FROZEN_PREDICTIONS_R2B.json') + '\n')
    dump(H / 'CURRENT_WORK_STATE.json', dict(lane='X34-stl-design-gate', status='R2B_PREDICTIONS_FROZEN', latest_gate='R2 open fixture retained; new bottom ring frozen', next_operation='Run same 3D acceptance and independent witness controls', updated_utc=f['frozen_utc']))
    print('R2B frozen')
if __name__ == '__main__':
    main()
