from pathlib import Path
import sys, json
R = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(R))
from occlusion_module import load_demo_input
from occlusion_module.api import clean
row = json.loads((R / 'inputs/cohort.json').read_text())[4]
(b, c) = load_demo_input(R, row)
npz = Path(row['package_input']).name
for (name, key) in [('xy_mm', 'xy'), ('faces', 'grid_faces'), ('antagonist_z_mm', 'ceiling'), ('reference_gap_mm', 'reference_gap')]:
    b[name] = {'npz': npz, 'key': key}
for (name, key) in [('vertices_mm', 'vertices'), ('faces', 'faces'), ('face_roles', 'face_roles')]:
    c[name] = {'npz': npz, 'key': key}
cc = c['local_boundary_certificate']
pn = Path(row['local_constraints']['path']).name
for (name, key) in [('data', 'A_data'), ('indices', 'A_indices'), ('indptr', 'A_indptr'), ('shape', 'A_shape'), ('b', 'b')]:
    cc[name] = {'npz': pn, 'key': key}
(R / 'inputs/example_packet.json').write_text(json.dumps(clean(dict(bite=b, crown=c, request={'local_boundary': True})), indent=2, allow_nan=False) + '\n')
