from common import *
import trimesh, argparse
p = argparse.ArgumentParser()
p.add_argument('uid')
p.add_argument('destination')
a = p.parse_args()
r = next((r for r in read(P / 'raw/PHYSICAL_ROWS.json') if r['uid'] == a.uid))
if r['status'] != 'ADAPTED':
    raise ValueError('No adapted geometry')
if sha(r['mesh_path']) != r['mesh_sha256']:
    raise ValueError('Frozen mesh drift')
rec = next((r for r in cohort() if r['uid'] == a.uid))
z = load_npz(r['mesh_path'])
pp = load_npz(rec['prep_path'])
dest = Path(a.destination)
if dest.exists():
    raise ValueError('Refuse overwrite of existing export')
trimesh.Trimesh(z['vertices'] @ pp['source_R'].T + pp['source_base'], z['faces'], process=False).export(dest)
print(json.dumps(dict(uid=a.uid, destination=str(dest), sha256=sha(dest), units='mm', geometry_accepted=r['geometry_pass'], manufactured_physical_validation='UNKNOWN')))
