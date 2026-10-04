from local import *
(a, b) = first()
m = npz(a['mesh_path'])
(pv, pf) = compact(m['vertices'], m['faces'][m['roles'] == 1][:, ::-1])
top = trimesh.Trimesh(pv, pf, process=False)
bd = loops(top)[0]
rim = pv[bd]
n = np.cross(pv[pf][:, 1] - pv[pf][:, 0], pv[pf][:, 2] - pv[pf][:, 0])
poly = Polygon(rim[:, :2])
center = npz(b['public_path'])['cavity_apex']
print('XY valid', poly.is_valid, 'nz negative', np.sum(n[:, 2] < 0), 'nz positive', np.sum(n[:, 2] > 0), 'z extents', pv[:, 2].min(), pv[:, 2].max(), 'apex', center, 'min normals', np.min(np.linalg.norm(n, axis=1)))
print('igl', igl.signed_distance.__doc__[:600])
print('modules')
import importlib.util
for s in ['pymeshfix', 'pymeshlab', 'manifold3d', 'triangle', 'vtk', 'pyvista']:
    print(s, importlib.util.find_spec(s))
