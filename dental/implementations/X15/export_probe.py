"""Voxel-boundary mesh + center-sampled SDF of one actual tooth, laboratory probe."""
import json
import numpy as np
from scipy import ndimage as ndi
from scipy.spatial.transform import Rotation
from measure_r1 import HERE, DATA, dump, sha

def mesh(mask, offset, spacing):
    vertices = []
    faces = []
    lookup = {}
    for axis in range(3):
        for sign in (-1, 1):
            neighbor = ndi.shift(mask.astype(np.uint8), shift=[-sign if d == axis else 0 for d in range(3)], order=0, mode='constant', cval=0)
            cells = np.argwhere(mask & ~neighbor.astype(bool))
            others = [d for d in range(3) if d != axis]
            for cell in cells:
                keys = []
                for (u, v) in ((-1, -1), (1, -1), (1, 1), (-1, 1)):
                    key = 2 * (cell + offset)
                    key[axis] += sign
                    key[others[0]] += u
                    key[others[1]] += v
                    keys.append(tuple(key[::-1]))
                ids = []
                for key in keys:
                    if key not in lookup:
                        lookup[key] = len(vertices)
                        vertices.append(np.array(key) * spacing[::-1] / 2)
                    ids.append(lookup[key])
                xyz = np.array([vertices[k] for k in ids])
                normal = np.cross(xyz[1] - xyz[0], xyz[2] - xyz[0])
                target = np.zeros(3)
                target[2 - axis] = sign
                if normal @ target < 0:
                    ids = ids[::-1]
                faces.extend(((ids[0], ids[1], ids[2]), (ids[0], ids[2], ids[3])))
    return (np.array(vertices), np.array(faces, dtype=int))

def ply(path, v, f):
    with path.open('w') as out:
        out.write('ply\nformat ascii 1.0\ncomment units mm; published TF2 voxel-boundary laboratory probe\n')
        out.write(f'element vertex {len(v)}\nproperty float x\nproperty float y\nproperty float z\nelement face {len(f)}\nproperty list uchar int vertex_indices\nend_header\n')
        for x in v:
            out.write(' '.join((f'{q:.9f}' for q in x)) + '\n')
        for tri in f:
            out.write('3 ' + ' '.join((str(q) for q in tri)) + '\n')

def main():
    model = next((m for m in json.load(open(HERE / 'MOTION_MODELS.json')) if m['case'] == 'ToothFairy2F_001' and m['tooth'] == 44 and (m['source'] == 'label')))
    p = DATA / 'ToothFairy2F_001' / 'tooth_44_labels.npz'
    d = np.load(p)
    mask = d['label'] == 44
    sp = d['spacing_zyx']
    offset = d['offset_zyx']
    sdf = (ndi.distance_transform_edt(~mask, sampling=sp) - ndi.distance_transform_edt(mask, sampling=sp)).astype(np.float32)
    out = DATA / 'laboratory_probe'
    out.mkdir(exist_ok=True)
    np.savez_compressed(out / 'tooth_sdf.npz', sdf_mm=sdf, spacing_zyx_mm=sp, offset_zyx=offset)
    (v, f) = mesh(mask, offset, sp)
    a = np.array(model['axis_xyz'])
    b = np.array(model['buccal_xyz'])
    base = np.array(model['base_xyz_mm'])
    pivot = base + model['pivot_height_mm'] * a
    rot = Rotation.from_rotvec(np.deg2rad(5) * np.cross(a, b))
    moved = pivot + rot.apply(v - pivot) + 0.5 * b
    ply(out / 'tooth_original.ply', v, f)
    ply(out / 'tooth_moved_0p5mm_5deg.ply', moved, f)
    edges = np.sort(np.vstack((f[:, [0, 1]], f[:, [1, 2]], f[:, [2, 0]])), axis=1)
    (_, counts) = np.unique(edges, axis=0, return_counts=True)
    volume = float(np.sum(np.einsum('ij,ij->i', v[f[:, 0]], np.cross(v[f[:, 1]], v[f[:, 2]]))) / 6)
    files = [{'path': str(q), 'sha256': sha(q), 'bytes': q.stat().st_size} for q in sorted(out.iterdir()) if q.is_file()]
    dump(HERE / 'PROBE_EXPORT.json', {'case': model['case'], 'tooth': 44, 'units': 'mm', 'vertex_count': len(v), 'triangle_count': len(f), 'voxel_union_volume_mm3': float(mask.sum() * np.prod(sp)), 'mesh_signed_volume_mm3': volume, 'non_two_incident_edges': int(sum(counts != 2)), 'mesh_role': 'Actual whole-tooth voxel union; motion is rigid kinematic lab probe, no restoration or treatment design', 'sdf_role': 'Center-sampled signed EDT, negative inside; no anatomical/subvoxel accuracy claim', 'manufacturing_status': 'MESH_TO_REVIEW' if sum(counts != 2) else 'CLOSED_VOXEL_MESH; anatomy and material UNKNOWN', 'physical_validation': 'UNKNOWN; no manufactured part or physical measurement', 'files': files})
    print(len(v), 'vertices', len(f), 'triangles', sum(counts != 2), 'nonmanifold edges')
if __name__ == '__main__':
    main()
