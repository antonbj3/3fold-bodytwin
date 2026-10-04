"""Exact X12 ASCII voxel export functions; source remains INTE_REDO. No source outputs are written."""
import numpy as np
import hashlib

def ascii_stl(triangles, path):
    with path.open('w') as f:
        f.write('solid X12_mm\n')
        for t in triangles:
            n = np.cross(t[1] - t[0], t[2] - t[0])
            n /= np.linalg.norm(n)
            f.write(' facet normal ' + ' '.join((f'{v:.15g}' for v in n)) + '\n  outer loop\n')
            for p in t:
                f.write('   vertex ' + ' '.join((f'{v:.15g}' for v in p)) + '\n')
            f.write('  endloop\n endfacet\n')
        f.write('endsolid X12_mm\n')

def reread_ascii(path):
    vertices = []
    with path.open() as f:
        for line in f:
            q = line.split()
            if q and q[0] == 'vertex':
                vertices.append([float(v) for v in q[1:]])
    t = np.asarray(vertices, dtype=np.float64).reshape(-1, 3, 3)
    v = float(np.sum(np.einsum('ij,ij->i', t[:, 0], np.cross(t[:, 1], t[:, 2]))) / 6)
    return (t, v)

def voxel_stl(mask, spacing, path):
    padded = np.pad(mask, 1)
    all_tri = []
    for axis in range(3):
        for sign in [-1, 1]:
            sl = [slice(1, -1)] * 3
            sl[axis] = slice(2, None) if sign > 0 else slice(None, -2)
            surface = mask & ~padded[tuple(sl)]
            centers = np.argwhere(surface).astype(float)
            if not len(centers):
                continue
            u = (axis + 1) % 3
            v = (axis + 2) % 3
            corners = []
            for (a, b) in [(-0.5, -0.5), (0.5, -0.5), (0.5, 0.5), (-0.5, 0.5)]:
                q = centers.copy()
                q[:, axis] += 0.5 * sign
                q[:, u] += a
                q[:, v] += b
                corners.append(q * spacing)
            quad = np.stack(corners, axis=1)
            if sign < 0:
                quad = quad[:, ::-1, :]
            triangles = np.concatenate([quad[:, [0, 1, 2], :], quad[:, [0, 2, 3], :]])
            triangles = triangles[:, [0, 2, 1], ::-1]
            all_tri.append(triangles)
    triangles = np.concatenate(all_tri)
    ascii_stl(triangles, path)
    (readback, volume) = reread_ascii(path)
    import trimesh
    independent = trimesh.load(str(path), file_type='stl', process=False)
    independent_volume = float(independent.volume)
    expected = float(mask.sum() * np.prod(spacing))
    return {'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest(), 'bytes': path.stat().st_size, 'triangles': len(triangles), 'mesh_volume_mm3': volume, 'voxel_volume_mm3': expected, 'abs_volume_error_mm3': abs(volume - expected), 'max_vertex_readback_error_mm': float(np.max(np.abs(readback - triangles))), 'independent_trimesh_volume_mm3': independent_volume, 'independent_abs_volume_error_mm3': abs(independent_volume - expected), 'format': 'ASCII STL decimal15significantdigits', 'trimesh_version': trimesh.__version__}
