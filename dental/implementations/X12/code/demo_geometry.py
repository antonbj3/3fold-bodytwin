from dental_release.paths import expand as _release_expand
import hashlib, json, time, os
from pathlib import Path
import numpy as np
from scipy import ndimage as ndi
ROOT = Path(os.environ.get('X12_RUN_ROOT', str(Path(__file__).resolve().parents[1])))
DATA = Path(_release_expand('@DENTAL_WORK_ROOT@/X12'))

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

def main():
    t0 = time.time()
    p = DATA / 'R8_P1_36_paired.npz'
    if not p.exists():
        p = DATA / 'R6_P1_36_paired.npz'
    q = np.load(p)
    t = q['tooth']
    pulp = q['pulp']
    sp = q['spacing']
    depth = 0.5
    phi = ndi.distance_transform_edt(~t, sampling=sp) - ndi.distance_transform_edt(t, sampling=sp)
    preparation = t & (ndi.distance_transform_edt(t, sampling=sp) >= depth + 0.5 * float(sp.min()))
    files = [voxel_stl(t, sp, DATA / 'DEMO_tooth_voxel.stl'), voxel_stl(pulp, sp, DATA / 'DEMO_pulp_voxel.stl'), voxel_stl(preparation, sp, DATA / 'DEMO_uniform_offset_0p5mm.stl')]
    sdf = DATA / 'DEMO_sdf_and_offset.npz'
    np.savez_compressed(sdf, tooth=t, pulp=pulp, phi_tooth_mm=phi.astype(np.float32), preparation=preparation, spacing_mm=sp, depth_mm=depth)
    files.append({'path': str(sdf), 'sha256': hashlib.sha256(sdf.read_bytes()).hexdigest(), 'bytes': sdf.stat().st_size})
    import trimesh
    original = trimesh.load(files[0]['path'], file_type='stl', process=False)
    source_triangles = np.asarray(original.triangles)
    poison = []
    for (name, tri) in [('scale_x3', source_triangles * 3), ('reversed_winding', source_triangles[:, [0, 2, 1], :])]:
        path = DATA / f'POISON_{name}.stl'
        ascii_stl(tri, path)
        (_, v) = reread_ascii(path)
        vm = float(trimesh.load(str(path), file_type='stl', process=False).volume)
        poison.append({'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest(), 'injected_error': name, 'readback_volume_mm3': v, 'independent_volume_mm3': vm, 'injected_rejected': abs(vm - files[0]['voxel_volume_mm3']) > 1e-06})
    boundary = preparation & ~ndi.binary_erosion(preparation, structure=ndi.generate_binary_structure(3, 1))
    dist = ndi.distance_transform_edt(~boundary, sampling=sp)
    ztop = int(np.argwhere(pulp)[:, 0].max())
    horn = pulp.copy()
    horn[:max(0, ztop - round(2 / sp[0]))] = False
    remaining = float(dist[horn].min())
    removed = float((t.sum() - preparation.sum()) * np.prod(sp))
    result = {'files': files, 'removed_mm3': removed, 'demo_uniform_depth_mm': depth, 'remaining_horn_to_prepared_boundary_centres_mm': remaining, 'two_surface_voxel_envelope_mm': float(np.linalg.norm(sp)), 'clinical_status': 'UNKNOWN', 'scope': 'Voxel mesh and approximate uniform erosion illustration; no cervical margin, regional clinical preparation, cement or fabrication error. Total hard tissue; enamel/dentin unresolved.', 'volume_gate_pass': all((r.get('abs_volume_error_mm3', 0) <= 1e-06 and r.get('independent_abs_volume_error_mm3', 0) <= 1e-06 and (r.get('max_vertex_readback_error_mm', 0) <= 1e-09) for r in files)), 'poison': poison, 'poison_volume_gate_rejects': all((p['injected_rejected'] for p in poison)), 'external_referent': {'kind': 'published_code', 'locator': 'https://trimesh.org/trimesh.html#trimesh.Trimesh.volume', 'compared_quantity': 'saved mesh enclosed volume vs external published voxel-mask volume', 'refutes_us': False}, 'wall_s': time.time() - t0}
    (ROOT / 'raw/DEMO_GEOMETRY.json').write_text(json.dumps(result, indent=2) + '\n')
    assert result['volume_gate_pass'] and result['poison_volume_gate_rejects']
    print(json.dumps(result))
if __name__ == '__main__':
    main()
