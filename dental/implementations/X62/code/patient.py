"""Use published crown axes and ray/surface intersection; no learned tooth labels."""
from dental_release.paths import expand as _release_expand
import io, json, zipfile
import numpy as np
import trimesh
from common import ZIP, ROOT, save
UNN = [22, 23, 24, 25, 26, 27]
FDI = [33, 32, 31, 41, 42, 43]

def ray_hit(triangles, o, d):
    e1 = triangles[:, 1] - triangles[:, 0]
    e2 = triangles[:, 2] - triangles[:, 0]
    h = np.cross(np.broadcast_to(d, e2.shape), e2)
    a = np.einsum('ij,ij->i', e1, h)
    valid = np.abs(a) > 1e-12
    inv = np.zeros_like(a)
    inv[valid] = 1 / a[valid]
    s = o - triangles[:, 0]
    u = inv * np.einsum('ij,ij->i', s, h)
    q = np.cross(s, e1)
    v = inv * np.einsum('j,ij->i', d, q)
    t = inv * np.einsum('ij,ij->i', e2, q)
    good = valid & (u >= -1e-12) & (v >= -1e-12) & (u + v <= 1 + 1e-12) & (t > 1e-08)
    if not good.any():
        raise ValueError('no lingual crown ray hit; no silent substitute')
    ids = np.flatnonzero(good)
    i = int(ids[np.argmin(t[good])])
    return (o + t[i] * d, i, [float(1 - u[i] - v[i]), float(u[i]), float(v[i])])

def load():
    meshes = []
    rows = []
    members = []
    with zipfile.ZipFile(ZIP) as z:
        ap = _release_expand('@DENTAL_CASE_ID@/input/mandible/teeth_axes_mandible.json')
        axes = json.loads(z.read(ap))
        members.append(dict(member=ap, crc32=z.getinfo(ap).CRC, bytes=z.getinfo(ap).file_size))
        for (unn, fdi) in zip(UNN, FDI):
            name = _release_expand(f'@DENTAL_CASE_ID@/input/mandible/teeth/tooth_{unn}.stl')
            m = trimesh.load(io.BytesIO(z.read(name)), file_type='stl', process=False)
            a = axes[str(unn)]
            c = np.array(a['c'])
            x = np.array(a['x']) - c
            y = np.array(a['y']) - c
            long = np.array(a['z']) - c
            long /= np.linalg.norm(long)
            y -= long * np.dot(y, long)
            y /= np.linalg.norm(y)
            x = np.cross(y, long)
            d = -y
            (anchor, face, bary) = ray_hit(m.triangles, c, d)
            rows.append(dict(UNN=unn, FDI=fdi, crown_center_mm=c.tolist(), lingual_surface_anchor_mm=anchor.tolist(), long_axis=long.tolist(), labial_axis=y.tolist(), mesial_axis=x.tolist(), source_face=face, barycentric=bary, beam_neutral_axis_mm=(anchor + 4 / (3 * np.pi) * d).tolist(), crown_axis_radius_mm=float(np.max(np.linalg.norm(m.vertices - c, axis=1))), point_on_triangle_error_mm=float(np.linalg.norm(np.array(bary) @ m.triangles[face] - anchor))))
            meshes.append(m)
            members.append(dict(member=name, crc32=z.getinfo(name).CRC, bytes=z.getinfo(name).file_size))
    vertical = np.mean([r['long_axis'] for r in rows], axis=0)
    vertical /= np.linalg.norm(vertical)
    points = np.array([r['beam_neutral_axis_mm'] for r in rows])
    save('raw/PATIENT_GEOMETRY.json', dict(dataset='Open-Full-Jaw', patient=_release_expand('@DENTAL_CASE_ID@'), unit='mm', axes_convention='published UNN; crown long=z, labial=y, mesial=x; frame made orthonormal', resolution='PER_POINT', rows=rows, members=members, missing_candidates={_release_expand('@DENTAL_CASE_ID@'): 'UNN26 missing; 1/2 screened cases rejected'}, tooth_mesh_region='entire tooth; crown center used; vertex radius includes roots, conservative for displacement envelope', assumption='surface anchors are proposed bonding sites, not clinician annotations'))
    return (points, vertical, rows, meshes)

def export_section(points, vertical, filename):
    dense = np.concatenate([np.linspace(p, q, 21)[:-1] for (p, q) in zip(points[:-1], points[1:])] + [points[-1:]])
    theta = np.linspace(0, np.pi, 25)
    sec = np.array([0.75 * np.cos(theta), np.sin(theta) - 4 / (3 * np.pi)]).T
    vv = []
    for (i, p) in enumerate(dense):
        tangent = dense[min(i + 1, len(dense) - 1)] - dense[max(i - 1, 0)]
        tangent /= np.linalg.norm(tangent)
        y = vertical - tangent * np.dot(tangent, vertical)
        y /= np.linalg.norm(y)
        z = np.cross(tangent, y)
        vv.extend(p + sec[:, 0, None] * y + sec[:, 1, None] * z)
    n = len(sec)
    faces = []
    for i in range(len(dense) - 1):
        for j in range(n):
            a = i * n + j
            b = i * n + (j + 1) % n
            c = (i + 1) * n + j
            d = (i + 1) * n + (j + 1) % n
            faces.extend([[a, b, c], [b, d, c]])
    for j in range(1, n - 1):
        faces.extend([[0, j + 1, j], [(len(dense) - 1) * n, (len(dense) - 1) * n + j, (len(dense) - 1) * n + j + 1]])
    mesh = trimesh.Trimesh(np.array(vv), np.array(faces), process=False)
    target = ROOT / filename
    target.parent.mkdir(exist_ok=True)
    mesh.export(target)
    return dict(file=filename, watertight=bool(mesh.is_watertight), bytes=target.stat().st_size, clearance='UNKNOWN: no whole-path tooth clearance certified; straight interbond segments', geometry_only=True, clinical_device=False)
