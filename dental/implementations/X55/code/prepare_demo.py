"""Freeze local design labels and generate explicitly simulated manufactured scans."""
from dental_release.paths import expand as _release_expand
import json
import shutil
import sys
from pathlib import Path
import numpy as np
import trimesh
from scipy.ndimage import map_coordinates
from scipy.spatial.transform import Rotation
from metrology import sha, write_json, transform
R = Path(__file__).resolve().parents[1]
SOURCE = R.parent / _release_expand('X1B')
SINTER = R.parent / _release_expand('X14')

def run():
    lock = R / 'DEMO_INPUT_LOCK.json'
    if lock.exists():
        for (relative, expected) in json.loads(lock.read_text())['files'].items():
            if sha(R / relative) != expected:
                raise ValueError('Demo input hash mismatch: ' + relative)
        return json.loads((R / 'raw/GEOMETRY_MANIFEST.json').read_text())
    sys.path.insert(0, str(SINTER))
    from sinter_position import deform
    frozen = json.loads((R / 'FROZEN_PREDICTIONS.json').read_text())
    manifest = []
    for name in ['D1', 'M2']:
        dest = R / 'raw' / name
        dest.mkdir(parents=True, exist_ok=True)
        src = SOURCE / f'exports/{name}/crown.stl'
        shutil.copyfile(src, dest / 'design.stl')
        mesh = trimesh.load_mesh(dest / 'design.stl', process=False)
        a = np.load(SOURCE / f'inputs/geometry/{name}_model.npz')
        grid = json.loads((SOURCE / f'inputs/geometry/{name}_grid.json').read_text())
        p = mesh.triangles_center
        idx = (p - np.asarray(grid['origin'])) / grid['h']
        inside = map_coordinates(a['cavity'], idx.T, order=1)
        outside = map_coordinates(a['outer'], idx.T, order=1)
        z0 = grid['z_m'] + grid.get('margin_shift_mm', 0)
        labels = np.full(len(p), 'axial', dtype='<U8')
        labels[abs(inside) < abs(outside)] = 'intaglio'
        labels[(labels == 'axial') & (mesh.face_normals[:, 2] > 0.35) & (p[:, 2] > z0 + 2)] = 'occlusal'
        labels[p[:, 2] < z0 + 0.45] = 'marginal'
        assert set(labels) == {'marginal', 'intaglio', 'occlusal', 'axial'}
        write_json(dest / 'regions.json', dict(design_sha256=sha(dest / 'design.stl'), face_regions=labels.tolist(), resolution='PER_POINT', definition={'marginal': 'All surfaces within 0.45 mm above original CAD margin', 'intaglio': 'Inner surface by closest original cavity/outer SDF, excluding marginal band', 'occlusal': 'External surface nz>0.35 and >2 mm above margin', 'axial': 'Remaining external surface'}, label_uncertainty='CAD-defined simulation regions; lab must verify its own region labels', source_sdf_sha256=sha(SOURCE / f'inputs/geometry/{name}_model.npz')))
        merged = trimesh.load_mesh(dest / 'design.stl', process=True)
        V = merged.vertices.copy()
        N = merged.vertex_normals.copy()
        (_, _, f) = trimesh.proximity.closest_point(mesh, V)
        vlabels = labels[f]
        active = (vlabels == 'occlusal').astype(float)
        q = V.copy()
        q[:, 2] += active * frozen['synthetic_inputs']['occlusal_lift_um'] / 1000
        pose = np.eye(4)
        pose[:3, :3] = Rotation.from_euler('xyz', [2, -3, 4], degrees=True).as_matrix()
        pose[:3, 3] = [0.35, -0.2, 0.25]
        merged.vertices = transform(V, pose)
        merged.export(dest / 'rigid.stl')
        meta = dict(reference_kind='design', repeat_conditions='same_scanner_operator_protocol', simulation=True, scans={})
        meta['scans']['rigid.stl'] = dict(specimen_id='SYNTHETIC_' + name, initial_scan_to_reference=np.linalg.inv(pose).tolist(), scan_to_reference=np.linalg.inv(pose).tolist())
        for j in range(frozen['synthetic_inputs']['repeats']):
            rng = np.random.default_rng(5501 + j)
            noise = rng.normal(0, 0.008, len(V))[:, None] * N
            merged.vertices = transform(q + noise, pose)
            fn = f'repeat_{j}.stl'
            merged.export(dest / fn)
            meta['scans'][fn] = dict(specimen_id='SYNTHETIC_' + name, initial_scan_to_reference=np.linalg.inv(pose).tolist(), scan_to_reference=np.linalg.inv(pose).tolist())
        warped = deform(V, frozen['synthetic_inputs']['residual_sinter_stretch_xyz'])
        rounded = V + 0.035 * N * (vlabels == 'marginal')[:, None]
        for (tag, vertices) in [('occlusal', q), ('sinter', warped), ('rounding', rounded)]:
            merged.vertices = transform(vertices, pose)
            merged.export(dest / (tag + '.stl'))
            meta['scans'][tag + '.stl'] = dict(specimen_id='SYNTHETIC_' + name, initial_scan_to_reference=np.linalg.inv(pose).tolist(), scan_to_reference=np.linalg.inv(pose).tolist())
        write_json(dest / 'metadata.json', meta)
        manifest.append(dict(design=name, source=str(src), source_sha256=sha(src), region_counts={k: int((labels == k).sum()) for k in set(labels)}, source_sdf_grid_h_mm=grid['h'], sim_noise='8 um chosen vertex-normal sigma; interpolated scan-surface SD differs', sinter_operator_locator=str(SINTER / 'sinter_position.py'), sinter_operator_sha256=sha(SINTER / 'sinter_position.py')))
    write_json(R / 'raw/GEOMETRY_MANIFEST.json', manifest)
    return manifest
if __name__ == '__main__':
    run()
