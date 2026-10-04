import json, time
from pathlib import Path
import numpy as np
import trimesh
from identity import ROOT, DATA, sha, write
from geometry import global_mm
from x12_voxel_export import voxel_stl, ascii_stl, reread_ascii

def world_export(mask, spacing, frame, path):
    local = path.with_name(path.stem + '_local.stl')
    receipt = voxel_stl(mask, spacing, local)
    (tri, _) = reread_ascii(local)
    q = tri.reshape(-1, 3)[:, ::-1] / spacing
    world = global_mm(q, frame['shape_zyx'], frame['crop_origin_native_zyx'], frame['crown_axis_sign'], frame['source_header']).reshape(-1, 3, 3)
    M = np.array([float(x) for x in frame['source_header']['TransformMatrix'].split()]).reshape(3, 3)
    if frame['crown_axis_sign'] * np.linalg.det(M) < 0:
        world = world[:, [0, 2, 1]]
    ascii_stl(world, path)
    (back, volume) = reread_ascii(path)
    independent = trimesh.load(str(path), file_type='stl', process=False)
    expected = float(mask.sum() * np.prod(spacing))
    bounds = np.stack([world.min(axis=(0, 1)), world.max(axis=(0, 1))])
    err_bounds = float(np.max(np.abs(independent.bounds - bounds)))
    world_receipt = {'path': str(path), 'sha256': sha(path), 'bytes': path.stat().st_size, 'unit': 'mm', 'frame': 'native published TF2 world XYZ, source Offset/TransformMatrix/0.3mm scale, clinical registration UNKNOWN', 'expected_voxel_volume_mm3': expected, 'manual_signed_volume_mm3': volume, 'independent_trimesh_volume_mm3': float(independent.volume), 'manual_volume_error_mm3': abs(volume - expected), 'independent_volume_error_mm3': abs(float(independent.volume) - expected), 'bounds_error_mm': err_bounds, 'readback_vertex_error_mm': float(np.max(np.abs(back - world)))}
    for r in [receipt, world_receipt]:
        if 'abs_volume_error_mm3' in r:
            assert r['abs_volume_error_mm3'] <= 1e-06 and r['independent_abs_volume_error_mm3'] <= 1e-06 and (r['max_vertex_readback_error_mm'] <= 1e-09)
        else:
            assert r['manual_volume_error_mm3'] <= 1e-06 and r['independent_volume_error_mm3'] <= 1e-06 and (r['bounds_error_mm'] <= 1e-09) and (r['readback_vertex_error_mm'] <= 1e-09)
    sidecar = path.with_suffix('.frame.json')
    write(sidecar, frame | {'STL_coordinates': 'native TF2 WORLD_XYZ_mm', 'stl_sha256': sha(path), 'local_STL_sha256': sha(local)})
    return (receipt, world_receipt, world)

def main():
    started = time.monotonic()
    w = json.loads((ROOT / 'raw/R1_WITNESSES.json').read_text())
    witness = next((r for r in w if r['downstream_difference_mm'] > 0), w[0])
    fdi = witness['fdi']
    rows = json.loads((ROOT / 'raw/R1_ROWS.json').read_text())
    receipt = []
    first = None
    for op in ['REGIONAL_A_LOW_X', 'REGIONAL_B_HIGH_X']:
        row = next((r for r in rows if r['fdi'] == fdi and r['operation'] == op))
        a = np.load(row['spatial_export'])
        frame = json.loads(str(a['frame_json']))
        h = a['spacing']
        masks = [('prepared_' + op, a['prepared'])]
        if first is None:
            masks.extend([('original_tooth', a['tooth']), ('full_annotated_pulp', a['pulp'])])
        for (name, mask) in masks:
            (local, world, tri) = world_export(mask, h, frame, DATA / f'EXPORT_P48_FDI{fdi}_{name}.stl')
            receipt.extend([local, world])
            if first is None:
                first = (world, tri)
    (world, tri) = first
    controls = []
    for (name, wrong) in [('scale_times3', tri * 3), ('winding_reversed', tri[:, [0, 2, 1]])]:
        path = DATA / f'POISON_EXPORT_{name}.stl'
        ascii_stl(wrong, path)
        independent = trimesh.load(str(path), file_type='stl', process=False)
        rejected = abs(float(independent.volume) - world['expected_voxel_volume_mm3']) > 1e-06
        controls.append({'check': name, 'valid_pass': world['independent_volume_error_mm3'] <= 1e-06, 'serialized_poison_path': str(path), 'poison_volume_mm3': float(independent.volume), 'injected_rejected': rejected})
    assert all((r['valid_pass'] and r['injected_rejected'] for r in controls))
    write(ROOT / 'raw/EXPORT_RECEIPTS.json', receipt)
    write(ROOT / 'raw/EXPORT_CONTROLS.json', controls)
    write(ROOT / 'raw/EXPORT_OUTCOME.json', {'witness_fdi': fdi, 'files': len(receipt), 'all_pass': True, 'all_poison_rejected': True, 'world_export_volume_max_error_mm3': max((r['independent_volume_error_mm3'] for r in receipt if 'independent_volume_error_mm3' in r)), 'world_export_bounds_max_error_mm': max((r['bounds_error_mm'] for r in receipt if 'bounds_error_mm' in r)), 'wall_s': time.monotonic() - started})
    print('ASCII STL world/frame export PASS, actual scale and winding poisons rejected; FDI', fdi)
if __name__ == '__main__':
    main()
