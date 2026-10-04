"""Build physical dies from the same saved preparation SDF used by the surface gate."""
import io
import numpy as np
import trimesh
from skimage.measure import marching_cubes
from common import R, dump
from surface_gate import source_fields, evaluate_vertices
from matched_contrasts import three_mf

def run():
    results = []
    for design in ['M1', 'M2']:
        (_, field, g) = source_fields(design)
        for level in [0.0, 1e-06, -1e-06, 1e-05, -1e-05, 0.0001, -0.0001, 0.001, -0.001]:
            (V, F, _, _) = marching_cubes(field, level=level, spacing=(g['h'],) * 3, allow_degenerate=False)
            mesh = trimesh.Trimesh(V + g['origin'], F, process=True)
            mesh.fix_normals()
            mesh = trimesh.load_mesh(io.BytesIO(mesh.export(file_type='stl')), file_type='stl')
            result = evaluate_vertices(mesh.vertices, field, g)
            if mesh.is_watertight and mesh.volume > 0 and result['pass_surface']:
                break
        else:
            raise ValueError('No watertight source-consistent matched die; keep failed grid')
        for part in ['die', 'preparation']:
            mesh.export(R / f'exports/{design}/{part}.stl')
            three_mf(mesh, R / f'exports/{design}/{part}.3mf')
        results.append(dict(design=design, contour_level_mm=level, watertight=True, source_SDF_check=result))
    dump('raw/MATCHED_EXPORTS.json', results)
    return results
if __name__ == '__main__':
    print(run())
