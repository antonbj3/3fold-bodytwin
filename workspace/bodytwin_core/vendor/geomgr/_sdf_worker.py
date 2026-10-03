"""Field-engine SDF worker (runs in ~/projects/3fold_staging/.venv-bodytwin, CPU warp).

Same seam as 3fold-bodytwin/examples/anatomy/compose_mesh_to_field.py: strict mesh ingest (surface_mesh,
exact weld) -> faltkarna_v1_mesh_to_sdf.mesh_to_sdf_del(..., 'cpu', grind='flagga', metod='raypar_vindning',
returnera_falt=True); volume of the solid voxels vs the mesh volume and the raster bound of that script.
usage: python _sdf_worker.py input.npz outdir pitch staging_root
"""
import json
import os
import sys
import time
from pathlib import Path

import numpy as np

inp, outdir, pitch, root = sys.argv[1], Path(sys.argv[2]), float(sys.argv[3]), Path(sys.argv[4])
os.environ['CUDA_VISIBLE_DEVICES'] = ''
sys.path.insert(0, str(root / '3fold-bodytwin' / 'src'))
sys.path.insert(0, str(root / '3fold-field-engine' / 'src' / 'field_engine'))
from bodytwin.geometry.mesh_ingest_v1 import surface_mesh  # noqa: E402
import faltkarna_v1_mesh_to_sdf as field  # noqa: E402
import warp as wp  # noqa: E402

wp.config.quiet = True
wp.init()


def boundary_bound(v, f, origin, shape, P):
    """compose_mesh_to_field.boundary_bound (copied): volume of voxels touched by triangle bounding boxes."""
    band = np.zeros(shape, dtype=bool)
    for tri in v[f]:
        low = np.ceil((tri.min(0) - origin) / P - .5001).astype(int)
        high = np.floor((tri.max(0) - origin) / P + .5001).astype(int)
        low = np.maximum(low, 0)
        high = np.minimum(high, np.asarray(shape) - 1)
        if np.all(high >= low):
            band[tuple(slice(a, b + 1) for a, b in zip(low, high))] = True
    return int(band.sum()) * P ** 3


d = np.load(inp)
names = [str(x) for x in d['names']]
rows = []
for k, nm in enumerate(names):
    t0 = time.time()
    c0 = time.process_time()
    V, F = d[f'V_{k}'], d[f'F_{k}']
    mesh = surface_mesh(V, F, units='mm', weld_exact=True)
    v, f = mesh.field_arrays()
    lo = v.min(0) - 4 * pitch
    r = field.mesh_to_sdf_del(wp, v, f, pitch, 0., lo, 'cpu', grind='flagga', metod='raypar_vindning', returnera_falt=True)
    P = float(r.get('pitch_effektiv', pitch))
    origin = lo + np.asarray(r['gmin']) * P
    solid, sd = r['solid_final'], r['sd']
    vol = float(solid.sum()) * P ** 3
    bound = boundary_bound(v, f, origin, solid.shape, P)
    gate = r['vattentathet']
    np.savez_compressed(outdir / f'{nm}_sdf.npz', sdf_mm=sd.astype(np.float32), solid=solid, origin_mm=origin,
                        pitch_mm=np.array(P))
    rows.append(dict(name=nm, shape=list(solid.shape), pitch_mm=P, mesh_volume_mm3=float(mesh.volume_mm3),
                     voxel_volume_mm3=vol, volume_relative_error=abs(vol - mesh.volume_mm3) / mesh.volume_mm3,
                     volume_abs_error_mm3=abs(vol - mesh.volume_mm3), raster_bound_mm3=bound,
                     within_bound=abs(vol - mesh.volume_mm3) <= bound, gate_ok=gate['status'] == 'OK',
                     closed=bool(gate['vattentat']), finite=bool(np.isfinite(sd).all()),
                     wall_s=time.time() - t0, cpu_s=time.process_time() - c0))
rep = dict(backend='CPU (warp)', n=len(rows), rows=rows,
           ok=all(x['gate_ok'] and x['closed'] and x['finite'] and x['volume_relative_error'] <= 0.05 and x['within_bound']
                  for x in rows))
(outdir / 'sdf_report.json').write_text(json.dumps(rep, indent=1))
print(json.dumps(dict(ok=rep['ok'], n=len(rows))))
