"""Recompute the surface claim from delivered vertices and the source SDF.
No reported surface-error value or path existence can substitute for this test.
"""
import sys, zipfile, xml.etree.ElementTree as ET
import numpy as np
from scipy.ndimage import map_coordinates
import trimesh
from functools import lru_cache
from common import R, read
sys.path.insert(0, str(R / 'vendor'))
from crown_fit_geometry import Grid, redistance_lowmem

def load_3mf(path):
    with zipfile.ZipFile(path) as z:
        r = ET.fromstring(z.read('3D/3dmodel.model'))
    if r.attrib.get('unit') != 'millimeter':
        raise ValueError('3MF unit must be millimeter')
    ns = {'c': 'http://schemas.microsoft.com/3dmanufacturing/core/2015/02'}
    v = np.array([[float(e.attrib[k]) for k in ['x', 'y', 'z']] for e in r.findall('.//c:vertex', ns)])
    f = np.array([[int(e.attrib[k]) for k in ['v1', 'v2', 'v3']] for e in r.findall('.//c:triangle', ns)])
    return trimesh.Trimesh(v, f, process=False)

def evaluate_vertices(vertices, field, grid):
    v = np.asarray(vertices, float)
    q = (v - np.array(grid['origin'])) / grid['h']
    in_grid = np.isfinite(q).all() and bool(((q >= 0) & (q <= np.array(field.shape) - 1)).all())
    values = map_coordinates(field, q.T, order=1, mode='constant', cval=np.inf, prefilter=False)
    g = Grid(grid['origin'], np.array(grid['origin']) + grid['h'] * (np.array(grid['shape']) - 1), grid['h'])
    control = g.sample(field, v, cval=np.inf)
    p99 = float(np.percentile(abs(values), 99) * 1000) if np.isfinite(values).all() else None
    err = float(np.max(abs(values - control))) if np.isfinite(values).all() and np.isfinite(control).all() else None
    return dict(vertices=len(v), within_grid=in_grid, sdf_abs_p99_um=p99, control_difference_mm_max=err, pass_surface=bool(in_grid and p99 is not None and (p99 <= 60) and (err is not None) and (err <= 1e-06)))

@lru_cache(maxsize=5)
def source_fields(design):
    a = np.load(R / f'inputs/geometry/{design}_model.npz')
    g = read(f'inputs/geometry/{design}_grid.json')
    if 'margin_shift_mm' in g:
        z_m = g['z_m'] + g['margin_shift_mm']
    else:
        old = next((d for d in read('history/X1/FROZEN_PREDICTIONS.json')['designs'] if d['design'] == design))
        z_m = g['z_m'] + old['geometry_parameters']['margin_shift']
    z = g['origin'][2] + g['h'] * np.arange(g['shape'][2])
    grid = Grid(g['origin'], np.array(g['origin']) + g['h'] * (np.array(g['shape']) - 1), g['h'])
    physical_die = redistance_lowmem(np.maximum(a['prep'], z_m - 3 - z[None, None, :]).astype(np.float32), grid)
    return (a['crown'], physical_die, g)

def check_export(design, part, mesh=None):
    (crown, die, g) = source_fields(design)
    field = die if part in ('die', 'preparation') else crown
    own = mesh is None
    if own:
        mesh = trimesh.load_mesh(R / f'exports/{design}/{part}.stl', process=True)
    result = evaluate_vertices(mesh.vertices, field, g)
    result.update(design=design, part=part, watertight=bool(mesh.is_watertight), finite=bool(np.isfinite(mesh.vertices).all()), volume_mm3=float(mesh.volume))
    if own:
        other = load_3mf(R / f'exports/{design}/{part}.3mf')
        r = evaluate_vertices(other.vertices, field, g)
        dv = abs(mesh.volume / other.volume - 1) if other.volume > 0 else float('inf')
        result.update(three_mf_surface=r, roundtrip_volume_relative_error=dv)
        result['pass'] = bool(result['pass_surface'] and result['watertight'] and (mesh.volume > 0) and r['pass_surface'] and other.is_watertight and (dv <= 1e-05))
    else:
        result['pass'] = bool(result['pass_surface'] and result['watertight'] and (mesh.volume > 0))
    return result

def main():
    return [check_export(d, p) for d in ['D1', 'D2', 'D3', 'M1', 'M2'] for p in ['crown', 'die', 'preparation']]
