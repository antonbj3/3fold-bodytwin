"""Lock local meshes and public parameter cards. Never modify parent projects."""
import sys, collections
from common import *
import trimesh

def prepare():
    DATA.mkdir(parents=True, exist_ok=True)
    fnd = DENTAL / 'results/PROOF_LANE_GENCAD_FOUNDATION'
    v2 = DENTAL / 'results/PROOF_LANE_GENCAD_V2'
    x1 = DENTAL / 'results/LANE_X1B_CROWN_LOOP'
    sys.path.insert(0, str(fnd))
    from gencad_bench.geometry import shell, cap
    rows = []
    excluded = []
    for name in ['D1', 'D2', 'D3', 'M1', 'M2']:
        p = x1 / f'exports/{name}/crown.stl'
        g = read(x1 / f'inputs/geometry/{name}_grid.json')
        mesh = trimesh.load_mesh(p, process=True)
        z = g['z_m'] + g.get('margin_shift_mm', 0.0)
        if 'margin_shift_mm' not in g:
            old = next((d for d in read(x1 / 'history/X1/FROZEN_PREDICTIONS.json')['designs'] if d['design'] == name))
            z += old['geometry_parameters']['margin_shift']
        a = np.load(x1 / f'inputs/geometry/{name}_model.npz', allow_pickle=False)
        from scipy.ndimage import map_coordinates
        ct = mesh.triangles_center
        values = map_coordinates(a['cavity'], ((ct - np.asarray(g['origin'])) / g['h']).T, order=1, mode='constant', cval=100, prefilter=False)
        inner = (abs(values) < 2 * g['h']) & (ct[:, 2] > z + 0.06)
        inner &= np.einsum('ij,ij->i', mesh.face_normals[:, :2], ct[:, :2]) < 0
        inner |= (abs(values) < g['h']) & (ct[:, 2] > z + 0.06) & (mesh.face_normals[:, 2] < -0.5)
        dst = DATA / f'crownloop_{name}.npz'
        np.savez_compressed(dst, vertices=mesh.vertices, faces=mesh.faces, inner_faces=np.flatnonzero(inner))
        rows.append(dict(id=f'crownloop_{name}', family='crown_loop', file=str(dst), sha256=sha(dst), source=str(p), source_sha256=sha(p), die=str(x1 / f'exports/{name}/preparation.stl'), die_sha256=sha(x1 / f'exports/{name}/preparation.stl'), dataset='STS-Tooth3D', license='CC BY 4.0', locator='https://doi.org/10.5281/zenodo.10597292', geometry_error_mm=g['h'], scale='sintered_design_mm', intaglio_classifier='source cavity SDF plus inward radial/downward normal; marginal strip excluded'))
    tasks = {t['geometry_id']: t for t in read(fnd / 'data/tasks.json')}
    for row in read(fnd / 'rounds/R3/results.json')['rows']:
        gid = row['geometry_id']
        method = row['method']
        t = tasks[gid]
        p = row['profile']
        outer = t['public'].get('outer_template')
        if outer is None:
            path = sorted((fnd / 'raw/runs').glob('*'))[-1] / (t['task_id'] + '__parametric_IFU_design.json')
            outer = read(path)['outer']
        mesh = shell(outer, p)
        mesh.fix_normals()
        (iv, inf) = cap(p)
        (ov, of) = cap(outer)
        ids = np.arange(len(of), len(of) + len(inf))
        (pv, pf) = cap(t['preparation']['profile'])
        dst = DATA / (gid + '_' + method + '.npz')
        np.savez_compressed(dst, vertices=mesh.vertices, faces=mesh.faces, inner_faces=ids, prep_vertices=pv, prep_faces=pf)
        rows.append(dict(id=gid + '_' + method, family='GenCAD_Foundation', file=str(dst), sha256=sha(dst), source=str(fnd / 'rounds/R3/results.json'), source_sha256=sha(fnd / 'rounds/R3/results.json'), dataset='Teeth3DS-derived synthetic preparation', license='Teeth3DS terms at source; CC BY-SA reported by parent', locator='https://osf.io/xctdy/', geometry_error_mm=1e-06, scale='design_mm', parent_milling_status=row['milling']['status']))
    for row in read(v2 / 'raw/PREDICTIONS_R3.json'):
        if row['status'] != 'DESIGNED':
            excluded.append(dict(id=row['task_id'], reason=row['status']))
            continue
        if row['arm'] != 'axial_flare':
            continue
        path = Path(row['file'])
        if sha(path) != row['sha256']:
            raise ValueError('V2 parent hash drift')
        a = np.load(path, allow_pickle=False)
        dst = DATA / (row['task_id'] + '.npz')
        ids = np.arange(len(a['outer_faces']), len(a['outer_faces']) + len(a['inner_faces']))
        np.savez_compressed(dst, vertices=a['vertices'], faces=a['faces'], inner_faces=ids, prep_vertices=a['preparation_vertices'], prep_faces=a['preparation_faces'])
        rows.append(dict(id=row['task_id'], family='GenCAD_V2', file=str(dst), sha256=sha(dst), source=str(path), source_sha256=sha(path), dataset='Bits2Bites/Bite2Text-derived full cap', license='CC BY-NC-SA; see parent license lock', locator='https://doi.org/10.1038/s41597-025-05580-6', geometry_error_mm=1e-06, scale='design_mm', parent_milling_status='UNKNOWN_SHAFT_NOT_CERTIFIED'))
    tools = []
    for (dia, reach) in [(0.6, 3.0), (1.0, 16.0), (2.0, 16.0)]:
        tools.append(dict(id=f'VHF_Z{int(dia * 100):03d}', library='vhf', diameter_mm=dia, shank_mm=3.0, neck_reach_mm=reach, total_length_mm=40.0, gauge_mm=22.0, holder_diameter_mm=10.0, holder_length_mm=12.0, source='sources/VHF_tools.pdf p19', locator=read(ROOT / 'sources/URLS.json')['VHF_tools.pdf'], published_fields=['diameter_mm', 'shank_mm', 'total_length_mm', 'neck_reach_mm'], unknown_fields=['neck_radius_profile', 'installed_gauge', 'holder_shape'], model_status='CONDITIONAL_SIMPLIFIED_NECK_AND_HOLDER'))
    for (dia, reach) in [(0.6, 8.0), (1.0, 12.0), (2.5, 15.0)]:
        tools.append(dict(id=f'AG_{dia:g}', library='ceramill', diameter_mm=dia, shank_mm=3.0, neck_reach_mm=reach, total_length_mm=47.0, gauge_mm=22.0, holder_diameter_mm=10.0, holder_length_mm=12.0, source='Ha & Cho 2016 methods; only cutter diameters published in this paper', locator='https://doi.org/10.4047/jap.2016.8.6.439', published_fields=['diameter_mm'], unknown_fields=['shank_mm', 'neck_reach_mm', 'total_length_mm', 'neck_radius_profile', 'installed_gauge', 'holder_shape'], model_status='CONDITIONAL_SCENARIO_ALL_NONDIAMETER_GEOMETRY'))
    dump(ROOT / 'TOOL_LIBRARY.json', dict(tools=tools, profiles=[dict(id='vhf_4axis', axes=4, A_deg=[0, 360], B_deg=[0, 0], locator='https://www.vhf.com/wp-content/uploads/HB1226-ES_K4edition-Folder.pdf', status='published axes; local x mounting inferred, not calibrated'), dict(id='vhf_5axis', axes=5, A_deg=[0, 360], B_deg=[-35, 35], locator='https://www.vhf.com/en-us/products/dental/dental-milling-machines/k5/', status='published rotary limits; local x/y mapping inferred'), dict(id='ceramill_5axis_scenario', axes=5, A_deg=[0, 360], B_deg=[-35, 35], locator='https://www.amanngirrbach.com/en-us/equipment/production-cam/ceramill-motion-2/', status='5 axes published; +/-35 is borrowed scenario, actual Motion2 tilt UNKNOWN')], resolution='PER_POINT', units='mm and degrees'))
    dump(ROOT / 'INPUT_MANIFEST.json', dict(rows=rows, excluded=excluded, requested=5 + 24 + 32, unique_accepted=len(rows), duplicate_unflared_excluded=31, dropout_fraction=len(excluded) / (5 + 24 + 32), no_CAM_logs_found_in_scoped_inputs=True))
    return rows
if __name__ == '__main__':
    if not (ROOT / 'INPUT_MANIFEST.json').exists():
        prepare()
    print('locked meshes', len(read(ROOT / 'INPUT_MANIFEST.json')['rows']), 'own bytes', budget())
