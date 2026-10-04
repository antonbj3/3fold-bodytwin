from dental_release.paths import expand as _release_expand
from columns import *
start = time.perf_counter()
r2 = json.loads((ROOT / 'raw/R2.json').read_text())
disk = [r for r in r2['rows'] if r['topology_gate']]
disk.sort(key=lambda r: min((q['max_lift_source_unit'] for q in r['directions'])))
pilots = [disk[0], disk[-1]]
rows = []
for source in pilots:
    case = source['case']
    z = np.load(DATA / (case + '_R1.npz'))
    v = z['v']
    f = z['f'][z['lab'] == 1]
    tri = v[f]
    chosen = min(source['directions'], key=lambda q: q['max_lift_source_unit'])
    d = np.array(chosen['direction'])
    (q, basis) = coords(v, d)
    qtri = q[f]
    sub = []
    for h in [0.12, 0.06, 0.03]:
        began = time.perf_counter()
        (height, offset, cover) = raster(qtri, h)
        (mesh, size) = mesh_columns(height, h, offset, basis)
        outcome = continuous_gap(mesh, tri)
        cent = tri.mean(axis=1)
        selected = cent[np.linspace(0, len(cent) - 1, 16, dtype=int)]
        points = np.vstack([selected, np.array(outcome['witness'])])
        dist = distance(mesh, points)
        control = np.array([trimesh.proximity.closest_point_naive(mesh, p[None, :])[1][0] for p in points])
        error = float(np.max(np.abs(dist - control)))
        path = DATA / f'{case}_CAVITY_h{h:.2f}.stl'
        mesh.export(path)
        loaded = trimesh.load(path, force='mesh')
        back = distance(loaded, points)
        roundtrip = float(np.max(np.abs(back - dist)))
        if not (mesh.is_watertight and mesh.is_winding_consistent and cover['cover_gate'] and (error <= 1e-05) and (roundtrip <= 1e-05)):
            outcome['qualified_gate'] = 'UNKNOWN_VALIDATION_OR_TOPOLOGY'
        else:
            outcome['qualified_gate'] = outcome['gate']
        np.savez_compressed(DATA / f'{case}_COLUMNS_h{h:.2f}.npz', height=height, h=h, offset=offset, basis=basis)
        result = {'h_source_unit': h, 'cover': cover, 'mesh': size, 'watertight': bool(mesh.is_watertight), 'winding_consistent': bool(mesh.is_winding_consistent), 'components': len(mesh.split(only_watertight=False)), 'gap': outcome, 'distance_control_error': error, 'roundtrip_distance_error': roundtrip, 'stl': str(path), 'stl_sha256': sha(path), 'elapsed_seconds': time.perf_counter() - began}
        sub.append(result)
        write(f'raw/R3_{case}_h{h:.2f}.json', result)
        print(case, h, result['gap']['qualified_gate'], outcome['lower_max'], outcome['upper_max'], cover, size, flush=True)
    rows.append({'case': case, 'selection': 'R2 selected extreme pilot, not independent validation', 'direction_name': chosen['name'], 'direction': d.tolist(), 'source_native_units': 'unit metadata UNKNOWN', 'grids': sub, 'physical_gate': 'UNKNOWN'})
out = {'id': 'XBREAK5-R3', 'claim_type': 'capability', 'rows': rows, 'primary_gate': 'CONDITIONAL_PASS' if any((r['grids'][-1]['gap']['qualified_gate'] == 'CONDITIONAL_PASS' for r in rows)) else 'FAIL_OR_UNKNOWN', 'physical_gate': 'UNKNOWN', 'wall_seconds': time.perf_counter() - start}
write('raw/R3.json', out)
print(json.dumps({k: v for (k, v) in out.items() if k != 'rows'}))
write('CURRENT_WORK_STATE.json', {'lane': _release_expand('XBREAK_HUNT_5'), 'status': 'R3_COMPLETE', 'latest_gate': out['primary_gate'], 'next_operation': 'R4 budget-aware cover refinement or measured unit/margin/cavity consumer contract, guided by stored failure', 'physical_gate': 'UNKNOWN'})
