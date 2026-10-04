from common import *
from attach_core import exact, dot, cross

def run():
    count = 0
    for row in read(R / 'FROZEN_PREDICTIONS_F.json')['rows']:
        part = row['parent_part']
        m = load(part['mesh']['path'])
        d = exact(m['M'][2])
        h = F(part['full_preparation_formula']['h_exact'])
        pcs = trimesh.Trimesh(m['prep_vertices'], m['prep_faces'], process=False).split(only_watertight=False)
        pc = read(part['path_certificate']['path'])
        saved = next((x for x in read(R / 'RESULTS_F.json')['rows'] if x['key'] == row['key'] and x['material'] == row['material']))
        for b in row['bridges']:
            (v, f) = meshread(b['path'])
            top = np.array(b['top_coordinates'])
            q = pcs[b['component']]
            tri = q.vertices[q.faces]
            same = bool(np.any(np.all(tri == top, axis=(1, 2))))
            assert same, 'Missing actual shared source facet'
            assert np.array_equal(v[:3], top) and np.array_equal(v[3:], top - 0.0001 * m['M'][2])
            assert any((set(face) == {0, 1, 2} for face in f))
            bridge = trimesh.Trimesh(v, f, process=False)
            assert bridge.is_watertight and bridge.is_winding_consistent and (bridge.volume > 0)
            assert intersections(v, f, Path(b['path']).stem + '_REPLAY')['count'] == 0
            t = [exact(x) for x in top]
            n = cross([x - y for (x, y) in zip(t[1], t[0])], [x - y for (x, y) in zip(t[2], t[0])])
            assert all((dot(exact(p), d) < h and dot(n, [x - y for (x, y) in zip(exact(p), t[0])]) > 0 for p in v[3:]))
            z = subprocess.run([str(D / 'exact_projection_inside'), pc['projected_boundary']['path'], b['path']], capture_output=True, text=True, check=True)
            out = json.loads(z.stdout)
            assert out == saved['components'][b['component']]['projection'] and out['all_pass']
            count += 1
    row = read(R / 'FROZEN_PREDICTIONS_F.json')['rows'][0]
    b = row['bridges'][0]
    part = row['parent_part']
    m = load(part['mesh']['path'])
    pc = read(part['path_certificate']['path'])
    (v, f) = meshread(b['path'])
    bad = D / 'fault_F_projected_outside.mesh'
    meshwrite(bad, v + 100 * m['M'][0], f)
    z = subprocess.run([str(D / 'exact_projection_inside'), pc['projected_boundary']['path'], str(bad)], capture_output=True, text=True, check=True)
    outside = not json.loads(z.stdout)['all_pass']
    top = np.array(b['top_coordinates'])
    changed = top.copy()
    changed[0, 0] += 1e-05
    shared_not_equal = not np.array_equal(changed, v[:3])
    d = exact(m['M'][2])
    h = F(part['full_preparation_formula']['h_exact'])
    badbottom = v[3:] + m['M'][2]
    above = not all((dot(exact(p), d) < h for p in badbottom))
    result = dict(status='PASS' if outside and shared_not_equal and above else 'FAIL', attachments_replayed=count, original_carrier_disconnection='8/8 exact positive slabs retained', fault_controls=dict(projected_outside_rejected=outside, changed_shared_facet_rejected=shared_not_equal, bottom_above_native_lower_rejected=above), scope='Each original core component gains a shared-facet bridge into T_low. T_low connectivity and equality with nominal film cap still separate UNKNOWNs.')
    assert result['status'] == 'PASS'
    dump(R / 'raw/ATTACHMENT_REPLAY.json', result)
    print('Exact attachments replayed', count, 'new fault controls3/3', flush=True)
    return result
if __name__ == '__main__':
    run()
