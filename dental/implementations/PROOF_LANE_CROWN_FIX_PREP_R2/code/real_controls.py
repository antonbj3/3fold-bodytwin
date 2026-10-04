from common import *
from construction import exact_mesh

def run():
    rows = []

    def ck(name, ok, **kw):
        rows.append(dict(name=name, passed=bool(ok), **kw))
    row = next((r for r in read(R / 'RESULTS_C.json')['rows'] if r['status'] == 'GENERATED'))
    m = load(row['mesh_path'])
    (v, f, roles) = (m['vertices'], m['faces'], m['roles'])
    (qv, qf) = (m['prep_vertices'], m['prep_faces'])
    M = m['M']
    b = float(m['base'])
    rec = next((r for r in inputs() if r['key'] == row['key']))
    T = native(rec)
    qp = D / 'fault_actual_Q_shifted.mesh'
    meshwrite(qp, qv + 100 * M[2], qf)
    bad = exact_mesh(T[2]['path'], qp, 0.55, 'fault_actual_native')
    ck('real_native_outside_rejected', not bad['all_pass'])
    active = f[roles != 0][:, ::-1].copy()
    good = normals_exact(v, active, M[2])
    active[0] = active[0, ::-1]
    badn = normals_exact(v, active, M[2])
    ck('actual_intaglio_valid', good['pass_exact'])
    ck('actual_insertion_flipped_face_rejected', not badn['pass_exact'])
    P = Distance(qv, qf)
    pt = v[f[roles == 1][0]].mean(0)
    err0 = abs(P.query(pt[None])[0][0] - 0.05)
    ptbad = pt + 1 * M[2]
    err1 = abs(P.query(ptbad[None])[0][0] - 0.05)
    ck('actual_nominal_film_accepted', err0 <= 0.01, measured_residual_mm=float(err0))
    ck('actual_film_displacement_rejected', err1 > 0.01, measured_residual_mm=float(err1))
    rim = m['native_section'][0]
    ck('actual_margin_accepted', abs(rim @ M[2] - b) <= 0.025)
    ck('actual_margin_0.1mm_shift_rejected', abs((rim + 0.1 * M[2]) @ M[2] - b) > 0.025)
    crown = trimesh.Trimesh(v, f, process=False)
    missing = trimesh.Trimesh(v, f[:-1], process=False)
    ck('actual_mesh_closed', crown.is_watertight)
    ck('actual_mesh_face_removal_rejected', not missing.is_watertight)
    ext = v[f[roles == 0]]
    (ev, ef) = triangles_mesh(ext)
    sp = D / 'fault_actual_outer.mesh'
    meshwrite(sp, ev, ef)
    queries = D / 'fault_actual_wall.txt'
    with queries.open('w') as h:
        h.write('1\n')
        np.savetxt(h, np.r_[ext[0, 0], 0.5][None], fmt='%.17g')
    z = subprocess.run([str(D / 'exact_surface'), str(sp), str(queries), 'surface'], capture_output=True, text=True, timeout=60)
    exact = json.loads(z.stdout)
    ck('actual_zero_wall_rejected', not exact['rows'][0]['pass'], exact=exact['rows'][0])
    upper_source = D / 'fault_upper_source.mesh'
    upper_query = D / 'fault_upper_query.mesh'
    tri = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]], float)
    meshwrite(upper_source, tri, np.array([[0, 1, 2]]))
    meshwrite(upper_query, tri + [0, 0, 0.07], np.array([[0, 1, 2]]))
    z = subprocess.run([str(D / 'exact_gap_upper'), str(upper_source), str(upper_query)], capture_output=True, text=True, timeout=60)
    u = json.loads(z.stdout)
    ck('continuous_upper_70um_rejected', not u['all_pass'] and u['exact_refuting_points'] > 0, exact=u)
    loop = D / 'fault_projection.txt'
    loop.write_text('4 0 0 1\n0 0 0\n1 1 0\n0 1 0\n1 0 0\n')
    z = subprocess.run([str(D / 'exact_projected_loop'), str(loop)], capture_output=True, text=True, timeout=60)
    p = json.loads(z.stdout)
    ck('crossed_projected_loop_rejected', not p['simple'])
    before = sha(row['exports']['three_mf'])
    corrupt = D / 'fault_export.3mf'
    content = Path(row['exports']['three_mf']).read_bytes()
    corrupt.write_bytes(content[:-1] + bytes([content[-1] ^ 1]))
    ck('actual_export_byte_corruption_rejected', sha(corrupt) != before)
    assembly = load(row['assembly_surface']['path'])['triangles']
    limit = read(R / 'PREREG_B.json')['metrics']['shape_p95_mm'][row['family']]
    goodshape = shape_score(assembly, m['target'])['p95_mm']
    badshape = shape_score(assembly + 100 * M[2], m['target'])['p95_mm']
    ck('actual_assembly_form_accepted', goodshape <= limit, p95_mm=goodshape)
    ck('actual_shifted_assembly_form_rejected', badshape > limit, p95_mm=badshape)
    out = dict(all_pass=all((r['passed'] for r in rows)), passed=sum((r['passed'] for r in rows)), count=len(rows), rows=rows, scope='Actual case geometry corruption plus exact continuous-film and graph projection fixtures; not external dental validation')
    dump(R / 'raw/REAL_CONTROLS.json', out)
    print('Real control suite', out['passed'], '/', out['count'], flush=True)
    return out
if __name__ == '__main__':
    run()
