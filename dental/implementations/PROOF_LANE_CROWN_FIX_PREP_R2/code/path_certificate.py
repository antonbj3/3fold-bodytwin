from common import *
from construction import exact_mesh
from loop_crown import loops

def check(row):
    m = load(row['mesh_path'])
    (v, f, roles) = (m['vertices'], m['faces'], m['roles'])
    (qv, qf) = (m['prep_vertices'], m['prep_faces'])
    M = m['M']
    d = M[2]
    stem = Path(row['mesh_path']).stem
    cp = D / (stem + '_path_crown.mesh')
    qp = D / (stem + '_path_Q.mesh')
    meshwrite(cp, v, f)
    meshwrite(qp, qv, qf)
    qcs = []
    cqs = []
    parts = trimesh.Trimesh(qv, qf, process=False).split(only_watertight=False)
    for (i, part) in enumerate(parts):
        partpath = D / (stem + '_path_component' + str(i) + '.mesh')
        meshwrite(partpath, part.vertices, part.faces)
        qcs.append(exact_mesh(cp, partpath, 0.039, stem + '_seated_Q_' + str(i) + '_outside_C'))
        cqs.append(exact_mesh(partpath, cp, 0.039, stem + '_seated_C_outside_Q_' + str(i)))
    disjoint = all((q['boundary_clearance_pass'] and c['boundary_clearance_pass'] and (not q['inside_witness']) and (not c['inside_witness']) and (q['query_boundary_components'] == 1) and (q['support_boundary_components'] == 1) for (q, c) in zip(qcs, cqs)))
    qc = qcs
    cq = cqs
    active = f[roles != 0][:, ::-1]
    norm = normals_exact(v, active, d)
    (av, af) = compact(v, active)
    active_mesh = trimesh.Trimesh(av, af, process=False)
    cc = len(trimesh.graph.connected_components(active_mesh.face_adjacency, nodes=np.arange(len(af))))
    ll = loops(v, active, M)
    lp = D / (stem + '_projection.txt')
    with lp.open('w') as h:
        h.write(str(len(ll[0])) + ' ' + ' '.join(('%.17g' % x for x in d)) + '\n')
        np.savetxt(h, v[ll[0]], fmt='%.17g')
    z = subprocess.run([str(D / 'exact_projected_loop'), str(lp)], capture_output=True, text=True, timeout=60)
    if z.returncode:
        raise ValueError('PROJECTED_LOOP_EXIT ' + str(z.returncode))
    projection = json.loads(z.stdout)
    projection.update(path=str(lp), sha256=sha(lp))
    dq = list(map(lambda x: F(float(x)), d))
    dots = [sum((F(float(a)) * b for (a, b) in zip(p, dq))) for p in v]
    floor = min(dots)
    disc = cc == 1 and active_mesh.euler_number == 1 and (len(ll) == 1) and projection['simple']
    path = R / 'raw' / (stem + '_PATH.json')
    out = dict(claim_type='capability', seated_disjoint_exact=bool(disjoint), Q_vs_C=qc, C_vs_Q=cq, actual_intaglio_direction=norm, intaglio_disk=dict(components=cc, euler=int(active_mesh.euler_number), boundary_loops=len(ll), pass_gate=disc), projected_boundary=projection, full_preparation=dict(definition='P=(T intersect {d dot x <= h}) union Q', d_exact=[str(x) for x in dq], h_exact=str(floor), subset_certificate=row['containment'].get('path', row['containment'].get('components')), subset_exact=row['containment']['all_pass'], lower_native_separation='Every crown facet is a convex hull of vertices with d dot v>=h; positive d dot velocity gives continuous separation from T_low.', scope='Exact implicit CSG preparation. No exact-containment assertion for a separately rounded full-P STL.'), isolated_path_pass=bool(disjoint and norm['pass_exact'] and disc), path_theorem='The labelled intaglio+shoulder is an embedded disk with positive exact projected facet Jacobian and simple projected boundary, hence a PL height graph along d. The generated nested core occupies its lower void (exact seated disjointness and labelled core/cavity construction). Translation of the crown by t*d, t>=0, increases graph height. Native lower halfspace remains below every crown vertex. Whole arch and physical tissues excluded.', peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024)
    dump(path, out)
    return dict(path=str(path), sha256=sha(path), isolated_path_pass=out['isolated_path_pass'], subset_exact=row['containment']['all_pass'], full_preparation_formula=out['full_preparation'])

def run(tag='C'):
    rows = []
    for r in read(R / f'RESULTS_{tag}.json')['rows']:
        if r['status'] != 'GENERATED':
            continue
        try:
            c = check(r)
        except Exception as e:
            c = dict(isolated_path_pass=False, error=repr(e))
            __import__('traceback').print_exc()
        rows.append(dict(key=r['key'], material=r['material'], certificate=c))
        dump(R / f'RESULTS_PATH_{tag}.json', dict(claim_type='capability', rows=rows))
        print(r['key'], r['material'], c, flush=True)
if __name__ == '__main__':
    run(sys.argv[1] if len(sys.argv) > 1 else 'C')
