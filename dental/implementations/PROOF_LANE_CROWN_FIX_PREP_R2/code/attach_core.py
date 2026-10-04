from common import *
from scipy.spatial import ConvexHull
from construction import exact_mesh

def dot(a, b):
    return sum((x * y for (x, y) in zip(a, b)))

def exact(p):
    return [F(float(x)) for x in p]

def cross(a, b):
    return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]

def generate():
    rows = []
    for part in read(R / 'EXPORTS.json')['rows']:
        m = load(part['mesh']['path'])
        d = m['M'][2]
        parts = trimesh.Trimesh(m['prep_vertices'], m['prep_faces'], process=False).split(only_watertight=False)
        bridges = []
        for (i, q) in enumerate(parts):
            tri = q.vertices[q.faces]
            eligible = np.flatnonzero(abs(tri @ d - float(m['base'])).max(1) < 1e-07)
            j = eligible[np.argmax(q.area_faces[eligible])]
            top = tri[j].copy()
            v = np.r_[top, top - 0.0001 * d]
            hull = ConvexHull(v)
            mesh = trimesh.Trimesh(v, hull.simplices, process=False)
            trimesh.repair.fix_normals(mesh)
            p = D / (Path(part['mesh']['path']).stem + '_F_bridge' + str(i) + '.mesh')
            meshwrite(p, mesh.vertices, mesh.faces)
            bridges.append(dict(component=i, path=str(p), sha256=sha(p), shared_source_face=q.faces[j].tolist(), top_coordinates=top.tolist(), vertices=mesh.vertices.tolist(), faces=mesh.faces.tolist()))
        rows.append(dict(key=part['key'], material=part['material'], round='F', parent_part=part, bridges=bridges))
    lock('FROZEN_PREDICTIONS_F.json', dict(claim_type='capability', rows=rows, prediction='Every prism keeps existing geometry gates and attaches Q to native T_low; one-carrier full-P film equality remains UNKNOWN.'))

def validate():
    rows = []
    for row in read(R / 'FROZEN_PREDICTIONS_F.json')['rows']:
        tick = time.perf_counter()
        part = row['parent_part']
        m = load(part['mesh']['path'])
        d = exact(m['M'][2])
        h = F(part['full_preparation_formula']['h_exact'])
        pc = read(part['path_certificate']['path'])
        rec = next((x for x in inputs() if x['key'] == part['key']))
        (tv, tf, ti) = native(rec)
        cp = D / (Path(part['mesh']['path']).stem + '_path_crown.mesh')
        checks = []
        for b in row['bridges']:
            bp = Path(b['path'])
            (v, f) = meshread(bp)
            t = [exact(p) for p in b['top_coordinates']]
            normal = cross([x - y for (x, y) in zip(t[1], t[0])], [x - y for (x, y) in zip(t[2], t[0])])
            bottom = [exact(p) for p in v[3:]]
            attachment = bool(any((set(face) == {0, 1, 2} for face in f)) and np.array_equal(v[:3], np.array(b['top_coordinates'])) and all((dot(p, d) < h for p in bottom)) and all((dot(normal, [x - y for (x, y) in zip(p, t[0])]) > 0 for p in bottom)))
            nativecert = exact_mesh(ti['path'], bp, part['wall_lower_bound_mm'] + 0.05, bp.stem + '_native')
            cq = exact_mesh(cp, bp, 0.039, bp.stem + '_seated')
            qc = exact_mesh(bp, cp, 0.039, bp.stem + '_reciprocal')
            z = subprocess.run([str(D / 'exact_projection_inside'), pc['projected_boundary']['path'], str(bp)], capture_output=True, text=True, timeout=30)
            projection = json.loads(z.stdout)
            si = intersections(v, f, bp.stem)
            mesh = trimesh.Trimesh(v, f, process=False)
            ok = attachment and nativecert['all_pass'] and cq['boundary_clearance_pass'] and qc['boundary_clearance_pass'] and (not cq['inside_witness']) and (not qc['inside_witness']) and projection['all_pass'] and (si['count'] == 0) and mesh.is_watertight and (mesh.volume > 0)
            checks.append(dict(component=b['component'], all_pass=bool(ok), exact_facet_attachment=attachment, below_h_rational_margins=[str(h - dot(p, d)) for p in bottom], native=nativecert, crown_vs_bridge=cq, bridge_vs_crown=qc, projection=projection, projection_command=[str(D / 'exact_projection_inside'), pc['projected_boundary']['path'], str(bp)], mesh=dict(watertight=bool(mesh.is_watertight), self_intersections=si['count'], volume_mm3=float(mesh.volume)), resolution='PER_POINT'))
        rows.append(dict(key=row['key'], material=row['material'], all_pass=all((x['all_pass'] for x in checks)), components=checks, formula='Pplus=(T intersect H_low) union Q union native-certified facet-sharing prisms', scope='Every Q component attached to retained lower native set, same isolated crown path. Does not certify connectivity of T_low itself, equality to nominal film cap, CEJ, material or physical qualification.', seconds=time.perf_counter() - tick))
        dump(R / 'RESULTS_F.json', dict(claim_type='capability', external_referent=read(R / 'PREREG_F.json')['external_referent'], rows=rows))
        print(row['key'], row['material'], rows[-1]['all_pass'], flush=True)
if __name__ == '__main__':
    generate() if '--generate' in sys.argv else validate()
