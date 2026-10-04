from common import *
from construction import exact_mesh
from evaluate import boundary_edges, tool_sample, mf
from scipy.optimize import linprog
import warnings

def section_points(v, f, z, M):
    t = v[f] @ M.T
    out = []
    for (a, b) in [(0, 1), (1, 2), (2, 0)]:
        u = t[:, a]
        v2 = t[:, b]
        sel = (u[:, 2] - z) * (v2[:, 2] - z) < 0
        if sel.any():
            uu = u[sel]
            vv = v2[sel]
            p = uu + (vv - uu) * ((z - uu[:, 2]) / (vv[:, 2] - uu[:, 2]))[:, None]
            out.append(p)
    return np.concatenate(out) @ M if out else np.empty((0, 3))

def assess(row):
    rec = next((x for x in inputs() if x['key'] == row['key']))
    m = load(row['mesh_path'])
    (v, f, roles) = (m['vertices'], m['faces'], m['roles'])
    (qv, qf) = (m['prep_vertices'], m['prep_faces'])
    M = m['M']
    b = float(m['base'])
    mat = row['material_contract']
    stem = Path(row['mesh_path']).stem
    mm = trimesh.Trimesh(v, f, process=False)
    P = trimesh.Trimesh(qv, qf, process=False)
    ext = v[f[roles == 0]]
    inn = v[f[roles == 1]]
    if len(inn) == 0 or len(ext) == 0:
        raise ValueError('EMPTY_CROWN_REGION')
    S = m['target']
    sd = Distance(*triangles_mesh(S))
    shape = shape_score(ext, S)
    probes = np.r_[np.unique(inn.reshape(-1, 3), axis=0), sample(inn, 4096)]
    gapv = m.get('prep_cap_vertices', qv)
    gapf = m.get('prep_cap_faces', qf)
    gap = Distance(gapv, gapf).query(probes)[0]
    qlocal = qv[qf] @ M.T
    basal = abs(qlocal[:, :, 2] - b).max(1) < 1e-07
    insertion = normals_exact(qv, qf[~basal], M[2])
    n = P.face_normals[~basal]
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        lp = linprog([0, 0, 0, -1], A_ub=np.c_[-n, np.ones(len(n))], b_ub=np.zeros(len(n)), bounds=[(-1, 1)] * 3 + [(None, None)], method='highs', options={'threads': 1})
    ep = D / (stem + '_ext.mesh')
    ip = D / (stem + '_inn.mesh')
    meshwrite(ep, *triangles_mesh(ext))
    meshwrite(ip, *triangles_mesh(inn))
    wall = exact_mesh(ep, ip, mat['wall_mm'], stem + '_wall', True)
    outeredges = boundary_edges(f[roles == 0])
    marginpts = v[np.unique(outeredges)]
    (tv, tf, _) = __import__('common').native(rec)
    native_sec = section_points(tv, tf, b, M)
    nd = sd.query(native_sec)[0] if len(native_sec) else np.array([np.inf])
    crown_z = float((marginpts @ M.T)[:, 2].min())
    native_margin_error = sd.query(marginpts)[0]
    actual_margin_shift = float(abs(crown_z - b))
    margin = dict(native_section_points=len(native_sec), native_section_original_surface_max_mm=float(nd.max()), native_section_on_observed=bool(nd.max() <= 0.025), crown_margin_to_observed_max_mm=float(native_margin_error.max()), axial_finish_to_crown_mm=actual_margin_shift, scope='Planar digital finish-line candidate; CEJ UNKNOWN. All section intersections tested, no full curve source-error bound.')
    shoulder = dict(rounded_radius_mm=0.25, minimum_certified_radial_width_mm=mat['wall_mm'] + G, required_width_mm=mat['shoulder_min_mm'], full_material_compliance='UNKNOWN', reason='Anatomical support/regular thickness, physical CEJ and full preparation conditions unavailable; material width alone not sufficient')
    qparts = P.split(only_watertight=False) if 'components' in row['containment'] else [P]
    qchecks = [intersections(x.vertices, x.faces, stem + '_Q_part' + str(i)) for (i, x) in enumerate(qparts)]
    qsi = dict(count=sum((x.get('count', 0) for x in qchecks)), components=qchecks, scope='Each native-certified Q component; basal zero-measure contacts in graph union are handled by full P cap, not forced into one Surface_mesh')
    si = intersections(v, f, stem + '_crown')
    nc = len(trimesh.graph.connected_components(mm.face_adjacency, nodes=np.arange(len(f))))
    parityp = sample(inn, 64)
    aa = Distance(gapv, gapf).query(parityp)[0]
    (_, bb, _) = trimesh.proximity.closest_point(trimesh.Trimesh(gapv, gapf, process=False), parityp)
    parity = float(abs(aa - bb).max())
    stl = D / (stem + '_RESEARCH.stl')
    out3 = D / (stem + '_RESEARCH.3mf')
    mm.export(stl)
    mf.write(out3, v, f, roles, dict(case_key=row['case_key'], fdi=row['fdi'], material=row['material'], status='RESEARCH_NOT_QUALIFIED', frame='R4 native local millimetres'))
    (vv, ff, rr) = mf.readback(out3)
    exportok = np.array_equal(vv, v) and np.array_equal(ff, f) and np.array_equal(rr, roles)
    loaded = trimesh.load_mesh(stl, process=False)
    err = float(np.linalg.norm(loaded.triangles - mm.triangles, axis=2).max())
    import evaluate
    force = evaluate.contracts.field_query(dict(case_key=row['case_key'], geometry_sha256=row['mesh_sha256'], frame='R4_tooth_local_mm', fdi=row['fdi']))
    tool = tool_sample(v, f, roles)
    g = dict(native_containment='PASS' if row['containment']['all_pass'] and qsi.get('count') == 0 else 'FAIL', form='PASS' if shape['p95_mm'] <= read(R / 'PREREG_A.json')['metrics']['shape_p95_mm'][row['family']] else 'FAIL', margin='PASS' if margin['native_section_on_observed'] and actual_margin_shift <= 0.025 and (native_margin_error.max() <= 0.025) else 'FAIL', wall='PASS' if wall['all_pass'] else 'FAIL', nominal_gap='PASS' if abs(gap - G).max() <= 0.01 else 'FAIL', insertion='PASS' if insertion['pass_exact'] and lp.success and (lp.x[3] > 0) else 'FAIL', mesh='PASS' if mm.is_watertight and mm.is_winding_consistent and (mm.volume > 0) and (nc == 1) and (si.get('count') == 0) else 'FAIL', material_preparation='UNKNOWN', milling='UNKNOWN', STL='PASS' if exportok and err < 1e-05 else 'FAIL', full_arch_insertion='UNKNOWN', force=force['status'], fracture_q05='UNKNOWN', cement='UNKNOWN', lab='UNKNOWN')
    out = dict(**row, gates=g, complete_chain=all((x == 'PASS' for x in g.values())), complete_digital_crown=all((g[k] == 'PASS' for k in ['native_containment', 'form', 'margin', 'wall', 'nominal_gap', 'insertion', 'mesh', 'material_preparation', 'milling', 'STL'])), geometric_conjunction=all((g[k] == 'PASS' for k in ['native_containment', 'form', 'margin', 'wall', 'nominal_gap', 'insertion', 'mesh', 'STL'])), shape=shape, margin=margin, wall_certificate=wall, gap=dict(min_mm=float(gap.min()), max_mm=float(gap.max()), max_residual_mm=float(abs(gap - G).max()), samples=len(gap), continuous_enclosure='MISSING', resolution='PER_POINT'), insertion=insertion, lp_control=dict(success=lp.success, normal_margin=float(lp.x[3]) if lp.success else None), shoulder=shoulder, mesh=dict(watertight=bool(mm.is_watertight), self_intersections=si.get('count'), prep_self_intersections=qsi.get('count'), components=nc, core_volume_mm3=float(P.volume), crown_volume_mm3=float(mm.volume)), parity=dict(max_error_mm=parity, pass_gate=parity <= 1e-08), tools=tool, force=force, exports=dict(stl=str(stl), stl_sha256=sha(stl), three_mf=str(out3), three_mf_sha256=sha(out3), exact_roundtrip=bool(exportok), STL_max_vertex_rounding_mm=err, status='FROZEN_RESEARCH_ONLY'))
    dump(R / 'exports' / (stem + '_CERTIFICATE.json'), out)
    return out

def run(tag='A_PILOT'):
    rows = []
    for row in read(R / f'FROZEN_PREDICTIONS_{tag}.json')['rows']:
        t = time.perf_counter()
        if row['status'] == 'GENERATED':
            try:
                row = assess(row)
            except Exception as e:
                row = dict(**row, score_error=repr(e), complete_chain=False, complete_digital_crown=False, geometric_conjunction=False)
        else:
            row = dict(**row, complete_chain=False, complete_digital_crown=False, geometric_conjunction=False)
        row['score_seconds'] = time.perf_counter() - t
        rows.append(row)
        dump(R / f'RESULTS_{tag}.json', dict(rows=rows, claim_type='capability', external_referent=read(R / 'PREREG_A.json')['external_referent']))
        print(row['key'], row['material'], row.get('gates'), row.get('score_error'), flush=True)
if __name__ == '__main__':
    run(sys.argv[1] if len(sys.argv) > 1 else 'A_PILOT')
