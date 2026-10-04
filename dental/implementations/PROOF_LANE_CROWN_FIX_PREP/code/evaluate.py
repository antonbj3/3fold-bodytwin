from core import *
from scipy.spatial import cKDTree
from scipy.optimize import linprog
import warnings
cadmargin = module('fixprep_cadmargin', B / 'PROOF_LANE_CAD_10X/code/margin.py')
sys.path.insert(0, str(B / 'PROOF_LANE_CAD_10X/code'))
mf = module('fixprep_3mf', B / 'PROOF_LANE_CAD_10X/code/three_mf.py')
sys.path.pop(0)
contracts = module('fixprep_contracts', B / 'PROOF_LANE_COMBINE_CROWN_CHAIN/code/contracts.py')
threadpool_limits(1)

def boundary_edges(f):
    e = np.concatenate([f[:, [0, 1]], f[:, [1, 2]], f[:, [2, 0]]])
    (u, ct) = np.unique(np.sort(e, axis=1), axis=0, return_counts=True)
    return u[ct == 1]

def assess(row):
    tick = time.perf_counter()
    m = load(row['mesh_path'])
    (v, f, roles) = (m['vertices'], m['faces'], m['roles'])
    (qv, qf) = (m['prep_vertices'], m['prep_faces'])
    M = m['M']
    mat = row['material_contract']
    mm = trimesh.Trimesh(v, f, process=False)
    P = trimesh.Trimesh(qv, qf, process=False)
    ext = v[f[roles == 0]]
    inn = v[f[roles == 1]]
    (rv, rf) = triangles_mesh(inn)
    dist = Distance(qv, qf)
    probes = np.r_[rv, sample(inn, 4096)]
    gap = dist.query(probes)[0]
    shape = shape_score(ext, m['target'])
    di = Distance(rv, rf)
    wd = di.query(sample(ext, 4096))[0]
    qtri = qv[qf]
    basal = np.max(np.abs(qtri @ M[2] - float(m['base'])), axis=1) < 1e-07
    nonbasal = qf[~basal]
    exact = normals_exact(qv, nonbasal, M[2])
    normals = P.face_normals[~basal]
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        lp = linprog([0, 0, 0, -1], A_ub=np.c_[-normals, np.ones(len(normals))], b_ub=np.zeros(len(normals)), bounds=[(-1, 1)] * 3 + [(None, None)], method='highs', options={'threads': 1})
    detected = cadmargin.detect(P)
    (outercycles, _) = cadmargin.cycles(boundary_edges(f[roles == 0]), v)
    margin_error = None
    if detected['status'] == 'PROPOSED' and len(outercycles) == 1:
        margin_error = cadmargin.errors(detected['points'], outercycles[0]['points'])
    shoulder = dict(status='FAIL', reason='Compact core basal rim is not a whole-tooth shoulder; no measured/constructed rounded transition to source tooth in this round', material_required_width_mm=mat['shoulder_min_mm'])
    meshsi = intersections(v, f, Path(row['mesh_path']).stem + '_final')
    nc = len(trimesh.graph.connected_components(mm.face_adjacency, nodes=np.arange(len(f))))
    stem = Path(row['mesh_path']).stem
    stl = D / (stem + '_crown_RESEARCH.stl')
    out3 = D / (stem + '_crown_RESEARCH.3mf')
    mm.export(stl)
    mf.write(out3, v, f, roles, dict(case_key=row['case_key'], fdi=row['fdi'], material=row['material'], geometry_sha256=row['mesh_sha256'], status='RESEARCH_NOT_QUALIFIED', frame='original R4 tooth-local millimetres', preparation=row['prep_stl']))
    (vv, ff, rr) = mf.readback(out3)
    exportok = np.array_equal(vv, v) and np.array_equal(ff, f) and np.array_equal(rr, roles)
    loaded = trimesh.load_mesh(stl, process=False)
    stlerr = float(np.linalg.norm(loaded.triangles - mm.triangles, axis=2).max())
    crownport = dict(case_key=row['case_key'], geometry_sha256=row['mesh_sha256'], frame='R4_tooth_local_mm', fdi=row['fdi'])
    force = contracts.field_query(crownport)
    tool = tool_sample(v, f, roles)
    from wall_certificate import prove
    wallproof = prove(ext, inn, row)
    wallpass = bool(wallproof['all_pass'] and wd.min() >= mat['wall_mm'])
    globalcert = row['certificates']['Q_dilated_inside_native']
    globalstatus = 'PASS' if globalcert['all_pass'] else 'UNKNOWN'
    obstruction = None
    if not globalcert['all_pass']:
        cr = read(globalcert['path'])
        outside = [i for (i, x) in enumerate(cr['rows']) if x.get('inside') is False]
        if outside:
            qi = outside[0]
            qq = np.loadtxt(cr['query_path'], skiprows=1)[qi].copy()
            qq[3] = 0.0
            qp = D / (stem + '_core_for_membership.mesh')
            meshwrite(qp, qv, qf)
            member = certify(qp, qq[None], stem + '_OUTSIDE_WITNESS')
            if member['all_pass']:
                globalstatus = 'FAIL'
                obstruction = dict(point_mm=qq[:3], in_Q_exact=True, outside_T_exact=True, Q_membership=member, T_membership_record=globalcert['path'], T_membership_row=qi, statement='Q itself has a point outside T, so Q+B_w is not a subset of T for any nonnegative w.')
    straight_faces = np.arange(128)
    sn = P.face_normals[straight_faces] @ M.T
    taper = 2 * np.degrees(np.arctan2(sn[:, 2], np.linalg.norm(sn[:, :2], axis=1)))
    convergence = dict(min_total_degrees=float(taper.min()), max_total_degrees=float(taper.max()), resolution='PER_SURFACE_REGION', scope='first straight circumferential facet strip; top fillet excluded')
    g = dict(native_containment=globalstatus, convergence='PASS' if taper.min() >= 10 and taper.max() <= 20 else 'FAIL', form='PASS' if shape['p95_mm'] <= read(R / 'PREREG_A.json')['metrics']['shape_p95_mm'][row['family']] else 'FAIL', margin='PASS' if margin_error and margin_error['max_um'] <= 25 else 'FAIL', wall='PASS' if wallpass else 'FAIL', nominal_gap='PASS' if abs(gap - G).max() <= 0.01 else 'FAIL', insertion='PASS' if exact['pass_exact'] and lp.success and (lp.x[3] > 0) else 'FAIL', mesh='PASS' if mm.is_watertight and mm.is_winding_consistent and (mm.volume > 0) and (nc == 1) and (meshsi.get('count') == 0) else 'FAIL', material_preparation='FAIL', milling='UNKNOWN', STL='PASS' if exportok and stlerr < 1e-05 else 'FAIL', full_arch_insertion='UNKNOWN', force=force['status'], fracture_q05='UNKNOWN', cement='UNKNOWN', lab='UNKNOWN')
    return dict(**row, gates=g, complete_digital_crown=all((g[k] == 'PASS' for k in ['native_containment', 'form', 'margin', 'wall', 'nominal_gap', 'insertion', 'mesh', 'material_preparation', 'milling', 'STL'])), complete_chain=all((x == 'PASS' for x in g.values())), actual_mesh_wall_proof=wallproof, native_containment_obstruction=obstruction, convergence=convergence, shape=shape, wall=dict(sampled_min_mm=float(wd.min()), certified_bound_mm=mat['wall_mm'] if wallpass else None, scope='PER_SURFACE_REGION exterior/intaglio; cervical rim excluded'), gap=dict(min_mm=float(gap.min()), max_mm=float(gap.max()), max_error_mm=float(abs(gap - G).max()), samples=len(gap), scope='nominal geometric film, sampled at all cavity vertices plus4096 area-stratified centroids; continuous film enclosure absent'), insertion=exact, lp_control=dict(success=lp.success, normal_margin=float(lp.x[3]) if lp.success else None), cad10x_margin=dict(status=detected['status'], diagnostics=detected['diagnostics'], error_to_crown_outer_margin=margin_error), shoulder=shoulder, mesh=dict(watertight=bool(mm.is_watertight), components=nc, self_intersections=meshsi.get('count'), volume_mm3=float(mm.volume)), tools=tool, force=force, physical_debts=['Same-case registered bite', 'Whole-arch insertion trajectory', 'CAM tool/shank/blank path', 'Same-specimen seated spatial cement film', 'Material-process-matched fracture test', 'Pulp and true finish line'], exports=dict(stl=str(stl), stl_sha256=sha(stl), three_mf=str(out3), three_mf_sha256=sha(out3), three_mf_exact_roundtrip=bool(exportok), STL_max_vertex_rounding_mm=stlerr, status='FROZEN_RESEARCH_ONLY'), seconds_scoring=time.perf_counter() - tick)

def tool_sample(v, f, roles):
    saved = {n: sys.modules.get(n) for n in ['common', 'geometry', 'collision', 'run_r1']}
    old = list(sys.path)
    try:
        for n in saved:
            sys.modules.pop(n, None)
        p = B / 'LANE_X53_MILLING_TOOLS/code'
        sys.path.insert(0, str(p))
        collision = module('collision', p / 'collision.py')
        run = module('run_r1', p / 'run_r1.py')
        valid = np.linalg.norm(np.cross(v[f[:, 1]] - v[f[:, 0]], v[f[:, 2]] - v[f[:, 0]]), axis=1) / 2 >= 1e-14
        scene = collision.Scene(v, f[valid])
        points = run.sample(scene, np.flatnonzero(roles[valid] == 1), seed=5341, n=4)
        tools = run.lib('vhf', 5)
        rows = []
        for pt in points:
            trials = [dict(tool_id=t['id'], **collision.search(scene, pt['point'], pt['normal'], t, 5, offsets=(0.0,))) for t in tools]
            rows.append(dict(point=pt, tools=trials, has_witness=any((x['status'] == 'FOUND' for x in trials))))
        return dict(status='UNKNOWN', dropped_degenerate_facets=int((~valid).sum()), samples=rows, found=sum((x['has_witness'] for x in rows)), count=len(rows), scope='Executed X53 local sampled tool envelope; no continuous path or actual OEM/shank/blank calibration')
    except Exception as e:
        return dict(status='UNKNOWN', error=repr(e), scope='Tool computation failed, never interpreted as access')
    finally:
        sys.path[:] = old
        for (n, x) in saved.items():
            if x is None:
                sys.modules.pop(n, None)
            else:
                sys.modules[n] = x
        threadpool_limits(1)

def run(tag):
    pr = read(R / f'FROZEN_PREDICTIONS_{tag}.json')
    rows = []
    for r in pr['rows']:
        if r['status'] == 'GENERATED':
            try:
                r = assess(r)
            except Exception as e:
                r = dict(**r, score_error=repr(e), complete_digital_crown=False, complete_chain=False)
        else:
            r = dict(**r, complete_digital_crown=False, complete_chain=False, gates={k: 'UNKNOWN_UPSTREAM_REJECTED' for k in ['form', 'margin', 'wall', 'nominal_gap', 'insertion', 'mesh', 'material_preparation', 'milling', 'STL', 'full_arch_insertion', 'force', 'fracture_q05', 'cement', 'lab']})
        rows.append(r)
        dump(R / f'RESULTS_{tag}.json', dict(claim_type='capability', rows=rows, external_referent=read(R / 'PREREG_A.json')['external_referent']))
        print('SCORE', r['key'], r['material'], r.get('gates'), r.get('score_error'), flush=True)
if __name__ == '__main__':
    run(sys.argv[1])
