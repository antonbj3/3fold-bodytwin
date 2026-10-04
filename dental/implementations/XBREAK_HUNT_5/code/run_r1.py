from dental_release.paths import expand as _release_expand
from common import *
from scipy.optimize import linprog
start = time.perf_counter()
inputs = json.loads((ROOT / 'sources/INPUTS.json').read_text())
rows = []
for source in inputs['test_files']:
    (v, f, lab) = vtp(source['file'])
    (a, e1, e2, n, area, ct) = frame(v, f, lab)
    (topo, b) = topology(v, f, lab)
    na = n @ a
    nu = n @ e1
    nv = n @ e2
    opt = linprog([0, 0, -1], A_ub=np.c_[-nu, -nv, np.ones(len(n))], b_ub=na, bounds=[(-1, 1), (-1, 1), (None, None)], method='highs')
    feasible = linprog([0.0001, -0.0001], A_ub=np.c_[-nu, -nv], b_ub=na + 1e-08, bounds=[(-1, 1), (-1, 1)], method='highs')
    d = a + opt.x[0] * e1 + opt.x[1] * e2
    d /= np.linalg.norm(d)
    minimum = float(np.min(n @ d))
    cone = bool(minimum >= -1e-08)
    seq = orientloop(b)
    from shapely.geometry import LinearRing
    helper = np.eye(3)[np.argmin(np.abs(d))]
    q1 = np.cross(d, helper)
    q1 /= np.linalg.norm(q1)
    q2 = np.cross(d, q1)
    simple = bool(seq is not None and LinearRing(np.c_[v[seq] @ q1, v[seq] @ q2]).is_simple)
    import vtk
    pts = vtk.vtkPoints()
    from vtk.util.numpy_support import numpy_to_vtk
    pts.SetData(numpy_to_vtk(v, deep=True))
    ca = vtk.vtkCellArray()
    for face in f[lab == 1]:
        ca.InsertNextCell(3)
        for q in face:
            ca.InsertCellPoint(int(q))
    poly = vtk.vtkPolyData()
    poly.SetPoints(pts)
    poly.SetPolys(ca)
    fe = vtk.vtkFeatureEdges()
    fe.SetInputData(poly)
    fe.BoundaryEdgesOn()
    fe.FeatureEdgesOff()
    fe.ManifoldEdgesOff()
    fe.NonManifoldEdgesOff()
    fe.Update()
    boundary_control = fe.GetOutput().GetNumberOfCells()
    r = {'case': source['name'].split('/')[-1].split('_')[0], 'source_sha256': source['sha256'], 'vertices': len(v), 'faces': len(f), 'zero_area_top_facets': int(np.sum(area == 0)), 'topology': topo, 'vtk_boundary_edge_count': boundary_control, 'boundary_control_error': int(boundary_control - len(b)), 'axis': a.tolist(), 'direction': d.tolist(), 'cone_optimum_unnormalized': float(opt.x[2]), 'min_normal_dot_unit_direction': minimum, 'negative_area_fraction': float(np.sum(area[n @ d < -1e-08]) / area.sum()), 'candidate_cone_pass': cone, 'control_status': int(feasible.status), 'control_cone_pass': bool(feasible.success), 'agreement': bool(cone == feasible.success), 'digital_gate': bool(cone and topo['disk_gate'] and (boundary_control == len(b))), 'physical_gate': 'UNKNOWN', 'label_semantics': 'TOP_1_BY_PUBLIC_BINARY_CONTEXT; scanner/metrology UNKNOWN'}
    r['projected_boundary_simple'] = simple
    r['digital_gate'] = bool(r['digital_gate'] and simple)
    rows.append(r)
    np.savez_compressed(DATA / (r['case'] + '_R1.npz'), v=v, f=f, lab=lab, axis=a, e1=e1, e2=e2, d=d, boundary=b)
out = {'id': 'XBREAK5-R1', 'claim_type': 'capability', 'rows': rows, 'available_n': len(rows), 'source_reported_n': 13, 'missing_n': 13 - len(rows), 'digital_pass_n': sum((r['digital_gate'] for r in rows)), 'cone_pass_n': sum((r['candidate_cone_pass'] for r in rows)), 'control_agreement_n': sum((r['agreement'] for r in rows)), 'numeric_gate': 'PASS' if all((r['agreement'] and r['boundary_control_error'] == 0 for r in rows)) else 'FAIL', 'primary_gate': 'PASS' if any((r['digital_gate'] for r in rows)) else 'FAIL', 'physical_gate': 'UNKNOWN', 'wall_seconds': time.perf_counter() - start}
write('raw/R1.json', out)
print(json.dumps({k: x for (k, x) in out.items() if k != 'rows'}))
write('CURRENT_WORK_STATE.json', {'lane': _release_expand('XBREAK_HUNT_5'), 'status': 'R1_COMPLETE', 'latest_gate': out['primary_gate'], 'next_operation': 'Freeze directed-sweep representation R2; rigid normal cone is not a manufacturability verdict', 'physical_gate': 'UNKNOWN'})
