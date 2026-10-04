from common import *
import sys, trimesh
from scipy.sparse import coo_matrix
sys.path.insert(0, str(X1 / 'vendor'))
import crown_design_fe as FE
from crown_fit_geometry import Grid

def model():
    a = np.load(X1 / 'inputs/geometry/D1_model.npz')
    g = read(X1 / 'inputs/geometry/D1_grid.json')
    grid = Grid(g['origin'], np.array(g['origin']) + g['h'] * (np.array(g['shape']) - 1), g['h'])
    return (a, g, FE.Model((a['cV'], a['cT']), (a['dV'], a['dT']), a['crown'], a['die'], grid, g['z_m'], g['z_m'] - 3))

def deck_contract():
    lines = (X1 / 'inputs/fe/L005_lo_M1_k6_thin_3Y.inp').read_text().splitlines()
    blocks = []
    cur = None
    for line in lines:
        if line.startswith('**') or not line.strip():
            continue
        if line.startswith('*'):
            cur = {'header': line.upper(), 'rows': []}
            blocks.append(cur)
        elif cur is not None:
            cur['rows'].append(line)

    def block(key):
        return next((b for b in blocks if key in b['header']))
    slave = np.array([int(x) - 1 for r in block('*NSET, NSET=NSLAVE')['rows'] for x in r.split(',')], int)
    masters = np.array([[int(r.split(',')[0]) - 1, int(r.split(',')[1][1:])] for r in block('*SURFACE, NAME=SMASTER')['rows']], int)
    fixes = np.array([int(x) - 1 for r in block('*NSET, NSET=NFIX')['rows'] for x in r.split(',')], int)
    nodes = np.array([[float(x) for x in row.split(',')[1:4]] for row in block('*NODE')['rows']])
    loads = []
    for b in blocks:
        if b['header'].startswith('*CLOAD'):
            f = np.zeros_like(nodes)
            for row in b['rows']:
                (n, k, v) = row.split(',')
                f[int(n) - 1, int(k) - 1] += float(v)
            loads.append(f)
    mats = [b['rows'][0] for b in blocks if b['header'].startswith('*ELASTIC')]
    return dict(slaves=slave, master_owners=masters, fixes=fixes, nodes=nodes, loads=loads, materials=mats, tie_header=block('*TIE')['header'])

def project(m):
    slaves = np.unique(m.Fc_tied[:, 2:8])
    mesh = trimesh.Trimesh(m.Vd, m.Fd_master[:, 2:5], process=False)
    points = []
    dist = []
    face = []
    for start in range(0, len(slaves), 128):
        (p, d, i) = trimesh.proximity.closest_point(mesh, m.Vc[slaves[start:start + 128]])
        points.append(p)
        dist.append(d)
        face.append(i)
    points = np.concatenate(points)
    dist = np.concatenate(dist)
    face = np.concatenate(face)
    tri = m.Fd_master[face, 2:8]
    bary = trimesh.triangles.points_to_barycentric(m.Vd[tri[:, :3]], points)
    return (slaves, points, dist, face, tri, bary)

def q_matrix(n, slaves, master_global, w):
    active = np.setdiff1d(np.arange(n), slaves)
    mapping = np.full(n, -1, int)
    mapping[active] = np.arange(len(active))
    rows = np.r_[active, np.repeat(slaves, w.shape[1])]
    cols = np.r_[mapping[active], mapping[master_global].ravel()]
    vals = np.r_[np.ones(len(active)), w.ravel()]
    assert (cols >= 0).all()
    return (coo_matrix((vals, (rows, cols)), shape=(n, len(active))).tocsr(), mapping)

def p2(b):
    return np.c_[b * (2 * b - 1), 4 * b[:, 0] * b[:, 1], 4 * b[:, 1] * b[:, 2], 4 * b[:, 2] * b[:, 0]]

def port_checks(m, slaves, face, tri, Q, w, deck):
    p = read(SRC / 'PREREG_R1.json')['metrics']
    native_master = m.Fd_master[:, :2].copy()
    native_master[:, 0] += len(m.Tc)
    parity = dict(slave_nodes=bool(np.array_equal(np.sort(slaves), np.sort(deck['slaves']))), master_owners=bool(np.array_equal(native_master, deck['master_owners'])), fixed_nodes=bool(np.array_equal(np.sort(m.nc + m.fix_nodes), np.sort(deck['fixes']))))
    owner = bool(np.array_equal(tri, m.Fd_master[face, 2:8]))
    unity = float(np.max(abs(np.asarray(Q.sum(1)).ravel() - 1)))
    rng = np.random.default_rng(20261004)
    v = rng.normal(size=(Q.shape[1], 3))
    f = rng.normal(size=(Q.shape[0], 3))
    a = float(np.sum(Q @ v * f))
    b = float(np.sum(v * (Q.T @ f)))
    work = abs(a - b) / max(abs(a), abs(b), 1.0)
    return dict(set_parity=parity, facet_owner_parity=owner, partition_of_unity_error=unity, virtual_work_relative_error=work, algebra_pass=bool(all(parity.values()) and owner and (unity <= p['partition_of_unity_max_error']) and (work <= p['virtual_work_relative_error'])))

def run():
    verify_inputs()
    t0 = time.perf_counter()
    (a, g, m) = model()
    deck = deck_contract()
    (slaves, points, dist, face, tri, bary) = project(m)
    (Q, mapping) = q_matrix(len(m.Vc) + len(m.Vd), slaves, m.nc + tri, p2(bary))
    checks = port_checks(m, slaves, face, tri, Q, p2(bary), deck)
    d0 = dist.max()
    idx = int(np.argmax(dist))
    p = read(SRC / 'PREREG_R1.json')['metrics']
    alltri = m.Vd[m.Fd_master[:, 2:5]]
    point = m.Vc[slaves[idx]]
    allclosest = trimesh.triangles.closest_point(alltri, np.repeat(point[None, :], len(alltri), axis=0))
    alld = np.linalg.norm(allclosest - point, axis=1)
    witness = int(np.argmin(alld))
    exact_error = float(abs(alld.min() - d0))
    nn = trimesh.proximity.ProximityQuery(trimesh.Trimesh(m.Vd, m.Fd_master[:, 2:5], process=False)).vertex(m.Vc[slaves])
    (nd, ni) = nn
    (QN, _) = q_matrix(Q.shape[0], slaves, m.nc + ni[:, None], np.ones((len(slaves), 1)))
    rng = np.random.default_rng(17)
    u = rng.normal(size=(Q.shape[1], 3))
    f = rng.normal(size=(Q.shape[0], 3))
    vn = float(np.sum(QN @ u * f))
    wn = float(np.sum(u * (QN.T @ f)))
    wrong = tri.copy()
    wrong[0] = m.Fd_master[(int(face[0]) + 1) % len(m.Fd_master), 2:8]
    faults = dict(nearest_node={'partition_error': float(np.max(abs(np.asarray(QN.sum(1)).ravel() - 1))), 'virtual_work_relative_error': abs(vn - wn) / max(abs(vn), 1.0), 'port_identity_rejected': bool((QN - Q).nnz > 0), 'projection_max_mm': float(nd.max()), 'frozen_projection_rejected': bool(nd.max() > 0.06)}, wrong_facet_owner={'rejected': not np.array_equal(wrong, m.Fd_master[face, 2:8])}, lost_partition={'rejected': bool(abs(p2(bary)[0].sum() + 0.1 - 1) > p['partition_of_unity_max_error'])})
    native_tolerance = float(deck['tie_header'].split('POSITION TOLERANCE=')[1].split(',')[0])
    projection_pass = bool(d0 <= p['projection_max_mm'] and np.quantile(dist, 0.95) <= p['projection_p95_mm'])
    out = dict(round='R1_NATIVE_PORT', claim_type='capability', case='D1', interface_projection_mm={'max': float(d0), 'p95': float(np.quantile(dist, 0.95)), 'median': float(np.median(dist)), 'min': float(dist.min()), 'n': len(dist), 'resolution': 'PER_POINT'}, counts={'crown_tet4_nodes': len(a['cV']), 'crown_tet10_nodes': len(m.Vc), 'die_tet4_nodes': len(a['dV']), 'die_tet10_nodes': len(m.Vd), 'crown_tetra': len(m.Tc), 'die_tetra': len(m.Td), 'native_slave_faces': len(m.Fc_tied), 'native_master_faces': len(m.Fd_master), 'native_slave_nodes': len(slaves), 'native_fixed_nodes': len(m.fix_nodes)}, native_tolerances_mm={'face_classification': 1.2 * g['h'], 'solver_position_tolerance': native_tolerance, 'frozen_projection_max': 0.06}, checks=checks, gate={'native_sets': all(checks['set_parity'].values()), 'projection_max': projection_pass, 'algebra': checks['algebra_pass'], 'optimization_allowed': bool(projection_pass and checks['algebra_pass'])}, worst_point={'slave_node_0based': int(slaves[idx]), 'coordinate_mm': point, 'closest_point_mm': points[idx], 'master_facet_0based': int(face[idx]), 'master_owner': m.Fd_master[face[idx], :2], 'exhaustive_master_count': len(alltri), 'exhaustive_distance_mm': float(alld.min()), 'closest_search_error_mm': exact_error, 'exhaustive_best_face': witness}, source_material_cards=deck['materials'], source_loads=[{'angle_deg': ang, 'total_force_N': l.sum(0), 'loaded_nodes': int(np.count_nonzero(np.linalg.norm(l, axis=1))), 'resolution': 'PER_POINT'} for (ang, l) in zip([0, 30], deck['loads'])], source_node_rounding_max_mm=float(np.max(abs(deck['nodes'] - np.vstack([m.Vc, m.Vd])))), faults=faults, native_support_area_mm2=float(FE.face_geom(m.Vc, m.Fc_tied)[2].sum()), dropout={'cases': {'considered': 1, 'excluded': 0, 'fraction': 0}, 'slave_nodes': {'considered': len(slaves), 'removed': 0, 'over_006': int((dist > 0.06).sum()), 'over_006_fraction': float((dist > 0.06).mean())}}, cost=cost(t0))
    DATA.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(DATA / 'NATIVE_PORT.npz', slaves=slaves, projection_points=points, distances_mm=dist, master_face=face, master_nodes=tri, barycentric=bary, slave_faces=m.Fc_tied, all_crown_faces=m.Fc_all, master_faces=m.Fd_master, fixed_nodes=m.fix_nodes, cV10=m.Vc, dV10=m.Vd, cT10=m.Tc, dT10=m.Td, loads=np.array(deck['loads']))
    dump(R / 'raw/R1_NATIVE_PORT.json', out)
    state('R1_NATIVE_PORT_MEASURED', out['gate'], 'Write R1 handoff; freeze changed construction after failed projection, or freeze minimax if PASS')
    print(json.dumps(clean({k: out[k] for k in ['interface_projection_mm', 'counts', 'gate', 'worst_point', 'cost']}), indent=2))
if __name__ == '__main__':
    run()
