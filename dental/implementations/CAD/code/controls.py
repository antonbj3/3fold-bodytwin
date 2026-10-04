import thread_guard
from cadlib import *
from readonly import parents, cone_module, force_refusal, scope
from margin import cycles, line_distance
import trimesh, three_mf, zipfile
from fractions import Fraction
from scipy.optimize import linprog, minimize_scalar
P = parents()
cone = cone_module()

def stepped(r):
    r = np.asarray(r, float)
    base = np.array([[-1, -1], [1, -1], [1, 1], [-1, 1]], float)
    v = np.concatenate([np.c_[base * rr, np.full(4, z)] for (z, rr) in enumerate(r)])
    f = []
    for k in range(len(r) - 1):
        for i in range(4):
            a = 4 * k + i
            b = 4 * k + (i + 1) % 4
            c = b + 4
            d = a + 4
            f.extend([[a, b, c], [a, c, d]])
    f.extend([[12, 13, 14], [12, 14, 15]])
    return trimesh.Trimesh(v, np.array(f), process=False)

def run():
    checks = []

    def record(name, value, detail=None):
        checks.append(dict(name=name, pass_=bool(value), detail=detail))
        if not value:
            raise AssertionError(name)
    r1 = np.array([3.0, 2.5, 2.0, 1.5])
    r2 = np.array([3.0, 2.0, 2.5, 1.5])
    summary = lambda r: np.array([r[0], r[-1], r.mean(), 3.0, (r[0] - r[-1]) / 3])
    A = stepped(r1)
    B = stepped(r2)
    ca = cone.classify(A.face_normals)
    cb = cone.classify(B.face_normals)
    identity = float(np.max(abs(summary(r1) - summary(r2))))
    record('exact_summary_identity', identity == 0)
    record('downstream_insertion_changes', ca['status'] == 'YES_STRICT' and cb['status'] == 'NO_NONZERO_DIRECTION', dict(A=ca['status'], B=cb['status']))
    record('Farkas_all6_verified', len(cb['certificates']) == 6 and cone.recheck(cb))
    bad = json.loads(json.dumps(cb))
    bad['certificates'][0]['y'][0] = str(Fraction(bad['certificates'][0]['y'][0]) + 1)
    record('injected_Farkas_multiplier_rejected', not cone.recheck(bad))
    d = np.array(ca['unit_direction'])
    record('injected_reversed_insertion_rejected', bool(np.min(A.face_normals @ -d) < 0))
    lp = []
    for m in [A, B]:
        feasible = []
        for j in range(3):
            for s in [-1, 1]:
                bounds = [(-1, 1)] * 3
                bounds[j] = (s, s)
                q = linprog([0, 0, 0], A_ub=-m.face_normals, b_ub=np.zeros(len(m.faces)), bounds=bounds, options={'threads': 1})
                feasible.append(bool(q.success))
        lp.append(feasible)
    record('same_information_HiGHS_control', any(lp[0]) and (not any(lp[1])))
    suff = dict(summary_names=['basal_halfwidth_mm', 'occlusal_halfwidth_mm', 'mean_halfwidth_mm', 'height_mm', 'endpoint_slope'], A_summary=summary(r1), B_summary=summary(r2), identity_error=identity, A_radii_mm=r1, B_radii_mm=r2, downstream=dict(A_status=ca['status'], B_status=cb['status'], insertion_feasibility_difference=1, waist_rebound_difference_mm=0.5), minimum_extension='For this witness, maximum positive successive radius increment distinguishes states; for arbitrary3D preps retain oriented facet-normal cone and surface incidence, not average taper.', resolution='PER_SURFACE_REGION->PER_TOOTH', external_referent={'kind': 'closed_form', 'locator': 'Nested square cross sections imply straight extraction; symmetric waist rebound excludes every nonzero normal-cone direction', 'compared_quantity': 'zero-clearance extraction', 'refutes_us': False}, certificates=[ca, cb], LP_control=lp)
    dump(ROOT / 'raw/SUFFICIENCY.json', suff)
    v = np.c_[np.cos(np.arange(8) * np.pi / 4), np.sin(np.arange(8) * np.pi / 4), np.zeros(8)]
    e = np.c_[np.arange(8), np.roll(np.arange(8), -1)]
    (cs, _) = cycles(e, v)
    record('closed_cycle_baseline', len(cs) == 1)
    record('injected_open_curve_rejected', len(cycles(e[:-1], v)[0]) == 0)
    record('injected_branch_rejected', len(cycles(np.r_[e, [[0, 4]]], v)[0]) == 0)
    curve = np.array([[0.0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]])
    points = np.array([[0.3, 0.2, 0], [2, 0.5, 0.1]])
    own = line_distance(points, curve)
    other = []
    for p in points:
        dd = []
        for (a, b) in zip(curve, np.roll(curve, -1, axis=0)):
            sol = minimize_scalar(lambda t: np.sum((a + t * (b - a) - p) ** 2), bounds=(0, 1), method='bounded')
            dd.append(min(sol.fun, np.sum((a - p) ** 2), np.sum((b - p) ** 2)))
        other.append(np.sqrt(min(dd)))
    record('independent_curve_distance', np.max(abs(own - other)) < 1e-08)
    record('injected_curve_error_rejected', not np.max(abs(1000 * own - other)) < 1e-08)
    tri = np.array([[0.0, 0, 0], [1, 0, 0], [0, 1, 0]])
    f = np.array([[0, 1, 2]])
    record('parallel_wall_analytic_positive', P['cap'].wall_check(tri + [0, 0, 2], f, tri, f, 1)['status'] == 'PASS')
    record('injected_thin_wall_rejected', P['cap'].wall_check(tri + [0, 0, 0.1], f, tri, f, 1)['status'] == 'FAIL')
    row = next((r for r in read(ROOT / 'raw/R4B_PREDICTIONS.json') if r['status'] == 'GENERATED'))
    a = np.load(row['mesh_path'])
    inp = next((r for r in read(ROOT / 'INPUT_LOCK_R3.json')['rows'] if r['key'] == row['key']))
    from pipeline_r4 import canonical_prep
    (prep, ser) = canonical_prep(inp)
    pp = P['cap'].cap_planes(prep.vertices, prep.faces[np.array(row['cap']['face_indices'])])
    from gencad_bench.checks import cement_max
    record('matched_max_film_witness', cement_max.verify_positive(pp, a['iv'].tolist(), prep.vertices[:len(a['iv'])].tolist(), 0.12))
    record('injected_max_film_limit_rejected', not cement_max.verify_positive(pp, a['iv'].tolist(), prep.vertices[:len(a['iv'])].tolist(), 0.001))
    record('old_mismatched_witness_rejected', not cement_max.verify_positive(pp, a['iv'].tolist(), prep.vertices.tolist(), 0.12))
    xy = np.array([[0.0, 0], [1, 0], [1, 1], [0, 1]])
    ff = np.array([[0, 1, 2], [0, 2, 3]])
    g = np.array([0.05, 0.05, 0.05, 0.05])
    cm = P['v6'].compare(xy, ff, g, g)
    q = P['v6'].Q.contact_map(xy, ff, g, g)
    pr = read(ROOT / 'PREREG_R3.json')['metrics']['V6']
    record('same_information_V6_contact_control', cm['symdiff_mm2'] == q['contact_symdiff_mm2'] == 0)
    record('injected_wrong_contact_observable_rejected', not P['v6'].gates(dict(cm, symdiff_mm2=10), pr)['pattern'])
    record('force_missing_calibration_refused', force_refusal('control')['status'] == 'ABSTAIN')
    x13 = module('x13_port', BASE / 'LANE_X13_CEMENT_GAP/cement_port.py')
    xr = x13.hydraulic_input_from_group_mean(dict(row_id='control', measured_mean_um=80, state='nominal_CAD_not_measured'))
    record('mean_film_not_promoted_to_physical_field', xr['status'] == 'UNKNOWN_MISSING_SPATIAL_FILM')
    (v, f, roles) = three_mf.readback(row['three_mf_path'])
    record('3MF_exact_roundtrip', np.array_equal(v, a['vertices']) and np.array_equal(f, a['faces']) and np.array_equal(roles, a['roles']))
    scene = trimesh.load(row['three_mf_path'])
    mm = next(iter(scene.geometry.values())) if isinstance(scene, trimesh.Scene) else scene
    record('independent_3MF_reader', len(mm.faces) == len(f) and abs(mm.volume - trimesh.Trimesh(v, f, process=False).volume) < 1e-09)
    for (name, old, new) in [('unit', 'unit="millimeter"', 'unit="meter"'), ('region', 'p1="0" p2="0" p3="0"', 'p1="7" p2="7" p3="7"')]:
        path = DATA / ('control_bad_' + name + '.3mf')
        with zipfile.ZipFile(row['three_mf_path']) as z, zipfile.ZipFile(path, 'w') as w:
            for k in z.namelist():
                data = z.read(k)
                if k.endswith('.model'):
                    txt = data.decode()
                    assert old in txt
                    data = txt.replace(old, new, 1).encode()
                w.writestr(k, data)
        try:
            three_mf.readback(path)
            reject = False
        except ValueError:
            reject = True
        record('injected_3MF_' + name + '_rejected', reject)
    with scope(BASE / 'LANE_X89_MILLING_CODESIGN/code', ['codesign', 'collision']):
        mg = module('codesign', BASE / 'LANE_X89_MILLING_CODESIGN/code/codesign.py')
        c = read(BASE / 'LANE_X89_MILLING_CODESIGN/contracts/crownloop_D1.json')
        tool = read(c['tool_card'])['tools'][-1]
        res = mg.mechanical(tool, dict(c['scenario'], force_N=1000000000.0))
        record('injected_milling_force_rejected', res['mechanical_status'] == 'SCENARIO_INFEASIBLE')
    record('runtime_LP_after_cone_works', linprog([1.0], bounds=[(1, None)], options={'threads': 1}).success)
    dump(ROOT / 'raw/CONTROLS.json', dict(checks=checks, all_pass=all((x['pass_'] for x in checks)), count=len(checks), resolution='PER_POINT/PER_TOOTH', physical_validation=False))
    print('CONTROLS', len(checks), 'PASS')
if __name__ == '__main__':
    run()
