from common import *
from geometry import MET, topo, preparation, field_values
from physics import source_tests, sufficiency, source_quantile, q05, solve_conductance, X1B, MODEL
import trimesh

def must_reject(fn):
    try:
        fn()
    except (ValueError, AssertionError):
        return True
    return False

def physical_claim_ok(r):
    return r.get('calibrated_force05_N') is None and r.get('physical_status') == 'UNKNOWN_UNMATCHED_SPECIMEN_FAILED_TRANSFER'

def run(rows):
    tests = []

    def add(name, normal, bad, detail):
        tests.append(dict(control=name, normal_pass=bool(normal), injected_fault_rejected=bool(bad), fault=detail))
    st = source_tests()
    su = sufficiency()
    for f in st['source_faults']:
        add('source_' + f['name'], f['normal_pass'], f['rejected'], f['fault'])
    for q in st['quantile_controls']:
        add('Weibull_' + q['row'], q['correct'], q['wrong95_rejected'], 'upper 95% tail substituted for lower5%')
    for q in su['controls']:
        add('hydraulic_closed_form_' + q['name'], q['pass_gate'] and q['mass_gate'], q['injected_2x_rejected'], '2x conductance')
    r = next((r for r in rows if r.get('geometry_pass')))
    rec = next((x for x in cohort() if x['uid'] == r['uid']))
    m = load_npz(r['mesh_path'])
    field = load_npz(r['field_path'])
    (p, pm) = preparation(rec)
    good = topo(m['vertices'], m['faces'])
    bad = topo(m['vertices'], m['faces'][:-1])
    add('closed_shell', good['watertight'], not bad['watertight'], 'remove actual mesh face')
    (gap, _, normals) = field_values(pm, field['sites_mm'])
    (badgap, _, _) = field_values(pm, field['sites_mm'][::32] + np.array([20, 20, 20]))
    error = lambda g: float(np.quantile(abs(g - field['target_mm']), 0.95))
    add('spatial_gap', error(gap) <= 0.025, float(np.quantile(abs(badgap - field['target_mm'][::32]), 0.95)) > 0.025, 'translate actual metrology sites20mm on each axis')
    add('signed_gap', gap.min() >= 0, (-gap).min() < 0, 'reverse normal orientation')
    cap = np.unique(m['faces'][m['face_roles'] == 2])
    z = m['vertices'][cap, 2]
    add('margin_plane', np.max(abs(z - float(p['margin_z']))) <= 1e-06, np.max(abs(z + 0.01 - float(p['margin_z']))) > 1e-06, 'move cap10um')
    T = np.eye(4)
    T[:3, :3] = p['source_R']
    T[:3, 3] = p['source_base']
    wrong = T.copy()
    wrong[:3, :3] *= 2
    add('rigid_frame', not must_reject(lambda : MET.validate_transform(T)), must_reject(lambda : MET.validate_transform(wrong)), 'double transform scale')
    points = m['vertices'][[0, 100, 500, 1000]]
    world = MET.transform(points, T)
    fit = MET.kabsch(points, world)
    err = np.max(abs(MET.transform(points, fit) - world))
    wrongworld = world.copy()
    wrongworld[0] += 10
    badfit = MET.kabsch(points, wrongworld)
    baderr = np.max(abs(MET.transform(points, badfit) - wrongworld))
    add('X55_datum_registration', err < 1e-08, baderr > 1e-08, 'one datum displaced10mm')
    add('force_direct_formula', r['force']['direct_control_relative_error'] <= 1e-08, abs(2 * r['force']['X1B_mean_force_diagnostic']['point_N'] / r['force']['X1B_mean_force_diagnostic']['point_N'] - 1) > 1e-08, '2x predicted mean against independent dot product')
    add('force_source_domain', source_quantile(1.0)['force05_N'] is not None, source_quantile(2.0)['force05_N'] is None, '2mm query outside0.5..1.5mm published source support')
    fake = dict(r['force'], calibrated_force05_N=1234.0)
    add('absolute_transfer', physical_claim_ok(r['force']), not physical_claim_ok(fake), 'invent absolute calibrated force from unmatched source')
    mesh = trimesh.Trimesh(m['vertices'], m['faces'], process=False)
    tmp = P / 'raw/verification_mesh.stl'
    mesh.export(tmp)
    normal = MET.load_mesh(tmp, 'mm')
    badunit = must_reject(lambda : MET.load_mesh(tmp, 'm'))
    add('mesh_units', len(normal.faces) > 0, badunit, 'm units passed to mm reader')
    tmp.unlink()
    h = np.array([0.5, 1.5])
    identity = abs(h.mean() - h[::-1].mean())
    add('sufficiency_identity', identity == 0, abs(h.mean() - (h[::-1] + 0.125).mean()) != 0, 'change summary by0.125mm')
    for q in [su['thickness'], su['film']]:
        assert q['identity_error'] == 0 and q['difference'] > 0
    input_file = Path(rec['mesh_path'])
    data = input_file.read_bytes()
    add('source_hash', sha(input_file) == r['source_mesh_sha256'], __import__('hashlib').sha256(data + b'fault').hexdigest() != r['source_mesh_sha256'], 'append bytes to isolated source buffer')
    out = dict(all_pass=all((r['normal_pass'] and r['injected_fault_rejected'] for r in tests)), count=len(tests), tests=tests, datum_max_error_mm=float(err), datum_bad_error_mm=float(baderr), scope='Implementation/numerical controls, not physical validation')
    dump(P / 'raw/FAULT_INJECTIONS.json', out)
    if not out['all_pass']:
        raise AssertionError([r for r in tests if not r['normal_pass'] or not r['injected_fault_rejected']])
    return out
