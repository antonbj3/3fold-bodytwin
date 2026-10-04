import json, math, time
import numpy as np
from common import ROOT, CORPUS, save, state, sha, utc, tables, article, text
from patient import load, export_section
from beam import moments, integrated_moments, stiffness, condense, direct, affine_box

def source_cells():
    number = lambda s: float(s.replace('−', '-'))
    src = tables('PMC10824581')
    out = []
    for (table, kind) in [('tab4', 'linear_mm'), ('tab5', 'angular_deg')]:
        rows = src[table]
        start = next((i for (i, r) in enumerate(rows) if r[0] == 'Mandible'))
        for r in rows[start + 1:]:
            typ = int(r[0])
            assert len(r) == 10
            names = ['mesiodistal', 'occlusogingival', 'buccolingual'] if kind == 'linear_mm' else ['torque', 'rotation', 'angulation']
            for (j, name) in enumerate(names):
                out.append(dict(table=table, tooth_type=typ, quantity=name, unit='mm' if kind == 'linear_mm' else 'degree', mean=number(r[1 + 3 * j]), sd=number(r[2 + 3 * j]), threshold=0.5 if kind == 'linear_mm' else 2.0, resolution='PER_TOOTH', population='published tooth-type group; not local patient'))
    save('raw/EXTERNAL_SOURCE_TABLES.json', dict(PMC10824581=src, PMC10548050=tables('PMC10548050'), PMC10587025=tables('PMC10587025')))
    return out

def moment_certificate(cells):
    six = []
    for typ in [3, 2, 1, 1, 2, 3]:
        row = []
        for c in cells:
            if c['tooth_type'] != typ:
                continue
            h = c['threshold']
            mu = c['mean']
            sd = c['sd']
            upper_m2 = (abs(mu) + 0.005) ** 2 + (sd + 0.005) ** 2
            row.append(dict(quantity=c['quantity'], source_mean=mu, source_sd=sd, tail_upper=min(1.0, upper_m2 / h ** 2), mean_pass=abs(mu) <= h))
        six.append(dict(tooth_type=typ, axes=row, joint_success_lower=max(0.0, 1 - sum((c['tail_upper'] for c in row)))))
    assembly = max(0.0, 1 - sum((c['tail_upper'] for r in six for c in r['axes'])))
    c = next((c for c in cells if c['tooth_type'] == 1 and c['quantity'] == 'occlusogingival'))
    a = -0.5001
    mu = c['mean']
    var = c['sd'] ** 2
    p = var / (var + (mu - a) ** 2)
    b = mu + var / (mu - a)
    mean = p * a + (1 - p) * b
    variance = p * (a - mean) ** 2 + (1 - p) * (b - mean) ** 2
    return dict(per_bracket=six, assembly_success_lower=assembly, release_threshold=0.95, decision='RELEASE' if assembly >= 0.95 else 'NOT_CERTIFIED_FROM_SUMMARIES', mean_based_comparator_pass=all((a['mean_pass'] for r in six for a in r['axes'])), moment_condition='conditional on printed empirical moments; sample size and pairing differ after bond losses; future patient moments UNKNOWN', counterexample=dict(source_quantity=c, states=[a, b], probabilities=[p, 1 - p], reproduced_mean=mean, reproduced_variance=variance, mean_error=mean - mu, variance_error=variance - var, possible_failure_probability=p, status='constructed moment-compatible distribution, NOT independently measured failures'))

def rodrigues(axis, angle):
    axis = np.array(axis)
    axis = axis / np.linalg.norm(axis)
    (x, y, z) = axis
    W = np.array([[0, -z, y], [z, 0, -x], [-y, x, 0]])
    return np.eye(3) + np.sin(angle) * W + (1 - np.cos(angle)) * (W @ W)

def run():
    start = time.perf_counter()
    (points, vertical, rows, meshes) = load()
    state('R1_RUNNING', 'six surface anchors constructed', 'external beam forces and summary sufficiency')
    src = article('PMC10548050')
    excerpt = next((text(p) for p in src.findall('.//body//p') if '1.27' in text(p) and '0.95' in text(p)))
    assert '0.06' in excerpt and '0.02' in excerpt
    (A, Iy, Iz, cz) = moments()
    ctrl = integrated_moments()
    tests = []
    for (direction, I, delta, measured) in [('AP', Iy, 0.06, 1.27), ('vertical', Iz, 0.02, 0.95)]:
        pred = 48 * 4800 * I / 10 ** 3 * delta
        ci = ctrl[1] if direction == 'AP' else ctrl[2]
        control = 48 * 4800 * ci / 10 ** 3 * delta
        tests.append(dict(direction=direction, predicted_N=pred, observed_N=measured, displacement_mm=delta, control_N=control, relative_error=abs(pred - measured) / measured, gate_pass=abs(pred - measured) / measured <= 0.15, resolution='PHENOMENOLOGICAL', external_measurement_resolution='PER_SURFACE_REGION: one rod fixture'))
    (K, lengths) = stiffness(points, vertical)
    (S, ui, ri) = condense(K)
    delta = 1 / 64
    ua = np.repeat((delta * vertical)[None, :], 6, axis=0)
    ub = ua * np.array([1, -1, 1, -1, 1, -1])[:, None]
    summary = lambda u: float(np.sqrt(np.mean(np.sum(u * u, axis=1))))
    (sa, sb) = (summary(ua), summary(ub))
    assert sa == sb
    fa = (S @ ua.ravel()).reshape(-1, 3)
    fb = (S @ ub.ravel()).reshape(-1, 3)
    (fb_ctrl, q) = direct(K, ub)
    pair = dict(states_displacement_mm=[ua.tolist(), ub.tolist()], summary_name='RMS bond displacement', summaries_mm=[sa, sb], summary_hex=[sa.hex(), sb.hex()], identity_error_mm=abs(sa - sb), per_bond_absolute_magnitude_identity_error=float(np.max(np.abs(np.linalg.norm(ua, axis=1) - np.linalg.norm(ub, axis=1)))), peak_reactions_N=[float(np.max(np.linalg.norm(fa, axis=1))), float(np.max(np.linalg.norm(fb, axis=1)))], downstream_difference_N=float(np.max(np.linalg.norm(fb, axis=1)) - np.max(np.linalg.norm(fa, axis=1))), minimal_extension='signed relative displacement field + directional stiffness + bonding constraints', witness_kind='our_own_fixture on published geometry; not empirical force truth', resolution='PER_TOOTH')
    omega = np.array([0.003, -0.002, 0.001])
    ur = np.cross(np.broadcast_to(omega, points.shape), points)
    qr = np.zeros(36)
    qr[ui] = ur.ravel()
    qr[ri] = np.tile(omega, 6)
    numeric = dict(direct_control_max_abs_N=float(np.max(np.abs(fb - fb_ctrl))), net_force_N=float(np.linalg.norm(fb.sum(axis=0))), net_moment_Nmm=float(np.linalg.norm(np.cross(points, fb).sum(axis=0))), rigid_translation_max_N=float(np.max(np.abs(fa))), rigid_rotation_max_generalized_force=float(np.max(np.abs(K @ qr))), section_control_max_relative=float(np.max(np.abs(np.array(ctrl) - np.array(moments())) / np.array(moments()))), condensed_symmetry_error=float(np.max(np.abs(S - S.T))))
    av = np.array([-0.25, 0.25])
    ap = np.array([0.5, 0.5])
    bv = np.array([-1.0, 0.0, 1.0])
    bp = np.array([1 / 32, 30 / 32, 1 / 32])
    am = float(ap @ av)
    bm = float(bp @ bv)
    avariance = float(ap @ (av - am) ** 2)
    bvariance = float(bp @ (bv - bm) ** 2)
    assert am == bm and avariance == bvariance
    tail_pair = dict(summary=['signed mean', 'variance'], summary_A=[am, avariance], summary_B=[bm, bvariance], identity_error=0.0, success_within_half_mm=[float(ap @ (np.abs(av) <= 0.5)), float(bp @ (np.abs(bv) <= 0.5))], downstream_success_fraction_difference=0.0625, values=[av.tolist(), bv.tolist()], weights=[ap.tolist(), bp.tolist()], minimal_extension='actual exceedance indicator/quantile at specified manufacturing tolerance; moments alone do not suffice', witness_kind='our_own_fixture; not measured bracket distribution')
    angle = np.deg2rad(4)
    i = 2
    axis = rows[i]['mesial_axis']
    c = np.array(rows[i]['crown_center_mm'])
    R = rodrigues(axis, angle)
    v = meshes[i].vertices
    changed = (v - c) @ R.T + c
    disp = float(np.max(np.linalg.norm(changed - v, axis=1)))
    radius = float(np.max(np.linalg.norm(v - c, axis=1)))
    bound = 2 * radius * np.sin(angle / 2)
    pose_pair = dict(bracket_center_translation_states_mm=[np.zeros((6, 3)).tolist(), np.zeros((6, 3)).tolist()], center_summary_identity_error_mm=0.0, orientation_states_degrees=[0.0, 4.0], angular_acceptance=[True, False], max_tooth_point_displacement_mm=[0.0, disp], nonlinear_displacement_bound_mm=float(np.nextafter(bound, np.inf)), minimal_extension='relative SO(3) orientation at each bond; six-axis per-tooth pose', witness_kind='our_own_fixture on independent mesh; clinical tolerance from external source', point_domain='full tooth including root, not a measured achieved clinical movement', resolution='PER_POINT')
    cells = source_cells()
    cert = moment_certificate(cells)
    source_manifest = [dict(path=str(CORPUS / (x + '.xml')), sha256=sha(CORPUS / (x + '.xml'))) for x in ['PMC10548050', 'PMC10824581', 'PMC10587025']]
    save('raw/SOURCE_MANIFEST.json', source_manifest)
    export = export_section(points, vertical, 'exports/nominal_retainer_metrology.stl')
    save('raw/SUFFICIENCY_R1.json', dict(O18=pair, O20_pose=pose_pair, O20_moments=tail_pair))
    save('raw/O20_SOURCE_CELLS.json', cells)
    save('raw/O18_FORCE_REFERENCE.json', dict(source='PMC10548050 Discussion; E4H nominal', excerpt=excerpt, caveat='reference values known before prereg; mechanistic consistency test, not blind validation', tests=tests))
    save('raw/R1.json', dict(O18=dict(claim_type='capability', external_gate='PASS' if all((t['gate_pass'] for t in tests)) else 'FAIL', tests=tests, sufficiency=pair, numerical=numeric, lengths_mm=lengths.tolist(), beam_physical_validity='UNKNOWN on patient; external fixture tested', external_referent=json.loads((ROOT / 'PREREG_O18_R1.json').read_text())['external_referent']), O20=dict(claim_type='information_link', certificate=cert, pose_sufficiency=pose_pair, moment_sufficiency=tail_pair, external_referent=json.loads((ROOT / 'PREREG_O20_R1.json').read_text())['external_referent']), export=export, wall_seconds=time.perf_counter() - start, finished_utc=utc()))
    state('R1_COMPLETE', 'O18 external beam closure ' + ('PASS' if all((t['gate_pass'] for t in tests)) else 'FAIL') + '; O20 ' + cert['decision'], 'write handoff; change stiffness representation to empirical port and acquire full pose rather than moments')
    print(json.dumps(dict(retainer_external=tests, bracket_decision=cert['decision'], assembly_lower=cert['assembly_success_lower'], sufficiency_downstream_N=pair['downstream_difference_N'], numerical=numeric), indent=2))
if __name__ == '__main__':
    run()
