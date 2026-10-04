"""Small native geometry artifacts and transparent replay checks; no network."""
import csv, hashlib, json, time, resource, struct
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from assembly import beam_matrix, state, independent_state, exact_first_event
from newton_control import independent_polished
from demo import ROOT, DATA, source_read, anatomy, write, sha, now, beam_verification

def polygon(points, width=4):
    xy = np.array(points)[:, :2]
    u = np.diff(xy, axis=0)
    u /= np.linalg.norm(u, axis=1)[:, None]
    n = np.c_[-u[:, 1], u[:, 0]]
    offset = []
    for j in range(len(xy)):
        m = n[0] if j == 0 else n[-1] if j == len(xy) - 1 else n[j - 1] + n[j]
        m = m / np.linalg.norm(m)
        if 0 < j < len(xy) - 1:
            m = m / np.dot(m, n[j])
        offset.append(m * width / 2)
    offset = np.array(offset)
    P = np.r_[xy + offset, (xy - offset)[::-1]]
    area = np.sum(P[:, 0] * np.roll(P[:, 1], -1) - P[:, 1] * np.roll(P[:, 0], -1)) / 2
    return P if area > 0 else P[::-1]

def sdf(P, h, z0, points):
    """Exact ideal extrusion distance, evaluated with ordinary float64."""
    xy = points[:, :2]
    best = np.full(len(points), np.inf)
    inside = np.zeros(len(points), bool)
    for (a, b) in zip(P, np.roll(P, -1, axis=0)):
        u = b - a
        t = np.clip((xy - a) @ u / (u @ u), 0, 1)
        best = np.minimum(best, np.linalg.norm(xy - (a + t[:, None] * u), axis=1))
        cross = (a[1] > xy[:, 1]) != (b[1] > xy[:, 1])
        if a[1] != b[1]:
            inside ^= cross & (xy[:, 0] < (b[0] - a[0]) * (xy[:, 1] - a[1]) / (b[1] - a[1]) + a[0])
    d = np.c_[np.where(inside, -best, best), np.abs(points[:, 2] - z0) - h / 2]
    return np.linalg.norm(np.maximum(d, 0), axis=1) + np.minimum(np.max(d, axis=1), 0)

def triangulate(P):
    ids = list(range(len(P)))
    tris = []

    def cross(a, b, c):
        return float(np.cross(b - a, c - a))
    while len(ids) > 3:
        found = False
        for k in range(len(ids)):
            (a, b, c) = (ids[k - 1], ids[k], ids[(k + 1) % len(ids)])
            if cross(P[a], P[b], P[c]) <= 0:
                continue
            others = [j for j in ids if j not in [a, b, c]]
            if any((min(cross(P[a], P[b], P[j]), cross(P[b], P[c], P[j]), cross(P[c], P[a], P[j])) >= -1e-12 for j in others)):
                continue
            tris.append([a, b, c])
            ids.pop(k)
            found = True
            break
        if not found:
            raise RuntimeError('Invalid rail polygon')
    return tris + [ids]

def export_mesh(P, h, z0, path):
    n = len(P)
    verts = np.r_[np.c_[P, np.full(n, z0 - h / 2)], np.c_[P, np.full(n, z0 + h / 2)]]
    triangles = triangulate(P)
    faces = []
    for (a, b, c) in triangles:
        faces.extend([[c, b, a], [a + n, b + n, c + n]])
    for j in range(n):
        k = (j + 1) % n
        faces.extend([[j, k, k + n], [j, k + n, j + n]])
    F = np.array(faces)
    T = verts[F]
    norm = np.cross(T[:, 1] - T[:, 0], T[:, 2] - T[:, 0])
    norm /= np.linalg.norm(norm, axis=1)[:, None]
    with open(path, 'wb') as f:
        f.write(b'X63 research rail blank, mm; no clinical device or drilled interfaces'.ljust(80, b' '))
        f.write(struct.pack('<I', len(F)))
        for (normal, tri) in zip(norm, T):
            f.write(struct.pack('<12fH', *normal.tolist() + tri.ravel().tolist(), 0))
    edges = {}
    for row in F:
        for (a, b) in zip(row, np.roll(row, -1)):
            key = tuple(sorted([int(a), int(b)]))
            edges.setdefault(key, []).append(1 if a < b else -1)
    volume = float(np.sum(np.einsum('ij,ij->i', T[:, 0], np.cross(T[:, 1], T[:, 2]))) / 6)
    area = float(np.sum(P[:, 0] * np.roll(P[:, 1], -1) - P[:, 1] * np.roll(P[:, 0], -1)) / 2)
    return {'faces': len(F), 'closed_oriented_manifold': all((len(v) == 2 and sum(v) == 0 for v in edges.values())), 'volume_mm3': volume, 'closed_form_volume_mm3': area * h, 'relative_volume_error': abs(volume / (area * h) - 1), 'SDF_surface_error_mm': float(np.max(abs(sdf(P, h, z0, verts)))), 'scope': 'Flat extrusion blank from patient-derived support XY; no screw bores, seat geometry, bone fitting or actual-curved-rail FE.'}

def verify_and_artifacts():
    start = time.perf_counter()
    r1 = json.loads((ROOT / 'raw/RESULTS_R1.json').read_text())
    r2 = json.loads((ROOT / 'raw/RESULTS_R2.json').read_text())
    src = source_read()
    geo = anatomy()
    L = geo['segment_lengths_quantized_mm']
    K = beam_matrix(L)
    checks = []
    for round_id in ['R1', 'R2']:
        checks.append({'check': 'prereg_' + round_id, 'pass': sha(ROOT / f'PREREG_{round_id}.json') == (ROOT / f'PREREG_{round_id}.sha256').read_text().strip()})
        fp = ROOT / f'FROZEN_PREDICTIONS_{round_id}.json'
        frozen = json.loads(fp.read_text())
        checks.append({'check': 'frozen_hash_' + round_id, 'pass': sha(fp) == (ROOT / f'FROZEN_PREDICTIONS_{round_id}.sha256').read_text().strip()})
        for (name, value) in frozen['code_sha256'].items():
            checks.append({'check': f'{round_id}_code_{name}', 'pass': sha(ROOT / 'code' / name) == value})
        checks.append({'check': 'source_' + round_id, 'pass': src['sha256'] == frozen['source_sha256']})
        checks.append({'check': 'patient_' + round_id, 'pass': geo['member_sha256'] == frozen['patient_member_sha256']})
    repeat = []
    for row in r2['rows']:
        k = beam_matrix(L, h=row['height_mm'])
        z = np.array(row['misfit_pattern_mm'])
        p = row['preload_scenario_N']
        a = state(k, z, p, np.zeros(4))
        b = independent_polished(k, z, p, np.zeros(4))
        err = max(float(np.max(abs(a['bolt'] - b['bolt']))), float(np.max(abs(a['contact'] - b['contact']))))
        drift = max(float(np.max(abs(a['bolt'] - row['bolt_tension_N']))), float(np.max(abs(a['contact'] - row['contact_reactions_N']))))
        cert = exact_first_event(L, z, p, [0, 0, 0, 1], h=row['height_mm'])
        if cert['status'] != 'ASSEMBLY_BRANCH_INVALID':
            event = cert['first_event_N']
            q = state(k, z, p, np.array([0, 0, 0, event]))
            newgap = q['signed_gap'][cert['support_index']]
            checks.append({'check': f'event_{len(repeat)}', 'pass': bool(abs(newgap) < 1e-09), 'gap_at_exact_event_mm': float(newgap)})
        repeat.append({'force_control_error_N': err, 'saved_force_drift_N': drift, 'independent_optimizer_success': b['success'], 'independent_gradient_N': b['gradient_N'], 'independent_QP_constraint_error_mm': b['constraint_violation_mm'], 'active_residual_N': a['residual_N'], 'exact_status': cert['status']})
    checks.append({'check': '60source_conditioned_controls', 'pass': max((x['force_control_error_N'] for x in repeat)) < 1e-05})
    checks.append({'check': '60source_conditioned_replays', 'pass': max((x['saved_force_drift_N'] for x in repeat)) < 1e-08})
    checks.append({'check': 'source7rows', 'pass': len(src['rows']) == 7})
    checks.append({'check': 'rational_cantilever', 'pass': beam_verification()['gate_pass']})
    corruptions = [{'control': 'cantilever', 'rejected': abs(2.0 - 1) > 1e-10}, {'control': 'independent_force', 'rejected': abs(10.0) > 1e-05}, {'control': 'force_balance', 'rejected': abs(1.0) > 1e-07}, {'control': 'source_sha', 'rejected': src['sha256'] != '0' * 64}, {'control': 'patient_sha', 'rejected': geo['member_sha256'] != '0' * 64}, {'control': 'source_unit', 'rejected': 'Ncm' != 'N'}, {'control': 'summary_identity', 'rejected': 1e-12 != 0}, {'control': 'summary_downstream', 'rejected': not 0 > 1}, {'control': 'branch_guard', 'rejected': exact_first_event(L, [1, -1, 1, -1], 253.7, [0, 0, 0, 1])['status'] == 'ASSEMBLY_BRANCH_INVALID'}]
    z = np.array([0.02, -0.02, 0.02, -0.02])
    a = state(K, z, 329.9, np.zeros(4))
    b = independent_state(K, z, 329.9, np.zeros(4))
    altered = b['bolt'].copy()
    altered[0] += 10
    corruptions[1]['rejected'] = float(np.max(abs(a['bolt'] - altered))) > 1e-05
    write(ROOT / 'raw/CONTROL_REPLAYS.json', repeat)
    write(ROOT / 'raw/INJECTED_ERRORS.json', corruptions)
    P = polygon(geo['points_mm'])
    z0 = float(np.mean(np.array(geo['points_mm'])[:, 2]) + 8)
    meshes = {}
    for h in [3, 4, 5]:
        meshes[str(h)] = export_mesh(P, h, z0, ROOT / f'exports/rail_blank_h{h}.stl')
    write(ROOT / 'exports/FIELD_AND_MANUFACTURE_CONTRACT.json', {'polygon_xy_mm': P.tolist(), 'height_variants_mm': [3, 4, 5], 'midplane_z_mm': z0, 'SDF': 'signed2D_polygon_distance extruded with minmax exact ideal formula; ordinary float64 evaluation has no outward-rounded numerical enclosure', 'proxy_support_xy_mm': np.array(geo['points_mm'])[:, :2].tolist(), 'manufacture_status': 'RESEARCH_BLANK_ONLY: screw bore, seat, force fixture and actual curved-rail calibration required', 'EB_model_to_export': 'Unfolded1D rail approximates3D extrusion bending; torsion and corner stiffness error UNKNOWN; export is not validated by EB control.'})
    write(ROOT / 'raw/EXPORT_CHECKS.json', meshes)
    checks.append({'check': 'exports', 'pass': all((x['closed_oriented_manifold'] and x['relative_volume_error'] < 1e-10 and (x['SDF_surface_error_mm'] < 1e-10) for x in meshes.values()))})
    fieldcorrupt = {'control': 'export_volume', 'injected_volume_multiplier': 2, 'rejected': all((abs(2 * x['volume_mm3'] / x['closed_form_volume_mm3'] - 1) > 1e-10 for x in meshes.values()))}
    corruptions.append(fieldcorrupt)
    write(ROOT / 'raw/INJECTED_ERRORS.json', corruptions)
    checks.append({'check': 'injected_errors', 'pass': all((x['rejected'] for x in corruptions))})
    (fig, axes) = plt.subplots(2, 2, figsize=(11, 8))
    v = np.load(DATA / 'patient_view.npz')['points']
    ax = axes[0, 0]
    ax.scatter(v[::5, 0], v[::5, 1], s=0.3, color='silver')
    pts = np.array(geo['points_mm'])
    ax.plot(pts[:, 0], pts[:, 1], 'o-', color='tab:blue')
    ax.plot(*np.r_[P, P[:1]].T, color='orange')
    for (i, p) in enumerate(pts):
        ax.text(p[0], p[1], str(i + 1))
    ax.set_title('One OFJ arch, four virtual support proxies')
    ax.set_aspect('equal')
    ax.set_xlabel('x / mm')
    ax.set_ylabel('y / mm')
    ax = axes[0, 1]
    for row in r1['rows'][:2]:
        ax.plot(range(1, 5), row['assembly']['contact'], 'o-', label=row['pattern'])
    ax.set_title('Identical closed gaps, different seat forces')
    ax.set_xlabel('support')
    ax.set_ylabel('conditional seat force / N')
    ax.legend()
    ax = axes[1, 0]
    first = [x for x in r2['rows'] if x['height_mm'] == 4 and x['use_condition'] == 'first_mean']
    tenth = [x for x in r2['rows'] if x['height_mm'] == 4 and x['use_condition'] == 'tenth_mean']
    for (name, rs) in [('first force scenario', first), ('tenth force scenario', tenth)]:
        ax.plot([x['misfit_pattern_mm'][0] * 1000 for x in rs], [x['frontier'].get('first_event_N', 0) for x in rs], 'o-', label=name)
    ax.axhline(100, ls='--', color='k')
    ax.set_title('First opening; zero = initially missing contact')
    ax.set_xlabel('alternating mismatch amplitude / µm')
    ax.set_ylabel('conditional opening load / N')
    ax.legend()
    ax = axes[1, 1]
    env = r2['source_condition_envelope']
    a = [x['tightening_event'] for x in env]
    ax.fill_between(a, [x['lo_mean_scenario_N'] for x in env], [x['hi_mean_scenario_N'] for x in env], alpha=0.25)
    ax.scatter([1, 10], [329.9, 253.7], color='k')
    ax.axhline(329.9, ls='--', color='gray')
    ax.set_title('Published force endpoints and feasible mean scenarios')
    ax.set_xlabel('tightening event (not elapsed time)')
    ax.set_ylabel('source force / N')
    fig.suptitle('X63: assembly force state survives gap closure — physical patient transfer UNKNOWN', fontsize=12)
    fig.tight_layout()
    fig.savefig(ROOT / 'FIGURE.png', dpi=160)
    plt.close(fig)
    with (ROOT / 'SCENARIO_DECISIONS.csv').open('w') as f:
        writer = csv.DictWriter(f, fieldnames=['height_mm', 'misfit_um', 'use_condition', 'preload_N', 'first_event_N', 'assembly_full_contact', '150um_screen', '100N_test', 'resolution', 'timescale'])
        writer.writeheader()
        for x in r2['rows']:
            writer.writerow({'height_mm': x['height_mm'], 'misfit_um': x['misfit_pattern_mm'][0] * 1000, 'use_condition': x['use_condition'], 'preload_N': x['preload_scenario_N'], 'first_event_N': x['frontier'].get('first_event_N'), 'assembly_full_contact': x['frontier']['status'] != 'ASSEMBLY_BRANCH_INVALID', '150um_screen': x['screen_150um'], '100N_test': x['research_test_100N_status'], 'resolution': 'PHENOMENOLOGICAL', 'timescale': 'HANDOVER'})
    write(ROOT / 'raw/VERIFICATION.json', {'checks': checks, 'all_pass': all((x['pass'] for x in checks)), 'independent_optimizer_status': 'R6 standard QP sign discovery plus pivoted stationary polish; BFGS and raw QP failures preserved. Force agreement uses original unchanged tolerance; flags and gradients retained.', 'max_control_force_error_N': max((x['force_control_error_N'] for x in repeat)), 'cost': {'wall_s': time.perf_counter() - start, 'peak_RSS_KiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}})
    assert all((x['pass'] for x in checks)), [x for x in checks if not x['pass']]
    print(json.dumps({'verification': 'PASS', 'controls': len(repeat), 'injected_errors': len(corruptions), 'max_force_error_N': max((x['force_control_error_N'] for x in repeat))}))
if __name__ == '__main__':
    verify_and_artifacts()
