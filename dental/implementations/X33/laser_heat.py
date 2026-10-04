"""Conservative pulsed heat probe. Fixed voxel geometry; no ablation closure fitted."""
from pathlib import Path
import hashlib, json, math, time, resource, sys, os
import numpy as np
from scipy import ndimage as ndi
from scipy.spatial import cKDTree
from scipy.integrate import quad
from scipy.optimize import lsq_linear
import numba as nb
ROOT = Path(os.environ.get('X33_RUN_ROOT', str(Path(__file__).resolve().parent)))
sys.path.insert(0, str(ROOT / 'inputs'))
from decidability import ScalarPort, interval_decision, refinement_sequence

def dump(name, obj):
    (ROOT / name).write_text(json.dumps(obj, indent=2, ensure_ascii=False, allow_nan=False) + '\n')

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

@nb.njit
def conductances(neighbors, k, h):
    n = len(k)
    g = np.zeros((n, 6))
    for i in range(n):
        for d in range(6):
            j = neighbors[i, d]
            if j >= 0:
                g[i, d] = 2 * k[i] * k[j] / (k[i] + k[j]) * h
    return g

def geometry(subdivision, pr):
    q = np.load(pr['sources']['geometry']['path'])
    tooth = q['tooth']
    pulp = q['pulp']
    prepared = tooth & (ndi.distance_transform_edt(tooth, sampling=0.3) >= 0.65)
    assert not np.any(pulp & ~prepared)
    s = subdivision
    mask = prepared.copy()
    pm = pulp.copy()
    for a in range(3):
        mask = np.repeat(mask, s, axis=a)
        pm = np.repeat(pm, s, axis=a)
    h_mm = 0.3 / s
    h = h_mm * 0.001
    ijk = np.argwhere(mask).astype(np.int32)
    n = len(ijk)
    ids = np.full(mask.shape, -1, np.int32)
    ids[mask] = np.arange(n, dtype=np.int32)
    neighbors = np.full((n, 6), -1, np.int32)
    faces = []
    face_nodes = []
    face_d = []
    centres = (ijk + 0.5) * h_mm - 0.15
    for a in range(3):
        for sign in [-1, 1]:
            d = 2 * a + (sign > 0)
            pos = ijk.copy()
            pos[:, a] += sign
            valid = (pos[:, a] >= 0) & (pos[:, a] < mask.shape[a])
            rows = np.flatnonzero(valid)
            neighbors[rows, d] = ids[tuple(pos[rows].T)]
            exposed = np.flatnonzero(neighbors[:, d] < 0)
            p = centres[exposed].copy()
            p[:, a] += 0.5 * sign * h_mm
            faces.append(p)
            face_nodes.append(exposed)
            face_d.append(np.full(len(exposed), d))
    faces = np.concatenate(faces)
    face_nodes = np.concatenate(face_nodes)
    face_d = np.concatenate(face_d)
    is_pulp = pm[mask]
    k = np.where(is_pulp, 0.5, 0.58)
    rc = np.where(is_pulp, 1100 * 3470.0, 1960 * 1590.0)
    C = rc * h ** 3
    g = conductances(neighbors, k, h)
    if s == 1:
        horn_idx = np.argwhere(pulp & (np.indices(pulp.shape)[0] >= np.argwhere(pulp)[:, 0].max() - round(2 / 0.3)))
        (dist, j) = cKDTree(faces).query(horn_idx * 0.3)
        a = int(np.argmin(dist))
        centre = faces[j[a]]
        normal_d = int(face_d[j[a]])
        axis = normal_d // 2
        tangent = (axis + 1) % 3
        positions = []
        for p in range(42):
            v = centre.copy()
            v[tangent] += -0.6 + 1.2 * (p % 7 / 6)
            positions.append(v.tolist())
        dump('raw/SCAN_PATH.json', {'centre_zyx_mm': centre.tolist(), 'outward_face_direction': normal_d, 'tangent_axis': tangent, 'positions_zyx_mm': positions, 'selected_horn_zyx_mm': (horn_idx[a] * 0.3).tolist(), 'scope': 'Deterministic nearest prepared boundary / proxy horn; no expert tooth axes or clinical path.'})
    path = json.loads((ROOT / 'raw/SCAN_PATH.json').read_text())
    positions = np.array(path['positions_zyx_mm'])
    d = path['outward_face_direction']
    eligible = face_d == d
    source_ids = []
    source_weights = []
    for p in positions:
        r2 = ((faces - p) ** 2).sum(1)
        sel = eligible & (r2 < (4 * 0.1575) ** 2)
        ids0 = face_nodes[sel]
        w = np.exp(-r2[sel] / (2 * 0.1575 ** 2)) * h ** 2
        if not len(w) or w.sum() == 0:
            raise ValueError('Gaussian support absent')
        (un, inv) = np.unique(ids0, return_inverse=True)
        weights = np.bincount(inv, weights=w)
        weights /= weights.sum()
        source_ids.append(un)
        source_weights.append(weights)
    length = max(map(len, source_ids))
    pi = np.full((42, length), -1, np.int32)
    pw = np.zeros((42, length))
    for (j, (ix, w)) in enumerate(zip(source_ids, source_weights)):
        pi[j, :len(ix)] = ix
        pw[j, :len(ix)] = w
    wetfaces = eligible & (((faces - np.array(path['centre_zyx_mm'])) ** 2).sum(1) <= 9.0)
    wet_count = np.bincount(face_nodes[wetfaces], minlength=n)
    return (neighbors, g, C, is_pulp, pi, pw, wet_count, k, h, prepared, pulp)

@nb.njit
def advance(T, Tn, neighbors, g, C, cool, Tw, dt):
    loss = 0.0
    for i in range(len(T)):
        flux = cool[i] * (Tw - T[i])
        loss -= flux * dt
        for d in range(6):
            j = neighbors[i, d]
            if j >= 0:
                flux += g[i, d] * (T[j] - T[i])
        Tn[i] = T[i] + dt * flux / C[i]
    return loss

@nb.njit
def simulate(neighbors, g, C, pm, pi, pw, cool, E, eta, Tw, end, dt):
    T = np.zeros(len(C))
    Tn = np.zeros(len(C))
    history = np.zeros((1200, 4))
    nr = 0
    t = 0.0
    loss = 0.0
    energy_in = 0.0
    p = 0
    next_record = 0.0
    maxp = 0.0
    maxt = 0.0
    cem = np.zeros(len(C))
    steps = 0
    while t < end - 1e-12:
        if p < 42 and t >= p / 6.0 - 1e-12:
            for j in range(pi.shape[1]):
                i = pi[p, j]
                if i >= 0:
                    T[i] += eta * E * pw[p, j] / C[i]
            energy_in += eta * E
            p += 1
        pmax = -1e+99
        pmean = 0.0
        count = 0
        for i in range(len(T)):
            maxt = max(maxt, T[i])
            if pm[i]:
                pmax = max(pmax, T[i])
                pmean += T[i]
                count += 1
        maxp = max(maxp, pmax)
        if t >= next_record - 1e-12:
            history[nr, 0] = t
            history[nr, 1] = pmax
            history[nr, 2] = pmean / count
            history[nr, 3] = T.max()
            nr += 1
            next_record += 0.05
        stop = end
        if p < 42:
            stop = min(stop, p / 6.0)
        delta = min(dt, stop - t, next_record - t)
        if delta < 1e-12:
            delta = min(dt, stop - t)
        loss += advance(T, Tn, neighbors, g, C, cool, Tw, delta)
        for i in range(len(T)):
            if pm[i]:
                ta = 37 + T[i]
                tb = 37 + Tn[i]
                ra = 0.25 if ta < 43 else 0.5
                rb = 0.25 if tb < 43 else 0.5
                cem[i] += 0.5 * (ra ** (43 - ta) + rb ** (43 - tb)) * delta / 60.0
        (T, Tn) = (Tn, T)
        t += delta
        steps += 1
    stored = np.dot(C, T)
    history[nr] = np.array([t, T[pm].max(), T[pm].mean(), T.max()])
    nr += 1
    maxp = max(maxp, T[pm].max())
    return (history[:nr], maxp, maxt, stored, loss, energy_in, steps, cem[pm].max())

@nb.njit
def slab_fv(n, tend, q, k, rc, L):
    dx = L / n
    dt = 0.8 * dx * dx / (2 * k / rc)
    T = np.zeros(n)
    Tn = np.zeros(n)
    t = 0.0
    while t < tend - 1e-12:
        d = min(dt, tend - t)
        for i in range(n):
            flux = 0.0
            if i > 0:
                flux += k / dx * (T[i - 1] - T[i])
            else:
                flux += q
            if i < n - 1:
                flux += k / dx * (T[i + 1] - T[i])
            Tn[i] = T[i] + d * flux / (rc * dx)
        (T, Tn) = (Tn, T)
        t += d
    return T

def slab_closed(n, t, q=1000.0, k=0.58, rc=1960 * 1590.0, L=0.002):
    x = (np.arange(n) + 0.5) * L / n
    alpha = k / rc
    j = np.arange(1, 801)
    series = (np.cos(np.outer(x / L, np.pi * j)) * np.exp(-alpha * (j * np.pi / L) ** 2 * t) / j ** 2).sum(1)
    return q * t / (rc * L) + q * L / k * (1 / 3 - x / L + 0.5 * (x / L) ** 2 - 2 / np.pi ** 2 * series)

def impulse(t, z=0.002, sigma=0.0001575, E=0.25, k=0.58, rc=1960 * 1590.0):
    t = np.asarray(t)
    a = k / rc
    t0 = sigma * sigma / (2 * a)
    safe = np.maximum(t, 1e-15)
    out = 2 * E / (rc * (4 * np.pi * a) ** 1.5 * np.sqrt(safe) * (safe + t0)) * np.exp(-z * z / (4 * a * safe))
    return np.where(t > 0, out, 0.0)

def green_response(f, E, N, step=0.01):
    ts = np.arange(0, (N - 1) / f + 8 + step / 2, step)
    pulses = np.arange(N) / f
    return (ts, impulse(ts[:, None] - pulses[None, :], E=E).sum(1))

def verify(pr):
    a = slab_fv(320, 3.0, 1000.0, 0.58, 1960 * 1590.0, 0.002)
    b = slab_closed(320, 3.0)
    norm = max(abs(b))
    err = float(np.max(abs(a - b)) / norm)
    poison = float(np.max(abs(slab_fv(320, 3.0, 1100.0, 0.58, 1960 * 1590.0, 0.002) - b)) / norm)
    tests = {'slab_normalized_Linf': err, 'slab_pass': err < 0.005, 'flux_x1p1_rejected': poison > 0.005}
    t = np.linspace(0.1, 4, 31)
    tau = 0.00025
    rect = np.array([quad(lambda u: float(impulse(tt - u)), 0, tau, epsabs=1e-12)[0] / tau for tt in t])
    kick = impulse(t - tau / 2)
    rel = float(np.max(abs(rect - kick)) / max(abs(rect)))
    tests.update({'rectangular_vs_midpoint_impulse_rel': rel, 'pulse_check_pass': rel < 0.001})
    tests['pulse_z_scale_error_rejected'] = float(np.max(abs(rect - impulse(t - tau / 2, z=0.001))) / max(abs(rect))) > 0.001
    dump('raw/NUMERICAL_VERIFICATION.json', tests)
    assert tests['slab_pass'] and tests['flux_x1p1_rejected'] and tests['pulse_check_pass'] and tests['pulse_z_scale_error_rejected']
    return tests

def external_probe(pr):
    rows = json.loads((ROOT / 'inputs/Geraldo2005_table2.json').read_text())
    dry = [r for r in rows if not r['water']]
    preds = []
    for r in dry:
        (_, v) = green_response(r['frequency_Hz'], r['energy_mJ'] * 0.001, r['pulses'])
        preds.append(float(v.max()))
    pred = np.array(preds)
    obs = np.array([r['reported_rise_C'] for r in dry])
    train = np.array([r['energy_mJ'] != 350 for r in dry])
    test = ~train
    fit = lsq_linear(pred[train, None], obs[train], bounds=(0, 1))
    eta = float(fit.x[0])
    out = []
    for (r, v, o) in zip(dry, pred, obs):
        out.append({**r, 'Green_unit_eta_peak_C': float(v), 'proxy_predicted_peak_C': eta * float(v), 'residual_C': eta * float(v) - float(o), 'split': 'TRAIN' if r['energy_mJ'] != 350 else 'HELD_ENERGY'})
    mae = float(np.mean(abs(eta * pred[test] - obs[test])))
    import xml.etree.ElementTree as ET
    root = ET.parse(ROOT / 'inputs/PMC4327696.xml').getroot()
    table = root.find('.//table-wrap')
    cells = [[''.join(c.itertext()).strip() for c in tr] for tr in table.findall('.//tr')]
    laser = [row for row in cells if row and 'Er:YAG' in row[0]][0]
    (m, sd) = map(float, laser[1].replace(' ', '').split('±'))
    ok = abs(m - 0.84) < 1e-12 and abs(sd - 0.55) < 1e-12

    def table_validator(value):
        return abs(value - m) < 1e-12
    scope = {'source': 'DOI:10.1590/S1678-77572008000300009 Table1 Er:YAG row', 'primary_xml_sha256': sha(ROOT / 'inputs/PMC4327696.xml'), 'measured_max_rise_C': m, 'sample_SD_C': sd, 'resolution_level': 'POPULATION', 'power_conflict_W': [3.5, 0.25 * 4], 'laser_duration_s': None, 'remaining_dentin_mm': 0.5, 'scope_outcome': 'REJECT_PREDICTION_MATCH: unknown irradiation duration/path; 3.5W nominal vs 1W pulse energy-frequency product; sensor-paste maximum != all-pulp point maximum', 'xml_table_match': ok, 'injected_plus_1C_rejected': not table_validator(m + 1)}
    result = {'eta_fit_proxy': eta, 'training_groups': int(train.sum()), 'held_groups': int(test.sum()), 'held_proxy_MAE_C': mae, 'frozen_proxy_MAE_gate_pass': mae <= 0.5, 'external_statistic_match': 'UNKNOWN for 2005 signed statistic; proxy fit is not validation', 'wet_negative_group_means': sum((r['water'] and r['reported_rise_C'] < 0 for r in rows)), 'wet_zero_drive_model_refuted': 'Nonnegative laser source with Tw=T0 cannot generate negative temperature change, regardless of mesh resolution (maximum principle)', 'independent_2008': scope, 'rows': out, 'physical_validation': 'UNKNOWN_UNMATCHED'}
    dump('raw/EXTERNAL_COMPARISON_R1.json', result)
    assert ok and scope['injected_plus_1C_rejected']
    return result

def main():
    start = time.time()
    pr = json.loads((ROOT / 'PREREG_R1.json').read_text())
    expected = (ROOT / 'PREREG_R1.sha256').read_text().split()[0]
    assert sha(ROOT / 'PREREG_R1.json') == expected
    for r in pr['sources'].values():
        assert sha(r['path']) == r['sha256']
    tests = verify(pr)
    external = external_probe(pr)
    summary = []
    for s in [1, 2, 4]:
        t0 = time.time()
        geo = geometry(s, pr)
        (neighbors, g, C, pm, pi, pw, wet, k, h, prep, pulp) = geo
        cool = np.zeros(len(C))
        dt = 0.8 / np.max(g.sum(1) / C)
        (hist, peak, hard, stored, loss, ein, steps, cem) = simulate(neighbors, g, C, pm, pi, pw, cool, 0.25, 1.0, 0.0, 10.0, dt)
        rel = abs(stored + loss - ein) / ein
        row = {'subdivision': s, 'pitch_mm': 0.3 / s, 'nodes': len(C), 'pulp_nodes': int(pm.sum()), 'dt_max_s': float(dt), 'peak_unit_eta_pulp_C': float(peak), 'peak_unit_eta_hard_C': float(hard), 'stored_J': float(stored), 'removed_J': float(loss), 'input_J': float(ein), 'relative_energy_balance': float(rel), 'energy_gate_pass': rel < 1e-07, 'steps': int(steps), 'wall_s': time.time() - t0, 'resolution_level': 'PER_POINT', 'scope': 'constant-property fixed voxel domain; eta=1 is an energy partition endpoint, not a physical ablation prediction'}
        np.savetxt(ROOT / f'raw/HISTORY_s{s}.csv', hist, delimiter=',', header='time_s,max_pulp_rise_C,mean_pulp_rise_C,max_hard_rise_C', comments='')
        summary.append(row)
        dump('raw/REFINEMENT_R1.json', summary)
        dump('CURRENT_WORK_STATE.json', {'lane': 'X33-laser-pulp', 'claim_type': 'capability', 'phase': 'R1_REFINEMENT', 'latest': row, 'last_gate': 'ENERGY_PASS' if row['energy_gate_pass'] else 'ENERGY_FAIL', 'next_operation': 'Complete frozen subdivisions and cooling scenarios; no physical certificate yet'})
        print(json.dumps(row), flush=True)
        assert row['energy_gate_pass']
    native = geometry(1, pr)
    (neighbors, g, C, pm, pi, pw, wet, k, h, prep, pulp) = native
    scenarios = []
    for hw in [0.0, 500.0, 2000.0]:
        cool = wet * h * h / (1 / max(hw, 1e-100) + h / (2 * k)) if hw else np.zeros(len(C))
        dt = 0.8 / np.max((g.sum(1) + cool) / C)
        (hist, peak, hard, stored, loss, ein, steps, cem) = simulate(neighbors, g, C, pm, pi, pw, cool, 0.25, 0.02, -17.0, 10.0, dt)
        scenarios.append({'h_W_m2K': hw, 'eta': 0.02, 'Twater_C': 20.0, 'max_pulp_rise_C': float(peak), 'final_max_pulp_rise_C': float(hist[-1, 1]), 'final_mean_pulp_rise_C': float(hist[-1, 2]), 'CEM43_max_pulp_min': float(cem), 'max_hard_rise_C': float(hard), 'relative_energy_balance': float(abs(stored + loss - ein) / ein), 'resolution_level': 'PER_POINT', 'empirical_status': 'PHENOMENOLOGICAL_SCENARIO_NOT_MEASUREMENT'})
        np.savetxt(ROOT / f'raw/COOLING_h{int(hw)}.csv', hist, delimiter=',', header='time_s,max_pulp_rise_C,mean_pulp_rise_C,max_hard_rise_C', comments='')
    unit = summary[-1]['peak_unit_eta_pulp_C']
    centre = unit / 2
    floor = unit / 2
    port = ScalarPort(5.5 - centre, 'degC', bias_bound=floor, standard_deviations=(), coefficient=0.0, physical_model_validated=False).evaluate(0.075)
    direct = interval_decision(0.0, unit, 5.5)
    diagnostic = refinement_sequence([r['peak_unit_eta_pulp_C'] for r in summary], [r['pitch_mm'] for r in summary])
    last = abs(summary[-1]['peak_unit_eta_pulp_C'] - summary[-2]['peak_unit_eta_pulp_C'])
    distribution = ndi.distance_transform_edt(~(prep & ~ndi.binary_erosion(prep)), sampling=0.3)[pulp]
    np.savez_compressed(ROOT / 'raw/GEOMETRY_SMALL.npz', prepared=prep, pulp=pulp, spacing_mm=0.3)
    out = {'claim_type': 'capability', 'round': 'R1', 'refinement': summary, 'scenarios': scenarios, 'numerical_tests': tests, 'external': external, 'remaining_hard_tissue_distribution_mm': {'min': float(distribution.min()), 'P5': float(np.quantile(distribution, 0.05)), 'median': float(np.median(distribution)), 'P95': float(np.quantile(distribution, 0.95)), 'resolution_level': 'PER_POINT', 'warning': 'distances from pulp voxels to prepared hard-tissue boundary centres, includes UNKNOWN enamel/dentin split; not clinical RDT'}, 'decision': {'threshold_C': 5.5, 'retained_energy_interval': [0, 1], 'peak_pulp_interval_C': [0, unit], 'X26': port, 'direct_interval': direct, 'dominant_sigma': None, 'binding_quantity': 'retained heat fraction eta; spray coverage/h/Tw, ablation evolution, properties and anatomy/observation UNKNOWN', 'physical_class': 'KLASS_2_INFORMATION_LIMITED; no empirical residual-heat/cooling bounds', 'physical_certification': 'UNKNOWN', 'class_scope': 'X26 class is conditional on the fixed constant-property model; finite refinement is diagnostic only'}, 'refinement_diagnostic': diagnostic, 'last_refinement_change_C': last, 'refinement_gate_pass': last < 0.15, 'mutation_tests': {'energy_extra_0p05J_rejected': abs(summary[-1]['stored_J'] + 0.05 - summary[-1]['input_J']) / summary[-1]['input_J'] > 1e-07, 'claiming_interval_below_rejected': direct != 'BELOW', 'pulse_error_rejected': tests['pulse_z_scale_error_rejected'], 'XML_plus_1C_rejected': external['independent_2008']['injected_plus_1C_rejected']}, 'cost': {'wall_s': time.time() - start, 'CPU_s': time.process_time(), 'maxRSS_KiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, 'threads': 1, 'GPU': False, 'preparation_and_total_agent_cost': 'UNKNOWN_UNINSTRUMENTED', 'fit_groups': 6, 'validation_groups': 3, 'queries': '42 impulses * three subdivisions + three scenarios; includes verification and JIT'}, 'external_referent': pr['external_referent'], 'dropout': {'literature_discovery_papers': 27, 'retained_primary_protocols': 2, 'excluded': 25, 'excluded_fraction': 25 / 27, 'reason': '27-paper modern LIT corpus supplies 2005 discovery citations only; primary 2005 author table and separate cached 2008 XML are outside this 27-paper set. 27-paper denominator is discovery scope, not thermal case eligibility.'}, 'physical_validation': False}
    dump('raw/R1_RESULTS.json', out)
    print(json.dumps({'finished': 'R1', 'unit_eta_peak_C': unit, 'X26_class': port['class'], 'last_change_C': last, 'wall_s': out['cost']['wall_s']}), flush=True)
if __name__ == '__main__':
    main()
