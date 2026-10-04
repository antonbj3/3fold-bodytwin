import os
for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[k] = '4'
import copy, json, math, resource, time
from pathlib import Path
import numpy as np
from skimage.color import deltaE_ciede2000
from common import ROOT, dump, sha, stamp, check_frozen, inputs, state
from sources import read_sources, source_error
from optics import de00

def freeze_predictions(name, payload):
    p = ROOT / name
    if p.exists():
        old = check_frozen(name)
        assert old['predictions'] == payload, 'Frozen prediction drift'
    else:
        dump(p, {'frozen_utc': stamp(), 'prospective_status': 'Retrospective source-held construction, not a blinded validation', 'prereg_sha256': sha(ROOT / 'PREREG_R1.json'), 'code_sha256': sha(__file__), 'predictions': payload})
        (ROOT / (name + '.sha256')).write_text(sha(p) + '\n')

def tensor_data():
    d = np.load(ROOT / 'raw/stress.npz')
    n = d['nf'].astype(float)
    n /= np.linalg.norm(n, axis=1)[:, None]
    q = {}
    controls = []
    for lc in ('axial', 'offaxis30'):
        a = d['Sf__' + lc].astype(float)
        S = np.zeros((len(a), 3, 3))
        S[:, 0, 0] = a[:, 0]
        S[:, 1, 1] = a[:, 1]
        S[:, 2, 2] = a[:, 2]
        S[:, 0, 1] = S[:, 1, 0] = a[:, 3]
        S[:, 1, 2] = S[:, 2, 1] = a[:, 4]
        S[:, 0, 2] = S[:, 2, 0] = a[:, 5]
        tn = np.einsum('nij,nj->ni', S, n)
        normal = np.einsum('ni,ni->n', tn, n)
        shear = tn - normal[:, None] * n
        q[lc] = np.linalg.norm(shear, axis=1)
        tx = a[:, 0] * n[:, 0] + a[:, 3] * n[:, 1] + a[:, 5] * n[:, 2]
        ty = a[:, 3] * n[:, 0] + a[:, 1] * n[:, 1] + a[:, 4] * n[:, 2]
        tz = a[:, 5] * n[:, 0] + a[:, 4] * n[:, 1] + a[:, 2] * n[:, 2]
        normal2 = tx * n[:, 0] + ty * n[:, 1] + tz * n[:, 2]
        control = np.sqrt((tx - normal2 * n[:, 0]) ** 2 + (ty - normal2 * n[:, 1]) ** 2 + (tz - normal2 * n[:, 2]) ** 2)
        controls.append({'loadcase': lc, 'max_error_MPa_per_N': float(abs(q[lc] - control).max()), 'injected_demand_x2_error_MPa_per_N': float(abs(2 * q[lc] - control).max())})
    return (d, n, q, controls)

def sufficiency(bond, optical, target):
    (d, n, q, control) = tensor_data()
    tied = d['tied']
    eligible = np.flatnonzero(tied & (q['offaxis30'] > 0))
    hi = int(eligible[np.argmax(q['offaxis30'][eligible])])
    lo = int(eligible[np.argmin(q['offaxis30'][eligible])])
    demands = [float(q['offaxis30'][hi]), float(q['offaxis30'][lo])]
    low = next((x['mean_MPa'] for x in bond if x['system'] == 'MEC' and x['with_cover']))
    high = next((x['mean_MPa'] for x in bond if x['system'] == 'RXU' and x['with_cover']))
    A = np.array([low, high])
    B = A[::-1].copy()
    ca = min(A / np.array(demands))
    cb = min(B / np.array(demands))
    rays = np.load(ROOT / 'raw/rays.npz')
    visible = np.flatnonzero(rays['visible'])
    hidden = np.flatnonzero(~rays['visible'])
    good = next((x['mean_Lab'] for x in optical if x['substrate'] == 'CR' and x['cement'] == 'CG'))
    bad = next((x['mean_Lab'] for x in optical if x['substrate'] == 'MT' and x['cement'] == 'TR'))
    OA = np.array([good, bad])
    OB = OA[::-1].copy()
    mean_a = OA.mean(axis=0)
    mean_b = OB.mean(axis=0)
    de_a = de00(OA[0], target)
    de_b = de00(OB[0], target)
    result = {'bond': {'summary': 'equal total patch area, mean strength and full strength histogram', 'equal_area_each_mm2': 1.0, 'identity_error': float(abs(A.mean() - B.mean())), 'states_MPa': [A.tolist(), B.tolist()], 'STS_facet_ids': [hi, lo], 'positions_mm': d['cf'][[hi, lo]].tolist(), 'demands_MPa_per_N': demands, 'conditional_force_screens_N': [float(ca), float(cb)], 'downstream_difference_N': float(abs(cb - ca)), 'ratio': float(max(ca, cb) / min(ca, cb)), 'minimum_extension': 'For this fixed tied linear field and shear query: min_i S_i/q_i. New load, debonding or stiffness changes require the spatially matched traction-separation state.', 'status': 'INSUFFICIENT', 'resolution_level': 'PER_SURFACE_REGION', 'time_scale': 'HANDOVER', 'physical_status': 'synthetic equal-area subpatch states with inherited simulated demand; not a measured bond map'}, 'optical': {'summary': 'same mean Lab, whole Lab histogram, material thickness and nominal material TP', 'identity_error': float(abs(mean_a - mean_b).max()), 'mean_Lab': mean_a.tolist(), 'states_Lab': [OA.tolist(), OB.tolist()], 'STS_ray_ids': [int(visible[0]), int(hidden[0])], 'positions_mm': rays['points'][[visible[0], hidden[0]]].tolist(), 'visible_de00': [de_a, de_b], 'downstream_difference_de00': abs(de_a - de_b), 'threshold': 3.6, 'decisions': [de_a <= 3.6, de_b <= 3.6], 'minimum_extension': 'For fixed visibility/target/threshold: max_visible DeltaE00. Changing illumination, target or view requires matched local Lab/spectral stack field.', 'status': 'INSUFFICIENT', 'resolution_level': 'PER_SURFACE_REGION', 'time_scale': 'SIMULTANEOUS', 'physical_status': 'synthetic equal-area visible/hidden patch assignment, not measured patient colours'}, 'traction_controls': control}
    dump(ROOT / 'raw/SUFFICIENCY_R1.json', result)
    return result

def main():
    started = time.perf_counter()
    p = check_frozen('PREREG_R1.json')
    inputs()
    check_frozen('DECOMPOSITION.json')
    (bond, optical, target) = read_sources()
    dump(ROOT / 'raw/BOND_TABLE.json', bond)
    dump(ROOT / 'raw/OPTICAL_TABLE.json', optical)
    bond_predictions = {r['system']: r['mean_MPa'] for r in bond if not r['with_cover']}
    rows = {(r['substrate'], r['cement']): np.array(r['mean_Lab']) for r in optical}
    opt_predictions = {c: (rows['MT', 'CG'] + sum((rows[s, c] - rows[s, 'CG'] for s in ['ZR', 'CR'])) / 2).tolist() for c in ['TR', 'BL', 'OP', 'A2']}
    freeze_predictions('FROZEN_PREDICTIONS_R1.json', {'bond_MPa': bond_predictions, 'metal_Lab': opt_predictions})
    suff = sufficiency(bond, optical, target)
    btests = []
    for r in bond:
        if r['with_cover']:
            pred = bond_predictions[r['system']]
            err = abs(pred - r['mean_MPa']) / r['mean_MPa']
            btests.append({'system': r['system'], 'predicted_MPa': pred, 'measured_MPa': r['mean_MPa'], 'relative_error': err, 'pass': err <= p['chains']['R05']['tolerance']})
    otests = []
    for (c, pred) in opt_predictions.items():
        observed = rows['MT', c].tolist()
        err = de00(pred, observed)
        otests.append({'cement': c, 'predicted_Lab': pred, 'measured_Lab': observed, 'error_de00': err, 'pass': err <= p['chains']['R16']['tolerance']})
    source_err = source_error(bond, optical, target)
    bad = copy.deepcopy(bond)
    bad[0]['mean_MPa'] += 1
    bad_opt = copy.deepcopy(optical)
    bad_opt[0]['mean_Lab'][0] += 20
    ciede_errors = []
    for r in optical:
        a = r['mean_Lab']
        ours = de00(a, target)
        ref = float(deltaE_ciede2000(np.array(a), np.array(target)))
        ciede_errors.append(abs(ours - ref))
    bad_ciede = de00([optical[0]['mean_Lab'][0] + 20, *optical[0]['mean_Lab'][1:]], target)
    correct_ref = float(deltaE_ciede2000(np.array(optical[0]['mean_Lab']), np.array(target)))
    checks = {'source_value_error': source_err, 'source_pass': source_err <= 1e-10, 'injected_bond_plus1_rejected': source_error(bad, optical, target) > 1e-10, 'injected_optical_L_plus20_rejected': source_error(bond, bad_opt, target) > 1e-10, 'ciede_independent_max_error': max(ciede_errors), 'ciede_pass': max(ciede_errors) <= 1e-10, 'injected_actual_Lab_value_rejected_by_ciede_control': abs(bad_ciede - correct_ref) > 1e-10, 'sufficiency_identity_pass': all((suff[k]['identity_error'] == 0 for k in ['bond', 'optical'])), 'traction_pass': all((c['max_error_MPa_per_N'] <= 1e-07 for c in suff['traction_controls'])), 'injected_actual_tensor_demand_x2_rejected': all((c['injected_demand_x2_error_MPa_per_N'] > 1e-07 for c in suff['traction_controls']))}
    assert all((v for (k, v) in checks.items() if k.endswith('pass') or 'rejected' in k)), checks
    out = {'round': 'R1', 'claim_type': 'capability', 'bond_transfer': {'rows': btests, 'max_relative_error': max((r['relative_error'] for r in btests)), 'gate': 'PASS' if all((r['pass'] for r in btests)) else 'FAIL'}, 'optical_transfer': {'rows': otests, 'max_error_de00': max((r['error_de00'] for r in otests)), 'gate': 'PASS' if all((r['pass'] for r in otests)) else 'FAIL'}, 'sufficiency': suff, 'checks': checks, 'external_referent': {'kind': 'independent_measurement', 'locator': ['doi:10.1016/j.jds.2023.05.011#Table5', 'doi:10.1016/j.heliyon.2023.e23046#Table2'], 'compared_quantity': 'covered SBS [MPa] and held MT mean Lab [CIEDE2000 error]', 'refutes_us': not all((r['pass'] for r in btests + otests))}, 'cost': {'elapsed_s': time.perf_counter() - started, 'max_rss_MiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 'threads_max': 4}}
    dump(ROOT / 'rounds/R1.json', out)
    (ROOT / 'HANDOFF_R1.md').write_text(f"R1 decided. Loss-free SBS- Transfer: {out['bond_transfer']['gate']}, max relativfel {out['bond_transfer']['max_relative_error']:.8f}Additive cement paint for metal: {out['optical_transfer']['gate']}, max ΔE00 {out['optical_transfer']['max_error_de00']:.8f}.\n\nBoth global summaries are insufficient in case of identity errors 0. R2 must change carrying representation to product-specifically covered state and substrate respectively×cement, local correspondence and strict support/UNKNOWNPrimary facits and failed validation gates persist.\n")
    state(status='R1_DECIDED', latest_gate={'bond': out['bond_transfer']['gate'], 'optics': out['optical_transfer']['gate']}, next_operation='Freeze R2: spatial protocol/stack states and rigorous conditional colour-box/force-screen enclosures; run measurement targets and coverage')
    print(json.dumps({'bond': out['bond_transfer'], 'optics': out['optical_transfer'], 'checks': checks}))
if __name__ == '__main__':
    main()
