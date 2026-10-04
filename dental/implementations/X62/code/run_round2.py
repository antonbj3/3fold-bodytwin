from dental_release.paths import expand as _release_expand
import csv, json, time
from fractions import Fraction
import numpy as np
from common import ROOT, DATA, save, state, utc, sha
from patient import load
from beam import stiffness, condense, direct, exact_binary_box
from measurement import assess_pose, assess_force

def run():
    start = time.perf_counter()
    (points, vertical, rows, meshes) = load()
    state('R2_RUNNING', 'R2 preregistered; no new physical measurement', 'freeze modal predictions before virtual controls')
    kAP = Fraction('1.27') / Fraction('.06')
    kV = Fraction('.95') / Fraction('.02')
    EIy = float(kAP * 1000 / 48)
    EIz = float(kV * 1000 / 48)
    (K, lengths) = stiffness(points, vertical, EI=(EIy, EIz))
    (S, ui, ri) = condense(K)
    (vals, vecs) = np.linalg.eigh(S)
    threshold = float(np.max(vals) * 1e-10)
    positive = np.flatnonzero(vals > threshold)
    rank = len(positive)
    assert rank == 12, 'unexpected additional constraints; no silent fixed rank'
    modes = []
    for (m, j) in enumerate(positive):
        v = vecs[:, j]
        v *= 1 if v[np.argmax(np.abs(v))] >= 0 else -1
        u = v.reshape(6, 3) / 1024
        f = (S @ u.ravel()).reshape(6, 3)
        modes.append(dict(mode=m + 1, displacement_mm=u.tolist(), force_N=f.tolist(), eigenvalue_N_per_mm=float(vals[j]), resolution='PER_TOOTH', timescale='SIMULTANEOUS'))
    nominal_pose = [dict(FDI=row['FDI'], **{k: 0.0 for k in ['mesiodistal_mm', 'occlusogingival_mm', 'buccolingual_mm', 'torque_deg', 'rotation_deg', 'angulation_deg']}) for row in rows]
    frozen = dict(schema='X62-frozen-lab-predictions-v1', created_utc=utc(), prereg_sha256=sha(ROOT / 'PREREG_R2.json'), patient=_release_expand('Open-Full-Jaw @DENTAL_CASE_ID@ mandible'), model='anisotropic source-calibrated Euler-Bernoulli; axial and torsional closures uncalibrated', force_mode_predictions=modes, nominal_pose_errors=nominal_pose, source_read_before_freeze=True, independent_physical_measurement='NOT_RUN', mechanical_uncertainty='UNKNOWN constitutive remainder; source rounding alone not a patient force guarantee', amplitude_mm=1 / 1024, claim_type='capability', measurement_domain='12 translation modes, six bonds held at prescribed translations with free rotations')
    frozen_path = ROOT / 'FROZEN_PREDICTIONS.json'
    if frozen_path.exists():
        old = json.loads(frozen_path.read_text())
        assert old['force_mode_predictions'] == modes and old['nominal_pose_errors'] == nominal_pose, 'prediction drift; create a new round instead of overwrite'
        frozen = old
    else:
        save('FROZEN_PREDICTIONS.json', frozen)
    h = sha(frozen_path)
    save('raw/FROZEN_RECEIPT.json', dict(sha256=h, frozen_utc=frozen['created_utc'], measurement_status='NOT_RUN'))
    virtual = []
    for mode in modes:
        for (fdi, force) in zip([33, 32, 31, 41, 42, 43], mode['force_N']):
            virtual.append(dict(mode=mode['mode'], FDI=fdi, **dict(zip(['Fx_N', 'Fy_N', 'Fz_N'], force))))
    for (name, items) in [('VIRTUAL_FORCE_IDENTITY.csv', virtual), ('VIRTUAL_POSE_IDENTITY.csv', nominal_pose)]:
        with open(ROOT / 'raw' / name, 'w') as f:
            w = csv.DictWriter(f, fieldnames=list(items[0]))
            w.writeheader()
            w.writerows(items)
    force_result = assess_force(virtual, h)
    pose_result = assess_pose(nominal_pose, rows, meshes)
    changed = [dict(x) for x in nominal_pose]
    changed[2]['torque_deg'] = 4.0
    pose_defect = assess_pose(changed, rows, meshes)
    v0 = vecs[:, positive[0]]
    v1 = vecs[:, positive[1]]
    S2 = S + vals[positive[1]] * np.outer(v1, v1)
    known_diff = float(np.linalg.norm((S2 - S) @ v0))
    unknown_diff = float(np.linalg.norm((S2 - S) @ v1))
    ca = np.diag([1.0, 3.0])
    cb = np.diag([3.0, 1.0])
    load_vector = np.array([1.0, 0.0])
    assert np.trace(ca) == np.trace(cb)
    stiffness_suff = dict(summary='trace of 2-direction stiffness', summary_A=4.0, summary_B=4.0, identity_error=0.0, response_A=(ca @ load_vector).tolist(), response_B=(cb @ load_vector).tolist(), downstream_difference=2.0, minimal_extension='two directional coefficients; full patient port additionally needs axial/torsional/interface data', witness_kind='our_own_fixture')
    sourceport = []
    for (label, f, d) in [('AP', Fraction('1.27'), Fraction('.06')), ('vertical', Fraction('.95'), Fraction('.02'))]:
        lo = (f - Fraction('.005')) / d
        hi = (f + Fraction('.005')) / d
        sourceport.append(dict(direction=label, stiffness_N_per_mm=float(f / d), printed_rounding_interval=[float(lo), float(hi)], observation_replay_N=float(f / d * d), heldout=False, resolution='PHENOMENOLOGICAL', debt='measure same manufactured retainer, same restraints, each direction and span', physical_interpolation_remainder='UNKNOWN; one endpoint is not a force law', exact_rationals=[str(lo), str(hi)]))
    center = np.array(modes[0]['displacement_mm'])
    radius = np.full((6, 3), 1 / 65536)
    (lo, hi) = exact_binary_box(S, center, radius)
    box = dict(center_mm=center.tolist(), radius_mm=radius.tolist(), lower_N=lo.tolist(), upper_N=hi.tolist(), enclosure='rigorous exact rational evaluation for frozen binary S coefficients; matrix construction error and physical model discrepancy UNKNOWN', matrix_binary_sha256=__import__('hashlib').sha256(S.tobytes()).hexdigest())
    DATA.mkdir(parents=True, exist_ok=True)
    np.savez(DATA / 'patient_port_R2.npz', K=K, S=S, points=points, vertical=vertical, positive_eigenvalues=vals[positive], positive_modes=vecs[:, positive], lengths=lengths)
    (ctrl, _) = direct(K, np.array(modes[0]['displacement_mm']))
    control_error = float(np.max(np.abs(ctrl - np.array(modes[0]['force_N']))))
    save('raw/R2.json', dict(claim_type='capability', source_port=sourceport, patient_port_nonrigid_rank=rank, rigid_modes=18 - rank, numerical_control_max_abs_N=control_error, one_mode_calibration_counterexample=dict(known_response_difference_N_per_mm=known_diff, unmeasured_response_difference_N_per_mm=unknown_diff, threshold_N_per_mm=threshold, interpretation='linear-algebra witness; data queries required to identify arbitrary symmetric linear port'), stiffness_sufficiency=stiffness_suff, force_box=box, pose_identity=pose_result, pose_4degree_defect=pose_defect, virtual_force_identity_pass=force_result['pass_gate'], frozen_prediction_sha256=h, data_file=dict(path=str(DATA / 'patient_port_R2.npz'), bytes=(DATA / 'patient_port_R2.npz').stat().st_size, sha256=sha(DATA / 'patient_port_R2.npz')), physical_validation='NOT_RUN; calibration replay is not external validation', wall_seconds=time.perf_counter() - start))
    state('R2_COMPLETE', 'frozen modal/pose interfaces constructed; physical measurements NOT_RUN', 'run independent source, geometry, interval and fault controls; package lab demo and graph feedback')
    print(json.dumps(dict(nonrigid_modes=rank, rigid_modes=18 - rank, source_port=sourceport, force_comparator_error_N=control_error, pose_defect_passes=pose_defect['all_six_pass'], frozen_sha256=h), indent=2))
if __name__ == '__main__':
    run()
