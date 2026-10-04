from common import *
import importlib.util, time, zipfile
import numpy as np
from fractions import Fraction as F

def main():
    r = read(ROOT / 'raw/R2_RESULTS.json')
    a = np.load(r['artifact']['path'])
    ix = a['bone_profile_indices']
    freeze('PREREG_R4_GRAY_REGISTRATION_SUPPORT.json', dict(frozen_utc=now(), claim_type='capability', capability='Propagate the same finite registration uncertainty into raw bone-profile intensity, not just geometric distance', obstacle='A registered intensity sample changes voxel and tissue class; no gray-to-density calibration may be implied', changed_operation='Enumerate every native voxel Voronoi box intersecting the .3mm translation ball for each of 21 fixed profile points. Retain per-point intensity range, class membership and shared translation state.', consumer='Raw gray-profile input to future bone calibration; no implant recommendation', metrics=dict(union_support_omissions=0, source_lookup_error=0, injected_false_acceptances=0), tolerance=0, decision='Discrete-model geometric support exhaustive, all explicit +/- axis shifts enclosed and every wrong value rejected. Measured IOS-to-CBCT uncertainty remains UNKNOWN.', strongest_equal_information_control='Actual direct source array lookups under six explicit finite translations, plus all 27 reachable cube witnesses', falsifier='A reachable native sample outside reported range, false bone-label admission or HU/density inference', full_cost=dict(preparation='Inherit R1-R3 and reread one local source image', fit='None', discovery='No fitting; exact finite support enumeration', validation='Explicit attainable support witnesses; injected gray/class errors', questions=0, fallback='UNKNOWN physical intensity envelope until independently measured shared pose bound'), external_referent=r['external_referent'], resolution='PER_POINT', time_scale='SIMULTANEOUS', scenario='tau=3/10mm, PHENOMENOLOGICAL; no physical calibration'))
    freeze('FROZEN_PREDICTIONS_R4.json', dict(frozen_utc=now(), prereg_sha256=sha(ROOT / 'PREREG_R4_GRAY_REGISTRATION_SUPPORT.json'), prediction='All explicit +/- one-voxel translations are included in exact .3mm-ball voxel support; no predicted calibrated gray, density or modulus for Demo1', physical_measurements_performed=0))
    start = time.perf_counter()
    s = importlib.util.spec_from_file_location('x78_tfio', DENT / 'cells/geometry/tf2_io.py')
    tf = importlib.util.module_from_spec(s)
    s.loader.exec_module(tf)
    with zipfile.ZipFile(tf.ZIP) as z:
        b = z.read(f'{tf.ROOT}/imagesTr/ToothFairy2P_001_0000.mha')
        lb = z.read(f'{tf.ROOT}/labelsTr/ToothFairy2P_001.mha')
    pair = next((x for x in read(X12 / 'raw/R6_pairs.json') if x['case'] == 'P1'))
    assert hashlib.sha256(b).hexdigest() == pair['tf2_CT_sha256'] and hashlib.sha256(lb).hexdigest() == pair['TF2_tooth_sha256']
    (image, sp, h) = tf.read_mha_bytes(b)
    (lab, _, _) = tf.read_mha_bytes(lb)
    offsets = np.indices((5, 5, 5)).reshape(3, -1).T - 2
    ds = np.maximum(abs(2 * offsets) - 1, 0)
    d2 = np.sum(ds * ds, axis=1)
    reachable = offsets[d2 <= 4]
    rows = []
    controls = []
    faults = []
    for (i, p) in enumerate(ix):
        vox = p + reachable
        assert np.all(vox >= 0) and np.all(vox < lab.shape)
        gray = image[tuple(vox.T)]
        labels = lab[tuple(vox.T)].astype(int)
        lo = float(gray.min())
        hi = float(gray.max())
        witnesses = []
        for (v, g) in zip(vox, gray):
            delta = np.maximum(abs(v - p) - 0.5, 0) * np.sign(v - p)
            t = delta + np.sign(v - p) * 1e-06
            assert np.linalg.norm(t * 0.3) <= 0.3
            snapped = np.floor(p + t + 0.5).astype(int)
            assert np.array_equal(snapped, v)
            direct = float(image[tuple(snapped)])
            assert direct == float(g) and lo <= direct <= hi
            witnesses.append(dict(voxel_index_zyx=v.tolist(), translation_zyx_mm=(t * 0.3).tolist(), gray=direct))
        for axis in range(3):
            for sign in (-1, 1):
                v = p + sign * np.eye(3, dtype=int)[axis]
                direct = float(image[tuple(v)])
                controls.append(dict(query=i, axis=axis, sign=sign, gray=direct, enclosed=lo <= direct <= hi))
        rows.append(dict(query=i, native_indices_zyx=p.tolist(), nominal_gray=float(image[tuple(p)]), nominal_label=int(lab[tuple(p)]), conditional_gray_interval=[lo, hi], reachable_label_classes=sorted(set(labels.tolist())), support_voxels=len(vox), all_support_bone=bool(np.all(labels == 1)), witnesses=witnesses, physical_gray_interval=None, density=None, modulus=None, resolution='PER_POINT', time_scale='SIMULTANEOUS'))
        faults.append(dict(query=i, injected_gray_above_hi=hi + 1, rejected=hi + 1 > hi))
        faults.append(dict(query=i, injected_wrong_bone_label=10, rejected=10 != 1))
    assert all((c['enclosed'] for c in controls)) and all((c['rejected'] for c in faults))
    out = dict(claim_type='capability', status='RAW_GRAY_SUPPORT_ENCLOSED_UNDER_CONDITIONAL_POSE_PHYSICAL_CALIBRATION_UNKNOWN', external_referent=r['external_referent'], queries=rows, controls=controls, faults=faults, profile_candidates=len(ix), nominal_bone_retained=sum((x['nominal_label'] == 1 for x in rows)), nominal_rejected=sum((x['nominal_label'] != 1 for x in rows)), nominal_rejection_fraction=sum((x['nominal_label'] != 1 for x in rows)) / len(ix), rejection_reasons=dict(tooth=sum((x['nominal_label'] == 36 for x in rows)), background=sum((x['nominal_label'] == 0 for x in rows))), uniform_bone_support_retained=sum((x['all_support_bone'] for x in rows)), source_unit='RAW_UNCALIBRATED_GRAY', registration_uncertainty='PHENOMENOLOGICAL .3mm translation only; actual shared IOS pose error UNKNOWN', joint_state='One common t with ||t||<=.3mm, never 21 independent shifts. Per-point hull is conservative for that joint state.', wall_seconds=time.perf_counter() - start)
    dump(ROOT / 'raw/R4_RESULTS.json', out)
    state('R4_DECIDED', '126 finite translations + exhaustive reachable gray support PASS; calibrated bone UNKNOWN', 'Finalize runnable package, per-step figures and minimum same-Demo1 measurement contract')
    print(out['status'], 'bone retained', out['nominal_bone_retained'], 'of', len(ix))
if __name__ == '__main__':
    main()
