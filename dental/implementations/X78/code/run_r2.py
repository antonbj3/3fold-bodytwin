from dental_release.paths import expand as _release_expand
from common import *
import csv, importlib.util, math, resource, sys, time, zipfile
import numpy as np
from fractions import Fraction as F
from scipy import ndimage as ndi
from tissue import exact_box_query, scipy_box_control

def module(name, p):
    s = importlib.util.spec_from_file_location(name, p)
    m = importlib.util.module_from_spec(s)
    s.loader.exec_module(m)
    return m

def main():
    start = time.perf_counter()
    cpu = time.process_time()
    diskcheck()
    tf = module('x78_tfio', DENT / 'cells/geometry/tf2_io.py')
    archive = module('x78_pulpy_io', X12 / 'code/archive_io.py')
    pair = next((x for x in read(X12 / 'raw/R6_pairs.json') if x['case'] == 'P1'))
    with zipfile.ZipFile(tf.ZIP) as z:
        lb = z.read(f'{tf.ROOT}/labelsTr/ToothFairy2P_001.mha')
        ib = z.read(f'{tf.ROOT}/imagesTr/ToothFairy2P_001_0000.mha')
    assert hashlib.sha256(lb).hexdigest() == pair['TF2_tooth_sha256']
    assert hashlib.sha256(ib).hexdigest() == pair['tf2_CT_sha256']
    (lab, sp, lh) = tf.read_mha_bytes(lb)
    (image, spi, ih) = tf.read_mha_bytes(ib)
    assert sp == spi == (0.3, 0.3, 0.3) and lh['TransformMatrix'] == '1 0 0 0 1 0 0 0 1' and (lh.get('Offset') == '0 0 0')
    idx = {r['name']: r for r in archive.index()}
    (rawim, rawsha) = archive.nifti(idx['Pulpy3D/P1/data.nii.gz'])
    (rawp, psha) = archive.nifti(idx['Pulpy3D/P1/gt_instance.nii.gz'])

    def pose(a):
        return np.flip(np.transpose(a, (2, 1, 0)), axis=1)
    assert rawsha == pair['pulpy_CT_sha256'] and psha == pair['pulpy_pulp_sha256']
    assert np.array_equal(pose(np.asanyarray(rawim.dataobj)), image), 'CT identity'
    pulpfull = pose(np.asanyarray(rawp.dataobj))
    crop = np.load(_release_expand('@DENTAL_WORK_ROOT@/X12/R8_P1_36_paired.npz'))
    origin = crop['crop_origin_zyx']
    sl = tuple((slice(int(v), int(v + n)) for (v, n) in zip(origin, crop['tooth'].shape)))
    crown = next((r for r in csv.DictReader((X12 / 'raw/R8_teeth.csv').open()) if r['case'] == 'P1' and r['fdi'] == '36'))
    flipped = int(crown['crown_axis_sign']) < 0
    tooth = crop['tooth'][::-1] if flipped else crop['tooth']
    pulp = crop['pulp'][::-1] if flipped else crop['pulp']
    assert np.array_equal(tooth, lab[sl] == 36)
    assert np.array_equal(pulp, pulpfull[sl] == 46), 'Source Pulp46 maps to canonical FDI36'
    prep = tooth & (ndi.distance_transform_edt(tooth, sampling=sp) >= 0.5 + 0.5 * 0.3)
    surface = prep & ~ndi.binary_erosion(prep, structure=ndi.generate_binary_structure(3, 1))
    allpts = np.argwhere(surface) + origin
    selected = np.linspace(0, len(allpts) - 1, min(96, len(allpts)), dtype=int)
    q = allpts[selected] * 2
    pc = np.argwhere(pulp) + origin
    pr = []
    oracle_rows = []
    for (i, qt) in enumerate(q):
        d = exact_box_query(qt, pc)
        (v, n) = scipy_box_control(qt, pc)
        err = max(abs(v - d['distance_interval_mm'][0]), abs(v - d['distance_interval_mm'][1]))
        assert err <= 1e-08
        pr.append(dict(query=i, query_ticks=qt.tolist(), **d, physical_distance_mm=None, registration_error_mm=None, anatomical_error_mm=None))
        oracle_rows.append(dict(stage='pulp', query=i, exact_mid_mm=sum(d['distance_interval_mm']) / 2, control_mm=v, error_mm=err, candidates=n))
    centre = np.argwhere(lab == 36).mean(0)
    target = centre + np.array([0, 5 / 0.3, 0])
    bone = np.argwhere(lab == 1)
    b = bone[np.argmin(np.sum((bone - target) ** 2, axis=1))]
    bi = b + np.arange(-10, 11)[:, None] * np.array([1, 0, 0])
    assert np.all(bi >= 0) and np.all(bi < lab.shape)
    gray = image[tuple(bi.T)].copy()
    gray_control = np.array([image[tuple(j)] for j in bi])
    assert np.array_equal(gray, gray_control)
    source = read(X73 / 'raw/patients/P1.json')
    cached = np.load(source['point_artifact']['path'])
    walls = {'tf2': np.argwhere(np.isin(lab, [3, 4]))}
    for src in source['source_bindings']:
        if src['dataset'] == 'tf2':
            continue
        with zipfile.ZipFile(src['archive']) as z:
            if z.getinfo(src['member']).compress_type == 9:
                import subprocess
                raw = subprocess.check_output(['7z', 'e', '-so', '-mmt=1', src['archive'], src['member']], stderr=subprocess.DEVNULL)
            else:
                raw = z.read(src['member'])
        assert hashlib.sha256(raw).hexdigest() == src['member_sha256']
        import io
        a = np.load(io.BytesIO(raw))
        assert a.shape == lab.shape
        walls[src['dataset']] = np.argwhere(a > 0)
    for (k, c) in walls.items():
        assert np.array_equal(c, cached[k + '_voxels_zyx'])
    canal = []
    for (i, qt) in enumerate(bi * 2):
        distances = {}
        for (k, c) in walls.items():
            d = exact_box_query(qt, c)
            (v, n) = scipy_box_control(qt, c)
            err = max(abs(v - d['distance_interval_mm'][0]), abs(v - d['distance_interval_mm'][1]))
            assert err <= 1e-08
            distances[k] = d
            oracle_rows.append(dict(stage='canal_' + k, query=i, exact_mid_mm=sum(d['distance_interval_mm']) / 2, control_mm=v, error_mm=err, candidates=n))
        canal.append(dict(query=i, query_ticks=qt.tolist(), distances=distances, revision_range_mm=[min((d['distance_interval_mm'][0] for d in distances.values())), max((d['distance_interval_mm'][1] for d in distances.values()))], union_distance_interval_mm=[min((d['distance_interval_mm'][0] for d in distances.values())), min((d['distance_interval_mm'][1] for d in distances.values()))], physical_clearance_mm=None, independent_annotator_spread_mm=None))
    faults = []
    for row in oracle_rows:
        faults.append(dict(stage=row['stage'], query=row['query'], injected_plus_mm=0.01, rejected=abs(row['exact_mid_mm'] + 0.01 - row['control_mm']) > 1e-08))
    faults.append(dict(stage='gray', injected_plus_gray=1, rejected=not np.array_equal(gray + 1, gray_control)))
    faults.append(dict(stage='crop_origin', injected_one_voxel=True, rejected=not np.array_equal(tooth, lab[tuple((slice(int(v + 1), int(v + 1 + n)) for (v, n) in zip(origin, tooth.shape)))] == 36)))
    faults.append(dict(stage='pulp_FDI', injected_identity_FDI=True, rejected=not np.array_equal(pulp, pulpfull[sl] == 36)))
    assert all((r['rejected'] for r in faults))
    asset = DATA / 'P1_operator_arrays.npz'
    np.savez_compressed(asset, prep_query_ticks=q, pulp_centres_indices=pc, bone_profile_indices=bi, gray_values=gray, profile_labels=lab[tuple(bi.T)], **{k + '_canal_indices': v for (k, v) in walls.items()}, tooth=tooth, pulp=pulp, preparation=prep, crop_origin=origin)
    result = dict(claim_type='capability', status='NATIVE_ANNOTATION_TISSUE_OPERATOR_DELIVERED_DEMO1_OPEN', validation_patient='P1 == ToothFairy2P_001; not Demo1', external_referent=read(ROOT / 'PREREG_R2_NATIVE_TISSUE_OPERATOR.json')['external_referent'], CT_voxels_compared=int(image.size), CT_identity_error=0, crop_identity_error=0, frame_repair=dict(X12_crop_z_reversed=flipped, source_crown_sign=int(crown['crown_axis_sign']), native_frame='z,y,x voxel centres; .3mm; zero offset; identity direction', local_crop_to_native='native_z=origin_z+size_z-1-local_z for reversed X12 crop'), exact_grid_tick_mm='3/20', resolution='PER_POINT', pulp_queries=pr, bone_profile=dict(indices_zyx=bi.tolist(), gray_values=gray.tolist(), labels=lab[tuple(bi.T)].astype(int).tolist(), probe_center_zyx_mm=(b * 0.3).tolist(), calibration='UNCALIBRATED_SOURCE_GRAY', implant_site='VIRTUAL_LABELLED_BONE_PROBE_NOT_REAL_OSTEOTOMY', density=None, modulus=None, gray_control_error=0), canal_queries=canal, control_max_error_mm=max((x['error_mm'] for x in oracle_rows)), controls_count=len(oracle_rows), faults=faults, source_bindings=source['source_bindings'], artifact=dict(path=str(asset), sha256=sha(asset), bytes=asset.stat().st_size), full_patient_gate=False, physical_calibration=False, independent_annotations=False, cost=dict(wall_seconds=time.perf_counter() - start, cpu_seconds=time.process_time() - cpu, peak_RSS_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, threads_max=4, gpu=False, data_bytes=diskcheck(), reasoning_preparation_fit_acquisition='UNKNOWN inherited costs; no speedup claim'))
    dump(ROOT / 'raw/R2_RESULTS.json', result)
    dump(ROOT / 'raw/CONTROLS_R2.json', oracle_rows)
    state('R2_DECIDED', 'Digital pointwise tissue geometry and controls PASS; Demo1 and physical calibration UNKNOWN', 'Freeze changed R3: common finite translation into actual fixed P1 tissue queries')
    (ROOT / 'HANDOFF_R2.md').write_text('# R2 handoff\n\nActual pulp, uncalibrated gray and canal revisions are now queried in verified P1 native CT. P1 remains separate from X77 Demo1. Exact box distances enclose integer-grid geometry; SciPy box projections are an executed independent control. Physical accuracy and independent annotator spread UNKNOWN. Next: replace scalar registration summary with directed finite displacement at the actual query points; freeze before running.\n')
    print(result['status'], 'controls', len(oracle_rows), 'max error', result['control_max_error_mm'])
if __name__ == '__main__':
    main()
