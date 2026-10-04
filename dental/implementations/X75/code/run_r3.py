"""Full streamed comparison of existing dataset versions, without voxel files on disk."""
from dental_release.paths import expand as _release_expand
import os
for k in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS']:
    os.environ[k] = '4'
import datetime, hashlib, json, pathlib, resource, time, zipfile
import numpy as np
from numpy.lib import format as fmt
P = pathlib.Path(__file__).resolve().parents[1]
TF2 = _release_expand('@DENTAL_CORPUS_ROOT@/geometry/ToothFairy2/ToothFairy2_Dataset.zip')

def run():
    start = time.perf_counter()
    pr = json.loads((P / 'PREREG_R3.json').read_text())
    rows = []
    assert hashlib.sha256((P / 'PREREG_R3.json').read_bytes()).hexdigest() == (P / 'PREREG_R3.sha256').read_text().strip()
    with zipfile.ZipFile(pr['external_referent']['locator']) as z1, zipfile.ZipFile(TF2) as z2:
        for v in pr['selection']:
            r = dict(v)
            n2 = 'Dataset112_ToothFairy2/imagesTr/' + r['case'] + '_0000.mha'
            with z1.open(v['original_member']) as original, z2.open(n2) as newer:
                version = fmt.read_magic(original)
                (shape, fortran, dt) = fmt._read_array_header(original, version)
                hdr = {}
                while True:
                    line = newer.readline().decode().strip()
                    if ' = ' in line:
                        (key, value) = line.split(' = ', 1)
                        hdr[key] = value
                    if line.startswith('ElementDataFile'):
                        break
                newshape = tuple((int(n) for n in hdr['DimSize'].split()[::-1]))
                newdt = {'MET_DOUBLE': np.dtype('<f8'), 'MET_INT': np.dtype('<i4')}[hdr['ElementType']]
                r.update(original_shape=list(shape), TF2_shape=list(newshape), original_dtype=str(dt), TF2_dtype=str(newdt), original_spacing_mm='UNKNOWN_FROM_NPY', TF2_spacing_mm=hdr['ElementSpacing'])
                if shape != newshape or fortran or dt != newdt:
                    r.update(status='REJECTED_SHAPE_TYPE_ORDER', independent_acquisition=False)
                    rows.append(r)
                    continue
                a = hashlib.sha256()
                b = hashlib.sha256()
                maxdiff = 0.0
                nvox = 0
                different = 0
                fault = None
                while True:
                    buf1 = original.read(1024 * 1024)
                    buf2 = newer.read(1024 * 1024)
                    if not buf1 and (not buf2):
                        break
                    if len(buf1) != len(buf2):
                        raise ValueError('Paired payload length mismatch')
                    a.update(buf1)
                    b.update(buf2)
                    x = np.frombuffer(buf1, dtype=dt)
                    y = np.frombuffer(buf2, dtype=newdt)
                    if not (np.isfinite(x).all() and np.isfinite(y).all()):
                        raise ValueError('Nonfinite original measurement')
                    diff = abs(x - y)
                    maxdiff = max(maxdiff, float(diff.max()))
                    different += int(np.count_nonzero(diff))
                    nvox += len(x)
                    if fault is None:
                        yy = y.copy()
                        yy[0] += 1
                        mutated = yy.tobytes()
                        fault = {'wrong_value_detected_by_voxel_route': not np.array_equal(x, yy), 'wrong_value_detected_by_hash_route': hashlib.sha256(buf1).hexdigest() != hashlib.sha256(mutated).hexdigest()}
                ha = a.hexdigest()
                hb = b.hexdigest()
                r.update(status='COMPARED_ALL_VOXELS', n_voxels=nvox, max_abs_difference_gray=maxdiff, different_voxels=different, original_payload_sha256=ha, TF2_payload_sha256=hb, exact_payload_identity=ha == hb, hash_voxel_routes_agree=(ha == hb) == (different == 0), new_independent_density_references=0, new_independent_torque_temperature_pairs=0, independent_acquisition_admissible=False, original_units='GRAY_NOT_CONFIRMED_HU', fault_injection=fault, original_archive_member_crc=z1.getinfo(v['original_member']).CRC, TF2_archive_member_crc=z2.getinfo(n2).CRC)
            rows.append(r)
            print(r['case'], 'n', r['n_voxels'], 'different', r['different_voxels'], 'max', r['max_abs_difference_gray'], flush=True)
            (P / 'CURRENT_WORK_STATE.json').write_text(json.dumps({'lane': 'X75-cbct-bone', 'phase': 'R3_COMPARING', 'updated_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'completed_original_pairs': len(rows), 'latest_gate': 'IN_PROGRESS', 'next_operation': 'finish8streamed counterpart comparisons; aggregate calibration admission'}, indent=2) + '\n')
    out = {'round': 'R3', 'claim_type': 'information_link', 'external_referent': pr['external_referent'], 'selected': len(rows), 'candidate_pairs': pr['candidate_count'], 'compared_all_voxel_pairs': sum((r['status'] == 'COMPARED_ALL_VOXELS' for r in rows)), 'exact_shared_payload_pairs': sum((r.get('exact_payload_identity', False) for r in rows)), 'new_independent_density_references': 0, 'new_independent_torque_temperature_pairs': 0, 'density_gain': 'UNKNOWN_NO_INDEPENDENT_REFERENCE', 'rejections': pr['rejected_before_selection'] + [r for r in rows if r['status'].startswith('REJECTED')], 'unselected_due_budget': pr['unselected_due_frozen_budget'], 'rows': rows, 'gates': {'full_stream_comparison': all((r['status'] == 'COMPARED_ALL_VOXELS' for r in rows)), 'hash_voxel_control': all((r.get('hash_voxel_routes_agree', False) for r in rows)), 'fault_injection': all((all(r.get('fault_injection', {}).values()) for r in rows)), 'independent_calibration_gain': False}, 'cost': {'wall_s': time.perf_counter() - start, 'peak_rss_kib': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, 'bytes_voxel_payload_compared': sum((r.get('n_voxels', 0) * 8 * 2 for r in rows)), 'fit_s': 0, 'threads': 4, 'physical_measurements': 0, 'large_intermediates_bytes': 0}}
    (P / 'results_R3.json').write_text(json.dumps(out, indent=2, allow_nan=False) + '\n')
    print(json.dumps({k: v for (k, v) in out.items() if k != 'rows'}, indent=2))
if __name__ == '__main__':
    run()
