import os
for k in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS']:
    os.environ[k] = '4'
import datetime, hashlib, json, pathlib, resource, sys, time, zipfile, traceback
import numpy as np
from site_measure import ray_grid, connected_run, slab_control, lateral_axis, shell_profile, core_samples, arch_centers, arch_normal
P = pathlib.Path(__file__).resolve().parents[1]
D = P.parent.parent
sys.path.insert(0, str(D / 'cells/geometry'))
import tf2_io

def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for b in iter(lambda : f.read(1024 * 1024), b''):
            h.update(b)
    return h.hexdigest()

def put(n, x):
    (P / n).write_text(json.dumps(x, indent=2, allow_nan=False) + '\n')

def run(round_name='R1'):
    start = time.perf_counter()
    pr = json.loads((P / f'PREREG_{round_name}.json').read_text())
    assert sha(P / f'PREREG_{round_name}.json') == (P / f'PREREG_{round_name}.sha256').read_text().strip()
    poses = json.loads((P / 'raw/FROZEN_POSES.json').read_text())
    data = pathlib.Path(pr['full_cost']['limits']['data_root'])
    data.mkdir(exist_ok=True, parents=True)
    out = P / f'raw/PER_SITE_{round_name}.jsonl'
    out.write_text('')
    manifest = []
    allrows = []
    controls = []
    groups = {}
    for r in poses:
        groups.setdefault(r['case'], []).append(r)
    with zipfile.ZipFile(tf2_io.ZIP) as z:
        for (case, rr) in groups.items():
            tc = time.perf_counter()
            member = f'{tf2_io.ROOT}/labelsTr/{case}.mha'
            buf = z.read(member)
            memberhash = hashlib.sha256(buf).hexdigest()
            (lab, sp, lh) = tf2_io.read_mha_bytes(buf)
            sp = np.array(sp)
            arrayhash = hashlib.sha256(memoryview(lab)).hexdigest()
            image_member = f'{tf2_io.ROOT}/imagesTr/{case}_0000.mha'
            imagebuf = z.read(image_member)
            ihash = hashlib.sha256(imagebuf).hexdigest()
            (img, isp, ih) = tf2_io.read_mha_bytes(imagebuf)
            assert np.array_equal(sp, isp) and img.shape == lab.shape
            assert lh['TransformMatrix'] == ih['TransformMatrix'] and lh['Offset'] == ih['Offset']
            cap = float(img.max())
            centers = arch_centers(lab, sp) if round_name == 'R2' else None
            manifest.append({'case': case, 'label_member': member, 'label_member_sha256': memberhash, 'label_array_sha256': arrayhash, 'image_member': image_member, 'image_member_sha256': ihash, 'header': {k: ih.get(k) for k in ['TransformMatrix', 'Offset', 'AnatomicalOrientation', 'ElementSpacing', 'DimSize']}, 'gray_cap': cap})
            for pose in rr:
                r = dict(pose)
                arrays = {}
                e = np.array(r['entry'])
                a = np.array(r['axis'])
                a /= np.linalg.norm(a)
                try:
                    sourcehash = memberhash if r['kind'].startswith('X8') else arrayhash
                    assert sourcehash == r['source_label_member_sha256'], 'POSE_SOURCE_HASH_MISMATCH'
                    if round_name == 'R2':
                        (b, r['orientation']) = arch_normal(centers, e, a, lab, sp)
                    else:
                        b = lateral_axis(lab, sp, e, a)
                    r['provisional_buccal_axis_zyx'] = b.tolist()
                    r['rays'] = []
                    for depth in [2.0, 4.0, 6.0]:
                        o = e + depth * a
                        (lo, hi, ix) = ray_grid(lab.shape, sp, o, b, -20, 20)
                        L = lab[tuple(ix.T)]
                        run = connected_run(lo, hi, L)
                        q = {'depth_mm': depth, 'envelope_chord_mm': None, 'buccal_extent_mm': None, 'lingual_extent_mm': None, 'profiles': {}}
                        if run:
                            (low, high) = run
                            q.update(envelope_chord_mm=high - low, buccal_extent_mm=high, lingual_extent_mm=-low)
                            for (side, t, u) in [('buccal', high, -b), ('lingual', low, b)]:
                                f = shell_profile(img, lab, sp, o + t * b, u, cap)
                                key = f'd{int(depth)}_{side}'
                                for k in ['gray_nearest', 'gray_interp', 'depth_mm', 'labels', 'ix']:
                                    arrays[key + '_' + k] = f.pop(k)
                                q['profiles'][side] = f
                                nearest = arrays[key + '_gray_nearest']
                                ids = arrays[key + '_ix']
                                valid = np.all((ids >= 0) & (ids < np.array(img.shape)), axis=1)
                                error = max([abs(float(img[tuple(idx)]) - float(v)) for (idx, v, ok) in zip(ids, nearest, valid) if ok], default=0.0)
                                controls.append({'case': case, 'site': r['site'], 'kind': 'gray', 'error': error, 'pass': error == 0})
                        r['rays'].append(q)
                        if depth == 4:
                            cr = slab_control(lab, sp, o, b, -20, 20)
                            error = 0.0 if run is None and cr is None else max((abs(x - y) for (x, y) in zip(run, cr))) if run is not None and cr is not None else 1000.0
                            controls.append({'case': case, 'site': r['site'], 'kind': 'ray_box', 'candidate': run, 'control': cr, 'error': error, 'pass': error <= 1e-08})
                    (lo, hi, ix) = ray_grid(lab.shape, sp, e, a, -5, 50)
                    ax = connected_run(lo, hi, lab[tuple(ix.T)], anchor=4.0)
                    r['axial_envelope_chord_mm'] = ax[1] - ax[0] if ax else None
                    r['axial_run_limits_from_pose_mm'] = ax
                    (ids, gv) = core_samples(img, lab, sp, e, a)
                    arrays['core_ix'] = ids
                    arrays['core_gray'] = gv
                    r['trabecular_candidate_roi'] = {'n_voxels': len(gv), 'gray_quantiles_05_25_50_75_95': np.quantile(gv, [0.05, 0.25, 0.5, 0.75, 0.95]).tolist() if len(gv) >= 30 else None, 'status': 'UNCALIBRATED_GRAY_OBSERVATION' if len(gv) >= 30 else 'INSUFFICIENT_CORE_VOXELS', 'density_g_cm3': None, 'density_calibration_error': 'UNBOUNDED', 'label': 'jaw interior; not independently validated trabecular tissue'}
                    r['status'] = 'MEASURED'
                    r['spacing_mm'] = sp.tolist()
                    r['anatomical_side_validation'] = False
                    r['resolution'] = {'ray': 'PER_POINT', 'site': 'PER_TOOTH', 'source_gray': 'PER_POINT', 'torque_temperature': 'UNKNOWN'}
                    fn = data / (case + '_' + r['site'] + ('_R2' if round_name == 'R2' else '') + '.npz')
                    np.savez_compressed(fn, **arrays)
                    r['raw_file'] = {'path': str(fn), 'sha256': sha(fn), 'bytes': fn.stat().st_size}
                except Exception as exc:
                    r.update(status='FAILED', reason=str(exc), traceback=traceback.format_exc())
                with out.open('a') as f:
                    f.write(json.dumps(r, allow_nan=False) + '\n')
                allrows.append(r)
            print(case, 'sites', len(rr), 'seconds', round(time.perf_counter() - tc, 2), flush=True)
            put('CURRENT_WORK_STATE.json', {'lane': 'X75-cbct-bone', 'phase': round_name + '_MEASURING', 'updated_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'completed_sites': len(allrows), 'latest_gate': 'IN_PROGRESS', 'next_operation': 'complete frozen cohort; then score coverage and missing truth'})
            del buf, lab, imagebuf, img
    put('raw/SOURCE_MANIFEST' + ('_R2' if round_name == 'R2' else '') + '.json', manifest)
    put(f'raw/CONTROLS_{round_name}.json', controls)
    train = []
    test = []
    profiles = []
    for r in allrows:
        if r['status'] != 'MEASURED':
            continue
        q = next((q for q in r['rays'] if q['depth_mm'] == 4))
        if q['envelope_chord_mm'] is not None:
            (train if r['case'] in pr['selection']['train_cases'] else test).append(q['envelope_chord_mm'])
        profiles.extend((v for q in r['rays'] for v in q['profiles'].values()))
    mean = float(np.mean(train)) if train else None
    rmse = float(np.sqrt(np.mean((np.array(test) - mean) ** 2))) if train and test else None
    usable = len(train) + len(test)
    retained = [v for v in profiles if v['rejection_reason'] is None]
    widths = [max(v['apparent_layer_mm_tau_04_05_06']) - min(v['apparent_layer_mm_tau_04_05_06']) for v in retained]
    result = {'round': round_name, 'claim_type': 'information_link', 'external_referent': pr['external_referent'], 'sites': len(allrows), 'usable_geometry_sites': usable, 'site_failures': len([r for r in allrows if r['status'] == 'FAILED']), 'annotation_geometry': {'train_count': len(train), 'test_count': len(test), 'training_population_mean_chord_mm': mean, 'test_RMSE_population_proxy_mm': rmse, 'test_direct_annotation_RMSE_mm': 0.0 if test else None, 'all_chord_quantiles_mm': np.quantile(train + test, [0, 0.25, 0.5, 0.75, 1]).tolist() if usable else None, 'physical_boundary_accuracy': 'UNKNOWN', 'resolution': 'PER_POINT'}, 'apparent_layer': {'candidate_profiles': len(profiles), 'retained_contrast_profiles': len(retained), 'rejected_profiles': len(profiles) - len(retained), 'rejection_fraction': (len(profiles) - len(retained)) / len(profiles) if profiles else None, 'median_threshold_span_mm': float(np.median(widths)) if widths else None, 'anatomical_cortex_accuracy': 'UNKNOWN', 'resolution': 'PER_POINT'}, 'density_accuracy': 'UNKNOWN_NO_VALIDATED_CALIBRATION', 'torque_variance_explained': None, 'temperature_variance_explained': None, 'underdrill_protocol_selection': 'UNKNOWN_NO_PAIRED_TRUTH', 'gates': {'source_control': all((q['pass'] for q in controls)), 'geometry_coverage': usable / len(allrows) >= 0.9, 'site_information_vs_population_proxy': rmse is not None and rmse >= 0.3, 'absolute_cortex_accuracy': 'UNKNOWN', 'absolute_density_accuracy': 'UNKNOWN', 'torque_temperature_gain': 'UNKNOWN'}, 'full_cost': {'wall_s': time.perf_counter() - start, 'peak_rss_kib': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, 'fit': 'training mean only', 'discovery': '24 frozen case IDs', 'validation_checks': len(controls), 'queries': len(allrows), 'fallback_sites': len(allrows), 'threads': 4, 'GPU': False, 'physical_measurements': 0, 'biological_patient_disjointness': 'case-disjoint; cross-dataset biological duplicate check absent'}, 'max_control_error_mm': max([q['error'] for q in controls if q['kind'] == 'ray_box'], default=0.0), 'data_files': [r['raw_file'] for r in allrows if 'raw_file' in r]}
    if round_name == 'R2':
        valid = [r for r in allrows if r.get('orientation', {}).get('status') == 'LOCAL_ARCH_PROXY']
        result['orientation'] = {'local_arch_sites': len(valid), 'fallback_sites': len(allrows) - len(valid), 'max_dot': max([max(abs(r['orientation']['normal_axis_dot']), abs(r['orientation']['normal_tangent_dot'])) for r in valid], default=1), 'median_tangent_residual_mm': float(np.median([r['orientation']['tangent_residual_rms_mm'] for r in valid])) if valid else None, 'physical_plane_accuracy': 'UNKNOWN'}
        result['gates']['local_orientation_coverage'] = len(valid) / len(allrows) >= 0.9
        result['gates']['basis_orthogonality'] = result['orientation']['max_dot'] <= 1e-12
    put(f'results_{round_name}.json', result)
    print(json.dumps({k: v for (k, v) in result.items() if k not in ['data_files']}, indent=2))
if __name__ == '__main__':
    run(sys.argv[1] if len(sys.argv) > 1 else 'R1')
