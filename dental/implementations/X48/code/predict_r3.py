from common import *
import time, sys, importlib.util
sys.dont_write_bytecode = True
spec = importlib.util.spec_from_file_location('x19_frozen_model', X19 / 'code/model.py')
model = importlib.util.module_from_spec(spec)
spec.loader.exec_module(model)

def predict():
    start = time.monotonic()
    (rows, ix) = pose_index()
    arches = []
    pred = []
    checks = []
    reg = read(ROOT / 'PREREG_R3.json')
    for p in reg.get('selection', {}).get('patients', ['3485', '6457']):
        for j in ['OK', 'UK']:
            ids = sorted({r['fdi'] for r in rows if r['patient'] == p and r['jaw'] == j})
            teeth = {}
            for f in ids:
                d = np.load(CACHE / f'{p}_Sirona_T0_Z{f}.npz')
                pts = d['points'].astype(float)
                h = -pts[:, 1]
                c = d['centroid']
                b = axes(ix, p, j, f)['buccal']
                tip = pts[h >= np.quantile(h, 0.97)].mean(0)
                teeth[str(f)] = dict(fdi=f, centroid_mm=c.tolist(), contact_mm=tip.tolist(), crown_height_mm=float(np.quantile(h, 0.95) - np.quantile(h, 0.05)), buccal_unit=b.tolist(), label_status='published dentist labels propagated by source', centroid_resolution='PER_TOOTH', contact_resolution='PER_SURFACE_REGION', unit='mm')
            a = dict(patient=p, jaw=j, teeth=teeth, dataset='Zenodo11280343', geometry='Source measured baseline crown surface; area samples; not manufactured aligner')
            arches.append(a)
            for f in ids:
                thick = reg['material']['h_anterior_mm' if f % 10 <= 3 else 'h_posterior_mm']
                E = reg['material']['E_MPa']
                base = np.array(ix[p, j, 'Bottmedical', 0, f]['H']['all_crowns'])
                c = np.array(teeth[str(f)]['centroid_mm'])
                b = np.array(teeth[str(f)]['buccal_unit'])
                fixed = model.response(a, f, 0.2, E, thick)
                for t in range(1, 10):
                    P = np.array(ix[p, j, 'Bottmedical', t, f]['H']['all_crowns']) @ np.linalg.inv(base)
                    Q = np.array(ix[p, j, 'Bottmedical', t - 1, f]['H']['all_crowns']) @ np.linalg.inv(base)
                    vector = apply(P, c) - apply(Q, c)
                    activation = float(vector @ b)
                    r = model.response(a, f, activation, E, thick)
                    ctrl = model.response(a, f, activation, E, thick, True)
                    if r['status'].startswith('REFUSED'):
                        pred.append(dict(patient=p, jaw=j, fdi=f, week=t, status=r['status']))
                        continue
                    err = max((np.max(np.abs(np.array(r['wrenches'][k]) - ctrl['wrenches'][k])) for k in r['wrenches']))
                    force = sum((np.array(v[:3]) for v in r['wrenches'].values()), np.zeros(3))
                    torque = sum((np.array(v[3:]) + np.cross(teeth[k]['centroid_mm'], v[:3]) for (k, v) in r['wrenches'].items()), np.zeros(3))
                    check = dict(patient=p, jaw=j, fdi=f, week=t, hermite_max_error=float(err), force_balance_N=float(np.linalg.norm(force)), torque_balance_Nmm=float(np.linalg.norm(torque)), passed=bool(err <= 1e-08 and np.linalg.norm(force) <= 1e-08 and (np.linalg.norm(torque) <= 1e-08)), wrong_force_rejected=bool(np.linalg.norm(force + np.array([1.0, 0.0, 0.0])) > 1e-08), wrong_torque_rejected=bool(np.linalg.norm(torque + np.array([1.0, 0.0, 0.0])) > 1e-08), wrong_hermite_rejected=bool(err + 1.0 > 1e-08))
                    checks.append(check)
                    pred.append(dict(patient=p, jaw=j, fdi=f, week=t, unit='N', moment_unit='Nmm', moment_origin='baseline crown centroid', resolution='PER_TOOTH', timescale='HANDOVER_WEEK', material_resolution='PHENOMENOLOGICAL', status='SCENARIO_ONLY_PHYSICAL_TRANSFER_FAILED', planned_buccal_activation_mm=activation, planned_centroid_vector_mm=vector.tolist(), force_N=r['active_force_norm_N'], moment_Nmm=r['active_moment_norm_Nmm'], geometry_stiffness_N_per_mm=fixed['active_force_norm_N'] / 0.2, force_balance_N=check['force_balance_N'], torque_balance_Nmm=check['torque_balance_Nmm'], wrenches=r['wrenches'], numerical_control=check, patient_physical_force='UNKNOWN_NATURALIGNER_MATERIAL_CONTACT_AND_EXPOSURE', planned_frame='all-crown best-fit, sensitivity-dependent activation'))
    out = dict(claim_type='information_link', rows=pred, arches=arches, prereg_sha256=sha(ROOT / 'PREREG_R3.json'), model_sha256=sha(X19 / 'code/model.py'), physical_validation=False, known_before_freeze='paper findings; X37 published numerical results; R1 signed component tables; no new R2 relative values yet', control=dict(all_pass=all((x['passed'] for x in checks)), max_error=max((x['hermite_max_error'] for x in checks))), material=reg['material'])
    dest = ROOT / 'raw/R3_FROZEN_FORCE_PREDICTIONS.json'
    if dest.exists():
        assert read(dest) == out, 'Frozen X19 prediction drift'
    else:
        write(dest, out)
        write(ROOT / 'FROZEN_PREDICTIONS_R3.json', dict(frozen_utc=now(), prereg_sha256=sha(ROOT / 'PREREG_R3.json'), prediction_file=str(dest.relative_to(ROOT)), sha256=sha(dest), blind=False, known_before_freeze=out['known_before_freeze'], future_local_measurement=False, material_mismatch='PETG proxy versus cellulose-coated Naturaligner; physically uncalibrated'))
    write(ROOT / 'raw/R3_PREDICTION_COST.json', dict(wall_s=time.monotonic() - start, fit_s=0.0, rows=len(pred)))
    print('R3 frozen/replayed force predictions', len(pred), out['control'])
if __name__ == '__main__':
    predict()
