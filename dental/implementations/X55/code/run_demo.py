import datetime
import json
import sys
import time
import shutil
from pathlib import Path
import numpy as np
import trimesh
from scipy.optimize import least_squares
from scipy.spatial.transform import Rotation
from metrology import analyze, sha, write_json, transform, sample_surface, closest, load_mesh, load_regions, REGIONS
from prepare_demo import run as prepare
from literature import run as literature
R = Path(__file__).resolve().parents[1]

def check_frozen(name):
    expected = (R / (name + '.sha256')).read_text().split()[0]
    if sha(R / name) != expected:
        raise ValueError('Frozen file hash mismatch: ' + name)

def known_deformation(path, out):
    d = load_mesh(path / 'design.stl', 'mm')
    labels = load_regions(path / 'regions.json', path / 'design.stl', d)
    (pts, n, faces, _, _) = sample_surface(d, 2400, 5501, labels)
    bary = trimesh.triangles.points_to_barycentric(d.triangles[faces], pts)
    meta = json.loads((path / 'metadata.json').read_text())
    results = {}
    for tag in ['sinter', 'rounding', 'occlusal']:
        rep = analyze(path / 'design.stl', [path / (tag + '.stl')], path / 'regions.json', out / tag, mode='datum', metadata=path / 'metadata.json')
        scan = load_mesh(path / (tag + '.stl'), 'mm')
        scan.apply_transform(meta['scans'][tag + '.stl']['scan_to_reference'])
        expected_pts = np.einsum('ij,ijk->ik', bary, scan.triangles[faces])
        truth = np.einsum('ij,ij->i', expected_pts - pts, n) * 1000
        measured = np.load(out / tag / 'point_fields.npz')['normal_deviation_mm'][0] * 1000
        results[tag] = dict(normal_field_rmse_um=float(np.sqrt(np.mean((truth - measured) ** 2))), true_field_rms_um=float(np.sqrt(np.mean(truth ** 2))), gate=bool(np.sqrt(np.mean((truth - measured) ** 2)) <= 8), reference='our_own_fixture withheld vertex correspondence; not external physical facit')
    return results

def control(path, candidate):
    tick = time.perf_counter()
    d = load_mesh(path / 'design.stl', 'mm')
    labels = load_regions(path / 'regions.json', path / 'design.stl', d)
    scan = load_mesh(path / 'occlusal.stl', 'mm')
    q = sample_surface(scan, 2400, 5502)[0]
    T = np.array(candidate['registration'][0]['scan_to_reference'])
    q = transform(q, T)

    def fun(v):
        q1 = q @ Rotation.from_rotvec(v[:3]).as_matrix().T + v[3:]
        (p, dist, f) = closest(d, q1)
        n = d.face_normals[f]
        groups = labels[f]
        e = np.einsum('ij,ij->i', q1 - p, n)
        w = np.zeros(len(e))
        for g in np.unique(groups):
            take = (groups == g) & (dist < 0.8)
            w[take] = 1 / max(1, take.sum())
        a = abs(e)
        rho = np.where(a <= 0.02, e * e, 0.04 * a - 0.0004)
        return np.sign(e) * np.sqrt(w * rho)
    before = float(fun(np.zeros(6)) @ fun(np.zeros(6)))
    fit = least_squares(fun, np.zeros(6), max_nfev=30, ftol=1e-09, xtol=1e-09, gtol=1e-09)
    return dict(control='SciPy full nearest-triangle robust nonlinear rigid objective, same scan sites/region labels', initialization='candidate final pose; this checks local objective stationarity, not global convergence', candidate_objective_mm2=before, control_objective_mm2=float(fit.fun @ fit.fun), objective_improvement_relative=float((before - fit.fun @ fit.fun) / max(before, 1e-15)), success=bool(fit.success), seconds=time.perf_counter() - tick, no_algorithm_superiority_claim=True)

def faults(path, out):
    cases = []

    def reject(name, call):
        try:
            call()
            cases.append(dict(name=name, rejected=False))
        except (ValueError, KeyError) as e:
            cases.append(dict(name=name, rejected=True, reason=str(e)))

    def call(**kw):
        args = dict(design_path=path / 'design.stl', scan_paths=[path / 'rigid.stl'], regions_path=path / 'regions.json', out=out / 'fault', samples=240)
        args.update(kw)
        return analyze(**args)
    reject('unsupported_unit', lambda : call(units='inch'))
    reject('duplicate_repeats', lambda : call(scan_paths=[path / 'rigid.stl'] * 2, metadata=path / 'metadata.json'))
    broken = json.loads((path / 'regions.json').read_text())
    broken['design_sha256'] = 'bad'
    write_json(out / 'bad_regions.json', broken)
    reject('region_hash_corruption', lambda : call(regions_path=out / 'bad_regions.json'))
    reject('different_or_unknown_specimens', lambda : call(scan_paths=[path / 'repeat_0.stl', path / 'repeat_1.stl']))
    m = load_mesh(path / 'rigid.stl', 'mm')
    m.vertices *= 1000
    m.export(out / 'wrong_units.stl')
    reject('stl_unit_scale_1000', lambda : call(scan_paths=[out / 'wrong_units.stl']))
    reject('missing_datum', lambda : call(mode='datum'))
    return cases

def sufficiency():
    a = np.array([40.0, 0.0, 0.0, 0.0])
    b = np.array([0.0, 40.0, 0.0, 0.0])
    ra = float(np.sqrt(a @ a / 4))
    rb = float(np.sqrt(b @ b / 4))
    da = float(100 - a[0])
    db = float(100 - b[0])
    return dict(summary='global RMS normal deviation', summary_a_um=ra, summary_b_um=rb, identity_error_um=abs(ra - rb), resolution='PER_POINT -> PER_SURFACE_REGION', downstream_quantity='specified local marginal clearance = 100 um - d_marginal', downstream_a_um=da, downstream_b_um=db, difference_um=abs(da - db), sufficient=False, minimal_extension='Region identity and signed local marginal deviation for this one-site consumer', external_referent=dict(kind='closed_form', locator='d = 100 um - u_n; exact arithmetic under stated local geometry', compared_quantity='local clearance', refutes_us=True), assumption='Chosen 100 um nominal clearance; not measured clinical fit')

def run():
    start = time.perf_counter()
    for name in ['PREREG_R1.json', 'DECOMPOSITION_R1.json', 'FROZEN_PREDICTIONS.json']:
        check_frozen(name)
    stamp = datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%S_%fZ')
    out = R / 'replays' / stamp
    out.mkdir(parents=True)
    (out / 'code').mkdir()
    for src in (R / 'code').glob('*.py'):
        shutil.copyfile(src, out / 'code' / src.name)
    code_hashes = {p.name: sha(p) for p in (out / 'code').glob('*.py')}
    prepare()
    lit = literature()
    write_json(R / 'CURRENT_WORK_STATE.json', dict(status='R1_RUNNING', latest_gate='frozen_hashes_PASS', output=str(out), next_operation='known rigid/deformation then registration and precision checks'))
    print('prepared sources and simulated scans', flush=True)
    measures = {}
    rigid = {}
    known = {}
    all_faults = []
    for name in ['D1', 'M2']:
        path = R / 'raw' / name
        rep = analyze(path / 'design.stl', [path / 'rigid.stl'], path / 'regions.json', out / name / 'rigid')
        rigid[name] = dict(max_regional_rms_um=max((r['pooled_normal_rms_um'] for r in rep['regions'].values())), converged=rep['registration'][0]['converged'], initialization='centroid translation only; injected rotation withheld')
        print(name, 'rigid', rigid[name], flush=True)
        known[name] = known_deformation(path, out / name / 'known')
        print(name, 'known deformations', known[name], flush=True)
        modes = {}
        for mode in ['balanced', 'marginal', 'intaglio']:
            modes[mode] = analyze(path / 'design.stl', [path / f'repeat_{j}.stl' for j in range(5)], path / 'regions.json', out / name / mode, mode=mode, metadata=path / 'metadata.json')
            print(name, mode, 'occlusal bias um', modes[mode]['regions']['occlusal']['signed_bias_um'], flush=True)
        measures[name] = modes
        all_faults.extend(faults(path, out / name))
    path = R / 'raw/D1'
    candidate = analyze(path / 'design.stl', [path / 'occlusal.stl'], path / 'regions.json', out / 'control_candidate', metadata=path / 'metadata.json')
    nonlinear = control(path, candidate)
    sf = sufficiency()
    differences = {name: abs(m['balanced']['regions']['occlusal']['signed_bias_um'] - m['intaglio']['regions']['occlusal']['signed_bias_um']) for (name, m) in measures.items()}
    precision_error = max((r['pairwise_identity_error_um'] for m in measures.values() for rep in m.values() for r in rep['regions'].values()))
    gates = dict(rigid_recovery=all((v['max_regional_rms_um'] <= 2 for v in rigid.values())), known_deformation=all((v['gate'] for n in known.values() for v in n.values())), registration_changes_answer=all((v >= 10 for v in differences.values())), precision_identity=precision_error <= 1e-09, global_summary_insufficient=sf['identity_error_um'] == 0 and sf['difference_um'] >= 30, fault_rejection=all((c['rejected'] for c in all_faults)))
    result = dict(claim_type='capability', status='COMPUTATIONAL_REGIONAL_CAPABILITY_PHYSICAL_VALIDATION_UNKNOWN', frozen_prereg_sha256=sha(R / 'PREREG_R1.json'), code_hashes=code_hashes, gates=gates, rigid=rigid, known_deformations=known, registration_difference_occlusal_um=differences, pairwise_identity_error_um_max=precision_error, measurement_reports={n: {m: str(out / n / m / 'report.json') for m in measures[n]} for n in measures}, regions={n: {m: rep['regions'] for (m, rep) in modes.items()} for (n, modes) in measures.items()}, sufficiency=sf, control=nonlinear, faults=all_faults, literature=lit, external_referent=dict(kind='independent_measurement', locator='DOI 10.7759/cureus.39819 Table 1; DOI 10.4047/jap.2023.15.3.155 Table 15', compared_quantity='Published between-crown manufacturing RMS; raw paired field validation unavailable', refutes_us=True), secondary_precision_referent=dict(kind='independent_measurement', locator='DOI 10.1055/s-0042-1758796 Table 5', compared_quantity='Between-manufactured-INLAY pairwise RMS, not crown repeat-scan precision', refutes_us=True), exclusions={'selected_table_candidates': lit['cells'] + lit['rejected'], 'scope': 'Exact selected tables; 68-title exploratory shortlist not a systematic review'}, uncertainty='Conditional repeat-only intervals in CLI; physical reference/pose/sampling uncertainty not enclosed', physical_validation='UNKNOWN_NO_ASBUILT_REFERENCE_SCAN_PAIR', timescale='HANDOVER', resolution=['PER_POINT', 'PER_SURFACE_REGION', 'PER_TOOTH', 'POPULATION'], full_cost={'executed_pipeline_seconds': time.perf_counter() - start, 'coding_discovery_and_acquisition_cost': 'UNKNOWN', 'physical_fallback_cost': 'UNKNOWN independent reference and external datum acquisition'}, replay_output=str(out))
    write_json(out / 'R1_RESULTS.json', result)
    write_json(R / 'R1_RESULTS.json', result)
    write_json(R / 'results.json', result)
    write_json(R / 'CURRENT_WORK_STATE.json', dict(status='R1_COMPLETE', latest_gate=gates, next_operation='R2: keep point sign/location, external datum and sensor/manufacture identifiability', output=str(out)))
    print('R1 gates', gates, flush=True)
    return result
if __name__ == '__main__':
    run()
