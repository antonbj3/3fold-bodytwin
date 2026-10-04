"""One-command frozen numerical evaluation. Failures are saved before exit."""
from dental_release.paths import expand as _release_expand
import copy
import datetime as dt
import json
import resource
import sys
import time
from pathlib import Path
import numpy as np
import trimesh
from parents import HERE, sha, transport as T, milling, inherited, verify_dependencies, geometry as G
from design_gate import check, cavity_planes, raw_normals
from export_bridge import export, check_native

def verify_freeze(f):
    active = json.loads((HERE / 'ACTIVE_ROUND.json').read_text())
    if sha(HERE / active['predictions_file']) != (HERE / active['predictions_hash_file']).read_text().strip():
        raise ValueError('Frozen predictions hash changed')
    if sha(HERE / active['prereg_file']) != f['prereg_sha256']:
        raise ValueError('PREREG changed')
    for (p, h) in f['code'].items():
        if sha(HERE / p) != h:
            raise ValueError('Frozen implementation changed: ' + p)
    for row in f['files']:
        if sha(row['path']) != row['sha256']:
            raise ValueError('Frozen input changed: ' + row['path'])
    verify_dependencies()

def independent_witness(r):
    """Published all-facet naive query receives identical target information."""
    w = r['witness']
    d = r.get('measurement', {})
    tri = np.array(d.get('target_triangles_mm', []))
    return (w, d, tri)

def main():
    start = time.perf_counter()
    cpu = time.process_time()
    active = json.loads((HERE / 'ACTIVE_ROUND.json').read_text())
    f = json.loads((HERE / active['predictions_file']).read_text())
    verify_freeze(f)
    tag = dt.datetime.now(dt.timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    run = HERE / 'raw/runs' / tag
    run.mkdir(parents=True)
    rows = []
    tests = []
    reports = {}
    exceptions = []

    def test(name, passed, **details):
        tests.append(dict(name=name, passed=bool(passed), **details))
    for (name, path) in f['cases'].items():
        try:
            paths = {k: str(Path(path) / (k + '.stl')) for k in ('prep', 'antagonist', 'crown')}
            r = check(paths, material='katana-ht', contract=Path(path) / 'contract.json', sinter_factor=1.2)
            reports[name] = r
            T.dump(run / (name + '.json'), r)
            statuses = {k: v['status'] for (k, v) in r['rules'].items()}
            for (rule, expected) in f['expected'][name].items():
                test(name + ':' + rule, statuses[rule] == expected, expected=expected, observed=statuses[rule])
            rows.append(dict(case=name, rule_status=statuses, wall_s=r['cost']['wall_s']))
            meta = json.loads((Path(path) / 'contract.json').read_text())
            unit = meta.get('units')
            scale = 0.001 if unit == 'um' else 1.0
            for key in ('material_wall', 'film_min', 'film_max'):
                rule = r['rules'][key]
                if 'robust_interval_mm' not in rule:
                    continue
                w = rule['witness']
                target = 'crown' if key == 'material_wall' else 'prep'
                region = 'exterior' if target == 'crown' else 'preparation'
                (tri, _) = T.load_mesh(paths[target], 'micron' if unit == 'um' else unit)
                mesh = T.physical_mesh(tri[meta['regions'][region]])
                pt = np.array([w['source_point_mm']])
                (q, dist, ids) = trimesh.proximity.closest_point_naive(mesh, pt)
                err = abs(float(dist[0]) - w['distance_mm'])
                test(name + ':' + key + ':external_naive', err <= 1e-07, parity_mm=err, external='published trimesh closest_point_naive')
                bad = abs(float(dist[0]) - (w['distance_mm'] + 0.01))
                test(name + ':' + key + ':corrupt_witness_rejected', bad > 1e-07, injected_error_mm=0.01)
        except Exception as e:
            exceptions.append(dict(case=name, error=repr(e)))
            test(name + ':exception', False, error=repr(e))
        T.dump(run / 'checkpoint.json', dict(rows=rows, tests=tests, exceptions=exceptions))
    valid = reports.get('valid')
    near = reports.get('near_limit')
    if valid and near:
        a = valid['rules']['film_min']
        b = near['rules']['film_min']
        test('uncertainty_changes_decision', a['status'] == 'PASS' and b['status'] == 'UNKNOWN')
        intervals = []
        for case in ('valid', 'near_limit'):
            r = reports[case]['rules']['film_min']
            e = sum(r['supplied_errors_mm'])
            interval = r['robust_interval_mm']
            base = Path(f['cases'][case])
            (ct, _) = T.load_mesh(base / 'crown.stl', 'mm')
            (pt, _) = T.load_mesh(base / 'prep.stl', 'mm')
            meta = json.loads((base / 'contract.json').read_text())
            inner = ct[meta['regions']['intaglio']]
            for delta in np.linspace(-e, e, 17):
                moved = inner.copy()
                moved[:, :, 2] += float(delta) * r['supplied_errors_mm'][0] / e
                target_tri = pt[meta['regions']['preparation']].copy()
                target_tri[:, :, 2] -= float(delta) * r['supplied_errors_mm'][1] / e
                mesh = T.physical_mesh(target_tri)
                (_, dist, _) = trimesh.proximity.closest_point_naive(mesh, moved.mean(1))
                value = float(dist.min())
                miss = not interval[0] - 1e-07 <= value <= interval[1] + 1e-07
                intervals.append(dict(case=case, delta_mm=float(delta), distance_mm=value, envelope_mm=interval, miss=miss))
        test('hausdorff_envelope_no_misses', not any((x['miss'] for x in intervals)))
        test('injected_envelope_rejected', not 0.1 <= 0.05 - 0.002 <= 0.11, injected_envelope_mm=[0.1, 0.11])
        T.dump(run / 'uncertainty.json', intervals)
    axial_rows = []
    for name in ('axial_valid', 'axial_penetration', 'axial_near'):
        if name not in reports:
            continue
        base = Path(f['cases'][name])
        meta = json.loads((base / 'contract.json').read_text())
        r = reports[name]['rules']['occlusal_contact']
        interval = r['robust_interval_mm']
        (a, _) = T.load_mesh(base / 'antagonist.stl', 'mm')
        (c, _) = T.load_mesh(base / 'crown.stl', 'mm')
        U = a[meta['regions']['antagonist']]
        L = c[meta['regions']['occlusal']]
        rng = np.random.default_rng(49)
        for i in range(11):
            u = U.copy()
            l = L.copy()
            eu = meta['axial_error_mm']['antagonist']
            el = meta['axial_error_mm']['crown']
            if i < 2:
                sign = 1 if i == 0 else -1
                u[:, :, 2] += sign * eu
                l[:, :, 2] -= sign * el
            else:
                for (tris, e) in ((u, eu), (l, el)):
                    (v, inv) = np.unique(tris.reshape(-1, 3), axis=0, return_inverse=True)
                    delta = rng.uniform(-e, e, len(v))
                    tris[:, :, 2] += delta[inv].reshape(-1, 3)
            candidate = G.continuous(u, l, [0, 0, 1], True)
            full = G.continuous(u, l, [0, 0, 1], False)
            x = candidate['minimum_gap_mm']
            parity = abs(x - full['minimum_gap_mm'])
            miss = not interval[0] <= x <= interval[1]
            axial_rows.append(dict(case=name, perturbation=i, minimum_mm=x, envelope_mm=interval, miss=miss, parity_mm=parity))
            test(name + ':axial_enclosure_' + str(i), not miss)
            test(name + ':full_pair_' + str(i), parity <= 1e-07, parity_mm=parity)
            test(name + ':injected_gap_' + str(i), abs(x + 0.01 - full['minimum_gap_mm']) > 1e-07)
        bad = [r['nominal_minimum_gap_mm'] + 0.2, r['nominal_minimum_gap_mm'] + 0.21]
        test(name + ':injected_envelope', not bad[0] <= r['nominal_minimum_gap_mm'] <= bad[1])
    if axial_rows:
        T.dump(run / 'axial_uncertainty.json', axial_rows)
    if 'valid' in reports:
        base = Path(f['cases']['valid'])
        paths = {k: str(base / (k + '.stl')) for k in ('prep', 'crown', 'antagonist')}
        original = json.loads((base / 'contract.json').read_text())
        bad = copy.deepcopy(original)
        bad['input_sha256']['crown'] = '0' * 64
        path = run / 'bad_hash_contract.json'
        T.dump(path, bad)
        rejected = False
        try:
            check(paths, contract=path)
        except ValueError:
            rejected = True
        test('corrupt_source_hash_rejected', rejected)
        r = check(paths)
        test('unitless_STL_unknown', r['rules']['scale']['status'] == 'UNKNOWN')
        classonly = copy.deepcopy(original)
        classonly.pop('ifu_profile')
        path = run / 'class_only.json'
        T.dump(path, classonly)
        r = check(paths, material='3Y', contract=path)
        test('generic_3Y_no_universal_IFU', r['rules']['material_wall']['status'] == 'UNKNOWN')
        threshold = reports['valid']['rules']['material_wall']['threshold_mm']
        test('product_IFU_external_threshold', threshold == 0.5, source='Manufacturer TI-010 ver.019, p2 posterior crown')
        test('injected_IFU_value_rejected', threshold + 0.1 != 0.5)
    sa = np.array(f['sufficiency']['a'])
    sb = np.array(f['sufficiency']['b'])
    identity = float(abs(sa.mean() - sb.mean()))
    downstream = float(abs(sa.min() - sb.min()))
    suff = dict(summary='equal-area mean wall thickness', summary_a_mm=float(sa.mean()), summary_b_mm=float(sb.mean()), identity_error_mm=identity, identity_to_machine_precision=identity == 0.0, downstream='minimum regional separation', a_min_mm=float(sa.min()), b_min_mm=float(sb.min()), downstream_difference_mm=downstream, minimum_sufficient_extension='minimum field for this wall rule; per-facet field and labels for location/other questions', mean_resolution='PER_TOOTH', downstream_resolution='PER_SURFACE_REGION', constructed_states=True, physical_measurements=False)
    inner = np.array([[[0, 0, 0], [1, 0, 0], [0, 1, 0]], [[10, 0, 0], [11, 0, 0], [10, 1, 0]]], float)
    quantities = []
    for (name, heights) in [('A', sa), ('B', sb)]:
        outer = inner.copy()
        outer[:, :, 2] = heights[:, None]
        mesh = T.physical_mesh(outer)
        (_, dist, _) = trimesh.proximity.closest_point_naive(mesh, inner.mean(1))
        measured = float(dist.min())
        test('sufficiency_surface_' + name, measured == float(heights.min()))
        quantities.append(measured)
        np.savez(run / ('sufficiency_' + name + '.npz'), intaglio_triangles_mm=inner, exterior_triangles_mm=outer)
    suff['external_published_code_minima_mm'] = quantities
    suff['external_quantity_difference_mm'] = abs(quantities[0] - quantities[1])
    test('sufficiency_external_difference', suff['external_quantity_difference_mm'] == downstream)
    test('exact_summary_identity', identity == f['sufficiency']['expected_summary_error'])
    test('summary_does_not_determine_minimum', downstream == f['sufficiency']['expected_min_difference_mm'])
    if 'convex_cup' in reports:
        r = reports['convex_cup']['rules']['insertion']
        cone = r.get('reused_finite_cone', {})
        test('cone_used_on_bound_facets', cone.get('status') == 'YES_STRICT')
        pth = Path(f['cases']['convex_cup'])
        meta = json.loads((pth / 'contract.json').read_text())
        data = G.load(pth / 'crown.stl', 1.0)
        (planes, vertices) = cavity_planes(data, meta['regions']['intaglio'])
        norm = [[float(__import__('fractions').Fraction(x)) for x in p[:3]] for p in planes]
        impossible = inherited.insertion_cone(np.r_[np.eye(3), -np.eye(3)])
        test('farkas_six_chart_refutation', inherited.verify_cone_negative(np.r_[np.eye(3), -np.eye(3)], impossible))
        bad = copy.deepcopy(impossible)
        if bad.get('certificates'):
            bad['certificates'][0]['b'][0] = '1000000'
        test('corrupt_farkas_rejected', not inherited.verify_cone_negative(np.r_[np.eye(3), -np.eye(3)], bad))
        mr = reports['convex_cup']['rules']['ball_milling']
        if mr.get('bur_reports') and mr['bur_reports'][0]['status'] == 'PASS':
            rr = mr['bur_reports'][0]
            cert = rr['certificate']
            wrong = copy.deepcopy(cert['centres'])
            wrong[0] = ['100', '100', '100']
            test('corrupt_positive_bur_rejected', not milling.verify_positive(planes, vertices, rr['radius_final_mm'], 0.3, cert['direction'], wrong))
        neg = reports.get('oversized_bur', {}).get('rules', {}).get('ball_milling', {})
        if neg.get('bur_reports') and neg['bur_reports'][0]['status'] == 'FAIL':
            rr = neg['bur_reports'][0]
            w = rr['witness']
            test('corrupt_negative_bur_rejected', not milling.verify_negative(planes, w['point'], rr['radius_final_mm'], 0.3, ['0'] * len(planes)))
    x1b = []
    actual_exports = []
    for row in f['x1b']:
        try:
            r = check(row['paths'], units='mm', sinter_factor=1.2)
            report_path = run / ('X1b_' + row['case'] + '.json')
            T.dump(report_path, r)
            statuses = {k: v['status'] for (k, v) in r['rules'].items()}
            x1b.append(dict(case=row['case'], rules=statuses, report=str(report_path), antagonist=row['registered_antagonist']))
            test('X1b_' + row['case'] + ':missing_frame_not_pass', all((statuses[k] == 'UNKNOWN' for k in ('material_wall', 'occlusal_contact'))))
            if f['round'] == 'R3':
                out = Path(_release_expand('@DENTAL_WORK_ROOT@/X49_design_gate/exports')) / tag / row['case']
                ex = export(row['paths'], r, out, None, 1.2)
                actual_exports.append(dict(case=row['case'], **ex))
                test('X1b_' + row['case'] + ':crown_transport', ex['status'] == 'PASS')
                if ex['status'] == 'PASS':
                    receipt = ex['receipt']
                    asset = Path(receipt['asset'])
                    pub = trimesh.load(asset, force='scene')
                    tri = np.concatenate([m.triangles for m in pub.geometry.values()])
                    (original, _) = T.load_mesh(row['paths']['crown'], 'mm')
                    (_, err) = T.correspondence(original, tri)
                    test('X1b_' + row['case'] + ':published_3MF', err <= 1e-05, parity_mm=err)
        except Exception as e:
            test('X1b_' + row['case'] + ':exception', False, error=repr(e))
    exports = []
    for name in ('valid', 'thin', 'convex_cup'):
        if name not in reports:
            continue
        base = Path(f['cases'][name])
        paths = {k: str(base / (k + '.stl')) for k in ('prep', 'crown', 'antagonist')}
        try:
            bad_report = copy.deepcopy(reports[name])
            bad_report['inputs']['crown']['sha256'] = '0' * 64
            rejected = False
            try:
                export(paths, bad_report, run / ('wrong_source_' + name), base / 'contract.json', 1.2)
            except ValueError:
                rejected = True
            test(name + ':wrong_report_crown_rejected', rejected)
            r = export(paths, reports[name], run / ('export_' + name), base / 'contract.json', 1.2)
            exports.append(dict(case=name, **r))
            test(name + ':export_transport', r['status'] == 'PASS')
            if r['status'] != 'PASS':
                continue
            receipt = r['receipt']
            asset = Path(receipt['asset'])
            side = asset.parent / 'model.3mf.json'
            pub = trimesh.load(asset, force='scene')
            meshes = list(pub.geometry.values())
            tri = np.concatenate([m.triangles for m in meshes])
            (original, _) = T.load_mesh(paths['crown'], 'mm')
            (_, err) = T.correspondence(original, tri)
            test(name + ':published_3mf_geometry', err <= 1e-05, maximum_error_mm=err)
            for (fault, source, report, factor, h) in [('source', '0' * 64, receipt['design_report_sha256'], 1.2, receipt['sidecar_sha256']), ('report', receipt['source_sha256'], '0' * 64, 1.2, receipt['sidecar_sha256']), ('sinter', receipt['source_sha256'], receipt['design_report_sha256'], 2.0, receipt['sidecar_sha256']), ('sidecar', receipt['source_sha256'], receipt['design_report_sha256'], 1.2, '0' * 64)]:
                rejected = False
                try:
                    check_native(asset, side, h, source, report, factor)
                except ValueError:
                    rejected = True
                test(name + ':corrupt_export_' + fault, rejected)
        except Exception as e:
            test(name + ':export_exception', False, error=repr(e))
    counts = {s: sum((v['status'] == s for r in reports.values() for v in r['rules'].values())) for s in ('PASS', 'FAIL', 'UNKNOWN')}
    allpassed = all((t['passed'] for t in tests))
    res = dict(round=f['round'], claim_type='capability', outcome='SCOPED_CAPABILITY_VERIFIED' if allpassed else 'FROZEN_GATE_FAILED', run_directory=str(run), tests=tests, passed_tests=sum((t['passed'] for t in tests)), total_tests=len(tests), cases=rows, x1b=x1b, exports=exports, x1b_exports=actual_exports, sufficiency=suff, rule_counts=counts, unknown_fraction=counts['UNKNOWN'] / sum(counts.values()) if sum(counts.values()) else None, exceptions=exceptions, external_referent=json.loads((HERE / 'PREREG_R1.json').read_text())['external_referent'], wall_s=time.perf_counter() - start, cpu_s=time.process_time() - cpu, peak_RSS_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, physical_measurements=0, review_state='PENDING_INDEPENDENT_REVIEW')
    T.dump(run / 'results.json', res)
    T.dump(HERE / ('rounds/' + f['round'] + '.json'), res)
    T.dump(HERE / 'CURRENT_WORK_STATE.json', dict(lane='X49-design-gate', status='R1_FINISHED', latest_gate=res['outcome'], next_operation='Diagnose all frozen failures; preregister changed construction, retain R1', updated_utc=dt.datetime.now(dt.timezone.utc).isoformat()))
    print(json.dumps(dict(outcome=res['outcome'], tests=f"{res['passed_tests']}/{res['total_tests']}", run=str(run))))
    return 0 if allpassed else 1
if __name__ == '__main__':
    sys.exit(main())
