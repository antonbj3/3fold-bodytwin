from dental_release.paths import expand as _release_expand
import thread_guard
from cadlib import *
from readonly import parents, force_refusal, milling
from pipeline import create
import pipeline
from pipeline_r4 import canonical_prep, offset_cap
pipeline.offset_cap = offset_cap
from margin import errors
import three_mf
import trimesh, resource
P = parents()
V2 = DATA.parent / _release_expand('GENCAD_V2')

def prepare():
    st = time.perf_counter()
    rows = []
    cohort = {r['case_key']: r for r in read(BASE / 'PROOF_LANE_GENCAD_V2/data/COHORT_LOCK.json')['payload']['cases']}
    refs = {r['task_id']: r['reference'] for r in read(V2 / 'private/references.json')}
    cached = None
    last = None
    for path in sorted((V2 / 'public').glob('*_crown_normal.json')):
        if not any((x in path.name for x in ['molar_crown', 'premolar_crown'])):
            continue
        t = P['tasks'].load_task(path)
        key = t['task_id']
        die = V2 / 'lab_pairs' / key / 'preparation_die.stl'
        row = dict(key=key, case_key=t['case_key'], family=t['family'], task_path=path, die_path=die, status='AVAILABLE' if die.exists() else 'MISSING_DIE')
        if die.exists():
            ref = np.load(V2 / 'private' / refs[key]['file'])
            (base, R) = (ref['base'], ref['R'])
            if t['case_key'] != last:
                (cached, src) = P['build'].pair(cohort[t['case_key']])
                last = t['case_key']
            lower = cached['lower']
            v = (lower['v'] - base) @ R
            f = lower['f']
            owner = lower['owner']
            opposite = t['source_fdi'] + 10
            donor = v[f[owner == opposite]]
            neighbors = v[f[(owner != t['source_fdi']) & (owner != 0)]]
            lim = max(np.linalg.norm(t['xy'], axis=1)) + 4
            keep = np.linalg.norm(neighbors[:, :, :2].mean(1), axis=1) < lim
            neighbors = neighbors[keep]
            p = DATA / 'R4B_inputs' / (key + '.npz')
            p.parent.mkdir(exist_ok=True)
            np.savez_compressed(p, donor=donor, neighbors=neighbors)
            row.update(public_geometry=p, public_geometry_sha256=sha(p), die_sha256=sha(die), donor_fdi=opposite, donor_faces=len(donor), neighbor_faces=len(neighbors), source_members=src, reference_path=V2 / 'private' / refs[key]['file'], reference_sha256=sha(V2 / 'private' / refs[key]['file']))
        rows.append(row)
        print('INPUT', key, row['status'], flush=True)
    freeze(ROOT / 'INPUT_LOCK_R4B.json', dict(rows=rows, seconds=time.perf_counter() - st, scope='Public donor/neighbour extraction only; base/R metadata copied from original site frame, target triangles not passed to generator'))
    return rows

def predict():
    st = time.perf_counter()
    inputs = read(ROOT / 'INPUT_LOCK_R3.json')['rows']
    out = []
    profiles = read(ROOT / 'MATERIAL_PROFILES.json')['profiles']
    for r in inputs:
        for mat in profiles:
            t0 = time.perf_counter()
            uid = r['key'] + '__' + mat['id']
            row = dict(uid=uid, key=r['key'], material=mat['id'], family=r['family'], case_key=r['case_key'])
            if r['status'] != 'AVAILABLE':
                out.append(dict(**row, status=r['status']))
                continue
            try:
                assert sha(r['public_geometry']) == r['public_geometry_sha256']
                assert sha(r['die_path']) == r['die_sha256']
                a = np.load(r['public_geometry'])
                t = P['tasks'].load_task(r['task_path'])
                (m, serialization) = canonical_prep(r)
                row['serialization'] = serialization
                if not len(a['donor']):
                    out.append(dict(**row, status='ABSTAIN_NO_CONTRALATERAL'))
                    continue
                d = create(m, t, mat, P, a['donor'])
                row.update({k: v for (k, v) in d.items() if k not in ['mesh', 'roles', 'ov', 'of', 'iv', 'inf', 'outer_roof', 'prior', 'undercut_dot', 'undercut_faces']})
                if d['status'] == 'GENERATED':
                    dest = DATA / 'R4B' / uid
                    dest.mkdir(parents=True, exist_ok=True)
                    npz = dest / 'crown.npz'
                    np.savez_compressed(npz, vertices=d['mesh'].vertices, faces=d['mesh'].faces, roles=d['roles'], ov=d['ov'], of=d['of'], iv=d['iv'], inf=d['inf'], roof=d['outer_roof'], prior=d['prior'], undercut_dot=d['undercut_dot'], undercut_faces=d['undercut_faces'])
                    stl = dest / 'crown.stl'
                    d['mesh'].export(stl)
                    mf = dest / 'crown.3mf'
                    three_mf.write(mf, d['mesh'].vertices, d['mesh'].faces, d['roles'], dict(case_id=r['key'], material=mat['id'], physical_release='UNKNOWN', source_die_sha256=r['die_sha256'], region_resolution='PER_SURFACE_REGION'))
                    row.update(mesh_path=npz, mesh_sha256=sha(npz), stl_path=stl, stl_sha256=sha(stl), three_mf_path=mf, three_mf_sha256=sha(mf))
                row['force_port'] = force_refusal(uid)
            except Exception as exc:
                import traceback
                row.update(status='CONSTRUCTION_ERROR', reason=repr(exc))
                (ROOT / 'raw' / ('R4B_ERROR_' + uid + '.txt')).write_text(traceback.format_exc())
            row['seconds'] = time.perf_counter() - t0
            out.append(row)
            dump(ROOT / 'raw/R4B_PARTIAL.json', out)
            print(uid, row['status'], row.get('reason', ''), flush=True)
        state('R4B_GENERATING', dict(requested_so_far=len(out), generated=sum((r['status'] == 'GENERATED' for r in out))), 'Freeze proposals before reference scoring')
    dump(ROOT / 'raw/R4B_PREDICTIONS.json', out)
    freeze(ROOT / 'FROZEN_PREDICTIONS_R4B.json', dict(prereg_sha256=sha(ROOT / 'PREREG_R4.json'), rows_sha256=sha(ROOT / 'raw/R4B_PREDICTIONS.json'), code={str(p): sha(p) for p in (ROOT / 'code').glob('*.py')}, seconds=time.perf_counter() - st, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, physical_measurement_performed=False, independent_reference_scored=False))
    state('R4B_FROZEN', dict(generated=sum((r['status'] == 'GENERATED' for r in out)), requested=len(out)), 'Independent surface/contact scoring and milling operator')

def milling_contract(row, a):
    stl = trimesh.load_mesh(row['stl_path'], process=False)
    tri = stl.triangles
    v = tri.reshape(-1, 3)
    f = np.arange(len(v)).reshape(-1, 3)
    p = Path(row['mesh_path']).parent / 'milling_geometry.npz'
    np.savez_compressed(p, vertices=v, faces=f)
    template = read(BASE / 'LANE_X89_MILLING_CODESIGN/contracts/crownloop_D1.json')
    c = dict(template)
    c.update(crown_sha256=sha(row['stl_path']), geometry_file=str(p), geometry_sha256=sha(p), region_faces={'intaglio': np.flatnonzero(a['roles'] == 1).tolist(), 'margin': np.flatnonzero(a['roles'] == 2).tolist()}, samples_per_region=2, region_status='Explicit construction roles, transferred through STL facet order', physical_scope='PHENOMENOLOGICAL inherited X89 scenario; not material-specific process qualification')
    c['scenario'] = dict(c['scenario'], n_per_region=2)
    p = Path(row['mesh_path']).parent / 'milling_contract.json'
    dump(p, c)
    return p

def score():
    st = time.perf_counter()
    pr = read(ROOT / 'PREREG_R4.json')
    inputs = {r['key']: r for r in read(ROOT / 'INPUT_LOCK_R3.json')['rows']}
    profiles = {r['id']: r for r in read(ROOT / 'MATERIAL_PROFILES.json')['profiles']}
    refs = np.load(DATA / 'R1_references.npz')
    rows = []
    for row in read(ROOT / 'raw/R4B_PREDICTIONS.json'):
        r = dict(row)
        if row['status'] != 'GENERATED':
            r['digital_cap_status'] = 'NOT_GENERATED'
            r['full_qualified_no_manual'] = 'UNKNOWN'
            rows.append(r)
            continue
        t0 = time.perf_counter()
        inp = inputs[row['key']]
        mat = profiles[row['material']]
        a = np.load(row['mesh_path'])
        t = P['tasks'].load_task(inp['task_path'])
        m = trimesh.Trimesh(a['vertices'], a['faces'], process=False)
        r['margin_error'] = errors(np.array(row['margin']['points']), refs[row['key']])
        r['topology'] = dict(watertight=bool(m.is_watertight), consistent=bool(m.is_winding_consistent), positive_volume=bool(m.volume > 0), volume_mm3=float(m.volume), self_intersection='UNKNOWN_NOT_FULLY_CHECKED')
        r['wall'] = P['cap'].wall_check(a['ov'], a['of'], a['iv'], a['inf'], mat['wall_mm'])
        (prep, serialization) = canonical_prep(inp)
        (_, gap, _) = trimesh.proximity.closest_point(prep, np.r_[a['iv'], a['iv'][a['inf']].mean(1)])
        r['film'] = dict(status='CONDITIONAL_SAMPLED', min_mm=float(gap.min()), max_mm=float(gap.max()), physical='UNKNOWN_X13_UNCALIBRATED', scope='intaglio vertex and facet-centroid nearest complete die; not full continuous film proof')
        pp = P['cap'].cap_planes(prep.vertices, prep.faces[np.asarray(row['cap']['face_indices'])])
        ip = P['cap'].cap_planes(a['iv'], a['inf'])
        from gencad_bench.checks import cement, cement_max
        r['min_film'] = cement.check(prep.vertices.tolist(), ip, 0.04)
        r['max_film'] = dict(status='PASS' if cement_max.verify_positive(pp, a['iv'].tolist(), prep.vertices.tolist(), 0.12) else 'UNKNOWN', scope='convex prep corresponding witness for complete cavity cap')
        ext = a['ov'][a['of']]
        ref = np.load(inp['reference_path'])
        target = ref['target_tri']
        target = target[target.mean(1)[:, 2] >= np.array(row['margin']['points'])[:, 2].min()]
        shape = P['r4'].metrics(ext, target, n=2048)
        r['shape'] = dict(**shape, resolution='PER_SURFACE_REGION', external='Independent original IOS, 2048 area-stratified facet-centroid samples/direction', gate='PASS' if shape['p95_mm'] <= pr['metrics']['shape_p95_mm_by_family'][r['family']] else 'FAIL')
        native = P['v2g'].height(target, t['xy'])
        pred = P['v2g'].height(ext, t['xy'])
        ceiling = t['ceiling']
        cm = P['v6'].compare(t['xy'], t['faces'], ceiling - pred, ceiling - native)
        r['V6_contact'] = cm
        r['V6_gates'] = P['v6'].gates(cm, pr['metrics']['V6'])
        r['V6_uncertainty'] = P['v6'].enclosure(t['xy'], t['faces'], ceiling - pred, ceiling - native, 0.05)
        from continuous import broad_phase, extremum
        O = ext[abs(np.cross(ext[:, 1, :2] - ext[:, 0, :2], ext[:, 2, :2] - ext[:, 0, :2])) > 1e-09]
        U = t['antagonist']
        pairs = broad_phase(U, O)
        occ = extremum(U, O, pairs, True) if len(pairs) else None
        r['occlusion'] = occ
        n = np.load(inp['public_geometry'])['neighbors']
        (q, dd, _) = P['r4'].closest(n, P['r4'].sample(ext, 512)) if len(n) else (None, np.array([np.nan]), None)
        r['neighbor'] = dict(min_sampled_unsigned_mm=float(dd.min()), resolution='PER_POINT', physical_collision='UNKNOWN_UNSIGNED_SAMPLE', source='same-subject lower arch excludes target label')
        (rv, rf, rr) = three_mf.readback(row['three_mf_path'])
        r['export'] = dict(coordinate_max_error_mm=float(np.abs(rv - a['vertices']).max()), faces_identical=bool(np.array_equal(rf, a['faces'])), regions_identical=bool(np.array_equal(rr, a['roles'])), unit='millimeter')
        mp = milling_contract(row, a)
        r['milling'] = milling(row['stl_path'], mp)
        dump(Path(row['mesh_path']).parent / 'milling_result.json', r['milling'])
        r['milling'] = {k: v for (k, v) in r['milling'].items() if k != 'records'}
        geom = r['margin_error']['max_um'] <= 25 and all((r['topology'][k] for k in ['watertight', 'consistent', 'positive_volume'])) and (r['wall']['status'] == 'PASS') and (r['min_film']['status'] == 'PASS') and (r['max_film']['status'] == 'PASS') and occ and (occ['minimum_gap_mm'] >= -1e-06) and all(r['V6_gates'].values()) and r['export']['regions_identical']
        r['digital_cap_status'] = 'PASS' if geom else 'FAIL'
        r['full_qualified_no_manual'] = 'UNKNOWN'
        r['full_qualification_missing'] = ['annotated clinical margin', 'calibrated spatial cement gap', 'same-specimen regional force response', 'complete CAM fixture/path and mounted-tool measurements', 'complete self-intersection/source uncertainty certificate']
        r['score_seconds'] = time.perf_counter() - t0
        rows.append(r)
        dump(ROOT / 'raw/R4B_SCORED_PARTIAL.json', rows)
        print('SCORE', row['uid'], 'wall', r['wall']['status'], 'shape', r['shape']['gate'], 'digital', r['digital_cap_status'], 'mill', r['milling']['found'], flush=True)
    dump(ROOT / 'raw/R4B_SCORED.json', rows)
    result = dict(claim_type='capability', requested=len(rows), generated=sum((r['status'] == 'GENERATED' for r in rows)), digital_cap_pass=sum((r.get('digital_cap_status') == 'PASS' for r in rows)), full_qualified_no_manual_pass=0, full_qualified_status='UNKNOWN', external_referent=pr['external_referent'], seconds=time.perf_counter() - st, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024)
    dump(ROOT / 'RESULTS_R4B.json', result)
    state('R4B_SCORED', result, 'Controls, independent rerun and next load-bearing construction')
    print(result)
if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('phase', choices=['predict', 'score'])
    args = ap.parse_args()
    if args.phase == 'predict':
        predict()
    else:
        score()
