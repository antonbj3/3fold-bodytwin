from common import *
import argparse, shutil

def check(name, ok):
    if not ok:
        raise ValueError('FAILED ' + name)

def verify():
    start = time.perf_counter()
    checks = 0
    for p in list(R.glob('PREREG*.json')) + list(R.glob('FROZEN_PREDICTIONS*.json')) + list(R.glob('DECOMPOSITION*.json')):
        check(p.name, sha(p) == p.with_suffix('.sha256').read_text().split()[0])
        checks += 1
    for rec in inputs():
        for name in ['private', 'public']:
            check(rec['key'] + '_' + name, sha(rec[name + '_path']) == rec[name + '_sha256'])
            checks += 1
        check(rec['key'] + '_outer', sha(rec['outer']['mesh_path']) == rec['outer']['mesh_sha256'])
        checks += 1
    man = read(R / 'MANIFEST.json')
    for row in man['files']:
        check(row['path'], sha(row['path']) == row['sha256'])
        checks += 1
    print('Input, frozen protocol and artifact hashes passed', checks, flush=True)
    certs = 0
    for p in sorted((R / 'raw').glob('*_CERT.json')):
        c = read(p)
        if 'CGAL EPECK exact triangle-triangle' not in c.get('backend', ''):
            continue
        check(p.name + '_support', sha(c['support_path']) == c['support_sha256'])
        check(p.name + '_query', sha(c['query_path']) == c['query_sha256'])
        z = subprocess.run(c['command'], capture_output=True, text=True, timeout=240)
        check(p.name + '_exit', z.returncode == 0)
        out = json.loads(z.stdout)
        check(p.name + '_exact_replay', all((c[k] == v for (k, v) in out.items())))
        certs += 1
        if certs % 10 == 0:
            print('Exact distance certificates replayed', certs, flush=True)
    for p in sorted((R / 'raw').glob('*_FILM_UPPER.json')):
        c = read(p)
        z = subprocess.run([str(D / 'exact_gap_upper'), c['support_path'], c['query_path']], capture_output=True, text=True, timeout=240)
        check(p.name + '_exit', z.returncode == 0)
        out = json.loads(z.stdout)
        check(p.name + '_exact_replay', all((c[k] == v for (k, v) in out.items())))
        certs += 1
        print('Continuous film replay', p.stem, flush=True)
    for tag in ['C', 'E']:
        for row in read(R / f'RESULTS_{tag}.json')['rows']:
            if row['status'] != 'GENERATED':
                continue
            m = load(row['mesh_path'])
            n = normals_exact(m['vertices'], m['faces'][m['roles'] != 0][:, ::-1], m['M'][2])
            fpath = R / 'raw' / (Path(row['mesh_path']).stem + '_PATH.json')
            if fpath.exists():
                saved = read(fpath)
                check(row['key'] + '_actual_normals', n == saved['actual_intaglio_direction'])
                proj = saved['projected_boundary']
                z = subprocess.run([str(D / 'exact_projected_loop'), proj['path']], capture_output=True, text=True, timeout=60)
                check(row['key'] + '_projected_loop', z.returncode == 0 and all((proj[k] == v for (k, v) in json.loads(z.stdout).items())))
    ctrl = module('r2_replay_controls', Path(__file__).with_name('controls.py')).run()
    check('fault_controls', ctrl['all_pass'])
    from real_controls import run as real_run
    real = real_run()
    check('actual_geometry_fault_controls', real['all_pass'])
    from evaluate import mf
    for part in read(R / 'EXPORTS.json')['rows']:
        m = load(part['mesh']['path'])
        (vv, ff, rr) = mf.readback(part['three_mf']['path'])
        check(part['key'] + '_3mf', np.array_equal(vv, m['vertices']) and np.array_equal(ff, m['faces']) and np.array_equal(rr, m['roles']))
        stl = trimesh.load_mesh(part['ascii_stl']['path'], file_type='stl', process=False)
        check(part['key'] + '_exact_text_STL', np.array_equal(stl.triangles, m['vertices'][m['faces']]))
        si = intersections(m['vertices'], m['faces'], Path(part['mesh']['path']).stem + '_REPLAY')
        check(part['key'] + '_complete_crown_no_intersection', si.get('count') == 0)
        print('Selected crown mesh and exact exports replayed', part['key'], part['material'], flush=True)
    from carrier_audit import audit
    check('original_exact_carrier_audit', audit() == read(R / 'raw/EXACT_CARRIER_AUDIT.json'))
    from attachment_replay import run as attach_run
    attachments = attach_run()
    check('new_facet_attachments', attachments['status'] == 'PASS')
    from lower_topology import run as lower_run
    expected = read(R / 'RESULTS_LOWER_TOPOLOGY.json')
    lower = lower_run()
    check('lower_native_exact_topology', lower == expected and lower['all_pass'])
    out = dict(status='PASS', hash_checks=checks, exact_certificates_replayed=certs, controls=ctrl['count'] + real['count'] + len(attachments['fault_controls']) + len(lower['controls']), facet_attachments=attachments['attachments_replayed'], native_lower_topology_cases=len(lower['rows']), seconds=time.perf_counter() - start, peak_process_RSS_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024)
    dump(R / 'raw/REPLAY_VALIDATION.json', out)
    print(json.dumps(out, indent=2), flush=True)

def recompute():
    size = sum((p.stat().st_size for p in D.rglob('*') if p.is_file() and (not p.is_symlink())))
    if size > 1.4 * 1024 ** 3:
        raise ValueError('Recompute would threaten3GiB lane intermediate limit')
    dest = D / ('recompute_' + datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%S'))
    dest.mkdir()
    art = dest / 'arrays'
    art.mkdir()
    (dest / 'raw').mkdir()
    for p in list(R.glob('PREREG*')) + list(R.glob('DECOMPOSITION*')):
        shutil.copy2(p, dest / p.name)
    for name in ['exact_mesh_clearance', 'exact_gap_upper', 'exact_projected_loop']:
        (art / name).symlink_to(D / name)
    env = dict(os.environ, DENT_R2_RUN_ROOT=str(dest), DENT_R2_DATA_ROOT=str(art))
    code = Path(__file__).resolve().parent
    for script in ['loop_crown.py', 'height_crown.py', 'offset_crown.py']:
        print('Fresh construction', script, flush=True)
        subprocess.run([sys.executable, '-u', str(code / script)], env=env, check=True)
    checks = []
    for tag in ['C', 'D', 'E']:
        old = read(R / f'FROZEN_PREDICTIONS_{tag}.json')['rows']
        new = read(dest / f'FROZEN_PREDICTIONS_{tag}.json')['rows']
        check(tag + '_same36', len(new) == len(old) == 36)
        for (a, b) in zip(old, new):
            check(a['key'] + '_status', a['status'] == b['status'])
            if a['status'] == 'GENERATED':
                x = load(a['mesh_path'])
                y = load(b['mesh_path'])
                fields = ['vertices', 'faces', 'roles', 'prep_vertices', 'prep_faces']
                same = all((np.array_equal(x[k], y[k]) for k in fields))
                check(a['key'] + '_fresh_geometry', same)
                checks.append(dict(round=tag, key=a['key'], material=a['material'], same_arrays=same))
    out = dict(status='PASS', fresh_directory=str(dest), rows_compared=108, generated_parts=len(checks), checks=checks)
    dump(R / 'raw/FRESH_RECOMPUTE.json', out)
    print(json.dumps(out, indent=2))
if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--recompute', action='store_true')
    a = p.parse_args()
    recompute() if a.recompute else verify()
