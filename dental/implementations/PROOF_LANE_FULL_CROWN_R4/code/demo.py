"""One-command reproduction of three selected generators and held-surface scores."""
from score import *
import argparse

def replay():
    st = time.perf_counter()
    rows = []
    inputs = read(ROOT / 'FROZEN_INPUTS_B.json')
    inp = {r['key']: r for r in inputs['records']}
    frozen = read(ROOT / 'FROZEN_EXPORTS.json')
    original = {r['key']: r for r in read(ROOT / 'rounds/C.json')['rows'] if r['method'] == 'exact_margin'}
    for r in inputs['records']:
        for kind in ['public', 'private']:
            assert sha(r[kind + '_path']) == r[kind + '_sha256']
    for r in inputs['templates']:
        if r['status'] == 'AVAILABLE':
            assert sha(r['path']) == r['sha256']
    from generate import run as unused
    generator = module('demo_generate', ROOT / 'code/generate.py')
    generator.DATA = DATA / 'demo_replay'
    generator.dump = lambda p, x: dump(ROOT / 'raw' / ('replay_' + Path(p).name), x)
    generator.state = lambda *a: None
    from generate_c import exact_trace, exact_shell
    for chosen in frozen['exports']:
        key = chosen['key']
        generator.run(only=key)
        br = next((r for r in read(ROOT / 'raw/replay_B_GENERATION.json') if r['method'] == 'rigid_control'))
        raw = npz(br['mesh_path'])
        ext = mesh(raw['vertices'][raw['faces'][raw['roles'] == 0]], True)
        p = npz(inp[key]['public_path'])
        (e, trace) = exact_trace(ext, p['margin_curve'])
        (m, roles) = exact_shell(e, p)
        expected = npz(chosen['mesh_path'])
        parity = dict(vertices=np.array_equal(m.vertices, expected['vertices']), faces=np.array_equal(m.faces, expected['faces']), roles=np.array_equal(roles, expected['roles']))
        assert all(parity.values()), (key, parity)
        target = npz(inp[key]['private_path'])['target']
        et = m.triangles[roles == 0]
        it = m.triangles[roles == 1]
        metric = metrics(et, target)
        ww = wall(et, it)
        (x, f) = grid(target)
        ceil = height(p['antagonist'], x, True)
        cm = contact.compare(x, f, ceil - height(et, x), ceil - height(target, x))
        old = original[key]
        discrepancy = max(abs(metric['p95_mm'] - old['in_situ']['p95_mm']), abs(ww['sampled_min_mm'] - old['wall']['sampled_min_mm']), abs(ww['continuous_lower_mm'] - old['wall']['continuous_lower_mm']))
        assert discrepancy <= 1e-10
        if cm['status'] == 'SCORED':
            assert abs(cm['symdiff_mm2'] - old['contact']['symdiff_mm2']) <= 1e-10
        rows.append(dict(key=key, family=chosen['family'], geometry_bitwise_parity=parity, metric_error_mm=discrepancy, fresh_p95_mm=metric['p95_mm'], fresh_wall=ww, fresh_v6_contact=cm, exact_trace=trace))
    result = dict(requested=3, reproduced=len(rows), rows=rows, seconds=time.perf_counter() - st, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, scope='Freshly regenerate all three selected C exact-margin exports from public collar/donor data; repeated nearest-template control generation also runs. Re-evaluate hidden target, wall and GenCAD v6 fields. Other cohort rows are preserved prior measurements unless --full is requested.')
    dump(ROOT / 'raw/DEMO_REPLAY.json', result)
    return result

def run(full=False):
    if (ROOT / 'RELEASE_CODE.json').exists():
        for (path, digest) in read(ROOT / 'RELEASE_CODE.json')['files'].items():
            assert sha(path) == digest, ('Code/dependency hash drift', path)
    from validate import run as validate
    v = validate()
    r = replay()
    if full:
        for tag in ['B', 'C', 'D2']:
            evaluate(tag)
    from report import run as report
    report()
    print('Verified:3/3 selected crowns regenerated identically; actual fault injections rejected. Full physical crown capability remains unestablished.', flush=True)
if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--full', action='store_true')
    args = p.parse_args()
    run(args.full)
