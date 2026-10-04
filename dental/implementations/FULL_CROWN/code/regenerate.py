"""Optional actual generation replay; source lanes and frozen predictions are never overwritten."""
from fc_common import *
import subprocess, importlib

def run():
    replay = DATA / 'replay'
    replay.mkdir(exist_ok=True)
    rows = []
    st = time.perf_counter()
    for tag in ['R2', 'R3']:
        target = replay / (tag + '_predictions')
        target.mkdir(exist_ok=True)
        cmd = read(ROOT / f'FROZEN_GENERATOR_{tag}.json')['argv']
        cmd = [str(target) if x == str(DATA / (tag + '_predictions')) else x for x in cmd]
        with (ROOT / 'raw' / f'replay_{tag}.log').open('w') as f:
            p = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT)
        if p.returncode:
            raise RuntimeError('generation replay failed ' + tag)
        old = read(DATA / (tag + '_predictions') / 'RECORDS.json')
        new = read(target / 'RECORDS.json')
        assert [(r['key'], r['participant'], r['status']) for r in old] == [(r['key'], r['participant'], r['status']) for r in new]
        for rec in old:
            if rec['status'] != 'EXPORTED':
                continue
            rel = rec['participant'] + '/' + rec['key'] + '/mesh.npz'
            a = npz(DATA / (tag + '_predictions') / rel)
            b = npz(target / rel)
            same = all((np.array_equal(a[k], b[k]) for k in a))
            rows.append(dict(round=tag, key=rec['key'], participant=rec['participant'], identical=same))
            assert same, 'prediction drift ' + rel
    for (tag, module_name) in [('R6', 'prescan_solid'), ('R7', 'prescan_local_holes')]:
        module = importlib.import_module(module_name)
        for rec in read(DATA / (tag + '_predictions') / 'RECORDS.json'):
            p = npz(V4 / 'payload/whole_inputs' / rec['key'] / 'preparation.npz')
            a = npz(V4 / 'payload/whole_private' / rec['key'] / 'reference.npz')
            try:
                (m, roles, die, cert) = module.build(a['source_triangles'], float(p['margin_z']))
                status = 'EXPORTED'
                reason = None
            except Exception as e:
                status = 'FAILED'
                reason = str(e)
            assert status == rec['status'], (tag, rec['key'], reason)
            if status == 'EXPORTED':
                q = npz(DATA / (tag + '_predictions') / rec['participant'] / rec['key'] / 'mesh.npz')
                same = np.array_equal(m.vertices, q['vertices']) and np.array_equal(m.faces, q['faces']) and np.array_equal(roles, q['face_roles'])
                rows.append(dict(round=tag, key=rec['key'], participant=rec['participant'], identical=same))
                assert same
            else:
                rows.append(dict(round=tag, key=rec['key'], status=status, reason=reason, identical=reason == rec['reason']))
                assert reason == rec['reason']
    dump(ROOT / 'raw/REGENERATION_RECEIPT.json', dict(all_pass=all((r['identical'] for r in rows)), rows=rows, seconds=time.perf_counter() - st, scope='Actual public-only generators R2/R3 rerun isolated; full-prescan builds R6/R7 rerun from original local scans; no frozen file replaced'))
if __name__ == '__main__':
    run()
