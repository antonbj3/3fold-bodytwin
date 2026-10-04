from common import *
import subprocess

def run():
    lock = read(ROOT / 'EXTERNAL_GENERATOR_R3.json')
    model = read(ROOT / 'raw/EXTERNAL_MODEL_LOCK.json')
    for rec in model['files']:
        if sha(rec['path']) != rec['sha256']:
            raise ValueError('Model drift ' + rec['path'])
    if sha(ROOT / 'gencad_bench_v5/toothcraft_plugin.py') != lock['plugin_sha256'] or sha(ROOT / 'gencad_bench_v5/external_worker.py') != lock['worker_sha256']:
        raise ValueError('Generator drift')
    dest = DATA / 'external_replay'
    dest.mkdir(exist_ok=True)
    cmd = lock['argv'].copy()
    old = str(DATA / 'external_predictions')
    cmd = [str(dest) if x == old else x for x in cmd]
    tic = time.perf_counter()
    with (ROOT / 'raw/external_replay.log').open('w') as f:
        r = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT)
    if r.returncode:
        raise RuntimeError('Replay failed; see raw/external_replay.log')
    frozen = read(ROOT / 'FROZEN_PREDICTIONS.json')
    rows = []
    for (name, info) in frozen['files'].items():
        if not name.endswith('.npz'):
            continue
        a = npz(dest / name)['tsdf']
        b = npz(DATA / 'external_predictions' / name)['tsdf']
        rows.append(dict(file=name, byte_identical=sha(dest / name) == info['sha256'], max_absolute_difference=float(np.max(abs(a - b)))))
    out = dict(status='PASS' if all((r['byte_identical'] for r in rows)) else 'FAIL', rows=rows, seconds=time.perf_counter() - tic, receipt=read(dest / 'RECEIPT.json'))
    dump(ROOT / 'raw/EXTERNAL_REPLAY.json', out)
    assert out['status'] == 'PASS'
    print(out)
if __name__ == '__main__':
    run()
