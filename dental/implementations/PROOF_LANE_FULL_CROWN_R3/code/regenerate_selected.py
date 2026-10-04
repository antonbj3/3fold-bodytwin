"""Actually rerun each selected generator in a fresh isolated output directory."""
from common_r3 import *
import subprocess, time

def run():
    start = time.perf_counter()
    rows = []
    for selected in read(ROOT / 'LAB_EXPORTS.json')['selected']:
        if 'key' not in selected:
            continue
        tag = selected['round']
        key = selected['key']
        method = selected['participant']
        out = DATA / 'demo_generation' / selected['family']
        out.mkdir(parents=True, exist_ok=True)
        if tag == 'A':
            cmd = read(ROOT / 'FROZEN_GENERATOR_A.json')['argv'][:-1] + ['/runner/resume_worker.py', key, method]
            cmd[cmd.index(str(DATA / 'A'))] = str(out)
        else:
            cmd = read(ROOT / f'FROZEN_GENERATOR_{tag}.json')['argv_prefix'] + [key]
            cmd[cmd.index(str(DATA / tag))] = str(out)
        full = ['/usr/bin/time', '-v', '-o', str(ROOT / 'raw' / f"REGENERATE_{selected['family']}_MEMORY.txt")] + cmd
        dump(ROOT / 'raw' / f"REGENERATE_CMD_{selected['family']}.json", dict(argv=full))
        with (ROOT / 'raw' / f"REGENERATE_{selected['family']}.log").open('w') as f:
            p = subprocess.run(full, stdout=f, stderr=subprocess.STDOUT)
        if p.returncode:
            raise RuntimeError('selected generator failed ' + key)
        a = npz(DATA / tag / method / key / 'mesh.npz')
        b = npz(out / method / key / 'mesh.npz')
        parity = {k: np.array_equal(a[k], b[k]) for k in ['vertices', 'faces', 'face_roles']}
        assert all(parity.values()), key
        rec = read(out / (key + '__' + method + '.json' if tag == 'A' else key + '.json'))
        rows.append(dict(round=tag, key=key, family=selected['family'], parity=parity, peak_worker_rss_MiB=rec['peak_rss_MiB'], seconds=rec['seconds']))
    result = dict(status='PASS', rows=rows, seconds=time.perf_counter() - start, source_boundaries='Fresh bwrap with same frozen inputs and no private target mount; outputs in own separate demo_generation directory')
    dump(ROOT / 'raw/GENERATOR_REPLAY.json', result)
    return result
if __name__ == '__main__':
    run()
