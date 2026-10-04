from common_r3 import *
import subprocess, time

def run(tag='B'):
    out = DATA / tag
    out.mkdir(exist_ok=True)
    old = read(ROOT / 'FROZEN_GENERATOR_A.json')
    start = time.perf_counter()
    cmd = old['argv'].copy()
    cmd[cmd.index(str(DATA / 'A'))] = str(out)
    cmd[-1] = '/runner/generate_b.py'
    if not (ROOT / 'FROZEN_GENERATOR_B.json').exists():
        freeze(ROOT / 'FROZEN_GENERATOR_B.json', dict(argv_prefix=cmd, files={p.name: sha(p) for p in [ROOT / 'code/generate_b.py', ROOT / 'code/memory_height.py', ROOT / 'code/run_b.py', ROOT / 'code/generate_a.py']}, prereg_sha256=sha(ROOT / 'PREREG_B.json'), acquisition_sha256=sha(ROOT / 'FROZEN_ACQUISITION_A.json'), whole_target_mounted=False))
    rows = []
    for rec in read(DATA / 'inputs/RECORDS.json'):
        key = rec['key']
        dest = out / (key + '.json')
        if not dest.exists():
            full = ['/usr/bin/time', '-v', '-o', str(ROOT / 'raw' / f'B_GENERATE_MEMORY_{key}.txt')] + cmd + [key]
            dump(ROOT / 'raw' / f'B_GENERATE_CMD_{key}.json', dict(argv=full))
            with (ROOT / 'raw' / f'B_generate_{key}.log').open('w') as f:
                r = subprocess.run(full, stdout=f, stderr=subprocess.STDOUT)
            if r.returncode:
                raise RuntimeError('B process failed ' + key)
        row = read(dest)
        rows.append(row)
        dump(out / 'RECORDS.json', rows)
        print(key, row['status'], row.get('reason'), flush=True)
        state('B_GENERATING', str(len(rows)) + '/18 processed', 'Preserve exterior prescription; then freeze all outputs')
    if not (ROOT / 'FROZEN_PREDICTIONS_B.json').exists():
        freeze(ROOT / 'FROZEN_PREDICTIONS_B.json', dict(files={str(p.relative_to(out)): dict(sha256=sha(p), bytes=p.stat().st_size) for p in out.rglob('*') if p.is_file()}, generator_sha256=sha(ROOT / 'FROZEN_GENERATOR_B.json'), target_mesh_queries_by_generator=0, physical_measurement='NOT_RUN'))
    dump(ROOT / 'raw/GENERATE_B_COST.json', dict(seconds=time.perf_counter() - start, peak_rss_MiB=max((r['peak_rss_MiB'] for r in rows))))
    state('B_FROZEN', 'All18 constructions recorded and outputs frozen', 'Score anatomy and function with unchanged gates')
if __name__ == '__main__':
    run()
