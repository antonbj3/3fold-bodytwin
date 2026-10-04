from common_r3 import *
import subprocess, time

def run():
    assert sha(ROOT / 'code/generate_b.py') == read(ROOT / 'FROZEN_GENERATOR_B.json')['files']['generate_b.py'], 'B geometry source changed'
    if not (ROOT / 'FROZEN_ACQUISITION_C.json').exists():
        subprocess.run([PYTHON, str(ROOT / 'code/acquire_c.py')], check=True)
    out = DATA / 'C'
    out.mkdir(exist_ok=True)
    old = read(ROOT / 'FROZEN_GENERATOR_A.json')
    cmd = old['argv'].copy()
    cmd[cmd.index(str(DATA / 'A'))] = str(out)
    cmd[cmd.index(str(DATA / 'inputs'))] = str(DATA / 'inputs_C')
    cmd[-1] = '/runner/generate_b.py'
    start = time.perf_counter()
    if not (ROOT / 'FROZEN_GENERATOR_C.json').exists():
        freeze(ROOT / 'FROZEN_GENERATOR_C.json', dict(argv_prefix=cmd, files={p.name: sha(p) for p in [ROOT / 'code/generate_b.py', ROOT / 'code/memory_height.py', ROOT / 'code/run_c.py', ROOT / 'code/acquire_c.py']}, prereg_sha256=sha(ROOT / 'PREREG_C.json'), acquisition_sha256=sha(ROOT / 'FROZEN_ACQUISITION_C.json'), whole_target_mounted=False, geometry_operator_identical_to_B=True))
    rows = []
    for rec in read(DATA / 'inputs_C/RECORDS.json'):
        key = rec['key']
        dest = out / (key + '.json')
        if not dest.exists():
            full = ['/usr/bin/time', '-v', '-o', str(ROOT / 'raw' / f'C_GENERATE_MEMORY_{key}.txt')] + cmd + [key]
            dump(ROOT / 'raw' / f'C_GENERATE_CMD_{key}.json', dict(argv=full))
            with (ROOT / 'raw' / f'C_generate_{key}.log').open('w') as f:
                p = subprocess.run(full, stdout=f, stderr=subprocess.STDOUT)
            if p.returncode:
                raise RuntimeError('C worker failed ' + key)
        row = read(dest)
        rows.append(row)
        dump(out / 'RECORDS.json', rows)
        print('C', key, row['status'], row.get('reason'), flush=True)
        state('C_GENERATING', str(len(rows)) + '/18 processed', 'Preserve source band-crossing information')
    if not (ROOT / 'FROZEN_PREDICTIONS_C.json').exists():
        freeze(ROOT / 'FROZEN_PREDICTIONS_C.json', dict(files={str(p.relative_to(out)): dict(sha256=sha(p), bytes=p.stat().st_size) for p in out.rglob('*') if p.is_file()}, generator_sha256=sha(ROOT / 'FROZEN_GENERATOR_C.json'), target_mesh_queries_by_generator=0, physical_measurement='NOT_RUN'))
    dump(ROOT / 'raw/GENERATE_C_COST.json', dict(seconds=time.perf_counter() - start, peak_rss_MiB=max((r['peak_rss_MiB'] for r in rows))))
    state('C_FROZEN', '18 constructions recorded', 'Score full anatomy and original function')
if __name__ == '__main__':
    run()
