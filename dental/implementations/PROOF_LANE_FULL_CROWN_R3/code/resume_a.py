from common_r3 import *
import subprocess, time
old = read(ROOT / 'FROZEN_GENERATOR_A.json')
for (name, h) in old['files'].items():
    assert sha(ROOT / 'code' / name) == h, name
for (name, r) in read(ROOT / 'FROZEN_ACQUISITION_A.json')['files'].items():
    assert sha(DATA / 'inputs' / name) == r['sha256'], name
out = DATA / 'A'
if not (ROOT / 'FROZEN_MEMORY_REPAIR.json').exists():
    freeze(ROOT / 'FROZEN_MEMORY_REPAIR.json', dict(parent_generator=sha(ROOT / 'FROZEN_GENERATOR_A.json'), change='Only chunk query_ball_point into16 queries; identical affine arithmetic and candidate order. Fresh process per tooth and method. No acquired inputs, thresholds or original generator modified.', files={p.name: sha(p) for p in [ROOT / 'code/memory_height.py', ROOT / 'code/resume_worker.py', ROOT / 'code/resume_a.py']}, memory_limit_MiB=3584, precision='Inherited small geometry and affine calculations kept float64 to preserve frozen numerics; volumetric fields float32; no bulk float64 case stacks.', already_completed={str(p.relative_to(out)): dict(sha256=sha(p), bytes=p.stat().st_size) for p in out.rglob('*') if p.is_file()}, previous_failure='User reports OOM at6GB; log ends after two first-molar outputs. Original complete runtime/RSS unavailable.'))
rows = read(out / 'RECORDS.json')
done = {(r['key'], r['participant']) for r in rows}
start = time.perf_counter()
for r in read(DATA / 'inputs/RECORDS.json'):
    for method in ['spatial_sheet', 'scalar_control']:
        if (r['key'], method) in done:
            continue
        ident = r['key'] + '__' + method
        cmd = old['argv'][:-1] + ['/runner/resume_worker.py', r['key'], method]
        timer = ROOT / 'raw' / ('MEMORY_' + ident + '.txt')
        full = ['/usr/bin/time', '-v', '-o', str(timer)] + cmd
        state('A_RESUMING', str(len(rows)) + '/36 completed; one process per tooth', 'Generate ' + ident)
        with (ROOT / 'raw' / ('resume_' + ident + '.log')).open('w') as f:
            proc = subprocess.run(full, stdout=f, stderr=subprocess.STDOUT)
        with (ROOT / 'COMMANDS.md').open('a') as f:
            f.write('\n- Resumed `' + ident + '`: `/usr/bin/time -v -o ' + str(timer) + ' <FROZEN_GENERATOR_A.argv with final script /runner/resume_worker.py> ' + r['key'] + ' ' + method + '`; exit ' + str(proc.returncode) + '.\n')
        if proc.returncode:
            raise RuntimeError('worker failed ' + ident)
        row = read(out / (ident + '.json'))
        rows.append(row)
        dump(out / 'RECORDS.json', rows)
        print(ident, row['status'], row['peak_rss_MiB'], flush=True)
for (rel, x) in read(ROOT / 'FROZEN_MEMORY_REPAIR.json')['already_completed'].items():
    if rel != 'RECORDS.json':
        assert sha(out / rel) == x['sha256'], rel
freeze(ROOT / 'FROZEN_PREDICTIONS_A.json', dict(files={str(p.relative_to(out)): dict(sha256=sha(p), bytes=p.stat().st_size) for p in out.rglob('*') if p.is_file()}, generator_sha256=sha(ROOT / 'FROZEN_GENERATOR_A.json'), memory_repair_sha256=sha(ROOT / 'FROZEN_MEMORY_REPAIR.json'), target_mesh_queries_by_generator=0, physical_measurement='NOT_RUN'))
dump(ROOT / 'raw/RESUME_A_COST.json', dict(seconds=time.perf_counter() - start, peak_rss_MiB=max((r.get('peak_rss_MiB', 0) for r in rows))))
state('A_FROZEN', '36 records and predictions frozen', 'Score unchanged18 sites')
