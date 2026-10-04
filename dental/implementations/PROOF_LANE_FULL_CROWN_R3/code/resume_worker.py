import resource
resource.setrlimit(resource.RLIMIT_AS, (3584 * 1024 ** 2, 3584 * 1024 ** 2))
import sys, time, json
from pathlib import Path
import generate_a as gen
from memory_height import height
import numpy as np
gen.height = height
(key, method) = sys.argv[1:3]
rec = next((r for r in json.loads(Path('/inputs/RECORDS.json').read_text()) if r['key'] == key))
out = Path('/output')
row = {k: rec[k] for k in ['key', 'case_key', 'family', 'dataset', 'source_fdi']}
row.update(participant=method, status='FAILED')
start = time.perf_counter()
try:
    with np.load(Path('/inputs') / (key + '.npz')) as z:
        p = {k: z[k] for k in z.files if k != 'contact_reference_gap_unclipped'}
    print('CLOSE', flush=True)
    outer = gen.close_outer(p)
    print('DEFORM', flush=True)
    (outer, di) = gen.deform(p, outer, method)
    print('SHELL', flush=True)
    (v, f, roles, gi) = gen.shell_from_outer(p, outer)
    print('EXPORT', flush=True)
    gen.export(out / method / key, v, f, roles, p, dict(deformation=di, **gi))
    row.update(status='EXPORTED', build=gi, deformation=di)
except Exception as e:
    row['reason'] = repr(e)
    import traceback
    traceback.print_exc()
row.update(seconds=time.perf_counter() - start, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024)
gen.dump(out / (key + '__' + method + '.json'), row)
print(json.dumps(row), flush=True)
