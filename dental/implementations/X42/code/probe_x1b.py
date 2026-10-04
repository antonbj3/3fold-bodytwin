from common import *
import collections, resource
rows = []
tic = time.perf_counter()
for orig in tasks(split()['dev_round1'][0]):
    if orig['status'] != 'READY':
        continue
    t = scene(orig)
    t['requirements'] = dict(t['requirements'], tool_radius_mm=[0.25, 0.5, 0.8, 1.0][['easy', 'normal', 'hard', 'boundary'].index(t['level'])])
    s = time.perf_counter()
    try:
        d = crown_loop(t)
        ch = all_checks(t, d)
        rows.append(dict(task_id=t['task_id'], radius_mm=t['requirements']['tool_radius_mm'], seconds=time.perf_counter() - s, verdict=ch['validity']))
    except Exception as e:
        rows.append(dict(task_id=t['task_id'], seconds=time.perf_counter() - s, error=repr(e)))
    print(rows[-1], flush=True)
dump(ROOT / 'raw/X1B_ADAPTER_PROBE.json', dict(rows=rows, seconds=time.perf_counter() - tic, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024))
