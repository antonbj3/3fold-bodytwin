"""Execute the existing GenCAD control unchanged, on its original information."""
import sys, time
from common import *
FND = DENTAL / 'results/PROOF_LANE_GENCAD_FOUNDATION'
sys.path.insert(0, str(FND))
from gencad_bench.geometry import planes_from_cap, cap
from gencad_bench.checks import milling

def run():
    start = time.perf_counter()
    rows = []
    for r in read(FND / 'rounds/R3/results.json')['rows']:
        (v, f) = cap(r['profile'])
        try:
            control = milling.check(planes_from_cap(r['profile']), v.tolist(), 0.5, 0.05)
        except ValueError as e:
            control = dict(status='UNKNOWN', reason=str(e))
        original = r['milling']['status']
        rows.append(dict(id=r['geometry_id'] + '_' + r['method'], original_status=original, replayed_status=control['status'], matches=original == control['status'], scope='unchanged ideal ball check; not a holder comparator', control=control))
    positive = next((r['control'] for r in rows if r['control']['status'] == 'PASS'))
    r = next((r for r in read(FND / 'rounds/R3/results.json')['rows'] if r['milling']['status'] == 'PASS'))
    planes = planes_from_cap(r['profile'])
    (v, _) = cap(r['profile'])
    cs = positive['certificate']['centres']
    corrupt = [list(x) for x in cs]
    corrupt[0][2] = '1000'
    rejected = not milling.verify_positive(planes, v.tolist(), 0.5, 0.05, [0, 0, -1], corrupt)
    assert rejected
    out = dict(rows=rows, all_match=all((r['matches'] for r in rows)), positive_corruption_rejected=rejected, seconds=time.perf_counter() - start, source=str(FND / 'gencad_bench/checks/milling.py'), source_sha256=sha(FND / 'gencad_bench/checks/milling.py'))
    dump(ROOT / 'raw/PARENT_REPLAY.json', out)
    print('parent replay', out['all_match'], len(rows))
if __name__ == '__main__':
    run()
