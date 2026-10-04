from common import *
from geometry import prepare
import time, resource

def run():
    check_lock()
    DATA.mkdir(parents=True, exist_ok=True)
    start = time.perf_counter()
    (contact, reuse) = parents()
    rows = []
    for r in read(ROOT / 'PREREG_R1.json')['cohort']:
        rows.append(prepare(r, contact, reuse))
        print('acquired', r['key'], rows[-1]['editable_vertices'], rows[-1]['legacy_gap_max_difference_mm'], flush=True)
    dump(ROOT / 'raw/ACQUISITION.json', dict(rows=rows, seconds=time.perf_counter() - start, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024))
    state('R1_GAPS_ACQUIRED', 'All10 source-native projected bands acquired on immutable grid', 'Run regional and ordinary mesh LP; freeze exact output arrays before contact validation')
if __name__ == '__main__':
    run()
