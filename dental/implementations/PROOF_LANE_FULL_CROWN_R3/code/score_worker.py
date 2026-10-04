import resource
resource.setrlimit(resource.RLIMIT_AS, (3584 * 1024 ** 2, 3584 * 1024 ** 2))
import sys, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'PROOF_LANE_FULL_CROWN_R2/code'))
import score as inherited
from common_r3 import *
from memory_height import height
inherited.VF.height = height

def run(tag, key, method, destination=None):
    rec = next((r for r in read(DATA / tag / 'RECORDS.json') if r['key'] == key and r['participant'] == method))
    r = dict(rec)
    start = time.perf_counter()
    if r['status'] == 'EXPORTED':
        p = DATA / tag / method / key / 'mesh.npz'
        fr = read(ROOT / f'FROZEN_PREDICTIONS_{tag}.json')
        assert sha(p) == fr['files'][method + '/' + key + '/mesh.npz']['sha256']
        try:
            r.update(inherited.score_row(rec, npz(p), read(ROOT / f'PREREG_{tag}.json')))
        except Exception as e:
            r.update(status='SCORE_FAILED', reason=repr(e))
            import traceback
            traceback.print_exc()
    r.update(score_seconds=time.perf_counter() - start, score_peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024)
    dump(Path(destination) if destination else ROOT / 'raw' / f'{tag}_SCORE_{key}__{method}.json', r)
    print(tag, key, method, r['status'], r.get('reconstruction_p95_mm'), r.get('gates'), flush=True)
if __name__ == '__main__':
    run(*sys.argv[1:])
