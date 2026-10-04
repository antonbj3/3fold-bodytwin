from common_r3 import *
import time
sys.path.insert(0, str(V6 / 'code'))
from contact import compare

def run():
    start = time.perf_counter()
    rows = []
    for rec in read(DATA / 'inputs/RECORDS.json'):
        p = npz(DATA / 'inputs' / (rec['key'] + '.npz'))
        q = npz(DATA / 'inputs_C' / (rec['key'] + '.npz'))
        g = p['contact_reference_gap_unclipped']
        xy = p['contact_xy']
        f = p['contact_faces']
        a = compare(xy, f, p['contact_gap_mm'], g)
        c = compare(xy, f, q['contact_gap_mm'], g)
        rows.append(dict(key=rec['key'], family=rec['family'], A_B_clipped_source_symdiff_mm2=a['symdiff_mm2'], C_crossing_source_symdiff_mm2=c['symdiff_mm2'], C_exact_tolerance_pass=c['symdiff_mm2'] <= 1e-07, additional_nodes=int(np.count_nonzero(np.nan_to_num(p['contact_gap_mm']) != np.nan_to_num(q['contact_gap_mm']))), finite_native_nodes=int(np.isfinite(g).sum()), resolution='PER_SURFACE_REGION; acquired nodes PER_POINT'))
    result = dict(rows=rows, seconds=time.perf_counter() - start, clipped_fail_exact_count=sum((r['A_B_clipped_source_symdiff_mm2'] > 1e-07 for r in rows)), augmented_fail_exact_count=sum((not r['C_exact_tolerance_pass'] for r in rows)), scope='acquisition-level PL gap only; separate from actual exported mesh scores and source physical uncertainty')
    dump(ROOT / 'raw/CLIPPING_DATA_DIAGNOSIS.json', result)
    return result
if __name__ == '__main__':
    run()
