import thread_guard
from cadlib import *
from margin import detect
from round2 import detect_cut
import trimesh

def run(all_cases=False):
    old1 = read(ROOT / 'raw/R1_PREDICTIONS.json')
    old2 = {r['key']: r for r in read(ROOT / 'raw/R2_PREDICTIONS.json')}
    selected = [r for r in old1 if r['status'] != 'UNAVAILABLE_SOURCE']
    if not all_cases:
        selected = [selected[0]] + [r for r in selected if r['key'] == 'D1_baseline']
    rows = []
    st = time.perf_counter()
    for r in selected:
        if sha(r['source']) != r['source_sha256']:
            raise ValueError('Margin source drift')
        if Path(r['source']).suffix == '.npz':
            a = np.load(r['source'])
            m = trimesh.Trimesh(a['vertices'], a['faces'], process=False)
        else:
            m = trimesh.load_mesh(r['source'], process=True)
        (a, b) = (detect(m), detect_cut(m))
        o = old2[r['key']]
        ca = a['status'] == r['status'] and (a['status'] != 'PROPOSED' or np.array_equal(a['points'], r['points']))
        cb = b['status'] == o['status'] and (b['status'] != 'PROPOSED' or np.array_equal(b['points'], o['points']))
        if not ca or not cb:
            raise AssertionError('Margin regeneration changed ' + r['key'])
        rows.append(dict(key=r['key'], R1_identical=bool(ca), R2_identical=bool(cb)))
    print('MARGIN_REPLAY', len(rows), 'PASS', flush=True)
    out = dict(requested=len(rows), all_identical=True, rows=rows, seconds=time.perf_counter() - st)
    dump(ROOT / 'raw/MARGIN_REPLAY_ALL.json' if all_cases else ROOT / 'raw/MARGIN_REPLAY_DEMO.json', out)
    return out
if __name__ == '__main__':
    run('--all' in sys.argv)
