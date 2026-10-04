import thread_guard
from cadlib import *
from cad_tool import generate_case
import resource

def run(all_cases=False):
    t0 = time.perf_counter()
    rows = []
    inputs = {r['key']: r for r in read(ROOT / 'INPUT_LOCK_R3.json')['rows']}
    profiles = {r['id']: r for r in read(ROOT / 'MATERIAL_PROFILES.json')['profiles']}
    selected = [r for r in read(ROOT / 'raw/R5_PREDICTIONS.json') if r['status'] == 'GENERATED']
    if not all_cases:
        selected = selected[:2]
    for old in selected:
        new = generate_case(inputs[old['key']], profiles[old['material']])
        a = np.load(old['mesh_path'])
        m = new['mesh']
        err = float(np.max(abs(a['vertices'] - m.vertices))) if a['vertices'].shape == m.vertices.shape else None
        equal = err == 0 and np.array_equal(a['faces'], m.faces) and np.array_equal(a['roles'], new['roles'])
        rows.append(dict(uid=old['uid'], coordinate_error_mm=err, faces_and_regions_identical=bool(equal)))
        if not equal:
            raise AssertionError('Regeneration differs: ' + old['uid'])
        print('REGENERATED', old['uid'], flush=True)
    result = dict(requested=len(selected), rows=rows, all_identical=True, seconds=time.perf_counter() - t0, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, scope='Fresh margin/axis/gap/mirror/roof/neighbour/ring operations from frozen public inputs; array equality to pre-score predictions; no target loaded during regeneration')
    dump(ROOT / 'raw/REGENERATION_ALL.json' if all_cases else ROOT / 'raw/REGENERATION_DEMO.json', result)
    return result
if __name__ == '__main__':
    run('--all' in sys.argv)
