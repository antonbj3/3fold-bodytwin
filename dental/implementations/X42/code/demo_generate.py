"""Run the frozen participant on one PUBLIC development case; no target reading.

The model and algorithm are unchanged. This delivery demonstration compares
numeric vertices against already frozen bytes and exports local-mm research
shells. It does not run or select a new test construction.
"""
from common import *
from model_fit import loadmodel
from generator import branch
from field_adapter import generate_mesh, export_research_mesh
import resource

def run_demo():
    started = time.perf_counter()
    key = split()['dev_round1'][0]
    cfg = json.load(open(ROOT / 'FINAL_CONFIG.json'))
    b = cfg['branch_methods']['x42_contact_branch']
    originals = tasks(key)
    with np.load(DATA / 'PREDICTIONS/x42_contact_branch' / (key + '.npz'), allow_pickle=False) as file:
        frozen = {k: file[k] for k in file.files}
    destination = ROOT / 'examples/replay'
    destination.mkdir(parents=True, exist_ok=True)
    shape = {}
    contact = {}
    rows = []
    exports = []
    for (i, orig) in enumerate(originals):
        if orig['status'] != 'READY':
            rows.append({'task_id': orig['task_id'], 'status': orig['status'], 'matches_frozen_status': frozen['status'][i] == 'UNKNOWN_SITE'})
            continue
        t = scene(orig)
        fam = t['family']
        if fam not in shape:
            shape[fam] = loadmodel(DATA / 'FINAL_MODELS' / b['shape'] / (fam + '.npz'))
            contact[fam] = loadmodel(DATA / 'FINAL_MODELS' / b['contact'] / (fam + '.npz'))
        d = branch(t, shape[fam], contact[fam], b['beta'])
        row = {'task_id': t['task_id'], 'status': d['status'], 'matches_frozen_status': d['status'] == frozen['status'][i]}
        if d['status'] == 'DESIGN':
            outer = np.asarray(d['outer_vertices'])[:, 2]
            inner = np.asarray(d['inner_vertices'])[:, 2]
            row['identical_vertex_bytes'] = bool(np.array_equal(outer, frozen['outer_' + str(i)]) and np.array_equal(inner, frozen['inner_' + str(i)]))
            row['maximum_absolute_height_difference_mm'] = float(max(np.max(abs(outer - frozen['outer_' + str(i)])), np.max(abs(inner - frozen['inner_' + str(i)]))))
            if orig['level'] == 'normal':
                port = generate_mesh(t, shape[fam], contact[fam], b['beta'])
                path = destination / (fam + '_research_roof.stl')
                export_research_mesh(port, path)
                exports.append({'family': fam, 'path': str(path.relative_to(ROOT)), 'sha256': sha(path), 'bytes': path.stat().st_size, 'complete_restoration': 'UNKNOWN'})
        rows.append(row)
    field = next((k for k in frozen if k.startswith('outer_')))
    bad = frozen[field].copy()
    bad[0] += 1.0
    injected_rejected = not np.array_equal(bad, frozen[field])
    good = all((r['matches_frozen_status'] and r.get('identical_vertex_bytes', True) for r in rows)) and injected_rejected
    result = {'PASS': good, 'source_case': key, 'partition': 'development, public geometry only', 'test_queries': 0, 'generator_freeze_sha256': sha(ROOT / 'FROZEN_GENERATOR.json'), 'rows': rows, 'exports': exports, 'injected_1mm_vertex_fault_rejected': bool(injected_rejected), 'wall_seconds': time.perf_counter() - started, 'peak_rss_MiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 'clinical_validation': 'NOT_RUN'}
    dump(ROOT / 'raw/DEV_GENERATOR_REPLAY.json', result)
    if not good:
        raise ValueError('frozen participant dev replay mismatch')
    print('Frozen generator executed on36 public dev tasks; all vertex bytes identical. Research STL exports: ' + str(len(exports)))
    return result
if __name__ == '__main__':
    run_demo()
