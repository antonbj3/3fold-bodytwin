from common import *
from model_fit import dataset, fit_arrays
from analysis_tools import decide
import resource

def run():
    r2 = json.load(open(ROOT / 'rounds/R2.json'))
    r3 = json.load(open(ROOT / 'rounds/R3.json'))
    (ds, rejected, access) = dataset(split()['dev_all'])
    base = DATA / 'FINAL_MODELS'
    base.mkdir(exist_ok=True)
    cost = []
    params = {'kernel': (1.0, 0.1), 'forest': 2, 'ridge': 100.0, 'mean': 0}
    for (fam, rr) in ds.items():
        for (kind, p) in params.items():
            tic = time.perf_counter()
            m = fit_arrays(rr, kind, p)
            folder = base / kind
            folder.mkdir(exist_ok=True)
            fn = folder / (fam + '.npz')
            np.savez_compressed(fn, **m)
            cost.append(dict(family=fam, method=kind, seconds=time.perf_counter() - tic, sha256=sha(fn), bytes=fn.stat().st_size, donors=len(rr)))
        for kind in ['kernel', 'forest']:
            cc = [dict(r, y=r['contact']) for r in rr]
            tic = time.perf_counter()
            m = fit_arrays(cc, kind, params[kind], output_modes=0)
            folder = base / (kind + '_contact')
            folder.mkdir(exist_ok=True)
            fn = folder / (fam + '.npz')
            np.savez_compressed(fn, **m)
            cost.append(dict(family=fam, method=kind + '_contact', seconds=time.perf_counter() - tic, sha256=sha(fn), bytes=fn.stat().st_size, donors=len(rr)))
    cfg = {'participants': ['x42_shape', 'x42_contact_branch', 'control_forest', 'control_contact_forest', 'control_ridge', 'control_mean', 'constraint_optimizer', 'population', 'parametric', 'x18_antagonist', 'x1b_v2_radius'], 'model_methods': {'x42_shape': 'kernel', 'control_forest': 'forest', 'control_ridge': 'ridge', 'control_mean': 'mean'}, 'branch_methods': {'x42_contact_branch': {'shape': 'kernel', 'contact': 'kernel_contact', 'beta': float(r3['selected']['kernel'].split('_')[-1])}, 'control_contact_forest': {'shape': 'forest', 'contact': 'forest_contact', 'beta': float(r3['selected']['forest'].split('_')[-1])}}, 'x1b_radius_mm': {'easy': 0.25, 'normal': 0.5, 'hard': 0.8, 'boundary': 1.0}, 'selected_candidate': 'x42_contact_branch' if float(r3['selected']['kernel'].split('_')[-1]) else 'x42_shape', 'selection_rule': 'R3 predeclared minimum projected dev loss including beta0; no test scores available', 'max_cpu_workers': 3, 'max_numeric_threads_per_worker': 1, 'fit_used_reference_keys': access, 'test_refit': 'NONE', 'mirror': 'NOT_AVAILABLE_PUBLIC_INPUT'}
    dump(ROOT / 'FINAL_CONFIG.json', cfg)
    dump(ROOT / 'raw/FINAL_FIT.json', dict(cost=cost, rejections=rejected, attempted=len(access) * 9, rejected=len(rejected), rejection_fraction=len(rejected) / (len(access) * 9), native_dev_cases=len(access), reference_access_keys=access, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024))
    print('Final refit cases', len(access), 'rejected', len(rejected), 'selected', cfg['selected_candidate'], 'beta', cfg['branch_methods'])
if __name__ == '__main__':
    run()
