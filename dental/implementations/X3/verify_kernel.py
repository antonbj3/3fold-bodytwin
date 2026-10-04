from crownbench import *
from exact_distance import exact_distances
files = [DATA / 'train_example.npz']
fr = json.load(open(P / 'FROZEN_PREDICTIONS_R1.json'))
files += [pathlib.Path(r['path']) for r in fr['prediction_files'] if r['method'] in ['population_mean_sdf', 'neighbor_sdf']][:2]
rows = []
for p in files:
    (v, f) = read_surface(p)
    pts = sampled(v, f, 1024, seed(str(p)))
    pts += np.random.default_rng(6103).uniform(-2, 2, pts.shape)
    ts = time.perf_counter()
    control = distances(pts, v, f)
    slow = time.perf_counter() - ts
    ts = time.perf_counter()
    fast = exact_distances(pts, v, f)
    elapsed = time.perf_counter() - ts
    err = float(np.max(abs(control - fast)))
    assert err <= 1e-10
    rows.append({'path': str(p), 'sha256': sha(p), 'query_points': len(pts), 'offsets_mm': [-2, 2], 'max_absolute_difference_mm': err, 'Rtree_seconds': slow, 'KD_enclosure_seconds': elapsed, 'speedup': slow / elapsed})
put(P / 'NUMERICAL_KERNEL_VERIFICATION.json', {'prereg_sha256': sha(P / 'PREREG_NUMERICAL_KERNEL.json'), 'rows': rows, 'all_exactness_pass': True, 'kind': 'published_code', 'locator': 'https://github.com/mikedh/trimesh/blob/main/trimesh/proximity.py', 'scientific_metric_changed': False})
print(json.dumps(rows, indent=2))
