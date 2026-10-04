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
    mm = mesh(v, f)
    ts = time.perf_counter()
    control = np.concatenate([trimesh.proximity.closest_point_naive(mm, pts[a:a + 16])[1] for a in range(0, len(pts), 16)])
    slow = time.perf_counter() - ts
    ts = time.perf_counter()
    fast = exact_distances(pts, v, f)
    elapsed = time.perf_counter() - ts
    err = float(np.max(abs(control - fast)))
    rtree = distances(pts, v, f)
    row = {'path': str(p), 'sha256': sha(p), 'query_points': len(pts), 'offsets_mm': [-2, 2], 'max_absolute_difference_mm': err, 'brute_force_seconds': slow, 'KD_enclosure_seconds': elapsed, 'Rtree_difference_mm': float(np.max(abs(rtree - fast))), 'exactness_pass': err <= 1e-10}
    rows.append(row)
    print(row, flush=True)
put(P / 'NUMERICAL_KERNEL_VERIFICATION_K2.json', {'prereg_sha256': sha(P / 'PREREG_NUMERICAL_KERNEL_K2.json'), 'rows': rows, 'all_exactness_pass': all((r['exactness_pass'] for r in rows)), 'kind': 'published_code', 'locator': 'https://github.com/mikedh/trimesh/blob/main/trimesh/proximity.py', 'scientific_metric_changed': False, 'R1_failed_gate_retained': True})
assert all((r['exactness_pass'] for r in rows))
