"""R2: fixed-region morphology sensitivity and conditional radius enclosures."""
import csv, json, resource, time
from pathlib import Path
import numpy as np
from scipy import ndimage as ndi
from operators import persistent_count, circumcircle, H
from atlas import ROOT, DATA, sha, dump, csvwrite

def counts_in_region(p, lo, hi):
    counts = []
    for z in range(lo, hi + 1):
        (labels, n) = ndi.label(p[z], structure=np.ones((3, 3)))
        counts.append(int(np.sum(np.bincount(labels.ravel())[1:] >= 2)))
    return (persistent_count(counts), counts)

def radius_interval(triplet, eps):
    (a, b, c) = np.asarray(triplet, float)
    (u, v) = (b - a, c - a)
    sides = np.array([np.linalg.norm(u), np.linalg.norm(c - b), np.linalg.norm(v)])
    cross = np.linalg.norm(np.cross(u, v))
    err = 2 * eps * (np.linalg.norm(u) + np.linalg.norm(v)) + 4 * eps ** 2
    lower = float(np.prod(np.maximum(0, sides - 2 * eps)) / (2 * (cross + err))) if cross + err > 0 else 0.0
    upper = float(np.prod(sides + 2 * eps) / (2 * (cross - err))) if cross > err else None
    return (lower, upper)

def topology_class_consistent(row):
    expected = row['nominal_count'] == row['eroded_count'] == row['dilated_count']
    return row['scenario_stable'] == expected

def main():
    start = time.monotonic()
    cpu = time.process_time()
    specs = json.loads((ROOT / 'PREREG_R2_TOPOLOGY_RADIUS_UNCERTAINTY.json').read_text())
    manifests = json.loads((ROOT / 'raw/R1_manifest.json').read_text())
    graphs = {(r['case'], r['fdi']): r for r in map(json.loads, (DATA / 'R1_graphs.jsonl').open())}
    topo = []
    for row in manifests:
        path = Path(row['path'])
        if sha(path) != row['sha256']:
            raise RuntimeError('Local crop hash drift')
        with np.load(path) as a:
            p = a['pulp']
        (lo, hi) = graphs[row['case'], row['fdi']]['z_range']
        eroded = np.zeros_like(p)
        dilated = np.zeros_like(p)
        for z in range(lo, hi + 1):
            eroded[z] = ndi.distance_transform_edt(p[z], sampling=H) > H + 1e-09
            dilated[z] = ndi.distance_transform_edt(~p[z], sampling=H) <= H + 1e-09
        (nominal, _) = counts_in_region(p, lo, hi)
        (er, _) = counts_in_region(eroded, lo, hi)
        (di, _) = counts_in_region(dilated, lo, hi)
        topo.append({'case': row['case'], 'fdi': row['fdi'], 'resolution_level': 'PER_TOOTH', 'nominal_count': nominal, 'eroded_count': er, 'dilated_count': di, 'scenario_min': min(nominal, er, di), 'scenario_max': max(nominal, er, di), 'scenario_stable': nominal == er == di, 'scenario_canal_vanished': er == 0, 'meaning': 'digital scenario envelope, not an anatomical confidence interval'})
    csvwrite(ROOT / 'raw/R2_topology.csv', topo)
    radius = []
    perturbation_triplets = []
    for row in map(json.loads, (DATA / 'R1_paths.jsonl').open()):
        if row['curvature_triplet_mm'] is None:
            continue
        perturbation_triplets.append(row['curvature_triplet_mm'])
        for eps in specs['frozen_parameters']['radius_point_epsilon_mm']:
            (low, up) = radius_interval(row['curvature_triplet_mm'], eps)
            radius.append({'case': row['case'], 'fdi': row['fdi'], 'branch_id': row['branch_id'], 'resolution_level': 'PER_SURFACE_REGION', 'point_epsilon_mm': eps, 'nominal_window_radius_mm': row['window_radius_mm'], 'radius_lower_mm': low, 'radius_upper_mm': up, 'status': 'CONDITIONAL_FINITE' if up is not None else 'UNIDENTIFIED_UPPER_RADIUS', 'closure': 'Smoothed centroid-position error <= epsilon; not validated by boundary error alone'})
    csvwrite(ROOT / 'raw/R2_radius.csv', radius)
    rng = np.random.default_rng(30)
    violations = 0
    poison_rejected = 0
    for i in range(10000):
        tr = np.asarray(perturbation_triplets[i % len(perturbation_triplets)])
        eps = specs['frozen_parameters']['radius_point_epsilon_mm'][i % 3]
        (lo, up) = radius_interval(tr, eps)
        delta = rng.normal(size=(3, 3))
        delta /= np.linalg.norm(delta, axis=1)[:, None]
        delta *= rng.random((3, 1)) * eps
        (r, _) = circumcircle(*tr + delta)
        violations += int(r < lo - 1e-08 or (up is not None and r > up + 1e-08))
        poison_rejected += int(r < max(1, lo) * 100000)
    unstable = [r for r in topo if not r['scenario_stable']]
    topology_poison = False
    if unstable:
        poisoned = dict(unstable[0], scenario_stable=True)
        dump(ROOT / 'raw/POISON_TOPOLOGY_ROW.json', poisoned)
        topology_poison = not topology_class_consistent(json.loads((ROOT / 'raw/POISON_TOPOLOGY_ROW.json').read_text()))
    stable = float(np.mean([r['scenario_stable'] for r in topo]))
    eps0 = specs['frozen_parameters']['radius_point_epsilon_mm'][0]
    subset = [r for r in radius if r['point_epsilon_mm'] == eps0]
    identified = float(np.mean([r['radius_upper_mm'] is not None for r in subset])) if subset else 0
    result = {'claim_type': 'capability', 'topology_stable_fraction': stable, 'topology_gate_pass': stable >= 0.8, 'finite_radius_upper_fraction': identified, 'radius_gate_pass': identified >= 0.8, 'topology_counts': {'n': len(topo), 'stable': sum((r['scenario_stable'] for r in topo)), 'vanished': sum((r['scenario_canal_vanished'] for r in topo))}, 'radius_n': len(subset), 'perturbation_control': {'n': 10000, 'violations': violations, 'gate_pass': violations == 0, 'poison_rejected_count': poison_rejected, 'all_poison_rejected': poison_rejected == 10000}, 'topology_poison_rejected': topology_poison, 'all_nominal_topology_classes_consistent': all((topology_class_consistent(r) for r in topo)), 'wall_s': time.monotonic() - start, 'cpu_s': time.process_time() - cpu, 'peak_RSS_MiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 'anatomical_topology': 'UNKNOWN', 'radius_interval_validity': 'CONDITIONAL_ON_POINT_ERROR_NOT_SEGMENTATION_ERROR'}
    dump(ROOT / 'raw/R2_result.json', result)
    print(json.dumps(result), flush=True)
if __name__ == '__main__':
    main()
