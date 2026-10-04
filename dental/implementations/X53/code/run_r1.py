import time, sys, collections, resource
from common import *
from collision import *

def sample(scene, inner_ids, seed=53, n=64):
    rng = np.random.default_rng(seed)
    sets = []
    mask = np.zeros(len(scene.f), bool)
    mask[inner_ids] = True
    for (region, selection) in [('intaglio', np.flatnonzero(mask)), ('exterior_and_rim', np.flatnonzero(~mask))]:
        if not len(selection):
            continue
        ids = rng.choice(selection, n, p=scene.area[selection] / scene.area[selection].sum(), replace=True)
        bary = rng.dirichlet([1, 1, 1], n)
        P = np.einsum('ij,ijk->ik', bary, scene.tri[ids])
        norm = np.cross(scene.tri[ids, 1] - scene.tri[ids, 0], scene.tri[ids, 2] - scene.tri[ids, 0])
        norm /= np.linalg.norm(norm, axis=1, keepdims=True)
        for i in range(n):
            sets.append(dict(point=P[i], normal=norm[i], face_id=int(ids[i]), region=region, area_weight_mm2=scene.area[selection].sum() / n))
    return sets

def lib(library, axes, scale=0.8):
    out = []
    for x in read(ROOT / 'TOOL_LIBRARY.json')['tools']:
        if x['library'] != library:
            continue
        t = dict(x)
        if library == 'vhf':
            t['total_length_mm'] = 35.0 if axes == 4 else 40.0
        for k in ['diameter_mm', 'neck_reach_mm', 'shank_mm', 'gauge_mm', 'holder_diameter_mm', 'holder_length_mm']:
            t[k] *= scale
        out.append(t)
    return sorted(out, key=lambda t: t['diameter_mm'])

def run():
    p = ROOT / 'PREREG_R1.json'
    pre = read(p)
    if sha(p) != p.with_suffix('.sha256').read_text().strip():
        raise ValueError('PREREG drift')
    if sha(ROOT / 'INPUT_MANIFEST.json') != pre['frozen_inputs']['manifest_sha256'] or sha(ROOT / 'TOOL_LIBRARY.json') != pre['frozen_inputs']['tool_library_sha256']:
        raise ValueError('input card drift')
    start = time.perf_counter()
    summary = []
    all_hash = {}
    manifest = read(ROOT / 'INPUT_MANIFEST.json')
    for (index, row) in enumerate(manifest['rows']):
        output = DATA / 'R1' / f"{row['id']}.json"
        if output.exists():
            r = read(output)
            if r['source_sha256'] != row['sha256'] or r['prereg_sha256'] != sha(p):
                raise ValueError('resume input drift')
            summary.append(r['summary'])
            all_hash[str(output)] = sha(output)
            continue
        if sha(row['file']) != row['sha256']:
            raise ValueError('mesh drift')
        a = np.load(row['file'], allow_pickle=False)
        scene = Scene(a['vertices'], a['faces'])
        pts = sample(scene, a['inner_faces'], 53 + index, pre['metrics']['sample_count_per_region'])
        variants = []
        for library in ['vhf', 'ceramill']:
            for axes in [4, 5]:
                tools = lib(library, axes)
                records = []
                seconds = time.perf_counter()
                for point in pts:
                    (P, N) = (np.asarray(point['point']), np.asarray(point['normal']))
                    access = []
                    for t in tools:
                        witness = search(scene, P, N, t, axes, offsets=(0.0,))
                        ball = scene.ball_clearance(P + (t['diameter_mm'] / 2 + 5 * EPS) * N, t['diameter_mm'] / 2) >= -EPS
                        access.append(dict(tool_id=t['id'], nominal_diameter_green_mm=t['diameter_mm'] / 0.8, diameter_design_mm=t['diameter_mm'], ball_only=bool(ball), **witness))
                    available = [x for x in access if x['status'] == 'FOUND']
                    records.append(dict(**point, tools=access, smallest_found_green_mm=min((x['nominal_diameter_green_mm'] for x in available)) if available else None, largest_found_green_mm=max((x['nominal_diameter_green_mm'] for x in available)) if available else None, status='CONDITIONAL_SCENARIO_WITNESS' if available else 'NOT_FOUND_ON_POSE_GRID', global_minimum_status='UNKNOWN_UNSAMPLED_CONTINUOUS_POSES'))
                by_region = {}
                for reg in ['intaglio', 'exterior_and_rim']:
                    rr = [r for r in records if r['region'] == reg]
                    if not rr:
                        continue
                    by_region[reg] = dict(n=len(rr), area_mm2=sum((r['area_weight_mm2'] for r in rr)), found_fraction=sum((r['smallest_found_green_mm'] is not None for r in rr)) / len(rr), ball_only_fraction=sum((any((t['ball_only'] for t in r['tools'])) for r in rr)) / len(rr), not_found_fraction=sum((r['smallest_found_green_mm'] is None for r in rr)) / len(rr), hoeffding_95_halfwidth=float(np.sqrt(np.log(40) / (2 * len(rr)))), unconditional_machine_status='UNKNOWN_MEASURE_HOLDER_NECK_MOUNT_FIXTURE')
                variants.append(dict(library=library, axes=axes, model='conditional simplified complete-tool geometry', records=records, by_region=by_region, seconds=time.perf_counter() - seconds))
        sums = dict(id=row['id'], family=row['family'], variants=[{k: v for (k, v) in vv.items() if k != 'records'} for vv in variants], empirical_CAM_status='UNKNOWN_NOT_MEASURED')
        dump(output, dict(id=row['id'], source_sha256=row['sha256'], prereg_sha256=sha(p), resolution='PER_POINT', points_and_tools=variants, summary=sums))
        summary.append(sums)
        all_hash[str(output)] = sha(output)
        dump(ROOT / 'raw/R1_PARTIAL.json', summary)
        state('R1_FIELDS_RUNNING', f"{index + 1}/{len(manifest['rows'])} designs saved", 'Continue full tool fields; do not infer no continuous solution from pose-grid absence')
        print(row['id'], [(v['library'], v['axes'], round(v['by_region']['intaglio']['found_fraction'], 3)) for v in variants], flush=True)
    result = dict(round='R1', claim_type='capability', prereg_sha256=sha(p), rows=summary, requested=manifest['requested'], accepted=len(summary), excluded=manifest['excluded'], wall_seconds_this_invocation=time.perf_counter() - start, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, own_bytes=budget(), raw_files=all_hash, physical_gate='UNKNOWN', scope='64 uniform-area samples per region, continuous straight sweeps at finite pose set, simplified measured/unknown geometry split')
    dump(ROOT / 'RESULTS_R1.json', result)
    freeze(ROOT / 'FROZEN_PREDICTIONS.json', dict(claim_type='capability', prereg_sha256=sha(p), result_sha256=sha(ROOT / 'RESULTS_R1.json'), point_files=all_hash, input_manifest_sha256=sha(ROOT / 'INPUT_MANIFEST.json'), tool_library_sha256=sha(ROOT / 'TOOL_LIBRARY.json'), code={str(f): sha(f) for f in (ROOT / 'code').glob('*.py')}, physical_measurement_performed=False))
    state('R1_FIELDS_FROZEN', 'Digital sampled-pose fields delivered; physical machine geometry UNKNOWN', 'Sufficiency counterexample, local stock-to-film successor and prospective 20-case CAM comparison')
    return result
if __name__ == '__main__':
    run()
