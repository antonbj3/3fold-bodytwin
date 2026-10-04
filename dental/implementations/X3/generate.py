from crownbench import *
import argparse, resource

def generate(group, round_name):
    start = time.perf_counter()
    records = []
    pub = json.load(open(DATA / 'public' / f'{group}_INDEX.json'))
    state(round_name + '_GENERATING', 'PUBLIC_CONTEXTS_READY', 'finish all predictions then freeze hashes')
    for r in pub:
        rr = {k: r[k] for k in ['key', 'case', 'jaw', 'fdi', 'eligible']}
        rr['methods'] = {}
        if not r['eligible']:
            rr['reason'] = r['reason']
            records.append(rr)
            continue
        try:
            assert sha(r['context']) == r['context_sha256'], 'PUBLIC_CONTEXT_HASH_DRIFT'
            with np.load(r['context']) as d:
                v = d['vertices']
                f = d['faces']
                labels = d['labels']
            k = r['fdi']
            assert k not in np.unique(labels), 'TARGET_LEAKAGE'
            teeth = tooth_dict(v, f, labels)
            q = 1 if r['jaw'] == 'upper' else 3
            n = axis(teeth)
            (rn, rc, pl) = reflection(teeth, q)
            (tv, tf) = teeth[swap(k)]
            mirror = reflect(tv, rn, rc)
            tf = tf[:, ::-1]
            F = local_frame(teeth, k, n)
            center = mirror.mean(0)
            (pv, pf) = read_surface(DATA / 'model' / f"{r['jaw']}_{k % 10}.npz")
            pv = pv @ F.T + center
            ts = time.perf_counter()
            (nv, nf, fit) = neighbor_fit(mirror, tf, teeth, k, F)
            fit['fit_seconds'] = time.perf_counter() - ts
            ts = time.perf_counter()
            (sv, sf, phi, sd) = sdf_roundtrip(nv, nf, F, center)
            sd['seconds'] = time.perf_counter() - ts
            methods = {'mirror': (mirror, tf, {'seconds': 0.0}), 'population_mean_sdf': (pv, pf, {'seconds': 0.0, 'train_only': True}), 'neighbor_mesh': (nv, nf, fit), 'neighbor_sdf': (sv, sf, dict(sd, shared_fit=fit))}
            if round_name == 'R2':
                from ring_method import ring_fit
                ts = time.perf_counter()
                (rv, rf, ring) = ring_fit(mirror, tf, v, f, labels, k, F, teeth)
                ring['seconds'] = time.perf_counter() - ts
                methods['ring_mesh'] = (rv, rf, ring)
                ts = time.perf_counter()
                (rsv, rsf, rphi, rsd) = sdf_roundtrip(rv, rf, F, center)
                rsd['seconds'] = time.perf_counter() - ts
                methods['ring_sdf'] = (rsv, rsf, dict(rsd, shared_fit=ring))
            for (name, (mv, mf, details)) in methods.items():
                path = DATA / 'predictions' / round_name / name / f"{r['key']}.npz"
                write_surface(path, mv, mf)
                rr['methods'][name] = {'path': str(path), 'sha256': sha(path), 'vertices': len(mv), 'faces': len(mf), 'details': details}
            if not (P / 'exports' / round_name).exists():
                exp = P / 'exports' / round_name
                exp.mkdir(parents=True, exist_ok=True)
                for (name, (mv, mf, _)) in methods.items():
                    mesh(mv, mf).export(exp / f'{name}.obj')
                put(exp / 'CONTEXT.json', {'key': r['key'], 'units': 'provisionally mm', 'surface': 'exposed natural crown reconstruction; open cervical rim; no intaglio', 'dataset_attribution': 'Teeth3DS + Bone -Hamadou et al; license conflict see README_DEMO'})
            rr.update(frame=F.tolist(), center=context_center(center), plane=pl, visible_labels=sorted(teeth), target_in_context=False)
        except Exception as e:
            rr['generation_failure'] = type(e).__name__ + ': ' + str(e)
        records.append(rr)
        put(P / f'PREDICTION_PROGRESS_{round_name}.json', {'completed': len(records), 'scheduled': len(pub), 'last': rr['key'], 'failures': sum(('generation_failure' in a for a in records))})
        if len(records) % 4 == 0:
            print('GENERATE', round_name, len(records), 'of', len(pub), 'seconds', round(time.perf_counter() - start, 1), flush=True)
        if len(records) % 16 == 0:
            disk_guard()
    idx = DATA / 'predictions' / round_name / 'INDEX.json'
    put(idx, records)
    freeze = {'schema': 'X3-frozen-predictions-v1', 'round': round_name, 'frozen_utc': now(), 'predictions_index': str(idx), 'index_sha256': sha(idx), 'prereg_sha256': sha(P / f'PREREG_{round_name}.json'), 'split_sha256': sha(P / 'SPLIT.json'), 'implementation_manifest': [{'path': str(p), 'sha256': sha(p)} for p in P.glob('*.py')], 'prediction_files': [{'key': r['key'], 'method': m, 'path': v['path'], 'sha256': v['sha256']} for r in records for (m, v) in r['methods'].items()], 'all_predictions_complete': len(records) == len(pub), 'hidden_evaluation_started': False, 'fit_and_generation_seconds': time.perf_counter() - start, 'peak_rss_kb': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, 'failure_count': sum(('generation_failure' in r for r in records)), 'data_bytes': disk_guard()}
    put(P / f'FROZEN_PREDICTIONS_{round_name}.json', freeze)
    if round_name == 'R1':
        put(P / 'FROZEN_PREDICTIONS.json', freeze)
    state(round_name + '_PREDICTIONS_FROZEN', {'complete': True, 'failures': freeze['failure_count']}, 'unblind external original tooth surfaces in evaluator only')
    print('FROZEN', round_name, 'seconds', freeze['fit_and_generation_seconds'], 'failures', freeze['failure_count'], flush=True)

def context_center(c):
    return c.tolist()
if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--round', default='R1', choices=['R1', 'R2'])
    a = ap.parse_args()
    generate('test_R1' if a.round == 'R1' else 'reserve_R2', a.round)
