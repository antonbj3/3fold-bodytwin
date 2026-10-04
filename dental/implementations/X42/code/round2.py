from common import *
from model_fit import dataset, fit_arrays, loadmodel
from generator import predict, conditional
import itertools, resource, collections
PARAMS = {'kernel': list(itertools.product([0.05, 0.2, 1.0], [0.1, 1.0, 10.0])), 'forest': [2, 5, 10], 'ridge': [1.0, 10.0, 100.0, 1000.0], 'mean': [0]}

def tag(kind, p):
    return kind + '_' + str(p).replace(' ', '').replace('(', '').replace(')', '').replace(',', '_').replace('.', 'p')

def run():
    start = time.perf_counter()
    sp = split()
    (ds, rejected, access) = dataset(sp['dev_fit'])
    base = DATA / 'R2_models'
    base.mkdir(exist_ok=True)
    fit_cost = []
    modelmap = {}
    for (fam, recs) in ds.items():
        for (kind, pp) in PARAMS.items():
            for p in pp:
                tic = time.perf_counter()
                m = fit_arrays(recs, kind, p)
                fn = base / (fam + '_' + tag(kind, p) + '.npz')
                np.savez_compressed(fn, **m)
                modelmap[fam, tag(kind, p)] = m
                fit_cost.append(dict(family=fam, kind=kind, param=p, seconds=time.perf_counter() - tic, donors=len(recs), bytes=fn.stat().st_size, sha256=sha(fn)))
        print('R2 fit', fam, len(recs), flush=True)
    dump(ROOT / 'raw/R2_FIT.json', dict(cost=fit_cost, rejected=rejected, dev_reference_keys=access))
    losses = collections.defaultdict(list)
    detail = []
    for key in sp['dev_validation']:
        refs = dev_reference(key)
        for t in tasks(key)[::4]:
            if t['status'] != 'READY':
                continue
            t = scene(t)
            ref = refs[t['family']]
            area = float(t['weights'].sum())
            for (kind, pp) in PARAMS.items():
                for p in pp:
                    name = tag(kind, p)
                    met = metric(t, predict(t, modelmap[t['family'], name]), ref)
                    if met['rmse_mm'] is not None:
                        loss = met['rmse_mm'] + (met['contact_error_mm2'] / area if met['contact_error_mm2'] is not None else 0.0)
                        losses[name].append(loss)
                        detail.append(dict(case_key=key, family=t['family'], model=name, loss=loss, **met))
    mean = {k: float(np.mean(v)) for (k, v) in losses.items()}
    selected = {kind: min([tag(kind, p) for p in pp], key=lambda k: mean[k]) for (kind, pp) in PARAMS.items()}
    control = min([selected['forest'], selected['ridge']], key=lambda k: mean[k])
    dump(ROOT / 'raw/R2_SELECTION.json', dict(selected=selected, strongest_conventional=control, unprojected_validation_mean_loss=mean, validation_cases=sp['dev_validation'], selection_metric='RMSE mm/1mm + absolute contact error / footprint mm2', no_test_query=True))
    dump(ROOT / 'raw/R2_VALIDATION_LOSSES.json', detail)
    rows = []
    cost = collections.Counter()
    methods = {'x42_shape': selected['kernel'], 'control_learned': control, 'control_mean': selected['mean']}
    for (ii, key) in enumerate(sp['dev_validation']):
        refs = dev_reference(key)
        ss = {}
        for orig in tasks(key):
            t = orig
            if t['status'] == 'READY':
                if t['geometry_file'] not in ss:
                    ss[t['geometry_file']] = scene(t)
                t = dict(ss[t['geometry_file']])
                t.update(orig)
                t['preparation_z'] = np.full(len(t['xy']), t['preparation_height_mm'])
            for method in [*methods, 'constraint_optimizer']:
                row = dict(case_key=key, family=t['family'], level=t['level'], method=method)
                if t['status'] != 'READY':
                    row['verdict'] = 'UNKNOWN_SITE'
                else:
                    tic = time.perf_counter()
                    d = standard_lp(t) if method == 'constraint_optimizer' else conditional(t, modelmap[t['family'], methods[method]])
                    cost[method] += time.perf_counter() - tic
                    ch = all_checks(t, d)
                    row.update(verdict=ch['validity'])
                    if d['status'] == 'DESIGN':
                        row.update(metric(t, np.asarray(d['outer_vertices'])[:, 2], refs[t['family']]))
                rows.append(row)
        if (ii + 1) % 10 == 0:
            dump(ROOT / 'raw/R2_ROWS.json', rows)
            state('R2_VALIDATING', {'cases': ii + 1}, 'Next frozen dev case', hidden_queries=0)
            print('R2 validate', ii + 1, 'seconds', round(time.perf_counter() - start), flush=True)
    dump(ROOT / 'raw/R2_ROWS.json', rows)
    out = dict(round='R2', claim_type='algorithm', summary=summarize(rows), selected=selected, strongest_conventional=control, unprojected_loss=mean, cost_seconds=dict(cost), fit_seconds=sum((r['seconds'] for r in fit_cost)), wall_seconds=time.perf_counter() - start, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, donor_rejected=len(rejected), attempted_donors=len(sp['dev_fit']) * 9, hidden_queries=0, prereg_sha256=sha(ROOT / 'PREREG_R2.json'), external_referent=json.loads((ROOT / 'PREREG_R2.json').read_text())['external_referent'])
    dump(ROOT / 'rounds/R2.json', out)
    print(json.dumps(out, indent=2))
if __name__ == '__main__':
    run()
