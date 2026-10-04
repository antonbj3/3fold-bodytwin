from common import *
from model_fit import dataset, fit_arrays, loadmodel
from generator import branch
from round2 import tag
from analysis_tools import paired, decide
import resource, collections, copy, concurrent.futures
MODELS = {}
BETA = [0.25, 0.5, 1.0]

def casework(key):
    refs = dev_reference(key)
    rows = []
    cost = collections.Counter()
    ss = {}
    for orig in tasks(key):
        t = orig
        if t['status'] == 'READY':
            if t['geometry_file'] not in ss:
                ss[t['geometry_file']] = scene(t)
            t = dict(ss[t['geometry_file']])
            t.update(orig)
            t['preparation_z'] = np.full(len(t['xy']), t['preparation_height_mm'])
        for kind in ['kernel', 'forest']:
            for beta in BETA:
                method = kind + '_branch_' + str(beta)
                row = dict(case_key=key, family=t['family'], level=t['level'], method=method)
                if t['status'] != 'READY':
                    row['verdict'] = 'UNKNOWN_SITE'
                else:
                    tic = time.perf_counter()
                    (sm, cm) = MODELS[t['family'], kind]
                    d = branch(t, sm, cm, beta)
                    cost[method] += time.perf_counter() - tic
                    row['verdict'] = all_checks(t, d)['validity']
                    if d['status'] == 'DESIGN':
                        row.update(metric(t, np.asarray(d['outer_vertices'])[:, 2], refs[t['family']]))
                        row['footprint_mm2'] = float(t['weights'].sum())
                rows.append(row)
    return (rows, dict(cost))

def run():
    start = time.perf_counter()
    sp = split()
    selection = json.load(open(ROOT / 'raw/R2_SELECTION.json'))
    (ds, rejected, access) = dataset(sp['dev_fit'])
    base = DATA / 'R3_models'
    base.mkdir(exist_ok=True)
    fcost = []
    for (fam, recs) in ds.items():
        for kind in ['kernel', 'forest']:
            shape = loadmodel(DATA / 'R2_models' / (fam + '_' + selection['selected'][kind] + '.npz'))
            rr = [dict(r, y=r['contact']) for r in recs]
            param = (1.0, 0.1) if kind == 'kernel' else 2
            tic = time.perf_counter()
            cm = fit_arrays(rr, kind, param, output_modes=0)
            fn = base / (fam + '_' + kind + '_contact.npz')
            np.savez_compressed(fn, **cm)
            MODELS[fam, kind] = (shape, cm)
            fcost.append(dict(family=fam, kind=kind, seconds=time.perf_counter() - tic, bytes=fn.stat().st_size, sha256=sha(fn)))
    dump(ROOT / 'raw/R3_FIT.json', dict(fits=fcost, rejected=rejected, reference_access_keys=access))
    rows = []
    cost = collections.Counter()
    with concurrent.futures.ProcessPoolExecutor(max_workers=2) as pool:
        for (ii, (rr, cc)) in enumerate(pool.map(casework, sp['dev_validation'])):
            rows.extend(rr)
            cost.update(cc)
            if (ii + 1) % 10 == 0:
                dump(ROOT / 'raw/R3_ROWS.json', rows)
                state('R3_VALIDATING', {'completed_dev_cases': ii + 1}, 'Next dev case; retain raw contact branch fields', hidden_queries=0)
                print('R3 cases', ii + 1, 'seconds', round(time.perf_counter() - start), flush=True)
    original = json.load(open(ROOT / 'raw/R2_ROWS.json'))
    for r in original:
        if r['method'] not in ['x42_shape', 'control_learned']:
            continue
        rr = dict(r)
        rr['method'] = ('kernel' if r['method'] == 'x42_shape' else 'forest') + '_branch_0.0'
        if 'rmse_mm' in rr:
            ot = [t for t in tasks(r['case_key']) if t['family'] == r['family'] and t['level'] == r['level']][0]
            rr['footprint_mm2'] = float(scene(ot)['weights'].sum())
        rows.append(rr)
    losses = {}
    chosen = {}
    for kind in ['kernel', 'forest']:
        names = [kind + '_branch_' + str(b) for b in [0.0, 0.25, 0.5, 1.0]]
        for name in names:
            rr = [r for r in rows if r['method'] == name and r['verdict'] == 'PASS' and (r.get('rmse_mm') is not None)]
            losses[name] = float(np.mean([r['rmse_mm'] + (r.get('contact_error_mm2') or 0) / r['footprint_mm2'] for r in rr]))
        chosen[kind] = min(names, key=lambda k: losses[k])
    dump(ROOT / 'raw/R3_ROWS.json', rows)
    cs = paired(rows, chosen['kernel'], chosen['forest'])
    r = dict(round='R3', claim_type='algorithm', summary=summarize(rows), selected=chosen, selection_loss=losses, contrasts=cs, gate_vs_learned=decide(cs[-1]), fit_seconds=sum((r['seconds'] for r in fcost)), cost_seconds=dict(cost), beta0_cost='raw/R2.json re-used, no free discovery', wall_seconds=time.perf_counter() - start, peak_rss_main_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, hidden_queries=0, prereg_sha256=sha(ROOT / 'PREREG_R3.json'), external_referent=json.load(open(ROOT / 'PREREG_R3.json'))['external_referent'])
    dump(ROOT / 'rounds/R3.json', r)
    print(json.dumps({k: v for (k, v) in r.items() if k != 'contrasts'}, indent=2))
if __name__ == '__main__':
    run()
