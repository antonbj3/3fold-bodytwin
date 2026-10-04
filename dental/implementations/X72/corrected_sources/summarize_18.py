"""Scores already frozen geometry-derived predictions against immutable anchors."""
import json, hashlib, datetime
from pathlib import Path
import numpy as np
from occlusion_operator import H, dump, sha, state, REGIONS

def agg(v):
    v = np.asarray(v)
    return v[:4] + v[4:8]

def boot_mean(x, seed=6102, alpha=0.05):
    x = np.asarray(x, float)
    rng = np.random.default_rng(seed)
    boot = x[rng.integers(0, len(x), (2000, len(x)))].mean(1)
    return dict(mean=float(x.mean()), interval=np.quantile(boot, [alpha / 2, 1 - alpha / 2]).tolist(), n=len(x), alpha=alpha)

def classify(r):
    l = r['annotation']['Left Class']
    rr = r['annotation']['Right Class']
    cl = lambda x: 'III' if 'III' in x else 'II' if 'II' in x else 'I'
    return cl(l) if cl(l) == cl(rr) else 'Mixed'

def main(round_name='R1'):
    directory = {'R1': 'cases', 'R2': 'cases_r2', 'R4': 'cases_3d'}[round_name]
    model_file = 'segmentation_r2_model.json' if round_name == 'R2' else 'segmentation_model.json'
    rows = [json.loads(p.read_text()) for p in sorted((H / 'raw' / directory).glob('*.json'))]
    if len(rows) != 200:
        raise RuntimeError('Requires all200 cases; found ' + str(len(rows)))
    prediction = dict(timestamp_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), prereg_sha256=sha(H / 'PREREG_R1.json'), cases_sha256={str(x['case']): sha(H / 'raw' / directory / f"{x['case']:03d}.json") for x in rows}, notes='Frozen registered-pose conditional predictions; no prospective force measurement has occurred', future_measurement='Acquire same-case T-Scan/Prescale without fitting; compare these unaltered values')
    fp = H / ('FROZEN_PREDICTIONS.json' if round_name == 'R1' else 'FROZEN_PREDICTIONS_' + round_name + '.json')
    if not fp.exists():
        dump(fp, prediction)
    ext = json.loads((H / 'raw/EXTERNAL_REFERENTS.json').read_text())
    ferr = ext['studies'][0]
    fv = np.asarray(ferr['mean_percent'])
    fdi = np.asarray(ferr['fdi']) % 10
    fmean = np.array([fv[np.isin(fdi, types)].sum() for types in [[1, 2], [3], [4, 5], [6, 7]]])
    fsd = np.array([np.asarray(ferr['sd_percent'])[np.isin(fdi, types)].sum() for types in [[1, 2], [3], [4, 5], [6, 7]]])
    ht = ext['studies'][1]
    hm = np.asarray(ht['single_side_mean_percent'])
    hs = np.asarray(ht['single_side_sd_percent'])
    hmean = 2 * np.array([hm[:2].sum(), hm[2], hm[3:5].sum(), hm[5:].sum()])
    hsd = 2 * np.array([hs[:2].sum(), hs[2], hs[3:5].sum(), hs[5:].sum()])
    uniform = np.array([4, 2, 4, 4]) / 14 * 100
    comparisons = []
    contrasts = []
    for (j, band) in enumerate([0.05, 0.1, 0.2]):
        eligible = [x for x in rows if x['metrics'][j]['mechanics'] is not None]
        invalid = [x['case'] for x in eligible if not x['metrics'][j]['mechanics']['valid']]
        for (jaw, off, ref, refsd, locator) in [('upper', 0, fmean, fsd, 'PMCID:PMC5726858 Table3; DOI:10.11138/ads/2017.8.2.089'), ('lower', 8, hmean, hsd, 'DOI:10.7144/sgf.2.111 Figure1 open bars, lower teeth')]:
            m = np.array([agg(np.asarray(x['metrics'][j]['mechanics']['shares_pp'])[off:off + 8]) for x in eligible])
            a = np.array([agg(np.asarray(x['metrics'][j]['area_shares_pp'])[off:off + 8]) for x in eligible])
            mm = m.mean(0)
            aa = a.mean(0)
            mae = float(np.abs(mm - ref).mean())
            aerr = float(np.abs(aa - ref).mean())
            uerr = float(np.abs(uniform - ref).mean())
            comparisons.append(dict(band_mm=band, jaw=jaw, n=len(eligible), predicted_mechanics_pp=mm.tolist(), predicted_area_pp=aa.tolist(), reference_pp=ref.tolist(), reference_sd_conservative_pp=refsd.tolist(), source_locator=locator, reference_sum_pp=float(ref.sum()), mechanics_MAE_pp=mae, area_MAE_pp=aerr, uniform_MAE_pp=uerr, gain_over_uniform_pp=uerr - mae, gates=dict(external_MAE=mae <= 5, area_superiority=bool(mae <= 0.8 * aerr and aerr - mae >= 2), uniform_improvement=uerr - mae >= 2), physical_validation='UNKNOWN: unmatched populations, source observation laws unknown, B2B M3 identity unverified', invalid_solver_cases=invalid))
        vals = {b: [] for b in ['Open Bite', 'Normal', 'Deep Bite', 'Inverted Bite']}
        for x in eligible:
            area = np.asarray(x['metrics'][j]['area_shares_pp'])
            vals[x['annotation']['Anterior Bite']].append(float((area[[0, 1, 4, 5]].sum() + area[[8, 9, 12, 13]].sum()) / 2))
        rng = np.random.default_rng(6102)

        def difference(x, y, alpha=0.05):
            if min(len(x), len(y)) < 10:
                return dict(status='UNKNOWN_SMALL_GROUP', n=[len(x), len(y)])
            x = np.asarray(x)
            y = np.asarray(y)
            v = x[rng.integers(0, len(x), (2000, len(x)))].mean(1) - y[rng.integers(0, len(y), (2000, len(y)))].mean(1)
            return dict(status='DESCRIPTIVE', mean=float(x.mean() - y.mean()), interval=np.quantile(v, [alpha / 2, 1 - alpha / 2]).tolist(), n=[len(x), len(y)])
        front = difference(vals['Open Bite'], vals['Normal'])
        front['band_mm'] = band
        front['quantity'] = 'Area anterior share (I+C), Open minus Normal; mean over two jaw-specific assignments'
        front['gate'] = front.get('mean', 0) <= -5 and front.get('interval', [0, 0])[1] < 0
        contrasts.append(front)
        for method in ['area', 'mechanics']:
            v = {c: [] for c in ['I', 'II', 'III']}
            for x in eligible:
                c = classify(x)
                if c not in v:
                    continue
                met = x['metrics'][j]
                w = np.asarray(met['area_shares_pp'] if method == 'area' else met['mechanics']['shares_pp'])
                v[c].append(float((w[[0, 1, 4, 5]].sum() + w[[8, 9, 12, 13]].sum()) / 2))
            for (c1, c2) in [('II', 'I'), ('III', 'I'), ('III', 'II')]:
                d = difference(v[c1], v[c2], 0.05 / 18)
                d.update(band_mm=band, method=method, contrast=c1 + '-' + c2, quantity='Anterior share I+C; pure bilateral Angle group')
                d['gate'] = bool(abs(d.get('mean', 0)) >= 5 and (d.get('interval', [0, 0])[0] > 0 or d.get('interval', [0, 0])[1] < 0))
                contrasts.append(d)
    val = json.loads((H / 'raw' / model_file).read_text())
    seg = all((v['validation']['gate'] for v in val.values()))
    core = [c for c in comparisons if c['band_mm'] == 0.1]
    num = all((x['metrics'][j]['mechanics'] is None or x['metrics'][j]['mechanics']['valid'] for x in rows for j in range(3)))
    result = dict(round=round_name, n_cases=200, segmentation_gate=seg, numerical_gate=num, comparisons=comparisons, contrasts=contrasts, anterior_bite_geometry_gate=all((x['gate'] for x in contrasts if 'Open minus Normal' in x['quantity'])), mechanical_superiority_gate=all((c['gates']['area_superiority'] for c in core)), patient_force_validity='UNKNOWN', external_referent=dict(kind='independent_measurement', locator='DOI:10.11138/ads/2017.8.2.089 Table3 upper; DOI:10.7144/sgf.2.111 Fig1 lower', compared_quantity='Population relative regional sensor/force shares, upper and lower scored separately', refutes_us=True), external_segmentation_referent=dict(kind='published_dataset', locator='https://osf.io/xctdy/; exact local heldout mesh+FDI hashes in raw/segmentation_validation.json', compared_quantity='Four-region heldout FDI recall', refutes_us=True), costs=dict(wall_s=sum((x['wall_s'] for x in rows)), cpu_s=sum((x['cpu_s'] for x in rows)), peak_rss_MiB=max((x['rss_MiB'] for x in rows)), shared_preparation_s=sum((x['decode_s'] + x['segmentation_s'] + x['map_s'] for x in rows)), mechanics_queries_s=sum((m.get('mechanics_query_s', 0) for x in rows for m in x['metrics'])), area_queries_s=sum((m.get('area_query_s', 0) for x in rows for m in x['metrics'])), fit_validation=json.loads((H / 'raw/segmentation_validation.json').read_text())['wall_s'], discovery='Reading, code development and source retrieval cost not measured as CPU solver cost; wall/task timestamps in CURRENT_WORK_STATE', future_measurement='NOT_RUN'))
    result['outcome'] = 'CONDITIONAL_CONTACT_OPERATOR; ' + ('REGION_TRANSFER_PASSED' if seg else 'REGION_TRANSFER_FAILED') + '; ' + ('MECHANICS_SUPERIORITY' if result['mechanical_superiority_gate'] else 'MECHANICS_NOT_SUPERIOR') + '; PATIENT_FORCE_UNKNOWN'
    outdir = H / {'R1': 'round1', 'R2': 'round2', 'R4': 'round4'}[round_name]
    dump(outdir / 'results.json', result)
    dump(H / 'results.json', result)
    state('R1_FINISHED', result['outcome'], 'Change from point-force prediction to force-feasibility and minimum measurement operator')
    (outdir / 'HANDOFF.md').write_text(round_name + ' finished. ' + result['outcome'] + '\nFrozen gates preserved in PREREG_R1.json. All200 paired meshes read without extraction. Hattori evaluated only against lower. No acquisition force/axial posterior support or sensor law identified. Continue with a preregistered graph-flow force feasible set and measurement rank, preserving R1 failures.\n')
    print(result['outcome'])
    print(json.dumps(core, indent=2))
    print('front geometry', result['anterior_bite_geometry_gate'])
if __name__ == '__main__':
    import sys
    main(sys.argv[1] if len(sys.argv) > 1 else 'R1')
