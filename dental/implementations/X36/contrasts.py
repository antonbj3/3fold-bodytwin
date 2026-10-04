import json, time, datetime
import numpy as np
from pathlib import Path
import xml.etree.ElementTree as ET
from scipy import stats
from analyze import P, read, dump, write_csv, sha, pairs, text, XML
from meta_engine import fit

def build_contrasts(gap, crown):
    out = []
    for (kind, rows) in [('cement', [r for r in gap if r['eligible'] and r['region'] == 'marginal']), ('crown', [r for r in crown if r['thickness_eligible']])]:
        groups = {}
        for r in rows:
            key = (r['study'],) if kind == 'cement' else tuple((str(r.get(k, 'UNKNOWN')) for k in ['study', 'material_class', 'cement', 'abrasion', 'margin']))
            groups.setdefault(key, []).append(r)
        for (key, a) in sorted(groups.items()):
            x = 'CAD_um' if kind == 'cement' else 'log_thickness'
            if len(set((r[x] for r in a))) < 2:
                continue
            lo = min(a, key=lambda r: r[x])
            hi = max(a, key=lambda r: r[x])
            dx = hi[x] - lo[x]
            factor = 10.0 if kind == 'cement' else 1.0
            b = (hi['outcome'] - lo['outcome']) / dx * factor
            var = (hi['variance'] + lo['variance']) / dx ** 2 * factor ** 2
            varnn = (hi['variance_no_n'] + lo['variance_no_n']) / dx ** 2 * factor ** 2
            proxy = 10 if kind == 'cement' else 2
            ci = [b - 1.95996398454 * np.sqrt(var), b + 1.95996398454 * np.sqrt(var)]
            out.append(dict(dataset=kind, cluster=lo['cluster'], stratum='|'.join(key), study=lo['study'], doi=lo['doi'], low_row_id=lo['row_id'], high_row_id=hi['row_id'], low_outcome=lo['outcome'], high_outcome=hi['outcome'], low_x=lo[x], high_x=hi[x], estimate=float(b), variance=float(var), variance_no_n=float(varnn), ci95_conditional=list(map(float, ci)), practice_proxy=proxy, local_proxy_excluded=bool(ci[1] < proxy or ci[0] > proxy), unit='um per 10um' if kind == 'cement' else 'dimensionless log ratio exponent', resolution_level=lo['resolution_level'], time_scale='SIMULTANEOUS', angle30=lo.get('angle30', 0), source_locator=lo['locator'] + '; ' + hi['locator'], precision_status='Conditional on independent arms and specimen-based SD'))
    return out

def main():
    start = time.perf_counter()
    d = read('NORMALIZED_DATA.json')
    rows = build_contrasts(d['cement'], d['crown'])
    dump('WITHIN_STUDY_CONTRASTS.json', rows)
    write_csv('WITHIN_STUDY_CONTRASTS.csv', rows)
    models = {}
    for kind in ['cement', 'crown']:
        a = [r for r in rows if r['dataset'] == kind]
        X = np.ones((len(a), 1))
        cl = [r['cluster'] for r in a]
        models[kind] = fit([r['estimate'] for r in a], X, np.diag([r['variance'] for r in a]), cl, ['intercept'], dense_control=True, cell_residual=kind != 'cement')
        models[kind]['transportable_support_pass'] = len(set(cl)) >= 3
        models[kind]['interpretation'] = 'Illustrative pooled contrasts; per-study conditional estimates are primary. No universal causal exponent.'
    g = [r for r in rows if r['dataset'] == 'cement']
    models['cement_no_n_division'] = fit([r['estimate'] for r in g], np.ones((len(g), 1)), np.diag([r['variance_no_n'] for r in g]), [r['cluster'] for r in g], ['intercept'], cell_residual=False)
    models['crown_angle_moderator'] = dict(status='UNKNOWN_INSUFFICIENT_REPLICATED_PROTOCOLS', independent_studies_per_angle={'0': 1, '30': 1}, required_studies_per_angle=3)
    src = XML / 'PMC10478297.xml'
    root = ET.parse(src).getroot()
    tb = next((t for t in root.findall('.//table-wrap') if t.get('id') == 'Tab3'))
    tr = next((tr for tr in tb.findall('.//tbody/tr') if text(tr[0]) == 'RelyX U200'))
    (mean, sd) = pairs(text(tr[1]))[0]
    old = next((r for r in d['crown'] if r['study'] == 'PMC10478297' and (not r['original_pair_verified'])))
    repair = dict(row_id=old['row_id'], old_mean_N=old['mean_N'], corrected_mean_N=mean, SD_N=sd, source_path=str(src), source_sha256=sha(src), cell_xpath=".//table-wrap[@id='Tab3']/table/tbody/tr[2]/td[2]", cell_text=text(tr[1]), cause='Flattened row label RelyX U200 concatenated with next numeric mean', applies_to='R2 corrected extraction only; R1 excludes corrupted row; frozen raw inputs unchanged', resolution_level='PER_TOOTH', time_scale='SIMULTANEOUS')
    corrected = [{**r} for r in d['crown']]
    cr = next((r for r in corrected if r['row_id'] == old['row_id']))
    cr.update(mean_N=mean, sd_N=sd, outcome=float(np.log(mean)), variance=sd ** 2 / (cr['n'] * mean ** 2), variance_no_n=sd ** 2 / mean ** 2, eligible=True, endpoint_exclusion='', original_pair_verified=True)
    write_csv('EXTRACTION_CROWN_CORRECTED_R2.csv', corrected)
    dump('SOURCE_CORRECTION_R2.json', repair)
    changedg = [{**r, 'outcome': r['outcome'] + 100} for r in d['cement']]
    changedc = [{**r, 'outcome': r['outcome'] + np.log(2)} for r in d['crown']]
    changed = build_contrasts(changedg, changedc)
    err = max((abs(a['estimate'] - b['estimate']) for (a, b) in zip(rows, changed)))
    verr = max((abs(a['variance'] - b['variance']) for (a, b) in zip(rows, changed)))
    checks = {'offset_scale_cancellation': dict(max_abs_error=err, variance_max_abs_error=verr, pass_gate=bool(max(err, verr) < 1e-09)), 'wrong_prefixed_cell_rejected': dict(injected_mean_N=old['mean_N'], source_mean_N=mean, abs_error=abs(old['mean_N'] - mean), pass_gate=abs(old['mean_N'] - mean) > 0.011), 'correct_source_cell_accepted': dict(pass_gate=mean == 2119.94 and sd == 205.93)}
    reconstruction = []
    for r in rows:
        if r['dataset'] == 'cement':
            direct = 10 * (r['high_outcome'] - r['low_outcome']) / (r['high_x'] - r['low_x'])
        else:
            direct = np.log(np.exp(r['high_outcome']) / np.exp(r['low_outcome'])) / np.log(np.exp(r['high_x']) / np.exp(r['low_x']))
        injected = direct + 1
        reconstruction.append(dict(stratum=r['stratum'], abs_error=float(abs(direct - r['estimate'])), correct_pass=bool(abs(direct - r['estimate']) < 1e-09), injected_wrong_estimate_rejected=bool(abs(injected - r['estimate']) > 1e-09)))
    checks['endpoint_control'] = reconstruction
    dump('R2_CONTROLS.json', checks)
    result = dict(round='R2', claim_type='information_link', review_state='PENDING_INDEPENDENT_REVIEW', verdict='LOCAL_CONTRASTS_IDENTIFIED_TRANSFER_REMAINS_UNSUPPORTED', models=models, contrasts=rows, source_correction=repair, controls=checks, external_referent=dict(kind='independent_measurement', locator='WITHIN_STUDY_CONTRASTS.csv DOI/table/cell; SOURCE_CORRECTION_R2.json XML cell', compared_quantity='marginal gap contrast per CAD dose and log crown-force ratio per thickness ratio', refutes_us=True), elapsed_seconds=time.perf_counter() - start, next_construction='Three-dose randomized spacer sweep measuring per-specimen marginal/axial/occlusal means/covariance; thickness x angle factorial on matched material/die/indenter with fractography')
    dump('RESULTS_R2.json', result)
    dump('CURRENT_WORK_STATE.json', dict(lane='X36-meta-regression', phase='R2_COMPLETE', latest_gate=result['verdict'], next_operation='Validate numerical implementations; deliver journal methods and transparent PRISMA reconstruction', updated_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
    (P / 'HANDOFF_R2.md').write_text('R2 cancels study offsets before synthesis. Two spacer-dose studies, both marginal; two crown thickness studies with six endpoint strata. Transfer support remains below frozen minimum. A source-cell-boundary extraction repairs the U200/mean concatenation in a separate version. See RESULTS_R2.json and raw inputs; no R1 numerical cell edited. Next is specimen-level covariance and a matched thickness-by-angle experiment.\n')
    print('R2 contrasts', [(r['study'], r['stratum'], round(r['estimate'], 4)) for r in rows])
if __name__ == '__main__':
    main()
