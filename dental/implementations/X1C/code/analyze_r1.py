import ast, csv, json, math, time, hashlib, datetime
from pathlib import Path
import numpy as np
from scipy.special import gamma
R = Path(__file__).resolve().parents[1]

def read(p):
    return json.loads((R / p).read_text())

def dump(p, x):
    (R / p).write_text(json.dumps(x, indent=2, allow_nan=False, default=lambda o: o.item() if isinstance(o, np.generic) else str(o)) + '\n')

def sha(p):
    return hashlib.sha256((R / p).read_bytes()).hexdigest()

def state(status, gate, nxt):
    dump('CURRENT_WORK_STATE.json', dict(lane='X1c-crown-literature', updated_at=datetime.datetime.now(datetime.timezone.utc).isoformat(), status=status, latest_gate=gate, next_operation=nxt, claim_type='information_link'))

def slope(rs):
    x = np.log([r['t'] for r in rs])
    y = np.log([r['mean'] for r in rs])
    v = np.array([r['v'] for r in rs])
    w = 1 / v
    X = np.c_[np.ones(len(rs)), x]
    coef = np.linalg.solve(X.T @ (w[:, None] * X), X.T @ (w * y))
    cov = np.linalg.inv(X.T @ (w[:, None] * X))
    return dict(n=float(coef[1]), variance=float(cov[1, 1]), intercept=float(coef[0]), groups=len(rs), residuals=(y - X @ coef).tolist())

def curated():
    primary = {r['pmcid']: r for r in read('inputs/primary_tables.json')}
    rows = []

    def add(st, table, product, cl, protocol, ts, means, sds, n, eligible=True, note='', origin='measured'):
        txt = next((t['text'] for t in primary[st]['tables'] if t['label'] == table))
        for (t, mu, sd) in zip(ts, means, sds):
            import re
            fmt = lambda z: '(?:' + '|'.join((re.escape(v) for v in sorted({str(float(z)), format(z, '.2f'), format(z, 'g')}, key=len, reverse=True))) + ')'
            sep = '\\s*' if st == 'PMC10286427' else '\\s*±\\s*'
            prefix = re.escape(str(t) + ' mm') if st == 'PMC10286427' else ''
            assert re.search(prefix + fmt(mu) + sep + fmt(sd), txt), (st, product, t, mu, sd)
            if st == 'PMC10489006':
                pattern = ''.join((fmt(a) + '\\s*±\\s*' + fmt(b) + '\\s*' for (a, b) in zip(means, sds)))
                assert re.search(pattern, txt), (product, 'complete row binding failed')
            rows.append(dict(study=st, doi=primary[st]['doi'], locator=table, product=product, material_class=cl, protocol=protocol, t=t, mean=mu, sd=sd, n=n, v=(sd / mu) ** 2 / n, eligible=eligible, note=note, origin=origin, resolution='PER_TOOTH' if 'crown' in protocol else 'PHENOMENOLOGICAL'))
    for (mat, cl, mu, sd) in [('ESS', 'leucite', [14.65, 48.65, 90.92, 154.02, 246.45], [2.25, 6.17, 11.32, 13.5, 31.35]), ('EMX', 'LS2', [41.29, 124.95, 209.29, 355.65, 709.0], [7.03, 15.6, 36.39, 100.76, 87.37]), ('LP', '3Y', [104.97, 348.76, 768.17, 1244.02, 1844.81], [19.71, 68.4, 140.31, 129.01, 216.77])]:
        add('PMC10004144', 'Table 2', mat, cl, 'free_disc_3ball_0.5', [0.4, 0.7, 1, 1.3, 1.6], mu, sd, 12)
    for (mat, cl, mu, sd) in [('IPS Empress Esthetic', 'leucite', [257.0, 424.3, 499.89], [52.6, 82.9, 73.89]), ('Vita Enamic', 'PICN', [449.7, 509.1, 576.6], [236.2, 42.35, 80.63]), ('IPS e.max Press', 'LS2', [456.1, 658.9, 1044.4], [67.79, 99.52, 111.2]), ('Ceramill Zolid HT', '4Y', [1086.1, 1640.0, 1569.0], [239.75, 200.33, 252.34])]:
        add('PMC10489006', 'Table 3', mat, cl, 'veneer_phantom_axial', [1, 1.5, 2], mu, sd, 10)
    for (mat, mu, sd) in [('Lava Ultimate', [1378.25, 2222.74, 1383.84], [232.76, 320.36, 208.54]), ('Vita Enamic', [1545.04, 1707.09, 1204.96], [331.74, 289.31, 130.5])]:
        add('PMC10286427', 'Table 2', mat, 'hybrid', 'implant_crown_45_confounded', [1, 2, 3], mu, sd, 8, False, '3mm switches custom to prefabricated abutment and proximal thickness; no clean exponent or plateau')
    for (table, protocol, vals, sds) in [('Table 2', 'natural_veneer_TC', [2131.0, 1919.0, 1413.0], [441.2, 306.2, 276.5]), ('Table 3', 'natural_veneer_TC_SAL', [1333.0, 1313.0, 1134.0], [44.6, 42.5, 115.4]), ('Table 4', 'natural_veneer_SAL', [1591.0, 1517.0, 1325.0], [82.1, 95.7, 163.6])]:
        for (mat, cl, mu, sd) in zip(['Cerasmart', 'Straumann Nice', 'Tetric CAD'], ['nanoceramic_resin', 'glass_ceramic', 'resin_composite'], vals, sds):
            add('PMC10020293', table, mat, cl, protocol, [0.5], [mu], [sd], 10)
    for (substrate, mu, sd) in [('dentin_analogue', [1393.0, 1789.0], [301.0, 233.0]), ('GIC_core', [428.0, 541.0], [93.0, 173.0]), ('RC_core', [1043.0, 1187.0], [339.0, 317.0])]:
        add('PMC10756807', 'Table 1', 'IPS e.max CAD', 'LS2', 'aged_LS2_disc_' + substrate, [0.5, 1], mu, sd, 10, note='not zirconia; crosshead 0.1 mm/min; 1e6 cycles plus thermocycling')
    for r in read('inputs/x1b_primary_groups.json'):
        if r['study'] != 'PMC10817558':
            continue
        abrasion = r['abrasion']
        cement = r['cement']
        prot = 'crown_30_Chen_' + r['material'] + '_' + abrasion + '_' + cement
        add(r['study'], 'Table 2', 'Katana HT 3Y' if r['material'] == '3Y' else 'Katana UTML 5Y', r['material'], prot, [r['thickness_mm']], [r['mean_N']], [r['sd_N']], r['n'], note='3.5mm indenter; NextDent C&B 2.1GPa; 0.5mm/min; 7 days water 37C; rubber interlayer')
    txt = (R / 'inputs/PMC8558575.txt').read_text()
    for r in read('inputs/x1b_matched.json')['groups']:
        if r['study'] != 'PMC8558575':
            continue
        assert str(r['characteristic_load_N']).rstrip('0').rstrip('.') in txt
        mu = r['characteristic_load_N'] * gamma(1 + 1 / r['Weibull_m'])
        assert abs(mu - r['mean_N']) < 1e-08
        rows.append(dict(study=r['study'], doi='10.4047/jap.2021.13.5.269', locator='Table 1 non-fatigued G' + str(r['thickness_mm']), product='inCoris TZI 3Y', material_class='3Y', protocol='crown_0_Prott_resin', t=r['thickness_mm'], mean=float(mu), sd=None, n=14, v=r['log_mean_variance'], eligible=True, note='6.36mm indenter; die modulus 8.6GPa manufacturer claim; 1.5mm/min; resin cement; mean derived via Weibull closure, zero parameter covariance', origin='Weibull_derived_mean', resolution='PER_TOOTH', F0=r['characteristic_load_N'], m=r['Weibull_m']))
    return rows

def main():
    t0 = time.perf_counter()
    rs = curated()
    dump('raw/CURATED_ROWS.json', rs)
    blocks = {}
    for r in rs:
        if r['eligible']:
            blocks.setdefault((r['study'], r['protocol'], r['product']), []).append(r)
    fitted = [dict(study=k[0], protocol=k[1], product=k[2], **slope(v)) for (k, v) in blocks.items() if len(set((r['t'] for r in v))) >= 2]
    dump('raw/SLOPES.json', fitted)
    fkeys = [k for (k, v) in blocks.items() if len(set((r['t'] for r in v))) >= 2]
    data = [r for k in fkeys for r in blocks[k]]
    X = np.zeros((len(data), 2 * len(fkeys)))
    for (i, r) in enumerate(data):
        j = fkeys.index((r['study'], r['protocol'], r['product']))
        X[i, j] = 1
        X[i, len(fkeys) + j] = math.log(r['t'])
    y = np.log([r['mean'] for r in data])
    w = 1 / np.array([r['v'] for r in data])
    beta = np.linalg.solve(X.T @ (w[:, None] * X), X.T @ (w * y))
    delta = max((abs(beta[len(fkeys) + j] - slope(blocks[k])['n']) for (j, k) in enumerate(fkeys)))
    folds = []
    h = read('PREREG_R1_PROTOCOL_TRANSFER.json')['metrics']['operational_halfwidth_log']
    for st in sorted({f['study'] for f in fitted}):
        train = [f for f in fitted if f['study'] != st]
        pts = []
        pool = sum((f['n'] / f['variance'] for f in train)) / sum((1 / f['variance'] for f in train))
        for f in fitted:
            if f['study'] != st:
                continue
            rows = sorted(blocks[st, f['protocol'], f['product']], key=lambda r: r['t'])
            ref = rows[0]
            matches = [q for q in train if q['protocol'] == f['protocol'] and q['product'] == f['product']]
            n = 2.0 if f['protocol'] == 'free_disc_3ball_0.5' else sum((q['n'] / q['variance'] for q in matches)) / sum((1 / q['variance'] for q in matches)) if matches else None
            for r in rows[1:]:
                x = math.log(r['t'] / ref['t'])
                actual = math.log(r['mean'] / ref['mean'])
                errs = {key: exp * x - actual for (key, exp) in [('practice_1.398', 1.398), ('inherited_1.5', 1.5), ('pooled_all_protocols', pool)]}
                pts.append(dict(product=r['product'], protocol=r['protocol'], locator='doi:' + r['doi'] + ' ' + r['locator'], reference_t_mm=ref['t'], heldout_t_mm=r['t'], actual_ratio=math.exp(actual), candidate_n=n, candidate_log_error=None if n is None else n * x - actual, candidate_covered=None if n is None else abs(n * x - actual) <= h, practice_log_errors=errs, heldout_calibration='ONE lowest-thickness mean, never target mean', resolution='PER_TOOTH' if 'crown' in r['protocol'] else 'PHENOMENOLOGICAL'))
        folds.append(dict(held_out_study=st, training_studies=sorted({f['study'] for f in train}), points=pts))
    free = [p for f in folds for p in f['points'] if p['protocol'] == 'free_disc_3ball_0.5']
    rms = lambda v: math.sqrt(sum((x * x for x in v)) / len(v))
    free_rmse = rms([p['candidate_log_error'] for p in free])
    bad_rmse = rms([math.log(p['heldout_t_mm'] / p['reference_t_mm']) - math.log(p['actual_ratio']) for p in free])
    import itertools
    ranking = []
    natural = [r for r in rs if r['study'] == 'PMC10020293']
    protocols = sorted({r['protocol'] for r in natural})
    products = sorted({r['product'] for r in natural})
    for (a, b) in itertools.combinations(products, 2):
        for (p, q) in itertools.combinations(protocols, 2):
            get = lambda mat, prot: next((r for r in natural if r['product'] == mat and r['protocol'] == prot))
            da = math.log(get(a, p)['mean'] / get(b, p)['mean'])
            db = math.log(get(a, q)['mean'] / get(b, q)['mean'])
            ranking.append(dict(products=[a, b], protocols=[p, q], log_ratios=[da, db], reversal=da * db < 0))
    pooled_resin = [r['mean'] for r in rs if r['material_class'] == 'resin_composite']
    pooled_glass = [r['mean'] for r in rs if r['material_class'] in ['LS2', 'glass_ceramic']]
    diag = dict(actual_TC_Nice_over_Tetric=1919 / 1413, geometric_pooled_glass_over_resin=math.exp(np.mean(np.log(pooled_glass)) - np.mean(np.log(pooled_resin))), status='CLASS_POOLING_DIAGNOSTIC; different products/geometries, not exact-product validation')
    predictions = []
    for (protocol, mat) in [('crown_0_Prott_resin', 'inCoris TZI 3Y'), ('crown_30_Chen_3Y_Yes_RMGI', 'Katana HT 3Y')]:
        rr = blocks[next((r['study'] for r in rs if r['protocol'] == protocol)), protocol, mat]
        sf = slope(rr)
        for (design, t) in [('D1', 0.8), ('M1', 1.0), ('M2', 1.5)]:
            supported = min((r['t'] for r in rr)) <= t <= max((r['t'] for r in rr))
            ratio = (t / 0.8) ** sf['n']
            predictions.append(dict(design=design, protocol=protocol, t_mm=t, ratio_point=ratio if supported else None, operational_ratio_window=[ratio * math.exp(-h), ratio * math.exp(h)] if supported and design != 'D1' else [1, 1] if design == 'D1' else None, validity='PROSPECTIVE_HYPOTHESIS_SINGLE_STUDY' if supported else 'UNKNOWN_OUTSIDE_SOURCE_SPAN', absolute_force_interval_N=None, exponent=sf['n'], exponent_kind='empirical protocol-specific closure', resolution='PER_TOOTH'))
    unknown = sum((p['candidate_n'] is None for f in folds for p in f['points']))
    total = sum((len(f['points']) for f in folds))
    out = dict(claim_type='information_link', free_disc_logRMSE=free_rmse, free_disc_wrong_n1_logRMSE=bad_rmse, folds=folds, rank_comparisons=ranking, source_corrected_ranking_comparison_count=len(ranking), rank_reversals=sum((p['reversal'] for p in ranking)), pool_diagnostic=diag, source_check='all curated mean/SD pairs adjacent in primary tables; four Weibull means independently reconstructed', controls=dict(block_GLS_max_slope_difference=delta, agreement=bool(delta < 1e-09), wrong_n1_rejected=bad_rmse > 0.3), LOSO_abstention=dict(unknown=unknown, total=total, fraction=unknown / total, reason='no same protocol and exact product survives held-out study'), predictions=predictions, gates=dict(free_disc_error=free_rmse <= 0.3, exact_protocol_replication=False, crown_LOSO_validation=False, narrow_statistical_crown_intervals=False), outcome='PARTIAL_FREE_DISC; CROWN_TRANSFER_UNKNOWN', elapsed_seconds=time.perf_counter() - t0)
    out.pop('elapsed_seconds')
    dump('raw/R1_RESULTS.json', out)
    state('R1_DECIDED', out['gates'], 'R2: endpoint calibration per exact lab protocol; hold interior design blind; prereg before computation')
    (R / 'HANDOFF_R1.md').write_text('R1: free-disc derived n=2 checked against primary Table 2; strict crown LOSO has no replicated protocol. No validated narrower absolute crown intervals. Implant plateau rejected as thickness/support confound. Next: freeze bracket operation using local endpoint forces and interior validation under the same product/support/contact/speed/cement/conditioning.\n')
    print(json.dumps({k: out[k] for k in ['outcome', 'free_disc_logRMSE', 'free_disc_wrong_n1_logRMSE', 'LOSO_abstention', 'controls', 'pool_diagnostic']}, indent=2))
if __name__ == '__main__':
    main()
