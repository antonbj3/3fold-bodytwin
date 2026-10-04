from common import *
import csv, collections, math
PROTOCOL_FIELDS = ['product', 'specimen_geometry', 'support', 'cement', 'angle_deg', 'indenter_diameter_mm', 'crosshead_mm_min', 'contact_layer', 'aging']

def calibrated_force(anchors, query):
    """A laboratory API. Anchors carry full protocols and actual measured endpoint forces."""
    if len(anchors) != 2:
        return {'status': 'UNKNOWN_NEEDS_TWO_ENDPOINTS'}
    for key in PROTOCOL_FIELDS:
        val = query.get(key)
        if val is None or val == 'UNKNOWN' or any((a.get(key) != val for a in anchors)):
            return {'status': 'UNKNOWN_SETUP_MISMATCH', 'field': key}
    (lo, hi) = sorted(anchors, key=lambda a: a['thickness_mm'])
    t = query.get('thickness_mm')
    if t is None or not 0 < lo['thickness_mm'] < t < hi['thickness_mm']:
        return {'status': 'UNKNOWN_OUTSIDE_CALIBRATED_INTERIOR'}
    if any((a.get('mean_force_N', 0) <= 0 for a in anchors)):
        return {'status': 'UNKNOWN_MEASUREMENT'}
    w = math.log(t / lo['thickness_mm']) / math.log(hi['thickness_mm'] / lo['thickness_mm'])
    f = lo['mean_force_N'] ** (1 - w) * hi['mean_force_N'] ** w
    return dict(status='CONDITIONAL_CALIBRATION', mean_force_N=f, operational_window_N=[f * math.exp(-0.2), f * math.exp(0.2)], rigorous_enclosure='MISSING', scope='Only the fully matched setup; no patient load or safety margin inferred')

def predict():
    start = time.perf_counter()
    cement = list(csv.DictReader((DATA / 'literature/CEMENT.csv').open()))
    crown = read(DATA / 'literature/CROWN_CURATED.json')
    out = []
    rejected = []
    regional = collections.defaultdict(list)
    for r in cement:
        if r['region'] in ['marginal', 'axial', 'occlusal']:
            regional[r['study'], r['region']].append(dict(value_um=float(r['measured_mean_um']), sd_um=float(r['reported_sd_um']) if r['reported_sd_um'] else None, doi=r['doi'], table=r['source_table'], row=r['source_row'], n=r['n_specimens'], arm=r['arm'], method=r['method'], state=r['state']))
    cb = collections.defaultdict(list)
    for r in cement:
        if r['region'] == 'marginal' and r['internal_spacer_um'] and (float(r['internal_spacer_um']) > 0):
            cb[r['study'], r['method'], r['state']].append(r)
    for (key, rs) in cb.items():
        rs = sorted(rs, key=lambda r: float(r['internal_spacer_um']))
        doses = [float(r['internal_spacer_um']) for r in rs]
        if len(set(doses)) < 3 or len(set(doses)) != len(doses):
            rejected.append(dict(kind='cement', key=key, rows=len(rs), reason='fewer than3 unique CAD doses or duplicate-arm protocol ambiguity'))
            continue
        (lo, hi) = rs[:2]
        (x1, x2) = map(float, [lo['internal_spacer_um'], hi['internal_spacer_um']])
        (y1, y2) = map(float, [lo['measured_mean_um'], hi['measured_mean_um']])
        b = (y1 - y2) / (1 / x1 - 1 / x2)
        a = y1 - b / x1
        for r in rs[2:]:
            x = float(r['internal_spacer_um'])
            pred = a + b / x
            direct = float(np.array([1, 1 / x]) @ np.linalg.solve([[1, 1 / x1], [1, 1 / x2]], [y1, y2]))
            out.append(dict(kind='cement', study=key[0], row_id=r['row_id'], x=x, predicted=pred, control=direct, observed=float(r['measured_mean_um']), units='um', region='marginal', resolution='PER_SURFACE_REGION', doi=r['doi'], locator=r['source_table'] + ' row ' + r['source_row'], calibration_rows=[lo['row_id'], hi['row_id']]))
    blocks = collections.defaultdict(list)
    for r in crown:
        if r['eligible'] and r['protocol'].startswith('crown'):
            blocks[r['study'], r['protocol'], r['product']].append(r)
    for (key, rs) in blocks.items():
        rs = sorted(rs, key=lambda r: r['t'])
        if len(set((r['t'] for r in rs))) < 3:
            rejected.append(dict(kind='crown', key=key, rows=len(rs), reason='fewer than3 distinct protocol-matched thicknesses'))
            continue
        (lo, hi) = (rs[0], rs[-1])
        for r in rs[1:-1]:
            w = math.log(r['t'] / lo['t']) / math.log(hi['t'] / lo['t'])
            pred = lo['mean'] ** (1 - w) * hi['mean'] ** w
            direct = math.exp(float(np.array([1, math.log(r['t'])]) @ np.linalg.solve([[1, math.log(lo['t'])], [1, math.log(hi['t'])]], np.log([lo['mean'], hi['mean']]))))
            out.append(dict(kind='crown', study=key[0], protocol=key[1], product=key[2], x=r['t'], predicted=pred, control=direct, observed=r['mean'], units='N', resolution='PER_TOOTH', doi=r['doi'], locator=r['locator'], calibration_thicknesses_mm=[lo['t'], hi['t']], quantity_origin=r['origin']))
    freeze(ROOT / 'FROZEN_LITERATURE_PREDICTIONS.json', dict(rows=out, rejected_blocks=rejected, regional_observation_sets=[dict(study=k[0], region=k[1], observations=v, resolution='PER_SURFACE_REGION', scope='Distribution of reported group means, NOT within-crown pointwise film distribution') for (k, v) in regional.items()], input_hashes={p.name: sha(p) for p in (DATA / 'literature').iterdir() if p.is_file()}, seconds=time.perf_counter() - start, retrospective=True))

def evaluate():
    from integrity import verify
    verify()
    p = read(ROOT / 'FROZEN_LITERATURE_PREDICTIONS.json')
    rows = []
    ctrl = []
    for r0 in p['rows']:
        r = dict(r0)
        error = abs(r['observed'] - r['predicted']) if r['kind'] == 'cement' else abs(math.log(r['observed'] / r['predicted']))
        threshold = 20 if r['kind'] == 'cement' else 0.2
        r.update(error=error, gate=error <= threshold)
        rows.append(r)
        bad = r['predicted'] + 100 if r['kind'] == 'cement' else 1.6 * r['predicted']
        e = abs(bad - r['predicted']) if r['kind'] == 'cement' else abs(math.log(bad / r['predicted']))
        ctrl.append(dict(study=r['study'], kind=r['kind'], wrong_value_rejected=e > threshold, numerical_control_matches=abs(r['control'] - r['predicted']) < 1e-08))
    protocol = dict(product='fixture_3Y', specimen_geometry='same_tooth', support='same_die', cement='same_cement', angle_deg=30, indenter_diameter_mm=3.5, crosshead_mm_min=0.5, contact_layer='same_rubber', aging='same_water7d')
    anchors = [dict(protocol, thickness_mm=0.5, mean_force_N=500), dict(protocol, thickness_mm=1.5, mean_force_N=1500)]
    query = dict(protocol, thickness_mm=1.0)
    valid = calibrated_force(anchors, query)
    controls = []
    for k in PROTOCOL_FIELDS:
        bad = dict(query)
        bad[k] = 'WRONG'
        result = calibrated_force(anchors, bad)
        controls.append(dict(field=k, rejected=result['status'] == 'UNKNOWN_SETUP_MISMATCH'))
    for k in PROTOCOL_FIELDS:
        bad = dict(query)
        bad.pop(k)
        result = calibrated_force(anchors, bad)
        controls.append(dict(field=k + '_missing', rejected=result['status'] == 'UNKNOWN_SETUP_MISMATCH'))
    generated = calibrated_force([], dict(thickness_mm=1.0, specimen_geometry='new_generated_crown'))
    out = dict(round='R4', claim_type='capability', decision='PROTOCOL_CONDITIONED_PORTS_WITH_PHYSICAL_UNKNOWN', rows=rows, regional_observation_sets=p['regional_observation_sets'], rejected_blocks=p['rejected_blocks'], controls=ctrl, protocol_controls=controls, positive_fixture_control=valid, generated_crown_result=generated, physical_candidate_cement_distribution='UNKNOWN: no regional matched process or spatial film measurements', fracture_rigorous_enclosure='MISSING: +/-0.2log is an operational test band, not rigorous enclosure', external_referent=read(ROOT / 'PREREG_R4.json')['external_referent'], dropout={'cement_input_rows': sum((1 for _ in csv.DictReader((DATA / 'literature/CEMENT.csv').open()))), 'primary_region_cells': sum((len(r['observations']) for r in p['regional_observation_sets'])), 'rejected_calibration_blocks': len(p['rejected_blocks']), 'no_assumption_of_independent_group_cells': True})
    dump(ROOT / 'rounds/R4.json', out)
    print('calibration targets', len(rows), 'wrong protocol controls', len(controls))
if __name__ == '__main__':
    if '--freeze' in sys.argv:
        predict()
    else:
        evaluate()
