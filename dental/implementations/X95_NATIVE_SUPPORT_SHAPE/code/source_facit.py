from common import *
import math, xml.etree.ElementTree as ET
X1C = D / 'results/LANE_X1C_CROWN_LITERATURE'

def run():
    t0 = time.perf_counter()
    xml = X1 / 'inputs/literature/PMC10817558.xml'
    root = ET.parse(xml).getroot()
    table = root.find(".//table-wrap[@id='materials-17-00365-t002']")
    carry = {}
    rows = []
    for tr in table.findall('.//tbody/tr'):
        vals = {i: v[0] for (i, v) in carry.items()}
        newcarry = {i: (v[0], v[1] - 1) for (i, v) in carry.items() if v[1] > 1}
        column = 0
        for td in tr.findall('td'):
            while column in vals:
                column += 1
            value = ''.join(td.itertext()).strip()
            vals[column] = value
            rs = int(td.get('rowspan', '1'))
            if rs > 1:
                newcarry[column] = (value, rs - 1)
            column += 1
        carry = newcarry
        (mat, t, abr, cement, force) = [vals[i] for i in range(5)]
        (mu, sd) = force.replace('*', '').split('±')
        rows.append(dict(material=mat, thickness_mm=float(t), abrasion=abr, cement=cement, mean_N=float(mu.strip()), sd_N=float(sd.strip()), locator='doi:10.3390/ma17020365 Table2 row' + str(len(rows) + 1)))
    keep = [r for r in rows if r['thickness_mm'] == 1.0 and r['abrasion'] == 'Yes' and (r['cement'] == 'RMGI')]
    assert len(keep) == 2
    by = {r['material']: r for r in keep}
    ratio = by['3Y-Z']['mean_N'] / by['5Y-Z']['mean_N']
    truth = read(X1C / 'raw/MATCHED_MATERIAL_PORT.json')
    err = abs(ratio - truth['ratio'])
    source_rows = read(X1C / 'raw/CURATED_ROWS.json')
    blocks = {}
    for r in source_rows:
        if r['eligible']:
            blocks.setdefault((r['study'], r['protocol'], r['product']), []).append(r)
    protocol_results = []
    rejections = []
    h = read(X1C / 'PREREG_R2_BRACKET_CALIBRATION.json')['metrics']['operational_halfwidth_log']
    for (key, rs) in sorted(blocks.items()):
        rr = sorted(rs, key=lambda r: r['t'])
        if len(set((r['t'] for r in rr))) < 3:
            rejections.append({'source_protocol': list(key), 'reason': 'fewer than3 source thicknesses', 'rows': len(rr)})
            continue
        (lo, hi) = (rr[0], rr[-1])
        for target in rr[1:-1]:
            w = math.log(target['t'] / lo['t']) / math.log(hi['t'] / lo['t'])
            pred = lo['mean'] ** (1 - w) * hi['mean'] ** w
            e = math.log(target['mean'] / pred)
            protocol_results.append(dict(study=key[0], protocol=key[1], product=key[2], thickness_mm=target['t'], observed_N=target['mean'], predicted_N=pred, log_error=e, covered=abs(e) <= h, locator='doi:' + target['doi'] + ' ' + target['locator'], origin=target['origin'], resolution=target['resolution']))
    studies = sorted(set((r['study'] for r in protocol_results)))
    rmse = math.sqrt(np.mean([np.mean([r['log_error'] ** 2 for r in protocol_results if r['study'] == st]) for st in studies]))
    coverage = float(np.mean([np.mean([r['covered'] for r in protocol_results if r['study'] == st]) for st in studies]))
    prior = read(X1C / 'raw/R2_RESULTS.json')['stats']['candidate']
    raw_error = max(abs(rmse - prior['study_balanced_logRMSE']), abs(coverage - prior['study_balanced_coverage']))
    facets = read(R / 'raw/R1_NATIVE_PORT.json')
    force_fault = by['3Y-Z']['mean_N'] * 2
    force_matches = lambda value: abs(value - by['3Y-Z']['mean_N']) <= 0.005
    source_protocol = {'geometry': 'premolar_coping', 'die_E_MPa': 2100.0, 'ball_diameter_mm': 3.5, 'angle_deg': 30.0, 'interlayer': 'rubber', 'speed_mm_min': 0.5, 'conditioning': 'water7days37C'}
    d1_protocol = {'geometry': 'STS_molar_D1', 'die_E_MPa': 18000.0, 'ball_diameter_mm': 5.0, 'angle_deg': 30.0, 'interlayer': None, 'speed_mm_min': None, 'conditioning': None}
    accepts_protocol = lambda protocol: all((protocol.get(k) == v for (k, v) in source_protocol.items()))
    out = {'claim_type': 'capability', 'external_referent': {'kind': 'independent_measurement', 'locator': 'doi:10.3390/ma17020365 Tables1-2; original X1C source-specific DOI/table rows retained below', 'compared_quantity': 'source mean fracture forces and exact-protocol contrasts, never D1 strength', 'refutes_us': True}, 'Chen_1mm_matched_source': {'source_rows': keep, 'ratio_3Y_5Y': ratio, 'X1C_ratio_error': err, 'source_snapshot_sha256': sha(xml), 'protocol': 'premolar coping; NextDent C&B2.1GPa; ball3.50mm;30deg;rubber interlayer;.5mm/min;7days water37C', 'origin': 'source independent group measurement', 'resolution': 'POPULATION', 'transport_to_D1': 'UNKNOWN: D1 molar/18GPa die/5mm ball/source force-only static port; no matched specimen fracture'}, 'original_X1C_bracket_recalculation': {'target_points': len(protocol_results), 'study_balanced_logRMSE': rmse, 'study_balanced_coverage': coverage, 'original_metrics_difference': raw_error, 'original_gate': {'logRMSE': rmse <= 0.2, 'coverage': coverage >= 0.9}, 'scope': 'Recalculation of already audited source rows; not new heldout validation and no pooling to infer D1', 'rows': protocol_results, 'dropout': {'blocks_considered': len(blocks), 'blocks_rejected': len(rejections), 'fraction': len(rejections) / len(blocks), 'reasons': rejections}}, 'faults': {'doubled_published_force_N': force_fault, 'doubled_force_rejected': not force_matches(force_fault), 'valid_force_accepted': force_matches(by['3Y-Z']['mean_N']), 'source_protocol_positive_accepted': accepts_protocol(source_protocol.copy()), 'source_protocol_transfer_rejected': not accepts_protocol(d1_protocol), 'source_protocol_fields': source_protocol, 'attempted_D1_fields': d1_protocol}, 'cost': cost(t0)}
    dump(R / 'raw/EXTERNAL_FACIT.json', out)
    print(json.dumps(clean({k: out[k] for k in ['Chen_1mm_matched_source', 'faults', 'cost']}), indent=2))
if __name__ == '__main__':
    run()
