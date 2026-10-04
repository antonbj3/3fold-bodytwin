"""Local primary-source extraction; neither ambient T nor motor limit is a measurement of drilling work."""
import re, csv, xml.etree.ElementTree as ET
from common import *
IDS = ['PMC10173627', 'PMC10906326', 'PMC10069805', 'PMC10130290', 'PMC10299697', 'PMC10889961', 'PMC10613809', 'PMC10110507', 'PMC10603942']

def txt(e):
    return ' '.join(''.join(e.itertext()).split())

def nums(s):
    s = s.replace('−', '-').replace('–', '-')
    s = re.sub('-\\s+(?=\\d)', '-', s)
    return [float(v) for v in re.findall('[-+]?\\d+(?:\\.\\d+)?', s)]

def main():
    verify('PREREG_R1.json')
    trees = {}
    tables = {}
    manifest = []
    records = []
    for sid in IDS:
        p = CORPUS / (sid + '.xml')
        r = ET.parse(p).getroot()
        trees[sid] = r
        ts = {t.get('id'): [[txt(c) for c in row] for row in t.findall('.//tr')] for t in r.findall('.//table-wrap')}
        tables[sid] = ts
        doi = next((a.text for a in r.findall('.//article-id') if a.get('pub-id-type') == 'doi'), None)
        manifest.append(dict(study=sid, path=str(p), sha256=sha(p), doi=doi, title=txt(r.find('.//article-title')), license=[txt(e) for e in r.findall('.//license')]))

    def add(sid, tid, row, col, group, quantity, value, **kw):
        locator = f"{sid}.xml//table-wrap[@id='{tid}']//tr[{row + 1}]/*[{col + 1}]"
        d = dict(id=f'{sid}:{tid}:{row}:{col}:{quantity}', study=sid, doi=next((s['doi'] for s in manifest if s['study'] == sid)), table_id=tid, table_locator=locator, raw_cell=tables[sid][tid][row][col], group=group, quantity=quantity, value_C=float(value), resolution='POPULATION', measurement_support='PER_POINT', full_temperature_history=False, calibration_use=False, rejection_reasons=[], **kw)
        records.append(d)
        return d
    sid = 'PMC10173627'
    rows = tables[sid]['Tab1']
    tissue = None
    for (i, row) in enumerate(rows):
        if row and row[0] == 'Overall':
            break
        if row and row[0] in ['Cortical', 'Cancellous']:
            tissue = row[0].lower()
        if 'Limiting' in row or 'Non-limiting' in row:
            off = 1 if row[0] in ['Cortical', 'Cancellous'] else 0
            guide = row[off]
            y = float(row[off + 1])
            sd = float(row[off + 2])
            d = add(sid, 'Tab1', i, off + 1, guide + '_' + tissue, 'mean_peak_rise', y, sd_C=sd, n=22, rpm=800, diameter_mm=2.2, feed_mm_s=None, cooling='external', flow_ml_min=26.67, guide=guide, bone_class=tissue, baseline_C=21.0, sensor_wall_distance_mm=1.0, sensor_depth_mm=1.5 if tissue == 'cortical' else 7.0, drill_depth_mm=9.0, wear_uses_max=10)
            d['calibration_use'] = True
    sid = 'PMC10906326'
    for (i, row) in enumerate(tables[sid]['tbl2'][2:], 2):
        feed = float(row[0]) / 60
        rpm = float(row[1])
        for (col, method) in [(2, 'conventional'), (3, 'ultrasonic')]:
            n = nums(row[col])
            force = nums(row[col + 2])
            d = add(sid, 'tbl2', i, col, f'{method}_{rpm:g}_{feed:g}', 'mean_peak_rise', n[0] - 22.0, temperature_peak_C=n[0], sd_C=n[1], n=3, force_N=force[0], force_sd_N=force[1], rpm=rpm, diameter_mm=4.0, feed_mm_s=feed, cooling=None, flow_ml_min=None, guide=None, bone_class='mixed_ovine_vertebra', baseline_C=22.0, sensor_wall_distance_mm=0.5, sensor_depth_mm=10.0, drill_depth_mm=None, wear_uses_max=None)
            d['calibration_use'] = method == 'conventional'
            if method == 'ultrasonic':
                d['rejection_reasons'].append('ultrasonic extra work not represented by rotary K3 source')
    sid = 'PMC10069805'
    for (i, row) in enumerate(tables[sid]['table2'][2:], 2):
        baseline = float(row[4])
        y = float(row[5])
        se = float(row[6])
        d = add(sid, 'table2', i, 5, row[0], 'mean_post_final_rise', y - baseline, temperature_post_C=y, baseline_C=baseline, se_C=se, n=10, rpm=800, diameter_mm=3.6 if row[0] != 'Alveolar expander' else 3.5, feed_mm_s=None, cooling='external', flow_ml_min=None, guide=None, bone_class='femoral_bone', sensor_wall_distance_mm=0.5, sensor_depth_mm=5.75, drill_depth_mm=11.5, wear_uses_max=None)
        d['rejection_reasons'] = ['post-final measurement is not continuous maximum; feed and prior sequence timing absent']
    sid = 'PMC10130290'
    rows = tables[sid]['Tab1']
    methods = ['standard', 'standard', 'standard', 'velodrill', 'velodrill', 'velodrill', 'piezo', 'piezo', 'piezo']
    for (row, depth) in [(3, 4.5), (4, 9.0)]:
        for j in range(1, 10):
            label = chr(64 + j)
            n = nums(rows[row][j])
            cold = j not in [3, 6, 9]
            d = add(sid, 'Tab1', row, j, label + f'_{depth:g}', 'mean_peak_rise', n[0], sd_C=n[1], n=35, rpm=800 if j <= 6 else None, diameter_mm=2.2 if j <= 6 else None, feed_mm_s=None, cooling='external', flow_ml_min=190 if j in [2, 5, 8] else 120, coolant_C=5.0 if cold else 20.0, guide=None, bone_class='PU_20pcf_plus_50pcf', baseline_C=None, baseline_nominal_C=22.0, sensor_wall_distance_mm=None, sensor_depth_mm=depth, drill_depth_mm=9.0, wear_uses_max=35, method=methods[j - 1])
            d['measurement_support'] = 'PER_SURFACE_REGION'
            if j <= 2:
                d['rejection_reasons'].append('duplicated Group B header: A/B column identities not uniquely verified')
            if j > 6:
                d['rejection_reasons'].append('piezoelectric extra power and different duty cycle')
            d['rejection_reasons'].append('PU and exposed IR ROI lack calibrated bone material and spatial observer')
    sid = 'PMC10889961'
    for (i, row) in enumerate(tables[sid]['medicina-60-00239-t001'][1:], 1):
        col = 2 if len(row) == 4 else 1
        n = nums(row[col])
        (region, guide) = row[0].split('-')
        d = add(sid, 'medicina-60-00239-t001', i, col, row[0], 'mean_peak_absolute', n[0], sd_C=n[1], n=10, rpm=1200, diameter_mm=3.5, feed_mm_s=None, cooling='none' if guide == 'NI' else 'external', flow_ml_min=None, guide=guide, bone_class='PU_20pcf_plus_40pcf', baseline_C=None, baseline_nominal_C=21.5, sensor_wall_distance_mm=None, sensor_depth_mm=None, drill_depth_mm=10.0, wear_uses_max=None, region=region)
        d['measurement_support'] = 'PER_SURFACE_REGION'
        d['rejection_reasons'] = ['PU material, missing feed, IR surface geometry and per-specimen baseline prohibit bone-class transfer']
    sid = 'PMC10299697'
    rows = tables[sid]['jcm-12-03944-t001']
    for (row, depth, kind) in [(1, 2.0, 'unspecified_mean_absolute'), (3, 2.0, 'group_extreme_absolute'), (4, 6.0, 'unspecified_mean_absolute'), (6, 6.0, 'group_extreme_absolute')]:
        offset = 2 if len(rows[row]) == 6 else 1
        for j in range(4):
            y = nums(rows[row][offset + j])[0]
            d = add(sid, 'jcm-12-03944-t001', row, offset + j, 'ABCD'[j] + f'_{depth:g}', kind, y, n=12, rpm=800, diameter_mm=2.2, feed_mm_s=None, cooling='none' if j == 3 else 'external', flow_ml_min=0 if j == 3 else 35, guide=['entry_exit', 'entry_only', 'handpiece', 'dry'][j], bone_class='bovine_rib_unspecified_region', baseline_C=None, sensor_wall_distance_mm=1.0, sensor_depth_mm=depth, drill_depth_mm=8.0, wear_uses_max=None, sensor_accuracy_C=0.8, sensor_logging_s=2.0)
            d['rejection_reasons'] = ['mean/max aggregation not same as mean individual peak rise; numeric baseline and feed absent']
    sid = 'PMC10613809'
    rows = tables[sid]['T1']
    deltas = tables[sid]['T2']
    consistency = []
    for (i, row) in enumerate(rows[2:], 2):
        for (j, cell) in enumerate(row):
            n = nums(cell)
            reported = float(deltas[i][j])
            valid = len(n) == 2 and '-' in cell
            n = [abs(v) for v in n]
            calculated = n[1] - n[0] if valid else None
            consistency.append(dict(study=sid, locator=f'T1/tr[{i + 1}]/td[{j + 1}] versus T2/tr[{i + 1}]/td[{j + 1}]', raw=cell, calculated_rise_C=calculated, reported_rise_C=reported, defect_C=None if calculated is None else calculated - reported, consistent=bool(valid and abs(calculated - reported) <= 0.15)))
    dump('raw/source_manifest.json', manifest)
    dump('raw/source_tables.json', tables)
    dump('raw/measurements.json', records)
    dump('raw/indian_table_consistency.json', consistency)
    fields = ['id', 'study', 'table_locator', 'group', 'quantity', 'value_C', 'calibration_use', 'bone_class', 'measurement_support']
    with (ROOT / 'raw/measurements.csv').open('w') as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction='ignore')
        w.writeheader()
        w.writerows(records)
    anchors = [('PMC10173627', 'Limiting_cortical', 0.82), ('PMC10173627', 'Non-limiting_cortical', 0.3), ('PMC10906326', 'conventional_500_0.5', 3.8), ('PMC10906326', 'ultrasonic_2000_1', 7.3), ('PMC10889961', 'Cor-Win', 21.451), ('PMC10889961', 'Cor-Sol', 25.005), ('PMC10299697', 'A_2', 22.1), ('PMC10299697', 'C_2', 24.9)]

    def parity(items):
        return all((abs(next((d['value_C'] for d in items if d['study'] == s and d['group'] == g)) - v) < 1e-10 for (s, g, v) in anchors))
    assert parity(records)
    corrupt = [dict(d) for d in records]
    corrupt[0]['value_C'] += 10
    assert not parity(corrupt)
    counts = {s: sum((d['study'] == s for d in records)) for s in IDS}
    counts['PMC10613809_internal_pairs'] = len(consistency)
    dump('raw/extraction_checks.json', dict(manual_parity=True, injected_plus10_rejected=True, counts=counts, strict_complete_rows=sum((all((d.get(k) is not None for k in ['feed_mm_s', 'cooling', 'drill_depth_mm', 'sensor_wall_distance_mm', 'baseline_C'])) and d['calibration_use'] for d in records)), numeric_table_candidates=len(records) + len(consistency), retained_for_nominal_temperature_fit=sum((d['calibration_use'] for d in records)), rejected_from_nominal_fit=sum((not d['calibration_use'] for d in records)) + len(consistency), indian_inconsistent_pairs=sum((not d['consistent'] for d in consistency))))
    print('Extracted', len(records), 'native measurement records; cross-table pairs', len(consistency))
if __name__ == '__main__':
    main()
