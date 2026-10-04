"""Primary measurements from local LIT_CROWN and already-local primary texts.
Every numerical group retains a table-cell or paragraph locator and raw source.
Missing restoration thickness or material class stays UNKNOWN.
"""
import re, json, csv, xml.etree.ElementTree as ET
from pathlib import Path
from common import R, sha, read, dump

def text(e):
    return ' '.join(''.join(e.itertext()).split())

def pairs(s):
    return [(float(a.replace(',', '')), float(b.replace(',', ''))) for (a, b) in re.findall('([\\d,]+(?:\\.\\d+)?)\\s*±\\s*([\\d,]+(?:\\.\\d+)?)', s)]

def table(root, label):
    return next((t for t in root.findall('.//table-wrap') if t.findtext('label') == label))

def expanded_rows(t):
    active = {}
    out = []
    for tr in t.findall('.//tbody/tr'):
        row = {}
        nxt = {}
        for (c, (value, left)) in active.items():
            row[c] = value
            if left > 1:
                nxt[c] = (value, left - 1)
        col = 0
        for cell in tr:
            while col in row:
                col += 1
            value = text(cell)
            row[col] = value
            span = int(cell.attrib.get('rowspan', '1'))
            if span > 1:
                nxt[col] = (value, span - 1)
            col += int(cell.attrib.get('colspan', '1'))
        active = nxt
        out.append([row[k] for k in sorted(row)])
    return out

def extract():
    corpus = list(map(json.loads, (R / 'inputs/literature/crown_fracture_papers.jsonl').read_text().splitlines()))
    byid = {r['pmcid']: r for r in corpus}
    rows = []
    audit = []

    def add(study, doi, locator, mean, sd, n, material, t, angle, cement, source, **extra):
        row = dict(row_id=f"{study}:{len([r for r in rows if r['study'] == study]) + 1}", study=study, doi=doi, locator=locator, mean_N=mean, sd_N=sd, n=n, material=material, thickness_mm=t, angle_deg=angle, cement=cement, source=source, source_sha256=sha(R / source), unit='N', quantity='group arithmetic mean fracture force', mechanism_origin='UNKNOWN', **extra)
        missing = [k for k in ['mean_N', 'sd_N', 'n', 'material', 'thickness_mm', 'angle_deg', 'cement'] if row[k] is None or row[k] == 'UNKNOWN']
        row['eligible_R1'] = not missing and material in ('3Y', '4Y', '5Y')
        row['exclusion_reason'] = None if row['eligible_R1'] else 'UNKNOWN ' + ','.join(missing) if missing else 'Material outside registered zirconia classes'
        rows.append(row)
    src = 'inputs/literature/PMC10817558.xml'
    root = ET.parse(R / src).getroot()
    for (i, c) in enumerate(expanded_rows(table(root, 'Table 2')), 1):
        (mean, sd) = pairs(c[4])[0]
        add('PMC10817558', '10.3390/ma17020365', f'Table 2 data row {i}; {c[0]}, {c[1]}mm, abrasion={c[2]}, {c[3]}', mean, sd, 10, c[0][:2], float(c[1]), 30, 'resin' if c[3] == 'Resin' else 'RMGI', src, abrasion=c[2], die_material='NextDent C&B', die_E_MPa=2100, indenter_diameter_mm=3.5, crosshead_mm_min=0.5, aging='24h water; no cyclic loading')
        if f'{mean:.2f}' not in byid['PMC10817558']['tables'][1]['text']:
            raise ValueError('XML/LIT measurement disagreement')
    src = 'inputs/literature/PMC10934854.xml'
    root = ET.parse(R / src).getroot()
    t2 = table(root, 'Table 2')
    rr = expanded_rows(t2)
    force = next((r for r in rr if 'Crown fracture force' in r[0]))
    mod = next((r for r in rr if 'Elastic modulus' in r[0]))
    for (i, name) in enumerate(['NextDent C&B', 'OnX', 'Lava Ultimate']):
        (mean, sd) = pairs(force[i + 1])[0]
        E = pairs(mod[i + 1])[0][0]
        add('PMC10934854', '10.3390/ma17051096', f'Table 2 crown fracture force row, column {name}', mean, sd, 12, '3Y', 0.8, 30, 'resin', src, die_material=name, die_E_MPa=E, indenter_diameter_mm=3.5, crosshead_mm_min=0.5, aging='24h water; no cyclic loading', crown_only_count=[8, 8, 10][i], crown_and_die_count=[4, 4, 2][i], mode_locator=f'Table 4 standardized {name}')
    rr = expanded_rows(table(root, 'Table 3'))
    force = next((r for r in rr if 'Crown fracture force' in r[0]))
    for (i, name) in enumerate(['Dentin', 'NextDent C&B replica']):
        (mean, sd) = pairs(force[i + 1])[0]
        add('PMC10934854', '10.3390/ma17051096', f'Table 3 crown fracture force row, column {name}', mean, sd, 12, '3Y', 0.8, 30, 'resin', src, die_material=name, die_E_MPa=None if i == 0 else 1951, indenter_diameter_mm=3.5, crosshead_mm_min=0.5, aging='24h water; no cyclic loading', crown_only_count=[8, 9][i], crown_and_die_count=[4, 3][i], mode_locator=f'Table 4 {name}')
    src = 'inputs/literature/PMC6642729.txt'
    s = (R / src).read_text()
    part = s.split('Table 1. Descriptive statistics')[1].split('Table 2.')[0]
    for (mat, label) in [('3Y', '3Y-PSZ'), ('5Y', '5Y-Z'), ('LS2', 'Lithium Disilicate')]:
        m = re.search(re.escape(label) + '\\s+(.+?)(?=5Y-Z|Lithium Disilicate|$)', part)
        vals = pairs(m.group(1))
        for (i, (mean, sd)) in enumerate(vals):
            cem = 'resin' if mat == 'LS2' or i % 2 == 0 else 'RMGI'
            abr = 'treated' if i < 2 else 'untreated'
            add('PMC6642729', '10.1111/jopr.13088', f'Table 1 {label}, {abr}, {cem}, numerical cell {i + 1}', mean, sd, 8, mat, 0.8, 0, cem, src, abrasion=abr, die_material='resin composite', die_E_MPa=None, indenter_diameter_mm=12, crosshead_mm_min=1, aging='100000 cycles at 100N', contact_layer='1.5mm rubber')
    src = 'inputs/literature/PMC4764450.txt'
    s = (R / src).read_text()
    part = s.split('Table 2 Critical splitting')[1]
    for (mat, label, t) in [('UNKNOWN', 'Zirconia', 1), ('LS2', 'Lithium disilicate', 2), ('composite', 'Dental composite', 2)]:
        vals = pairs(part.split(label, 1)[1])
        (mean, sd) = vals[0]
        add('PMC4764450', '10.1016/j.dental.2015.12.010', f'Table 2 {label}, EXPERIMENT P_S (kN); converted x1000', mean * 1000, sd * 1000, 6, mat, t, 0, 'resin', src, die_material='Z100', die_E_MPa=None, indenter_diameter_mm=6.4, crosshead_mm_min=0.1, aging='no cyclic aging', material_system=label, endpoint='critical splitting load; yttria class not explicit in read primary text')
    src = 'inputs/literature/PMC7274823.txt'
    s = (R / src).read_text()
    part = s.split('Table 2 Data for the failure load')[1].split('Table 3')[0]
    vals = [(float(a.replace(',', '')), float(b)) for (a, b) in re.findall('([\\d,]+(?:\\.\\d+)?)\\s*\\(([\\d.]+)\\)', part)]
    for (i, (mean, sd)) in enumerate(vals):
        add('PMC7274823', '10.1055/s-0040-1709342', f"Table 2 subgroup {['A1', 'A2', 'B1', 'B2'][i]}", mean, sd, 10, 'UNKNOWN', [1, 0.5, 1, 0.5][i], 0, 'resin', src, material_system='IPS e.max ZirCAD MT', die_material='NiCr', indenter_diameter_mm=4, crosshead_mm_min=0.5, aging='500 thermocycles', margin='chamfer' if i < 2 else 'shoulderless', mode_locator='Table 4; occlusal origin in fractography', mechanism_comment='Fractographic text reports occlusal origin; intaglio origin not established')
    src = 'inputs/literature/PMC9081244.txt'
    s = (R / src).read_text()
    part = s.split('Table 2. Mean and standard deviation')[1]
    for (mat, label) in [('3Y', '3Y-PSZ'), ('4Y', '4Y-PSZ')]:
        vals = pairs(part.split(label, 1)[1])[:2]
        for (i, (mean, sd)) in enumerate(vals):
            add('PMC9081244', '10.1111/jerd.12907', f'Table 2 {label}, aging column {i + 1}', mean, sd, 10, mat, None, 0, 'resin', src, die_material='Trilor', indenter_diameter_mm=4, crosshead_mm_min=None, aging='thermocycling only' if i == 0 else 'thermomechanical', thickness_comment='Primary text states 1.5mm preparation reduction; actual restoration thickness UNKNOWN')
    r = byid['PMC10478297']
    src = 'inputs/literature/crown_fracture_papers.jsonl'
    tt = next((t for t in r['tables'] if t['label'] == 'Table 3'))['text']
    for (i, (mean, sd)) in enumerate(pairs(tt)):
        add(r['pmcid'], r['doi'], f'Table 3 numerical mean/SD cell {i + 1}', mean, sd, 10, '5Y', None, 0, 'resin', src, material_system='Katana HTML; Table 1 says 5mol%', aging='10000 thermocycles + 240000 cycles at50N', indenter_diameter_mm=5, crosshead_mm_min=0.5)
    r = byid['PMC10560138']
    p = next((p for p in r['key_paragraphs'] if p.startswith('Descriptive statistics for fracture')))
    for (i, (mean, sd)) in enumerate(pairs(p)):
        add(r['pmcid'], r['doi'], f'Figure 4 and Results paragraph starting Descriptive statistics; mean/SD pair {i + 1}', mean, sd, 9 if i % 2 == 0 else 6, '3Y', None, 0, 'resin', src, fabrication='milled' if i < 2 else 'printed', aging='no cycling' if i % 2 == 0 else '1200000 cycles', thickness_comment='1.5mm preparation reduction is not a measured restoration thickness')
    audit.append(dict(study='PMC8558575', locator='Table 1', decision='EXCLUDED_R1', reason='Characteristic Weibull strength is not an arithmetic mean; no measured mean/SD in the local table; G2 control reused from earlier work', source_sha256=sha(R / 'inputs/literature/PMC8558575.txt')))
    audit.append(dict(study='PMC10264088', locator='Table 1', decision='EXCLUDED_R1', reason='0.5mm copings are not full-contour crowns'))
    audit.append(dict(study='BT-DW48-AUTO-98ae75d8c2a13a', decision='READ_EXISTING_WORK', reason='Completed LIT_CROWN census, not a primary extracted calibration; rate coefficient remains confounded'))
    dump('raw/LITERATURE_GROUPS.json', rows)
    dump('raw/EXTRACTION_AUDIT.json', dict(corpus_records=len(corpus), rows=len(rows), eligible_rows=sum((r['eligible_R1'] for r in rows)), eligible_studies=sorted({r['study'] for r in rows if r['eligible_R1']}), excluded_studies=audit))
    keys = ['row_id', 'study', 'doi', 'locator', 'material', 'thickness_mm', 'angle_deg', 'cement', 'mean_N', 'sd_N', 'n', 'eligible_R1', 'exclusion_reason']
    with (R / 'raw/LITERATURE_GROUPS.csv').open('w') as f:
        w = csv.DictWriter(f, keys, extrasaction='ignore')
        w.writeheader()
        w.writerows(rows)
    return rows
if __name__ == '__main__':
    rows = extract()
    print(read('raw/EXTRACTION_AUDIT.json'))
