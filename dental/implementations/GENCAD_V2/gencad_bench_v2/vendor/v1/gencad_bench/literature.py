"""Reproducible extraction from original local XML; no guessed merged columns."""
from dental_release.paths import expand as _release_expand
import re, json
from pathlib import Path
import xml.etree.ElementTree as ET
from .io import ROOT, sha, dump
XML = Path(_release_expand('@DENTAL_CORPUS_ROOT@/europepmc/fulltext'))
STUDIES = {'PMC10817558': '10.3390/ma17020365', 'PMC10478297': '10.1186/s12903-023-03365-5', 'PMC10413631': '10.1186/s12903-023-03247-w'}

def text(x):
    return ' '.join(''.join(x.itertext()).split())

def expanded_rows(table):
    carry = {}
    out = []
    for tr in table.findall('.//tr'):
        cells = list(tr)
        row = []
        col = 0
        ci = 0
        new = {}
        while ci < len(cells) or any((k >= col for k in carry)):
            if col in carry:
                (value, remaining) = carry[col]
                row.append(value)
                if remaining > 1:
                    new[col] = (value, remaining - 1)
                col += 1
                continue
            if ci >= len(cells):
                break
            cell = cells[ci]
            ci += 1
            value = text(cell)
            rs = int(cell.get('rowspan', '1'))
            cs = int(cell.get('colspan', '1'))
            for _ in range(cs):
                row.append(value)
                if rs > 1:
                    new[col] = (value, rs - 1)
                col += 1
        out.append(row)
        carry = new
    return out

def mean_sd(value):
    m = re.search('(\\d+(?:\\.\\d+)?)\\s*±\\s*(\\d+(?:\\.\\d+)?)', value)
    if not m:
        raise ValueError('Ambiguous measurement: ' + value)
    return (float(m[1]), float(m[2]))

def table(tree, label):
    t = next((t for t in tree.findall('.//table-wrap') if t.findtext('label') == label))
    return (t.get('id'), expanded_rows(t))

def extract():
    trees = {p: ET.parse(XML / (p + '.xml')) for p in STUDIES}
    records = []
    tables = {}
    sources = []
    for (pid, tr) in trees.items():
        sources.append(dict(pmcid=pid, path=str(XML / (pid + '.xml')), sha256=sha(XML / (pid + '.xml')), doi=STUDIES[pid], licence='See original article licence element'))

    def row(pid, table_id, label, values, **kw):
        (mean, sd) = mean_sd(values)
        return dict(row_id=f'{pid}_{len(records):02d}', study=pid, doi=STUDIES[pid], locator=f'https://pmc.ncbi.nlm.nih.gov/articles/{pid}/#{table_id}', table_id=table_id, table_row=label, mean=mean, sd=sd, quantity='crown_fracture_force', unit='N', provenance='published independent in-vitro group measurement', **kw)
    pid = 'PMC10817558'
    (tid, rows) = table(trees[pid], 'Table 2')
    tables[pid] = rows
    for r in rows[1:]:
        (material, thick, abrasion, cement, value) = r
        records.append(row(pid, tid, ' | '.join(r[:4]), value, n=10, material_product='Katana HT' if material == '3Y-Z' else 'Katana UTML', material_reported=material, thickness_mm=float(thick), specimen='monolithic_premolar_crown', support='3D-printed resin die', cement=cement, particle_abrasion=abrasion, load_angle_deg=30, ageing='no ageing reported in extracted protocol', loading='compressive load to fracture; 30 degrees', crosshead_mm_min=None, protocol_locator=f'https://pmc.ncbi.nlm.nih.gov/articles/{pid}/#sec2-materials-17-00365'))
    pid = 'PMC10478297'
    (tid, rows) = table(trees[pid], 'Table 3')
    tables[pid] = rows
    for r in rows[1:]:
        for (j, primer) in [(1, False), (2, True)]:
            records.append(row(pid, tid, r[0] + (' with primer' if primer else ' without primer'), r[j], n=10, material_product='Katana Zirconia HTML (as reported)', material_reported='reported 5 mol%; product-family conflict unresolved', thickness_mm=None, specimen='monolithic_premolar_crown', support='natural human maxillary premolar', cement=r[0], primer=primer, load_angle_deg=0, ageing='10000 thermal cycles 5-55 C; 240000 cycles at 50 N', loading='5 mm steel ball, axial static fracture', crosshead_mm_min=0.5, protocol_locator=f'https://pmc.ncbi.nlm.nih.gov/articles/{pid}/#Sec2'))
    pid = 'PMC10413631'
    (tid, rows) = table(trees[pid], 'Table 3')
    tables[pid] = rows
    for (j, material) in enumerate(['Zircon.x', 'breCAM BioHPP milled', 'BioHPP pressed']):
        records.append(row(pid, tid, rows[1][j], rows[2][j], n=6, material_product=material, material_reported=rows[1][j], thickness_mm=1.5, coping_thickness_mm=0.5, specimen='HIPC_veneered_molar_crown', support='cobalt-chromium die', cement=None, load_angle_deg=0, ageing='none reported', loading='5 mm semispherical head; central fossa', crosshead_mm_min=1.0, protocol_locator=f'https://pmc.ncbi.nlm.nih.gov/articles/{pid}/#Sec2'))
    verified = [(639.3, 111.77), (1378.1, 143.22), (522.67, 108.57), (516.86, 63.32), (635.89, 78.0), (865.3, 116.39), (663.78, 106.8), (980.1, 123.5), (1002.56, 149.58), (696.3, 165.72), (1272.44, 97.78)]
    assert [(r['mean'], r['sd']) for r in records[:11]] == verified
    (tid, rr) = table(trees['PMC10413631'], 'Table 2')
    gaps = []
    for j in range(3):
        (m, s) = mean_sd(rr[2][j])
        gaps.append(dict(study='PMC10413631', doi=STUDIES['PMC10413631'], locator='https://pmc.ncbi.nlm.nih.gov/articles/PMC10413631/#Tab2', group=rr[1][j], n=6, quantity='vertical_marginal_gap_before_cementation', unit='um', mean=m, sd=s, method='40x microscope; 20 locations per crown; sample unit crown', setting='50 um CAD spacer starting 1 mm from margin'))
    out = dict(rows=records, sources=sources, manual_checked_rows=11, xml_expanded_tables=tables, marginal_gap_rows=gaps, requested_gap_review=dict(doi='10.1016/j.prosdent.2025.01.007', status='UNKNOWN', reason='No local verified full-text table found; no numerical row fabricated'))
    dump(ROOT / 'data/published_measurements.json', out)
    return out
