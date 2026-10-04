"""Preserve primary cells, correct article identity, and separate endpoint kinds."""
import csv, hashlib, json, re, time, xml.etree.ElementTree as ET
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]

def tx(e):
    return re.sub('\\s+', ' ', ''.join(e.itertext())).strip()

def cell_rows(table):
    return [[tx(c) for c in list(tr) if c.tag in ['td', 'th']] for tr in table.findall('.//tr')]

def num(s):
    return float(re.search('[-+]?\\d+(?:\\.\\d+)?', s.replace(',', '')).group())

def main():
    t = time.perf_counter()
    sources = {}
    tables = {}
    for p in sorted((ROOT / 'inputs/primary').glob('*.xml')):
        r = ET.parse(p).getroot()
        doi = next((tx(e) for e in r.iter('article-id') if e.get('pub-id-type') == 'doi'), None)
        licence = next((tx(e) for e in r.iter() if 'license' in e.tag), 'UNKNOWN: local XML has no license statement')
        ss = {'pmcid': p.stem, 'doi': doi, 'locator': f'https://doi.org/{doi}' if doi else f'https://europepmc.org/articles/{p.stem}', 'local_path': str(p.relative_to(ROOT)), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest(), 'title': next((tx(e) for e in r.iter('article-title')), ''), 'license_statement': licence, 'tables': []}
        for table in r.iter('table-wrap'):
            tt = {'id': table.get('id'), 'label': tx(table.find('label')) if table.find('label') is not None else None, 'cells': cell_rows(table)}
            tables[p.stem, tt['id']] = tt
            ss['tables'].append(tt)
        ss['paragraphs'] = [{'paragraph_index': i, 'text': tx(e)} for (i, e) in enumerate(r.iter('p')) if re.search('fatigue limit|fatigue load limit|run.?out|residual stress|F = 200|F.*FL', tx(e), re.I)]
        sources[p.stem] = ss
    (ROOT / 'raw/primary_tables.json').write_text(json.dumps(sources, indent=2, ensure_ascii=False) + '\n')
    specs = [('K1-A3a', 'PMC7730231', 'ijerph-17-08988-t002', 1, 5, 190), ('K1-A3b', 'PMC7730231', 'ijerph-17-08988-t002', 2, 5, 150), ('K1-A3c', 'PMC7730231', 'ijerph-17-08988-t002', 3, 5, 325), ('N01', 'PMC6803993', 'materials-12-03264-t003', 5, 1, 285), ('N02', 'PMC6803993', 'materials-12-03264-t003', 11, 1, 316), ('N03', 'PMC9680401', 'jfb-13-00239-t003', 4, 1, 405), ('N04', 'PMC9680401', 'jfb-13-00239-t003', 13, 1, 397), ('N05', 'PMC9680401', 'jfb-13-00239-t003', 20, 1, 431), ('N10', 'PMC10971635', 'jcm-13-01632-t001', 4, 1, 118), ('N11', 'PMC10971635', 'jcm-13-01632-t001', 8, 1, 230)]
    audit = []
    direct = {}
    for (aid, pmc, table, ri, ci, expected) in specs:
        cells = tables[pmc, table]['cells'][ri]
        raw = cells[ci]
        value = num(raw.split('–')[-1])
        loc = f'{pmc} table {table}, zero-based tr {ri}, cell {ci}'
        audit.append({'id': aid, 'locator': loc, 'cell_text': raw, 'inspected_expected_N': expected, 'extracted_N': value, 'equal': value == expected, 'method': 'direct primary cell extraction checked against source-reading transcription'})
        direct[aid] = (value, loc)
    if not all((a['equal'] for a in audit)):
        raise ValueError('primary cell audit fails')
    records = []
    corrections = []
    old = json.loads((ROOT / 'inputs/k1b_eval_legacy.json').read_text())['per_anchor']
    for r in csv.DictReader((ROOT / 'inputs/anchors_legacy.csv').open()):
        s = sources[r['pmcid']]
        doi = s['doi']
        orig = r['doi']
        if doi != orig:
            corrections.append({'id': r['id'], 'pmcid': r['pmcid'], 'old_doi': orig, 'primary_doi': doi})
        reasons = []
        strict = []
        if r['role'].startswith('excluded'):
            reasons.append('implantoplasty arm not in legacy forward domain')
        if r['id'] == 'K1-A1':
            reasons.append('inconsistent percent/static force and nonstandard embedding')
        if r['id'] in ['N06', 'N07', 'N08', 'N19', 'N20']:
            reasons.append('asymptotic fitted force or pooled/unknown cycle endpoint')
        if r['id'] == 'N16':
            reasons.append('porous AM geometry and diameter unknown')
        (val, loc) = direct.get(r['id'], (float(r['F_lim_N']), r['source_loc']))
        if r['id'] not in direct:
            strict.append('not independently cell-extracted; legacy transcription with primary paragraph/figure locator')
        try:
            horizon = float(r['N_runout'])
        except ValueError:
            horizon = None
        if horizon != 5000000.0:
            strict.append('cycle horizon not explicit 5e6')
        if float(r['angle_deg']) != 30:
            strict.append('angle not 30 degrees')
        if r['R'] != '0.1':
            strict.append('R not explicit 0.1')
        try:
            lever = float(r['lever_l_mm'])
            lever_status = 'reported_or_derived_from_published_moment'
        except ValueError:
            lever = None
            lever_status = 'UNKNOWN'
            strict.append('lever not reported; default forbidden for strict cohort')
        if r['id'] in ['N13', 'N14', 'N15', 'N12']:
            strict.append('failure mode not modeled implant body')
        strict += reasons
        records.append({'id': r['id'], 'pmcid': r['pmcid'], 'doi': doi, 'locator': s['locator'], 'table_locator': loc, 'source_loc_legacy': r['source_loc'], 'force_N': val, 'force_kind': 'maximum_cyclic_force', 'force_sd_N': None, 'diameter_mm': float(r['D_mm']) if r['D_mm'] != 'UNKNOWN' else None, 'connection': r['connection'] if 'UNKNOWN' not in r['connection'] else None, 'material': r['material'], 'manufacturing': r['manufacturing'], 'angle_deg': float(r['angle_deg']), 'R': 0.1 if r['R'] == '0.1' else None, 'horizon_cycles': horizon, 'frequency_Hz': float(r['freq_Hz']), 'medium': '0.9% NaCl' if r['id'] == 'N09' else None, 'lever_mm': lever, 'lever_status': lever_status, 'criterion': r['runout_criterion'], 'failure_location': r['failure_location'], 'legacy_prediction_N': old.get(r['id'], {}).get('M_asm_P10_P50_P90', [None, None, None])[1], 'material_class': old.get(r['id'], {}).get('mat'), 'strict_eligible': len(strict) == 0, 'strict_reasons': strict, 'endpoint_eligible': len(reasons) == 0, 'endpoint_rejection_reasons': reasons, 'resolution': 'PER_TOOTH', 'source_sha256': s['sha256']})
    data = {'records': records, 'manual_cell_audit': audit, 'identity_corrections': corrections, 'screened': len(records), 'strict_kept': sum((r['strict_eligible'] for r in records)), 'strict_rejected': sum((not r['strict_eligible'] for r in records)), 'strict_rejection_fraction': sum((not r['strict_eligible'] for r in records)) / len(records), 'wall_s': time.perf_counter() - t, 'primary_articles': len(sources), 'claim_type': 'information_link'}
    (ROOT / 'raw/published_endpoints.json').write_text(json.dumps(data, indent=2, ensure_ascii=False) + '\n')
    print(json.dumps({k: v for (k, v) in data.items() if k not in ['records', 'manual_cell_audit', 'identity_corrections']}))
if __name__ == '__main__':
    main()
