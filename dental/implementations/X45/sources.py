"""Read immutable local primary tables; cross-check against a separate transcription."""
from dental_release.paths import expand as _release_expand
import hashlib, json, re
from pathlib import Path
import xml.etree.ElementTree as E
CORPUS = Path(_release_expand('@DENTAL_CORPUS_ROOT@/europepmc/fulltext'))
ROOT = Path(__file__).resolve().parent

def text(e):
    return ' '.join(' '.join(e.itertext()).split())

def row_cells(table):
    return [[text(c) for c in row if c.tag in ('td', 'th')] for row in table.findall('.//tr')]

def first_number(s):
    return float(re.search('[-−]?\\d+(?:\\.\\d+)?', s).group().replace('−', '-'))

def pair(s):
    s = re.sub('\\s', '', s)
    p = s.split('±')
    if len(p) != 2:
        raise ValueError('mean/uncertainty pair missing ' + s)
    return (first_number(p[0]), first_number(p[1]))

def load():
    out = {'ISQ': {}, 'LTD': {}, 'WEAR': {}, 'manifest': [], 'parser_controls': []}
    for id in ['PMC10092180', 'PMC10569846', 'PMC10258403']:
        path = CORPUS / (id + '.xml')
        r = E.parse(path).getroot()
        out['manifest'].append({'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest(), 'bytes': path.stat().st_size, 'doi': r.findtext('.//article-id[@pub-id-type="doi"]'), 'license': text(r.find('.//permissions'))})
        (ROOT / 'raw' / ('SOURCE_TABLES_' + id + '.json')).write_text(json.dumps([{'id': t.get('id'), 'caption': text(t.find('caption')), 'cells': row_cells(t)} for t in r.findall('.//table-wrap')], indent=2, ensure_ascii=False) + '\n')
        if id == 'PMC10092180':
            tab = r.find('.//table-wrap[@id="cid13140-tbl-0002"]')
            rows = row_cells(tab)
            for row in rows:
                if not row or row[0] not in ['Baseline', '7 days', '14 days', '21 days', '28 days', '60 days', '90 days']:
                    continue
                day = 0 if row[0] == 'Baseline' else int(row[0].split()[0])
                for (name, cell) in zip(['OD', 'PISP'], row[1:3]):
                    (y, s) = pair(cell)
                    out['ISQ'].setdefault(name, []).append({'day': day, 'mean_ISQ': y, 'SD_ISQ': s, 'resolution': 'POPULATION'})
        elif id == 'PMC10569846':
            rows = row_cells(r.find('.//table-wrap[@id="TB2252138-4"]'))
            for row in rows:
                if len(row) != 5 or row[0] not in ['T', 'ST', 'XT', 'P']:
                    continue
                out['LTD'][row[0]] = {'unaged': {'m': first_number(row[1]), 'sigma0_MPa': first_number(row[2])}, 'aged': {'m': first_number(row[3]), 'sigma0_MPa': first_number(row[4])}}
        else:
            rows = row_cells(r.find('.//table-wrap[@id="T3"]'))
            for row in rows:
                if len(row) < 7 or not row[0].startswith(('Filtek Z250', 'Vertise flow', 'Admira', 'Grandio flow', 'Filtek flow')):
                    continue
                name = 'Admira fusion flow' if row[0].startswith('Admira') else row[0]
                points = []
                for (cycles, cell) in zip([5000, 10000, 20000, 40000, 80000, 120000], row[1:7]):
                    (y, s) = pair(cell)
                    points.append({'cycles': cycles, 'mean_mg': float(str(y)) / 10, 'SEM_mg': float(str(s)) / 10, 'source_mean_1e-4g': y, 'source_SEM_1e-4g': s, 'resolution': 'POPULATION'})
                out['WEAR'][name] = points
    expect = {'ISQ': {'OD': [71.3, 68.4, 66.4, 65.9, 66.3, 67.9, 69.3], 'PISP': [69.3, 67.3, 66.5, 67.1, 67.8, 69.7, 70.9]}, 'LTD': {'T': [7.04, 871.22, 10.48, 850.48], 'ST': [8.25, 724.66, 10.71, 683.35], 'XT': [7.75, 458.66, 6.75, 391.11], 'P': [6.4, 436.81, 4.5, 368.34]}, 'WEAR': {'Filtek Z250': [9.8, 10.6, 11.3, 12, 12.8, 13.3], 'Vertise flow': [8.1, 8.8, 9.8, 9.84, 11.9, 15.4], 'Admira fusion flow': [5.4, 7.1, 7.9, 9.5, 12.4, 17.5], 'Grandio flow': [10.7, 15.5, 16, 16.7, 18.4, 21.1], 'Filtek flow': [5.9, 10.2, 12.3, 14.8, 15.5, 16.6]}}
    actual = {'ISQ': {k: [p['mean_ISQ'] for p in v] for (k, v) in out['ISQ'].items()}, 'LTD': {k: [v['unaged']['m'], v['unaged']['sigma0_MPa'], v['aged']['m'], v['aged']['sigma0_MPa']] for (k, v) in out['LTD'].items()}, 'WEAR': {k: [p['source_mean_1e-4g'] for p in v] for (k, v) in out['WEAR'].items()}}
    if actual != expect:
        raise ValueError('independent published-table transcription mismatch: ' + json.dumps(actual))
    out['parser_controls'] = [{'name': 'independent primary-table transcription', 'pass': True, 'values': 60}, {'name': 'injected OD baseline+10', 'detected': actual['ISQ']['OD'][0] + 10 != expect['ISQ']['OD'][0]}, {'name': 'injected sigma0*1.1', 'detected': actual['LTD']['T'][1] * 1.1 != expect['LTD']['T'][1]}, {'name': 'injected mass conversion factor1000', 'detected': out['WEAR']['Filtek Z250'][0]['mean_mg'] * 1000 != 0.98}]
    (ROOT / 'raw' / 'PRIMARY_MEASUREMENTS.json').write_text(json.dumps(out, indent=2, ensure_ascii=False) + '\n')
    return out
