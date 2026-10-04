"""Auditable table extraction from LOCAL primary XML only. No downloads."""
from dental_release.paths import expand as _release_expand
from pathlib import Path
import xml.etree.ElementTree as ET
import re, json, csv, hashlib, time
P = Path(__file__).resolve().parent
ROOT = Path(_release_expand('@DENTAL_CORPUS_ROOT@/europepmc/fulltext'))
LIT = Path(_release_expand('@DENTAL_INPUT_ROOT@/workspace/tasks/swarm48/sources/LIT/marginal_fit_papers.jsonl'))

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def text(e):
    return ' '.join(' '.join(e.itertext()).split())

def number(s):
    return float(re.search('[-−]?\\d+(?:\\.\\d+)?', s).group().replace('−', '-'))

def pair(s):
    m = re.search('(\\d+(?:\\.\\d+)?)\\s*±\\s*(\\d+(?:\\.\\d+)?)', s)
    if not m:
        raise ValueError('No mean/SD: ' + s)
    return tuple(map(float, m.groups()))

def grid(tb):
    """Expand row/col spans and retain original numeric cells."""
    occupied = {}
    out = []
    for (i, tr) in enumerate(tb.findall('.//tr')):
        j = 0
        for e in tr:
            while (i, j) in occupied:
                j += 1
            for a in range(int(e.get('rowspan', '1'))):
                for b in range(int(e.get('colspan', '1'))):
                    occupied[i + a, j + b] = text(e)
            j += int(e.get('colspan', '1'))
        out.append([occupied.get((i, c), '') for c in range(max((c for (a, c) in occupied if a == i)) + 1)])
    return out
META = {'PMC10246932': dict(cad='CEREC InLab 5.2.8', machine='MCXL', material='lithium_disilicate', process='milled', restoration='endocrown', method='direct_stereomicroscopy', state='dry', n=10, margin=None), 'PMC10721348': dict(cad='CORiTEC SmartControl', machine='CORiTEC 250i Loader Pro', material='zirconia', process='milled', restoration='crown', method='direct_microscopy', state='dry', n=10, margin=25), 'PMC10333096': dict(cad='', machine='', material='zirconia', process='milled', restoration='crown', method='silicone_replica', state='replica', n=40, margin=40), 'PMC10582242': dict(cad='UNKNOWN', machine='Digident plus', material='resin', process='printed', restoration='crown', method='silicone_replica', state='replica', n=12, margin=None), 'PMC10557992': dict(cad='Exocad 2016', machine='', material='resin', process='printed', restoration='crown', method='silicone_replica', state='replica', n=10, margin=None), 'PMC10828905': dict(cad='3Shape', machine='Ceramill', material='', process='milled', restoration='endocrown', method='silicone_replica', state='replica', n=10, margin=None), 'PMC10971874': dict(cad='3Shape 1.7.1.4', machine='Ceramill', material='ZLS', process='milled', restoration='endocrown', method='silicone_replica', state='replica', n=10, margin=None), 'PMC10912410': dict(cad='3Shape', machine='Meg-printer 2', material='resin', process='printed', restoration='crown', method='silicone_replica', state='replica', n=10, margin=None), 'PMC11010717': dict(cad='UNKNOWN', machine='UNKNOWN', material='zirconia', process='milled', restoration='bridge', method='replica_digital_superposition', state='replica', n=10, margin=40), 'PMC10082505': dict(cad='Exocad DentalCAD 3.0', machine='CEREC MCXL', material='', process='milled', restoration='endocrown', method='stereomicroscopy', state='cemented', n=8, margin=0), 'PMC10096728': dict(cad='Exocad', machine='Imes Icore 350 iPRO+', material='PEEK_PEKK', process='milled', restoration='onlay', method='microCT', state='cemented', n=8, margin=0)}

def extract():
    started = time.perf_counter()
    rows = []
    sources = []
    raw = {}

    def add(pmc, tid, i, col, arm, region, spacer, mean, sd, **extra):
        meta = META[pmc].copy()
        meta.update(extra)
        d = dict(row_id=f'{pmc}:{tid}:{arm}:{region}', study=pmc, source_family='ElGhoul_fit_family' if pmc in ['PMC10828905', 'PMC10971874'] else pmc, arm=arm, region=region, internal_spacer_um=spacer, marginal_spacer_um=meta.pop('margin'), measured_mean_um=mean, reported_sd_um=sd, n_specimens=meta.pop('n'), source_table=tid, source_row=i, source_column=col, source_cell=raw[pmc][tid][i][col], **meta)
        source = next((x for x in sources if x['study'] == pmc))
        d.update(doi=source['doi'], source_path=source['path'], source_sha256=source['sha256'], license=source['license'])
        rows.append(d)
    for pmc in META:
        path = ROOT / (pmc + '.xml')
        r = ET.parse(path).getroot()
        permission = text(r.find('.//permissions'))
        lic = 'CC BY-NC 4.0' if 'by-nc/4.0' in permission else 'CC BY 3.0' if 'by/3.0' in permission else 'CC BY 4.0' if 'by/4.0' in permission else 'CC BY (version unspecified)'
        sources.append(dict(study=pmc, title=text(r.find('.//article-title')), doi=next((x.text for x in r.findall('.//article-id') if x.get('pub-id-type') == 'doi'), None), path=str(path), sha256=sha(path), license=lic))
        raw[pmc] = {t.get('id'): grid(t) for t in r.findall('.//table-wrap')}
    for (i, row) in enumerate(raw['PMC10246932']['TAB2'][1:], 1):
        add('PMC10246932', 'TAB2', i, 1, row[0], 'marginal', number(row[0]), number(row[1]), number(row[2]))
    g = raw['PMC10721348']['tab1']
    spacer = None
    for (i, row) in enumerate(g):
        if row[0].startswith('Group'):
            spacer = number(row[0])
        if row[0].strip().lower() == 'mean':
            sd = next((number(z[1]) for z in g[i + 1:] if 'standard deviation' in z[0].lower() or z[0].strip() == 'SD'))
            add('PMC10721348', 'tab1', i, 1, f'G{int(spacer)}', 'marginal', spacer, number(row[1]), sd)
    for (tid, region) in [('T5', 'marginal'), ('T9', 'axial'), ('T13', 'occlusal')]:
        g = raw['PMC10333096'][tid]
        for (i, row) in enumerate(g[2:], 2):
            for (j, machine) in enumerate(['Aegis HM', 'Trione Z', 'Motion 2'], 2):
                (mu, sd) = pair(row[j])
                cad = row[1]
                add('PMC10333096', tid, i, j, cad + '/' + machine, region, 40, mu, sd, cad=cad, machine=machine)
    pmc = 'PMC10582242'
    tid = 'cre2758-tbl-0001'
    region = None
    for (i, row) in enumerate(raw[pmc][tid]):
        if row[0] in ['Marginal', 'Axial', 'Axio‐occlusal', 'Occlusal']:
            region = {'Axio‐occlusal': 'transition'}.get(row[0], row[0].lower())
        if re.fullmatch('25|50|100', row[0]):
            (mu, sd) = pair(row[1])
            add(pmc, tid, i, 1, 'layer_' + row[0], region, 80, mu, sd, layer_um=int(row[0]))
    pmc = 'PMC10557992'
    tid = 'T1'
    region = None
    for (i, row) in enumerate(raw[pmc][tid]):
        if row[0]:
            region = {'Gap marginal': 'marginal', 'Axial gap': 'axial', 'Occluso-axial gap': 'transition', 'Occlusal gap': 'occlusal'}.get(row[0], region)
        if len(row) > 2 and row[1] in ['Asiga', 'Digident']:
            (mu, sd) = pair(row[2])
            add(pmc, tid, i, 2, row[1], region, 50, mu, sd, machine=row[1], layer_um=50)
    pmc = 'PMC10828905'
    tid = 'cre2843-tbl-0002'
    for (i, row) in enumerate(raw[pmc][tid][1:], 1):
        region = {'Marginal': 'marginal', 'Axial': 'axial', 'Cervical': 'cervical', 'Pulpl': 'pulpal', 'Internal': 'internal_overall'}[row[0]]
        add(pmc, tid, i, 3, row[1] + '/' + row[2], region, 60, number(row[3]), number(row[4]), material='lithium_disilicate' if 'max' in row[2] else 'zirconia', n=int(row[5]))
    pmc = 'PMC10971874'
    tid = 'materials-17-01411-t002'
    for (i, row) in enumerate(raw[pmc][tid][2:], 2):
        for (j, region) in [(2, 'marginal'), (4, 'cervical'), (6, 'axial'), (8, 'pulpal'), (10, 'internal_overall')]:
            add(pmc, tid, i, j, row[0] + '/' + row[1], region, 60, number(row[j]), number(row[j + 1]))
    pmc = 'PMC10912410'
    tid = 'Tab2'
    tooth = None
    region = None
    for (i, row) in enumerate(raw[pmc][tid]):
        if row[0] in ['Posterior teeth', 'Anterior teeth']:
            tooth = row[0]
        if len(row) > 2 and row[2] == 'Optical microscope':
            region = {'MG': 'marginal', 'Axial': 'axial', 'Occlusal': 'occlusal', 'Incisal': 'incisal', 'AMD': 'absolute_marginal', 'Chamfer': 'cervical', 'Angle': 'transition'}[row[1]]
            add(pmc, tid, i, 3, tooth, region, 60, number(row[3]), number(row[4]))
    pmc = 'PMC11010717'
    tid = 'tbl4'
    for (i, row) in enumerate(raw[pmc][tid][3:], 3):
        region = {'Marginal gaps': 'marginal', 'Internal gaps at occlusal surface': 'occlusal', 'Internal gaps at axial wall': 'axial'}[row[0]]
        for (j, arm) in enumerate(['CI-LI-X-X-P2', 'CI-X-C-X-P2', 'CI-X-X-P1-P2'], 1):
            (mu, sd) = pair(row[j].split(',')[1])
            add(pmc, tid, i, j, arm, region, 50, mu, sd, sd_basis='measurement_locations', n_locations=int(row[j].split('/')[0]))
    pmc = 'PMC10082505'
    tid = 'Tab1'
    for (i, row) in enumerate(raw[pmc][tid][2:], 2):
        region = 'internal_overall' if row[0] == 'Internal fit' else 'marginal'
        (mu, sd) = pair(row[2])
        add(pmc, tid, i, 2, row[1], region, 50, mu, sd, material='PEEK_PEKK' if row[1].startswith('P') else 'lithium_disilicate')
    pmc = 'PMC10096728'
    tid = 'polymers-15-01715-t003'
    for (i, row) in enumerate(raw[pmc][tid][1:], 1):
        for (j, arm) in enumerate(['Pekkton', 'Juvora', 'BioHPP'], 1):
            (mu, sd) = pair(row[j])
            add(pmc, tid, i, j, arm, 'marginal' if row[0] == 'MG' else 'internal_overall', 30, mu, sd)
    assert len(set((d['row_id'] for d in rows))) == len(rows)
    keys = list(dict.fromkeys((k for d in rows for k in d)))
    with (P / 'measurements.csv').open('w') as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)
    (P / 'measurements.json').write_text(json.dumps(rows, indent=2, ensure_ascii=False) + '\n')
    (P / 'RAW_TABLES.json').write_text(json.dumps(raw, indent=2, ensure_ascii=False) + '\n')
    lit = {d['pmcid']: d for d in map(json.loads, LIT.read_text().splitlines())}
    for s in sources:
        s['lit_package_same_doi'] = lit.get(s['study'], {}).get('doi') == s['doi']
    report = dict(rows=len(rows), studies=len(sources), primary_rows=sum((d['region'] in ['marginal', 'axial', 'occlusal'] for d in rows)), source_families=len(set((d['source_family'] for d in rows))), sources=sources, measurements_sha256=sha(P / 'measurements.json'), csv_sha256=sha(P / 'measurements.csv'), lit_path=str(LIT), lit_sha256=sha(LIT), seconds=time.perf_counter() - started, exclusions=[{'source': 'PMC12237415', 'reason': 'Referenced in R04 but absent from prescribed local fulltext and LIT package; not downloaded, not treated as admitted numeric data'}, {'source': 'PMC10154143:T3', 'reason': 'Signed CAD surface deviations, not seated measured gap; negative values cannot be interpreted as film'}, {'source': 'PMC10318718', 'reason': 'Only direct marginal table; regional internal values shown in figures. Not digitized or guessed.'}], conflicts=[{'source': 'PMC10721348:tab1:G50', 'table_mean_um': 23.22, 'prose_mean_um': 22.22, 'resolution': 'table preserved'}, {'source': 'PMC10721348:tab1', 'reason': 'Reported SEM appears to use measurement locations; n_specimens remains 10. No inferred precision weights.'}, {'source': ['PMC10828905', 'PMC10971874'], 'reason': 'Overlapping authors/method and repeated values/SD; conservatively joint holdout. Common samples not established.'}, {'source': 'PMC10582242:occlusal:layer_100', 'reason': 'Published minimum/maximum and SD appear inconsistent; retain published mean/SD, flag source; no correction or specimen-level claim.'}])
    manifest = P / 'SOURCE_MANIFEST.json'
    if manifest.exists():
        old = json.loads(manifest.read_text())
        assert old['measurements_sha256'] == report['measurements_sha256'], 'Frozen extracted data changed'
        assert old['sources'] == report['sources'], 'Frozen source metadata changed'
        report = old
    else:
        manifest.write_text(json.dumps(report, indent=2, ensure_ascii=False) + '\n')
    print(json.dumps({k: report[k] for k in ['rows', 'studies', 'primary_rows', 'source_families', 'seconds']}))
    return (rows, report)
if __name__ == '__main__':
    extract()
