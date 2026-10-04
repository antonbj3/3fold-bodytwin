"""Small source-specific extraction. Preserves independent article/table locators."""
from dental_release.paths import expand as _release_expand
from pathlib import Path
import csv, hashlib, json, re, subprocess, xml.etree.ElementTree as ET
ROOT = Path(__file__).resolve().parents[1]
CORPUS = Path(_release_expand('@DENTAL_CORPUS_ROOT@/europepmc/fulltext'))
X13 = Path(_release_expand('@DENTAL_IMPLEMENTATIONS@/X13'))

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def save(name, value):
    (ROOT / name).write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + '\n')

def text(e):
    return ' '.join(' '.join(e.itertext()).split())

def grid(tb):
    (occupied, out) = ({}, [])
    for (i, tr) in enumerate(tb.findall('.//tr')):
        j = 0
        for e in tr:
            while (i, j) in occupied:
                j += 1
            (rs, cs) = (int(e.get('rowspan', '1')), int(e.get('colspan', '1')))
            for a in range(rs):
                for b in range(cs):
                    occupied[i + a, j + b] = text(e)
            j += cs
        out.append([occupied.get((i, c), '') for c in range(1 + max((c for (a, c) in occupied if a == i)))])
    return out

def extract_cunali():
    pdf = ROOT / 'raw/cunali2017.pdf'
    result = subprocess.run(['pdftotext', '-layout', str(pdf), '-'], check=True, capture_output=True, text=True)
    t = result.stdout
    lines = [l for l in t.splitlines() if re.search('(Replica|Micro-CT)\\s+\\d+\\.\\d+±', l)]
    if len(lines) != 4:
        raise ValueError('Expected four source method rows in Table1')
    vals = [re.findall('(\\d+\\.\\d+)±(\\d+\\.\\d+)', l) for l in lines]
    if any((len(v) != 4 for v in vals)):
        raise ValueError('Expected four regions per method row')
    regions = ['marginal', 'axial', 'transition', 'occlusal']
    rows = []
    for (k, system) in enumerate(['Amann', 'Dentsply_Sirona']):
        for (j, region) in enumerate(regions):
            (rm, rs) = vals[2 * k][j]
            (cm, cs) = vals[2 * k + 1][j]
            rows.append(dict(source='Cunali2017', system=system, region=region, replica_mean_um=rm, replica_sd_um=rs, ct_mean_um=cm, ct_sd_um=cs, n_specimens=10, material='zirconia', restoration='coping', physical_state='PVS_finger_pressure_5min_vs_clean_dry_reseating', pair_level='SAME_SPECIMEN_COHORT_REGION_SUMMARY_ONLY', covariance=None, resolution='PER_SURFACE_REGION', statistical_scope='mean and sample SD across ten copings; four points averaged per coping', source_locator='https://doi.org/10.1590/0103-6440201601531#Table1', source_page=470, source_row_replica=2 * k + 1, source_row_ct=2 * k + 2, source_column=j + 1, source_path=str(pdf), source_sha256=sha(pdf)))
    save('raw/CUNALI_REGION_SUMMARIES.json', rows)
    return rows

def extract_pasha():
    path = CORPUS / 'PMC10323296.xml'
    root = ET.parse(path).getroot()
    tables = {tw.get('id'): grid(tw) for tw in root.findall('.//table-wrap')}
    save('raw/PASHA_TABLES.json', tables)
    out = []
    for (tid, method) in [('TAB1', 'replica'), ('TAB2', 'ct')]:
        region = None
        for (i, row) in enumerate(tables[tid]):
            if row[0] in ['Occlusal cavosurface area', 'Proximal', 'Axial']:
                region = row[0]
            if not region:
                continue
            if len(row) < 5 or row[1].lower() not in ['cad-cam', '3d printing']:
                continue
            out.append(dict(source='PMC10323296', system='CAD-CAM' if row[1].lower() == 'cad-cam' else '3D printing', region=region, method=method, n_specimens=int(row[2]), mean_um=row[3], sd_um=row[4], source_table=tid, source_row=i, mean_column=3, sd_column=4, source_locator='https://doi.org/10.7759/cureus.40020#' + tid, source_path=str(path), source_sha256=sha(path), resolution='PER_SURFACE_REGION'))
    if len(out) != 10:
        raise ValueError(f'Expected ten present Pasha table rows, got {len(out)}')
    save('raw/PASHA_EXCLUSIONS.json', [dict(system='3D printing', region='Axial', method=m, reason='Axial 3D-printing row absent from local XML table; prose-only numbers not promoted to table data', table=t, source_path=str(path), source_sha256=sha(path)) for (m, t) in [('replica', 'TAB1'), ('ct', 'TAB2')]])
    pairs = []
    for r in [r for r in out if r['method'] == 'replica']:
        c = next((c for c in out if c['method'] == 'ct' and c['region'] == r['region'] and (c['system'] == r['system'])))
        pairs.append(dict(source='PMC10323296', system=r['system'], region=r['region'], replica_mean_um=r['mean_um'], replica_sd_um=r['sd_um'], ct_mean_um=c['mean_um'], ct_sd_um=c['sd_um'], n_specimens=r['n_specimens'], covariance=None, pair_level='SAME_SAMPLE_COHORT_SUMMARY_ONLY', material='zirconia' if r['system'].lower() == 'cad-cam' else 'resin', restoration='onlay', physical_state='replica and CT state linkage not sufficiently described', resolution='PER_SURFACE_REGION', source_locator=r['source_locator'] + ' + TAB2', source_path=str(path), source_sha256=sha(path), cells=[r, c]))
    save('raw/PASHA_REGION_SUMMARIES.json', pairs)
    return pairs

def extract_smartphone():
    path = CORPUS / 'PMC10912410.xml'
    root = ET.parse(path).getroot()
    tb = next((x for x in root.findall('.//table-wrap') if x.get('id') == 'Tab2'))
    g = grid(tb)
    save('raw/SMARTPHONE_TABLE2.json', g)
    out = []
    for (i, row) in enumerate(g):
        if len(row) > 4 and row[2] == 'Optical microscope':
            c = g[i + 1]
            if c[2] != 'Smartphone camera':
                raise ValueError('Device row pair changed')
            out.append(dict(source='PMC10912410', tooth_type=row[0], region=row[1], microscope_mean_um=row[3], microscope_sd_um=row[4], phone_mean_um=c[3], phone_sd_um=c[4], covariance=None, reported_ICC=row[8], physical_state='SAME_SECTIONED_SILICONE_REPLICA', pair_level='SAME_FILM_GROUP_SUMMARY_ONLY', resolution='PER_SURFACE_REGION', source_locator='https://doi.org/10.1038/s41598-024-55711-4#Tab2', source_row_microscope=i, source_row_phone=i + 1, source_path=str(path), source_sha256=sha(path)))
    if len(out) != 11:
        raise ValueError(f'Expected eleven device-region pairs, got {len(out)}')
    save('raw/SMARTPHONE_REGION_SUMMARIES.json', out)
    return out

def x13_rows():
    rows = json.loads((X13 / 'measurements.json').read_text())
    if len(rows) != 148:
        raise ValueError('X13 input count changed; require new prereg')
    for row in rows:
        if row['source_sha256'] != sha(row['source_path']):
            raise ValueError('X13 source drift: ' + row['row_id'])
    save('raw/X13_INPUT_SNAPSHOT.json', rows)
    return rows

def source_manifest():
    sources = [ROOT / 'raw/cunali2017.pdf', CORPUS / 'PMC10323296.xml', CORPUS / 'PMC10912410.xml', X13 / 'measurements.json']
    save('SOURCE_MANIFEST.json', [dict(path=str(p), bytes=p.stat().st_size, sha256=sha(p)) for p in sources])
if __name__ == '__main__':
    a = extract_cunali()
    b = extract_pasha()
    c = extract_smartphone()
    x = x13_rows()
    source_manifest()
    print(json.dumps(dict(cunali_region_pairs=len(a), pasha_region_pairs=len(b), device_region_pairs=len(c), x13_rows=len(x))))
