"""Re-extract primary cells; no fabricated individual specimens or geometry."""
import csv, re, shutil, subprocess, xml.etree.ElementTree as ET
from pathlib import Path
from common import R, CORPUS, dump, sha, now

def text(e):
    return ' '.join(' '.join(e.itertext()).split()) if e is not None else ''

def extract():
    proc = subprocess.run(['rg', '-l', '-i', 'connector', str(CORPUS), '-g', '*.xml'], capture_output=True, text=True, check=True)
    paths = sorted(proc.stdout.splitlines())
    candidates = []
    screening = []
    for p in paths:
        root = ET.parse(p).getroot()
        title = text(root.find('.//article-title'))
        body = text(root).lower()
        relevant = 'fracture' in body and ('zirconia' in body or 'disilicate' in body) and ('prosthes' in body or 'bridge' in body)
        eligible = Path(p).stem == 'PMC10471503'
        reason = 'RETAINED_SYSTEM_EVENT_NONCONNECTOR_ORIGIN' if eligible else 'NOT_TARGET_DENTAL_FRACTURE_QUERY' if not relevant else 'CANTILEVER_LOAD_PATH' if 'cantilever' in title.lower() else 'REVIEW_NO_PRIMARY_TARGET_MEASUREMENT' if any((x in title.lower() for x in ['review', 'overview'])) else 'SIMULATION_NO_EXTERNAL_BREAK_LOAD' if any((x in title.lower() for x in ['finite element', 'remodelling of connectors'])) else 'OTHER_MATERIAL_RESIN' if any((x in title.lower() for x in ['interim', 'provisional resin', 'polylactic'])) else 'WRONG_PROSTHESIS_OR_ENDPOINT_OR_NO_CONNECTOR_GEOMETRY'
        screening.append(dict(pmc=Path(p).stem, title=title, path=p, relevant=relevant, retained_system=eligible, retained_uncensored_connector=False, reason=reason))
        if relevant:
            candidates.append(p)
    additional = read_additional()
    known = {s['pmc'] for s in screening}
    for p in additional:
        if Path(p).stem in known:
            continue
        root = ET.parse(p).getroot()
        screening.append(dict(pmc=Path(p).stem, title=text(root.find('.//article-title')), path=p, relevant=True, retained_system=False, retained_uncensored_connector=False, reason='MATERIAL_COUPON_OR_REVIEW_NOT_3UNIT_BRIDGE'))
    src = CORPUS / 'PMC10471503.xml'
    target = R / 'raw/sources/PMC10471503.xml'
    target.parent.mkdir(exist_ok=True, parents=True)
    if target.exists() and sha(target) != sha(src):
        raise ValueError('Source snapshot drift')
    shutil.copyfile(src, target)
    root = ET.parse(target).getroot()
    table = root.find('.//table-wrap[@id="T2"]')
    rows = table.findall('.//tr')[2:]
    doi = text(root.find('.//article-id[@pub-id-type="doi"]'))
    data = []
    for tr in rows:
        cells = list(tr)
        product = text(cells[0])
        material = re.search('([345]Y)-', product)[1]
        for (j, (h, b)) in enumerate([(4.0, 2.25), (3.0, 3.0)], 1):
            raw = text(cells[j])
            m = re.search('([\\d.]+)\\s*±\\s*([\\d.]+)', raw)
            data.append(dict(id=f'PMC10471503_{material}_{j}', study='PMC10471503', doi=doi, material=material, product=product, n=6, mean_N=float(m[1]), sd_N=float(m[2]), area_mm2=9.0, height_mm=h, width_mm=b, geometry_order='height x width, methods explicitly decreasing height 4 to 3 mm', span_mm=None, gingival_radius_mm=None, setup='end-supported posterior 3-unit; CoCr dies; U200 cement; 5 mm steel ball at pontic center; axial; 1 mm/min; unaged', endpoint='catastrophic_system_fracture', origin='pontic_occlusal_contact', connector_observation='right_censored_at_system_break; NOT uncensored connector failure', mean_resolution='POPULATION', specimen_resolution='PER_TOOTH', geometry_resolution='PER_SURFACE_REGION', force_unit='N', dimension_unit='mm', area_unit='mm2', time_scale='SIMULTANEOUS', table_locator=f'/article/body//table-wrap[@id="T2"]/table//tr[{rows.index(tr) + 3}]/td[{j + 1}]', n_locator='Materials and methods: six groups, total 36; abstract n=6', geometry_locator='T1; connector configurations; Introduction height reduction', origin_locator='RESULT; DISCUSSION; Figures5-7: origin occlusal pontic, no connector fracture origin', raw_cell=raw, source_sha256=sha(target), source_path=str(src)))
    dump('raw/MEASUREMENTS.json', data)
    dump('raw/SCREENING.json', screening)
    keys = list(data[0])
    with (R / 'raw/MEASUREMENTS.csv').open('w') as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(data)
    dump('SOURCE_MANIFEST.json', dict(created_utc=now(), primary_sources=[dict(path=str(src), snapshot='raw/sources/PMC10471503.xml', sha256=sha(src), bytes=src.stat().st_size, license='CC BY-NC 4.0', doi=doi)], brief_missing_local_sources=['PMC7867680', 'PMC12697478'], search=dict(connector_fulltexts=len(paths), query_relevant=sum((s['relevant'] for s in screening)), additional_unique_sources=len(screening) - len(paths), retained_system_studies=1, uncensored_connector_studies=0), no_network_dataset_downloads=True))
    return (data, screening)

def read_additional():
    import json
    p = R / 'raw/additional_candidates.json'
    return [x['path'] for x in json.loads(p.read_text())] if p.exists() else []
if __name__ == '__main__':
    (d, s) = extract()
    print(f'{len(d)} system groups; {len(s)} screened sources; 0 uncensored connector studies')
