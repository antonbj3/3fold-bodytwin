from dental_release.paths import expand as _release_expand
from pathlib import Path
import urllib.request, urllib.parse, concurrent.futures, json, re, hashlib, time, datetime, xml.etree.ElementTree as ET
P = Path(__file__).resolve().parent
S = Path(_release_expand('@DENTAL_WORK_ROOT@/XFACIT-harvest/sources'))
S.mkdir(parents=True, exist_ok=True)
rows = json.loads((P / 'SELECTED_40.json').read_text())
tasks = {}
for r in rows:
    loc = r['external_referent'].get('locator', '')
    for pmc in set(re.findall('PMC\\d+', loc)):
        tasks[pmc + '.xml'] = 'https://www.ebi.ac.uk/europepmc/webservices/rest/' + pmc + '/fullTextXML'
    for pmid in set(re.findall('(?:PMID\\s*:?[ /]|EXT_ID:)(\\d+)', loc, re.I)):
        tasks['PMID' + pmid + '.json'] = 'https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=' + urllib.parse.quote('EXT_ID:' + pmid + ' AND SRC:MED') + '&resultType=core&format=json'
for pmid in ['25242248', '27994456', '33601367', '27272434', '34547819', '37367272', '34713996', '36170311', '26904493', '18548932', '20367100', '16775711', '11568312', '11334760', '19885413', '20528965', '29310875', '6576145', '26920510']:
    tasks['PMID' + pmid + '.json'] = 'https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=' + urllib.parse.quote('EXT_ID:' + pmid + ' AND SRC:MED') + '&resultType=core&format=json'
for doi in ['10.3390/ma16051997', '10.3390/ma16062413', '10.1038/s41598-021-90142-5', '10.1055/s-0042-1757910']:
    tasks['DOI' + doi.replace('/', '_') + '.json'] = 'https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=' + urllib.parse.quote('DOI:"' + doi + '"') + '&resultType=core&format=json'
tasks['frontiers.html'] = 'https://www.frontiersin.org/journals/mechanical-engineering/articles/10.3389/fmech.2026.1788600/full'
tasks['hattori.pdf'] = 'https://www.jstage.jst.go.jp/article/sgf1994/2/2/2_2_111/_pdf'
tasks['wear_corrigendum.pdf'] = 'https://iopscience.iop.org/article/10.1088/1402-4896/ad59d7/pdf'

def fetch(k, u):
    st = time.monotonic()
    p = S / k
    try:
        if p.exists():
            b = p.read_bytes()
            cached = True
        else:
            req = urllib.request.Request(u, headers={'User-Agent': 'DentalEvidenceReview/1.0 (public article numeric verification)'})
            with urllib.request.urlopen(req, timeout=22) as f:
                b = f.read(8000000)
            p.write_bytes(b)
            cached = False
        return {'file': str(p), 'url': u, 'sha256': hashlib.sha256(b).hexdigest(), 'bytes': len(b), 'wall_s': time.monotonic() - st, 'cached': cached, 'status': 'OK'}
    except Exception as e:
        return {'file': str(p), 'url': u, 'status': 'FAILED', 'error': str(e), 'wall_s': time.monotonic() - st}
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:
    rec = list(ex.map(lambda kv: fetch(*kv), tasks.items()))
(P / 'SOURCE_FETCH_LOG.json').write_text(json.dumps({'at': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'requests': rec}, indent=2))
print('fetch outcomes', len(rec), sum((x['status'] == 'OK' for x in rec)))
for x in rec:
    if x['status'] != 'OK':
        print(x)
for f in S.glob('PMC*.xml'):
    try:
        root = ET.parse(f).getroot()
    except Exception:
        continue
    out = []
    for el in root.findall('.//article-title')[:1]:
        out.append('TITLE ' + ''.join(el.itertext()))
    for el in root.findall('.//article-id'):
        out.append('ID ' + str(el.attrib) + ' ' + ''.join(el.itertext()))
    for el in root.findall('.//license'):
        out.append('LICENSE ' + ''.join(el.itertext()))
    for el in root.findall('.//abstract'):
        out.append('ABSTRACT ' + ''.join(el.itertext()))
    for el in root.findall('.//table-wrap'):
        out.append('TABLE ' + str(el.attrib) + ' ' + ''.join(el.find('label').itertext()) if el.find('label') is not None else 'TABLE ' + str(el.attrib))
        cap = el.find('caption')
        if cap is not None:
            out.append('CAPTION ' + ''.join(cap.itertext()))
        for row in el.findall('.//tr'):
            out.append(' | '.join((''.join(c.itertext()) for c in row)))
    for el in root.findall('.//body//p'):
        t = ''.join(el.itertext())
        if any((c.isdigit() for c in t)):
            out.append('PARAGRAPH ' + str(el.attrib) + ' ' + t)
    (S / (f.stem + '.extracted.txt')).write_text('\n'.join(out))
for f in S.glob('PMID*.json'):
    try:
        d = json.loads(f.read_text())
        recs = d['resultList']['result']
    except Exception:
        continue
    out = []
    for r in recs:
        out.append(json.dumps({k: r.get(k) for k in ['id', 'pmcid', 'doi', 'title', 'authorString', 'journalInfo', 'abstractText']}, ensure_ascii=False, indent=2))
    (S / (f.stem + '.extracted.txt')).write_text('\n'.join(out))
