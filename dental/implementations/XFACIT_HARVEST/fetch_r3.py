from dental_release.paths import expand as _release_expand
from pathlib import Path
import json, re, urllib.parse, urllib.request, concurrent.futures, hashlib, time, datetime
P = Path(__file__).resolve().parent
S = Path(_release_expand('@DENTAL_WORK_ROOT@/XFACIT-harvest/sources'))
tasks = {}
for r in json.loads((P / 'SELECTED_R3.json').read_text()):
    l = r['external_referent']['locator']
    for pmc in set(re.findall('PMC\\d+', l)):
        tasks[pmc + '.xml'] = 'https://www.ebi.ac.uk/europepmc/webservices/rest/' + pmc + '/fullTextXML'
    for pmid in set(re.findall('PMID\\s*(\\d+)', l)):
        tasks['PMID' + pmid + '.json'] = 'https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=' + urllib.parse.quote('EXT_ID:' + pmid + ' AND SRC:MED') + '&resultType=core&format=json'
    for doi in set(re.findall('10\\.\\d{4,9}/[A-Za-z0-9._\\-/]+', l)):
        tasks['DOI' + doi.replace('/', '_') + '.json'] = 'https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=' + urllib.parse.quote('DOI:"' + doi + '"') + '&resultType=core&format=json'

def f(kv):
    (k, u) = kv
    st = time.monotonic()
    p = S / k
    try:
        if not p.exists():
            with urllib.request.urlopen(u, timeout=18) as r:
                b = r.read(5000000)
            p.write_bytes(b)
        b = p.read_bytes()
        return {'path': str(p), 'url': u, 'bytes': len(b), 'sha256': hashlib.sha256(b).hexdigest(), 'status': 'OK', 'wall_s': time.monotonic() - st}
    except Exception as e:
        return {'path': str(p), 'url': u, 'status': 'FAILED', 'error': str(e), 'wall_s': time.monotonic() - st}
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:
    rec = list(ex.map(f, tasks.items()))
second = {}
for x in rec:
    p = Path(x['path'])
    if x['status'] == 'OK' and p.name.startswith('DOI'):
        try:
            for r in json.loads(p.read_text())['resultList']['result']:
                if r.get('pmcid'):
                    second[r['pmcid'] + '.xml'] = 'https://www.ebi.ac.uk/europepmc/webservices/rest/' + r['pmcid'] + '/fullTextXML'
        except Exception:
            pass
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:
    rec += list(ex.map(f, second.items()))
(P / 'SOURCE_FETCH_R3_LOG.json').write_text(json.dumps(rec, indent=2))
print('fetched', len(rec), 'ok', sum((x['status'] == 'OK' for x in rec)))
for x in rec:
    if x['status'] != 'OK':
        print(x)
for x in rec:
    p = Path(x['path'])
    if x['status'] == 'OK' and p.name.startswith('DOI'):
        d = json.loads(p.read_text())
        print(p.name, [(r.get('pmcid'), r.get('id'), r.get('title')) for r in d.get('resultList', {}).get('result', [])])
