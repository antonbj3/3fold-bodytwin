from dental_release.paths import expand as _release_expand
import concurrent.futures, datetime, hashlib, json, re, time, urllib.request, urllib.parse
from pathlib import Path
ROOT = Path(__file__).resolve().parent
DENTAL = ROOT.parent.parent
OUT = ROOT / 'sources'
OUT.mkdir(exist_ok=True)
LOCAL = Path(_release_expand('@DENTAL_CORPUS_ROOT@/europepmc/fulltext'))
API = 'https://www.ebi.ac.uk/europepmc/webservices/rest/'
queries = set()
inventory = ROOT / ('INPUT_INVENTORY_R4.json' if (ROOT / 'INPUT_INVENTORY_R4.json').exists() else 'INPUT_INVENTORY.json')
for row in json.loads(inventory.read_text()):
    if row['exists']:
        s = (ROOT / row['frozen_report_snapshot']).read_text() if row.get('frozen_report_snapshot') else Path(row['path']).read_text()
        queries.update(('EXT_ID:' + p + ' AND SRC:MED' for p in re.findall('PMID[:\\s*·]*(\\d+)', s)))
        queries.update(('DOI:' + d.rstrip('.,;') for d in re.findall('10\\.\\d{4,9}/[A-Za-z0-9._;()/:+-]+', s)))
queries.update(('DOI:' + d for d in ['10.1080/03036758.2019.1691612', '10.1111/joor.70173', '10.1038/s41598-017-05788-x', '10.1186/s12903-022-02484-9', '10.1016/j.jbiomech.2009.03.040', '10.7144/sgf.2.111']))
if (ROOT / 'EXTRA_QUERIES.json').exists():
    queries.update(json.loads((ROOT / 'EXTRA_QUERIES.json').read_text()))

def request(url):
    req = urllib.request.Request(url, headers={'User-Agent': 'Dental-independent-source-audit/1.0'})
    with urllib.request.urlopen(req, timeout=25) as r:
        return r.read()

def work(q):
    key = hashlib.sha256(q.encode()).hexdigest()[:14]
    dest = OUT / ('metadata_' + key + '.json')
    log = {'query': q, 'timestamp_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'metadata_path': str(dest)}
    try:
        if not dest.exists():
            b = request(API + 'search?' + urllib.parse.urlencode({'query': q, 'format': 'json', 'resultType': 'core', 'pageSize': '10'}))
            dest.write_bytes(b)
        data = json.loads(dest.read_text())
        rows = data['resultList']['result']
        log['hitCount'] = data['hitCount']
        log['records'] = []
        for row in rows:
            rec = {k: row.get(k) for k in ['id', 'source', 'title', 'doi', 'pmcid', 'firstPublicationDate', 'isOpenAccess']}
            pmc = row.get('pmcid')
            if pmc:
                p = OUT / (pmc + '.xml')
                src = LOCAL / (pmc + '.xml')
                if not p.exists():
                    try:
                        b = src.read_bytes() if src.exists() else request(API + pmc + '/fullTextXML')
                        if b.lstrip().startswith(b'<'):
                            p.write_bytes(b)
                    except Exception as e:
                        rec['fulltext_error'] = str(e)
                if p.exists():
                    rec.update(fulltext_path=str(p), fulltext_sha256=hashlib.sha256(p.read_bytes()).hexdigest(), origin='local_corpus' if src.exists() else 'Europe_PMC_API')
            log['records'].append(rec)
    except Exception as e:
        log['error'] = str(e)
    return log
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
    logs = list(pool.map(work, sorted(queries)))
(ROOT / 'FETCH_LOG.json').write_text(json.dumps(logs, ensure_ascii=False, indent=2) + '\n')
print('queries', len(logs), 'query failures', sum(('error' in x for x in logs)), 'XML', len(list(OUT.glob('PMC*.xml'))))
