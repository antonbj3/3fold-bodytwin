"""Freeze DOI identities from the primary registration agency; no datasets fetched."""
from pathlib import Path
import concurrent.futures
import datetime
import hashlib
import json
import urllib.parse
import urllib.request
ROOT = Path(__file__).resolve().parent
REFS = [('10.1148/ryai.240300', 'Checklist for Artificial Intelligence in Medical Imaging'), ('10.1038/s41592-023-02151-z', 'Metrics reloaded'), ('10.1038/s41592-023-02150-0', 'Understanding metric-related pitfalls'), ('10.1016/j.media.2026.104095', 'Multi-structure segmentation in CBCT volumes'), ('10.1109/CVPR52734.2025.00494', 'Segmenting Maxillofacial Structures'), ('10.1109/ACCESS.2024.3408629', 'Enhancing Patch-Based Learning'), ('10.1038/s41598-022-20605-w', 'Comparison of deep learning segmentation'), ('10.1111/clr.13578', 'Guidance means accuracy'), ('10.1186/s40729-020-00272-0', 'Accuracy of dynamic navigation')]

def fetch(entry):
    (doi, expected) = entry
    url = 'https://api.crossref.org/works/' + urllib.parse.quote(doi, safe='')
    req = urllib.request.Request(url, headers={'User-Agent': 'Dental-X41-local-reference-check/1.0'})
    try:
        data = urllib.request.urlopen(req, timeout=30).read()
        d = json.loads(data)['message']
        name = hashlib.sha256(doi.encode()).hexdigest()[:16] + '.json'
        (ROOT / 'references/metadata' / name).write_bytes(data)
        title = d['title'][0]
        ok = d['DOI'].lower() == doi.lower() and expected.lower() in title.lower()
        return {'doi': doi, 'title': title, 'container': d.get('container-title'), 'year': d.get('published', {}).get('date-parts'), 'authors': d.get('author'), 'volume': d.get('volume'), 'issue': d.get('issue'), 'page': d.get('page'), 'article_number': d.get('article-number'), 'url': d.get('URL'), 'identity_verified': ok, 'primary_metadata_url': url, 'metadata_file': 'references/metadata/' + name, 'sha256': hashlib.sha256(data).hexdigest()}
    except Exception as exc:
        return {'doi': doi, 'identity_verified': False, 'error': str(exc)}

def main():
    if (ROOT / 'REFERENCES_VERIFIED.json').exists():
        raise SystemExit('Reference snapshot already exists; do not silently refresh')
    (ROOT / 'references/metadata').mkdir(parents=True, exist_ok=True)
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
        rows = list(pool.map(fetch, REFS))
    d = {'checked_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'method': 'Live Crossref registered DOI and expected title; quantitative locators checked separately', 'references': rows}
    (ROOT / 'REFERENCES_VERIFIED.json').write_text(json.dumps(d, indent=2, ensure_ascii=False) + '\n')
    print(json.dumps([{'doi': r['doi'], 'title': r.get('title'), 'pass': r['identity_verified'], 'error': r.get('error')} for r in rows], ensure_ascii=False))
    if not all((r['identity_verified'] for r in rows)):
        raise SystemExit('Some references unresolved; preserve this attempt and repair in a new record')
if __name__ == '__main__':
    main()
