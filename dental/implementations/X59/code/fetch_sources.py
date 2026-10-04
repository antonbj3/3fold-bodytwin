"""Optional original acquisition replay. run_all.sh uses preserved local files only."""
from datetime import datetime, timezone
from pathlib import Path
import hashlib, json, urllib.request, urllib.error
ROOT = Path(__file__).resolve().parents[1]

def get(url):
    with urllib.request.urlopen(url, timeout=40) as r:
        return (r.read(), r.headers.get('Content-Type'))

def main():
    stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S')
    target = ROOT / 'raw' / ('acquisition_' + stamp)
    target.mkdir()
    (metadata, mime) = get('https://api.figshare.com/v2/articles/5669308')
    j = json.loads(metadata)
    (target / 'metadata.json').write_bytes(metadata)
    total = sum((f['size'] for f in j['files']))
    if total >= 2000000000:
        raise ValueError('Source exceeds frozen download size gate')
    if j['license']['name'] != 'CC BY 4.0':
        raise ValueError('Source license changed; do not silently proceed')
    outcomes = []
    for f in j['files']:
        try:
            (data, mime) = get(f['download_url'])
            if len(data) != f['size'] or hashlib.md5(data).hexdigest() != f['computed_md5']:
                raise ValueError('File size/checksum differs from source metadata')
            (target / f['name']).write_bytes(data)
            outcomes.append(dict(file=f['name'], status='DOWNLOADED', bytes=len(data), sha256=hashlib.sha256(data).hexdigest()))
        except urllib.error.HTTPError as e:
            outcomes.append(dict(file=f['name'], status='HTTP_FAILURE', http_status=e.code, url=f['download_url']))
    pdfurl = 'https://www.scielo.br/j/bdj/a/SqXx7S5vR7yn3gSCdTkMrrK/?format=pdf&lang=en'
    (pdf, mime) = get(pdfurl)
    if not pdf.startswith(b'%PDF'):
        raise ValueError('Original publication response is not PDF')
    if len(pdf) > 50000000:
        raise ValueError('Unexpected publication size; preserve metadata and stop')
    (target / 'cunali2017.pdf').write_bytes(pdf)
    (target / 'manifest.json').write_text(json.dumps(dict(created_utc=stamp, declared_total_bytes=total, license=j['license'], files=outcomes, original_pdf=dict(url=pdfurl, bytes=len(pdf), sha256=hashlib.sha256(pdf).hexdigest())), indent=2) + '\n')
    print(json.dumps(dict(acquisition=str(target), declared_bytes=total, downloaded=sum((o['status'] == 'DOWNLOADED' for o in outcomes)), original_pdf_bytes=len(pdf))))
if __name__ == '__main__':
    main()
