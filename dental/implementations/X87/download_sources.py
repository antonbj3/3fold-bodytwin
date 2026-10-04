"""One-time publication retrieval; demo replay never needs network."""
from pathlib import Path
import urllib.request, json, hashlib, shutil, subprocess
from datetime import datetime, timezone
R = Path(__file__).resolve().parent
S = R / 'sources'
sources = [('meta2024.xml', 'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC10956322/fullTextXML'), ('multigrader2022.xml', 'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC9633839/fullTextXML'), ('Gerlach2010.json', 'https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=EXT_ID:20036043%20AND%20SRC:MED&format=json&resultType=core'), ('IACAT2.pdf', 'https://federicobolelli.it/media/publications/pdfs/2023iciap_iacat2.pdf')]
manifest = []
for (name, url) in sources:
    try:
        req = urllib.request.Request(url, headers={'User-Agent': 'Dental research source verification'})
        with urllib.request.urlopen(req, timeout=30) as response:
            data = response.read(4000001)
        if len(data) > 4000000:
            raise ValueError('Source exceeds 4MB cap')
        (S / name).write_bytes(data)
        manifest.append({'file': name, 'url': url, 'sha256': hashlib.sha256(data).hexdigest(), 'bytes': len(data), 'status': 'RETRIEVED', 'retrieved_utc': datetime.now(timezone.utc).isoformat()})
    except Exception as exc:
        manifest.append({'file': name, 'url': url, 'status': 'FAILED', 'error': str(exc)})
for n in ['Varga2020.txt', 'PMC7683639.xml']:
    p = R.parent / 'LANE_X8_GUIDE_NERVE_RISK/sources' / n
    shutil.copyfile(p, S / n)
    manifest.append({'file': n, 'source': str(p), 'sha256': hashlib.sha256((S / n).read_bytes()).hexdigest(), 'bytes': p.stat().st_size, 'status': 'REUSED_LOCAL'})
if (S / 'IACAT2.pdf').exists():
    subprocess.run(['pdftotext', '-layout', str(S / 'IACAT2.pdf'), str(S / 'IACAT2.txt')], check=True)
(S / 'MANIFEST.json').write_text(json.dumps(manifest, indent=2) + '\n')
print(json.dumps([{'file': m['file'], 'status': m['status'], 'bytes': m.get('bytes')} for m in manifest]))
