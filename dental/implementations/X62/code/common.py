"""Small immutable-source readers. Patient datasets are read locally, in ZIPs."""
from dental_release.paths import expand as _release_expand
from pathlib import Path
import hashlib, json, datetime, xml.etree.ElementTree as ET
ROOT = Path(__file__).resolve().parents[1]
CORPUS = Path(_release_expand('@DENTAL_CORPUS_ROOT@/europepmc/fulltext'))
ZIP = Path(_release_expand('@DENTAL_DATA_ROOT@/geometry/Open-Full-Jaw/dataset/@DENTAL_CASE_ID@.zip'))
DATA = Path(_release_expand('@DENTAL_WORK_ROOT@/X62_ortho'))

def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda : f.read(1024 * 1024), b''):
            h.update(b)
    return h.hexdigest()

def save(name, obj):
    p = ROOT / name
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(obj, ensure_ascii=False, indent=2, allow_nan=False) + '\n')

def utc():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()

def state(status, gate, next_operation):
    p = ROOT / 'CURRENT_WORK_STATE.json'
    d = json.loads(p.read_text())
    d.update(status=status, latest_gate=gate, next_operation=next_operation, updated_utc=utc())
    save('CURRENT_WORK_STATE.json', d)

def text(e):
    return ' '.join(' '.join(e.itertext()).split()) if e is not None else ''

def article(pmc):
    return ET.parse(CORPUS / (pmc + '.xml')).getroot()

def tables(pmc):
    r = article(pmc)
    return {t.get('id'): [[text(c) for c in row] for row in t.findall('.//tr')] for t in r.findall('.//table-wrap')}
