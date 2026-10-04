from dental_release.paths import expand as _release_expand
from pathlib import Path
import datetime, hashlib, json, time, os, xml.etree.ElementTree as ET
PACKAGE = Path(__file__).resolve().parents[1]
ROOT = Path(os.environ.get('X10_OUTPUT_DIR', str(PACKAGE))).resolve()
PROJECT = Path(os.environ.get('DENTAL_PROJECT_ROOT', str(PACKAGE.parents[1]))).resolve()
CORPUS = Path(_release_expand('@DENTAL_CORPUS_ROOT@/europepmc'))

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def write(name, value):
    path = ROOT / name
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + '\n')
    return value

def load(name):
    return json.loads((ROOT / name).read_text())

def frozen(name):
    path = ROOT / name
    expected = path.with_suffix('.sha256').read_text().split()[0]
    if sha(path) != expected:
        raise ValueError('Frozen preregistration hash drift: ' + name)
    return load(name)

def freeze_once(name, value):
    path = ROOT / name
    if path.exists():
        old = load(name)
        if old['predictions'] != value['predictions'] or old['inputs'] != value['inputs']:
            raise ValueError('Frozen predictions changed: ' + name)
        if sha(path) != path.with_suffix('.sha256').read_text().split()[0]:
            raise ValueError('Frozen prediction hash drift')
        return old
    value['frozen_utc'] = datetime.datetime.now(datetime.timezone.utc).isoformat()
    write(name, value)
    path.with_suffix('.sha256').write_text(sha(path) + '  ' + path.name + '\n')
    return value

def state(stage, last_gate, next_operation):
    s = load('CURRENT_WORK_STATE.json')
    s.update(stage=stage, last_gate=last_gate, next_operation=next_operation, updated_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
    write('CURRENT_WORK_STATE.json', s)

def flatten(e):
    return ' '.join(''.join(e.itertext()).split())

def tables(pmc):
    path = CORPUS / 'fulltext' / f'{pmc}.xml'
    tree = ET.parse(path).getroot()
    out = []
    for t in tree.findall('.//table-wrap'):
        rows = []
        for row in t.findall('.//tr'):
            cells = [flatten(c) for c in row if c.tag in ('td', 'th')]
            rows.append(cells)
        out.append({'id': t.attrib.get('id'), 'caption': flatten(t.find('caption')) if t.find('caption') is not None else '', 'rows': rows})
    return (path, out, [flatten(e) for e in tree.findall('.//license')])

def clinical_primary_record(pmid):
    for name in ['marginal_fit', 'cadcam', 'prosthodontics']:
        for line in (CORPUS / 'meta' / f'{name}.jsonl').open():
            if f'"{pmid}"' in line:
                d = json.loads(line)
                if d.get('pmid') == pmid:
                    return (d, hashlib.sha256(line.encode()).hexdigest())
    raise ValueError('Missing local primary metadata ' + pmid)
