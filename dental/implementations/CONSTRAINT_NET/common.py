from pathlib import Path
import hashlib, json
from functools import lru_cache
HERE = Path(__file__).resolve().parent
DENT = HERE.parents[1]
R = DENT / 'results'

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def dump(p, obj):
    p = HERE / p
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(obj, ensure_ascii=False, indent=2, allow_nan=False) + '\n')

def read(p):
    return json.loads(Path(p).read_text())
MANIFEST = read(HERE / 'SOURCE_MANIFEST.json')
SOURCES = {x['source_file']: x for x in MANIFEST['files']}

def path(name):
    p = Path(name)
    return p if p.is_absolute() else DENT / p

def frozen(name):
    return HERE / SOURCES[str(path(name))]['snapshot']

@lru_cache(maxsize=256)
def src(name):
    return read(frozen(name))

def pointer(obj, key):
    if key == '':
        return obj
    if not key.startswith('/'):
        raise ValueError('RFC6901 pointer required')
    for part in key[1:].split('/'):
        part = part.replace('~1', '/').replace('~0', '~')
        obj = obj[int(part)] if isinstance(obj, list) else obj[part]
    return obj

def esc(s):
    return str(s).replace('~', '~0').replace('/', '~1')

def evidence(name, key, level, unit='1', role='source_value'):
    p = str(path(name))
    rec = SOURCES[p]
    return {'source_file': p, 'key': key, 'value': pointer(src(name), key), 'source_sha256': rec['sha256'], 'snapshot': rec['snapshot'], 'resolution_level': level, 'unit': unit, 'role': role}

def working_binding(nid):
    name = 'results/GRAPH_WORKING_VIEW_20260923/WORKING_VIEW.json'
    d = src(name)
    for branch in ['nodes', 'planning_nodes']:
        if nid in d[branch]:
            return {'id': nid, 'binding_kind': 'quantity_context_not_scientific_dependency', 'working_state': 'DRAFT' if branch == 'nodes' else 'PLANNING_ONLY', 'source_claim': evidence(name, '/' + branch + '/' + esc(nid) + '/claim', 'PHENOMENOLOGICAL', role='working_node_scope')}
    raise KeyError(nid)
