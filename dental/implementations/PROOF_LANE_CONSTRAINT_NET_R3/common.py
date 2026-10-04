from dental_release.paths import expand as _release_expand
import hashlib, json
from pathlib import Path
from functools import lru_cache
HERE = Path(__file__).resolve().parent
DENT = Path(_release_expand('@DENTAL_INPUT_ROOT@/workspace'))

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def read(path):
    return json.loads(Path(path).read_text())

def dump(name, value):
    (HERE / name).write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + '\n')

def pointer(value, path):
    if not path:
        return value
    for part in path.lstrip('/').split('/'):
        part = part.replace('~1', '/').replace('~0', '~')
        value = value[int(part)] if isinstance(value, list) else value[part]
    return value

@lru_cache(None)
def lock():
    return read(HERE / 'INPUT_LOCK.json')['files'] | read(HERE / 'AUX_INPUT_LOCK.json')['files']

def logical(name):
    return str(Path(name) if Path(name).is_absolute() else DENT / name)

@lru_cache(None)
def source(name):
    row = lock()[logical(name)]
    p = HERE / row['snapshot']
    if sha(p) != row['sha256']:
        raise ValueError('SOURCE_HASH:' + name)
    return read(p)

def ev_value(row):
    p = HERE / row['snapshot']
    expected = lock().get(row['source_file'])
    if not expected or expected['sha256'] != row['source_sha256'] or expected['snapshot'] != row['snapshot']:
        raise ValueError('SOURCE_HASH')
    return pointer(source(row['source_file']), row['key'])

def evidence(name, key, unit, level, role='residual_operand'):
    name = logical(name)
    row = lock()[name]
    return {'source_file': name, 'key': key, 'value': pointer(source(name), key), 'source_sha256': row['sha256'], 'snapshot': row['snapshot'], 'unit': unit, 'resolution_level': level, 'role': role}

def verify_frozen():
    for n in ['PREREG_R3A.json', 'PREREG_R3B.json', 'FROZEN_PREDICTIONS.json', 'SEMANTIC_RULES.json', 'MEASUREMENT_CONTRACTS.json']:
        assert sha(HERE / n) == (HERE / (n + '.sha256')).read_text().split()[0], n
    for (f, r) in lock().items():
        assert sha(HERE / r['snapshot']) == r['sha256'], f
