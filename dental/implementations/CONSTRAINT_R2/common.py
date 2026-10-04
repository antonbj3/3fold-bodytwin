"""Small immutable source reader. No third-party packages or external writes."""
from dental_release.paths import expand as _release_expand
import hashlib
import json
from pathlib import Path
HERE = Path(__file__).resolve().parent
DENT = HERE.parents[1]
R1 = HERE.parent / _release_expand('CONSTRAINT_NET')

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

def source_path(name):
    p = Path(name)
    return str(p if p.is_absolute() else DENT / p)

def source(name):
    manifest = read(HERE / 'INPUT_LOCK.json')
    row = manifest['files'][source_path(name)]
    p = HERE / row['snapshot']
    if sha(p) != row['sha256']:
        raise ValueError('SOURCE_HASH: ' + name)
    return read(p)

def evidence(name, key, unit, level, role='residual_operand'):
    row = read(HERE / 'INPUT_LOCK.json')['files'][source_path(name)]
    return {'source_file': source_path(name), 'key': key, 'value': pointer(source(name), key), 'source_sha256': row['sha256'], 'snapshot': row['snapshot'], 'unit': unit, 'resolution_level': level, 'role': role}

def ev_value(row):
    p = HERE / row['snapshot']
    if not p.exists():
        p = R1 / row['snapshot']
    if sha(p) != row['source_sha256']:
        raise ValueError('SOURCE_HASH: ' + str(p))
    return pointer(read(p), row['key'])
