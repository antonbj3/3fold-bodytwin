"""Trusted evaluator input boundary. All paths inside a manifest are relative.

The release anchor and evaluator are trusted. This cannot defend against an owner
replacing both the trusted evaluator and its independently recorded release hash.
Every consumption verifies the same bytes that are parsed (not hash-then-reopen).
"""
import os, stat, json, hashlib, io, zipfile, math, re
from pathlib import Path, PurePosixPath

class IntegrityError(ValueError):
    pass

def strict_json(blob):

    def pairs(items):
        d = {}
        for (k, v) in items:
            if k in d:
                raise IntegrityError('duplicate JSON key: ' + k)
            d[k] = v
        return d

    def bad(x):
        raise IntegrityError('nonfinite JSON: ' + x)

    def finite_float(x):
        n = float(x)
        if not math.isfinite(n):
            raise IntegrityError('JSON numeric overflow: ' + x)
        return n
    return json.loads(blob, object_pairs_hook=pairs, parse_constant=bad, parse_float=finite_float)

def h(blob):
    return hashlib.sha256(blob).hexdigest()

def relative(name):
    if not isinstance(name, str) or not name or '\\' in name or ('\x00' in name):
        raise IntegrityError('invalid relative path')
    p = PurePosixPath(name)
    if p.is_absolute() or '..' in p.parts or str(p) != name or (name == '.'):
        raise IntegrityError('noncanonical relative path: ' + name)
    return p.parts

def safe_bytes(root, name, limit=100000000):
    parts = relative(name)
    fd = os.open(root, os.O_RDONLY | os.O_DIRECTORY)
    try:
        for part in parts[:-1]:
            n = os.open(part, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW, dir_fd=fd)
            os.close(fd)
            fd = n
        f = os.open(parts[-1], os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK, dir_fd=fd)
        try:
            st = os.fstat(f)
            if not stat.S_ISREG(st.st_mode) or st.st_size > limit:
                raise IntegrityError('not a bounded regular file: ' + name)
            with os.fdopen(f, 'rb', closefd=False) as stream:
                blob = stream.read(limit + 1)
            if len(blob) > limit:
                raise IntegrityError('input grew beyond bound')
            return blob
        finally:
            os.close(f)
    except OSError as e:
        raise IntegrityError('unsafe/missing input ' + name + ': ' + str(e)) from e
    finally:
        os.close(fd)

def file_set(root):
    out = set()
    for (base, dirs, files) in os.walk(root, followlinks=False):
        for name in dirs + files:
            p = Path(base) / name
            if p.is_symlink():
                raise IntegrityError('symlink forbidden: ' + str(p.relative_to(root)))
        for name in files:
            out.add((Path(base) / name).relative_to(root).as_posix())
    return out

class Bundle:

    def __init__(self, root=None):
        self.root = Path(root or Path(__file__).resolve().parents[1]).resolve()
        from release_anchor import BENCHMARK_SHA256
        b = safe_bytes(self.root, 'BENCHMARK_LOCK.json', 50000000)
        if h(b) != BENCHMARK_SHA256:
            raise IntegrityError('untrusted benchmark manifest')
        self.lock = strict_json(b)
        if self.lock.get('schema') != 'DentalGenCAD-Bench-v3':
            raise IntegrityError('benchmark identity')
        self.payload = (self.root / 'payload').resolve()
        self.files = self.lock['files']
        for name in self.files:
            relative(name)
        self.verify()

    def location(self, name):
        relative(name)
        if name.startswith('payload/'):
            return (self.payload, name[len('payload/'):])
        return (self.root, name)

    def bytes(self, name):
        if name not in self.files:
            raise IntegrityError('unlisted input: ' + name)
        (root, rel) = self.location(name)
        entry = self.files[name]
        blob = safe_bytes(root, rel, entry['bytes'])
        if len(blob) != entry['bytes'] or h(blob) != entry['sha256']:
            raise IntegrityError('locked input drift: ' + name)
        return blob

    def json(self, name):
        return strict_json(self.bytes(name))

    def npz(self, name):
        return load_npz(self.bytes(name))

    def verify(self):
        for prefix in self.lock['closed_directories']:
            (root, rel) = self.location(prefix)
            base = root / rel
            actual = {prefix + '/' + p for p in file_set(base)} - {'code/release_anchor.py'}
            expected = {p for p in self.files if p.startswith(prefix + '/')}
            if actual != expected:
                raise IntegrityError('file set changed: ' + prefix)
        for name in self.files:
            self.bytes(name)

    def tasks(self):
        for name in sorted(self.files):
            if name.startswith('payload/public/tasks/') and name.endswith('.json'):
                yield (name, self.json(name))

def load_npz(blob, max_expanded=100000000):
    import numpy as np
    try:
        with zipfile.ZipFile(io.BytesIO(blob)) as z:
            names = z.namelist()
            if len(names) > 512 or len(set(names)) != len(names) or sum((i.file_size for i in z.infolist())) > max_expanded:
                raise IntegrityError('NPZ duplicate names or expanded size')
            if any(('/' in n or '\\' in n or len(n) > 128 or (not n.endswith('.npy')) for n in names)):
                raise IntegrityError('NPZ key path')
            out = {}
            total = 0
            for name in names:
                data = z.read(name)
                total += len(data)
                if total > max_expanded:
                    raise IntegrityError('actual NPZ expanded size')
                stream = io.BytesIO(data)
                version = np.lib.format.read_magic(stream)
                if version == (1, 0):
                    (shape, fortran, dtype) = np.lib.format.read_array_header_1_0(stream)
                elif version == (2, 0):
                    (shape, fortran, dtype) = np.lib.format.read_array_header_2_0(stream)
                else:
                    raise IntegrityError('unsupported NPY header version')
                if dtype.kind not in 'fiubU' or len(shape) > 3 or any((not isinstance(n, int) or n < 0 for n in shape)):
                    raise IntegrityError('unsupported NPY type/shape')
                size = math.prod(shape) * dtype.itemsize
                if size > max_expanded or stream.tell() + size != len(data):
                    raise IntegrityError('NPY shape allocation/length mismatch')
                out[name[:-4]] = np.load(io.BytesIO(data), allow_pickle=False)
            return out
    except (ValueError, OSError, EOFError, zipfile.BadZipFile) as e:
        if isinstance(e, IntegrityError):
            raise
        raise IntegrityError('malformed NPZ: ' + str(e)) from e

def predictions(root, expected_digest, expected_cases, benchmark_digest):
    root = Path(root)
    blob = safe_bytes(root, 'FROZEN_PREDICTIONS.json', 2000000)
    if h(blob) != expected_digest:
        raise IntegrityError('prediction freeze anchor mismatch')
    manifest = strict_json(blob)
    if manifest['benchmark_sha256'] != benchmark_digest:
        raise IntegrityError('wrong prediction benchmark')
    names = manifest.get('participants')
    if not isinstance(names, list) or not names or any((not isinstance(n, str) or not re.fullmatch('[a-z][a-z0-9_]{0,63}', n) for n in names)) or (len(set(names)) != len(names)):
        raise IntegrityError('nonempty unique participant domain required')
    files = manifest['files']
    wanted = {name + '/' + case + '.npz' for name in manifest['participants'] for case in expected_cases}
    if set(files) != wanted:
        raise IntegrityError('prediction task coverage')
    actual = file_set(root) - {'FROZEN_PREDICTIONS.json', 'COST.json'}
    if actual != wanted:
        raise IntegrityError('extra/missing prediction files')
    for (name, entry) in files.items():
        b = safe_bytes(root, name, 2000000)
        if h(b) != entry['sha256'] or len(b) != entry['bytes']:
            raise IntegrityError('prediction changed: ' + name)
    return manifest
