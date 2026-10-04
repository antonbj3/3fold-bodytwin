"""Read complete local ZIP records without modifying an incomplete source archive.

Only sizes recorded in local headers are accepted; data descriptors unsupported.
CRC verification is mandatory for every full member read. Headers may be probed
without allocating or extracting the full member, and are marked unverified.
"""
from dental_release.paths import expand as _release_expand
import gzip
import hashlib
import io
import os
import struct
import zlib
from pathlib import Path
ARCHIVE = Path(_release_expand('@DENTAL_DATA_ROOT@/geometry/Pulpy3D/Pulpy3D_Dataset.zip'))

def index(path=ARCHIVE):
    rows = []
    with open(path, 'rb') as f:
        size = os.fstat(f.fileno()).st_size
        off = 0
        while off + 30 <= size:
            f.seek(off)
            h = f.read(30)
            if h[:4] != b'PK\x03\x04':
                break
            (_, ver, flags, method, t, d, crc, cs, us, nl, el) = struct.unpack('<4s5H3I2H', h)
            name = f.read(nl).decode('utf8')
            extra = f.read(el)
            start = off + 30 + nl + el
            if cs == 4294967295 or us == 4294967295:
                pos = 0
                while pos + 4 <= len(extra):
                    (tag, n) = struct.unpack_from('<HH', extra, pos)
                    v = extra[pos + 4:pos + 4 + n]
                    pos += 4 + n
                    if tag == 1:
                        j = 0
                        if us == 4294967295:
                            us = struct.unpack_from('<Q', v, j)[0]
                            j += 8
                        if cs == 4294967295:
                            cs = struct.unpack_from('<Q', v, j)[0]
            if flags & 8:
                raise ValueError('Untrusted data-descriptor record: ' + name)
            rows.append(dict(name=name, start=start, compressed=cs, uncompressed=us, crc32=crc, method=method, flags=flags, complete=start + cs <= size))
            if start + cs > size:
                break
            off = start + cs
    return rows

def read_member(row, path=ARCHIVE):
    if not row['complete']:
        raise ValueError('Truncated member rejected: ' + row['name'])
    with open(path, 'rb') as f:
        f.seek(row['start'])
        compressed = f.read(row['compressed'])
    if row['method'] == 8:
        raw = zlib.decompress(compressed, -15)
    elif row['method'] == 0:
        raw = compressed
    else:
        raise ValueError('Unsupported compression')
    if len(raw) != row['uncompressed'] or zlib.crc32(raw) & 4294967295 != row['crc32']:
        raise ValueError('CRC/length failure: ' + row['name'])
    return raw

def nifti(row, path=ARCHIVE):
    import nibabel as nib
    raw = read_member(row, path)
    body = gzip.decompress(raw)
    img = nib.Nifti1Image.from_bytes(body)
    return (img, hashlib.sha256(raw).hexdigest())

def nifti_header(row, path=ARCHIVE):
    """Metadata-only probe; full member CRC is not yet validated."""
    import nibabel as nib
    with open(path, 'rb') as f:
        f.seek(row['start'])
        head = f.read(min(row['compressed'], 65536))
    if row['method'] == 8:
        head = zlib.decompressobj(-15).decompress(head)
    with gzip.GzipFile(fileobj=io.BytesIO(head)) as g:
        body = g.read(348)
    h = nib.Nifti1Header.from_fileobj(io.BytesIO(body))
    return dict(shape=list(h.get_data_shape()), spacing=[float(v) for v in h.get_zooms()], dtype=str(h.get_data_dtype()), affine=h.get_best_affine().tolist(), crc_verified=False)
