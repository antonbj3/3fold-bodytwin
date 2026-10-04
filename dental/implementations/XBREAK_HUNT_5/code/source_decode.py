"""Decode the published VTK compressed appended arrays without changing payload."""
import xml.etree.ElementTree as ET, numpy as np, base64, zlib, math

def decode(path):
    root = ET.parse(path).getroot()
    piece = root.find('PolyData/Piece')
    stream = root.find('AppendedData').text.strip()[1:]
    types = {'Float32': '<f4', 'Float64': '<f8', 'Int64': '<i8', 'UInt32': '<u4', 'Int32': '<i4'}

    def array(a):
        off = int(a.attrib['offset'])
        s = stream[off:]
        first = np.frombuffer(base64.b64decode(s[:16]), dtype='<u4')
        blocks = int(first[0])
        size = int(first[1])
        last = int(first[2])
        hl = 4 * math.ceil(4 * (3 + blocks) / 3)
        header = np.frombuffer(base64.b64decode(s[:hl]), dtype='<u4')
        lens = header[3:]
        cl = int(lens.sum())
        bl = 4 * math.ceil(cl / 3)
        compressed = base64.b64decode(s[hl:hl + bl])
        chunks = []
        pos = 0
        for (k, cs) in enumerate(lens):
            raw = zlib.decompress(compressed[pos:pos + int(cs)])
            pos += int(cs)
            expected = last if k == blocks - 1 else size
            if len(raw) != expected:
                raise ValueError('decompression byte count')
            chunks.append(raw)
        data = np.frombuffer(b''.join(chunks), dtype=types[a.attrib['type']]).copy()
        return data
    v = array(piece.find('Points/DataArray')).reshape(-1, 3)
    f = array(piece.find('Polys/DataArray[@Name="connectivity"]')).reshape(-1, 3)
    offsets = array(piece.find('Polys/DataArray[@Name="offsets"]'))
    lab = array(piece.find('CellData/DataArray[@Name="Label"]'))
    if not np.array_equal(offsets, np.arange(1, len(f) + 1) * 3):
        raise ValueError('polygon offsets not triples')
    if len(v) != int(piece.attrib['NumberOfPoints']) or len(f) != int(piece.attrib['NumberOfPolys']) or len(lab) != len(f):
        raise ValueError('XML counts disagree')
    if np.max(f) >= len(v) or np.min(f) < 0 or (not np.all(np.isfinite(v))):
        raise ValueError('invalid geometry')
    if not set(np.unique(lab)).issubset({0.0, 1.0}):
        raise ValueError('invalid labels')
    return (v.astype(np.float64), f.astype(np.int64), lab.astype(np.int64))
