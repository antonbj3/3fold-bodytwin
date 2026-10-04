from common import *

def write(path, tri):
    with open(path, 'w') as h:
        h.write('solid research_geometry\n')
        for t in tri:
            n = np.cross(t[1] - t[0], t[2] - t[0])
            nn = np.linalg.norm(n)
            n = n / nn if nn else n
            h.write(' facet normal ' + ' '.join(('%.17g' % x for x in n)) + '\n  outer loop\n')
            for p in t:
                h.write('   vertex ' + ' '.join(('%.17g' % x for x in p)) + '\n')
            h.write('  endloop\n endfacet\n')
        h.write('endsolid research_geometry\n')
    loaded = trimesh.load_mesh(path, file_type='stl', process=False)
    return dict(path=str(path), sha256=sha(path), exact_triangle_roundtrip=bool(np.array_equal(loaded.triangles, tri)), bytes=Path(path).stat().st_size, scope='17-digit ASCII STL parsed as binary64 by supplied pinned reader. Downstream importer quantization is outside this certificate.')
