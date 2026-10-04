"""Export small synthetic calibration witnesses, not patient crown designs."""
from common import *
import math
out = ROOT / 'exports'
out.mkdir(exist_ok=True)

def write_stl(name, triangles):
    path = out / name
    with path.open('w') as f:
        f.write('solid calibration_witness\n')
        for tri in triangles:
            f.write('facet normal 0 0 0\nouter loop\n')
            for v in tri:
                f.write('vertex %.9f %.9f %.9f\n' % tuple(v))
            f.write('endloop\nendfacet\n')
        f.write('endsolid calibration_witness\n')
    return {'path': str(path), 'sha256': sha(path), 'bytes': path.stat().st_size, 'unit': 'mm'}
N = 96
alpha = math.radians(3)

def ring(z, r):
    return [(r * math.cos(2 * math.pi * k / N), r * math.sin(2 * math.pi * k / N), z) for k in range(N)]

def side(p, q, reverse=False):
    ts = []
    for k in range(N):
        l = (k + 1) % N
        ts += [(p[k], p[l], q[l]), (p[k], q[l], q[k])]
    return [tuple(reversed(t)) for t in ts] if reverse else ts
base = ring(0, 5)
top = ring(5, 5 - 5 * math.tan(alpha))
tri = side(base, top)
for k in range(N):
    j = (k + 1) % N
    tri += [((0, 0, 0), base[j], base[k]), ((0, 0, 5), top[k], top[j])]
manifest = [write_stl('die_half_taper_3deg.stl', tri)]
for (e, label) in [(57.5, 'expanded'), (-57.5, 'contracted')]:
    b = (40 + e) / 1000 / math.cos(alpha)
    z0 = 0.04
    z1 = 5.04
    i0 = ring(z0, 5 - z0 * math.tan(alpha) + b)
    i1 = ring(z1, 5 - z1 * math.tan(alpha) + b)
    o0 = ring(z0, 6 - z0 * math.tan(alpha))
    o1 = ring(z1, 6 - z1 * math.tan(alpha))
    tri = side(o0, o1) + side(i0, i1, True) + side(i0, o0) + side(i1, o1, True)
    manifest.append(write_stl('ring_' + label + '_signed57p5um.stl', tri))
write('EXPORT_MANIFEST.json', {'files': manifest, 'kind': 'our_own_fixture', 'purpose': 'analytic cone-ring witnesses for independent lab scanner/metrology test', 'projected_gap_definition': 'vertical margin-plane separation, not absolute marginal discrepancy or Euclidean closest-edge SDF gap', 'limitations': 'faceted planar-cone approximation; unseated exported coordinates; no tool compensation or sintering calibration; not a clinical device', 'expected_axial_translation_um': {'expanded': 0.0, 'contracted': (57.5 - 40) / math.sin(alpha)}})
print('Synthetic mm STL witnesses exported.')
