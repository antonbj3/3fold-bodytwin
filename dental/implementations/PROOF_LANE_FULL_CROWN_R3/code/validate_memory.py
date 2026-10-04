from common_r3 import *
from memory_height import height
old = module('old_geometry', OLD / 'code/vendor/geometry.py')
p = npz(DATA / 'inputs/8017101de9a27cf6_premolar_crown.npz')
tri = p['prior_vertices'][p['prior_faces'][p['prior_roles'] == 0]]
xy = p['contact_xy'][::7]
rows = []
for upper in [False, True]:
    a = old.height(tri, xy, upper)
    b = height(tri, xy, upper)
    rows.append(dict(upper=upper, bitidentical=np.array_equal(a, b, equal_nan=True), finite_max_abs_error_mm=float(abs(a[np.isfinite(a)] - b[np.isfinite(a)]).max()), queries=len(xy)))
    assert np.array_equal(a, b, equal_nan=True)
dump(ROOT / 'raw/MEMORY_EQUIVALENCE.json', dict(rows=rows, scope='same affine arithmetic and global-radius search, representative measured mesh; no precision reduction in frozen operations', old_source_sha256=sha(OLD / 'code/vendor/geometry.py'), new_source_sha256=sha(ROOT / 'code/memory_height.py')))
