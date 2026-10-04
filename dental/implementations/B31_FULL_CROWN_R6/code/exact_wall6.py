from local import *
from score_a import exact_failure

def run():
    (a, b) = first()
    d = npz(D6 / 'D_exterior.npz')
    old = npz(D / 'B/distance_local_thickening' / (a['key'] + '.npz'))
    ext = d['vertices'][d['faces']]
    inner = old['vertices'][old['faces'][old['roles'] == 1]]
    (q, dd, j) = fast_nearest(inner, ext.mean(1))
    i = int(dd.argmin())
    r = exact_failure(ext[i], inner[j[i]])
    r.update(outer_triangle=ext[i], inner_triangle=inner[j[i]], outer_face_id=i, inner_face_id=int(j[i]), nearest_centroid_distance_mm=float(dd[i]), resolution='PER_POINT', scope='Exact rational binary64-geometry witness of a wall failure; no scan/process enclosure', outer_source=str(D6 / 'D_exterior.npz'), inner_source=str(D / 'B/distance_local_thickening' / (a['key'] + '.npz')))
    save(R6 / 'raw/EXACT_WALL_WITNESS.json', r)
    assert r['fails']
    print('Exact rational wall failure verified', r['distance_squared_mm2'])
if __name__ == '__main__':
    run()
