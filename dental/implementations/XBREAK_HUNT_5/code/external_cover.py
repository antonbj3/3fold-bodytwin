from common import *
from shapely.geometry import Polygon, box
start = time.perf_counter()
rows = []
for case in ['6158-11', '6138-21']:
    z = np.load(DATA / (case + '_R1.npz'))
    v = z['v']
    f = z['f'][z['lab'] == 1]
    col = np.load(DATA / f'{case}_COLUMNS_h0.03.npz')
    h = float(col['h'])
    height = col['height']
    qtri = (v @ col['basis'] - col['offset'])[f]
    min_clearance = float('inf')
    max_area_error = 0.0
    degenerate = 0
    checked = 0
    for face in qtri:
        polygon = Polygon(face[:, :2])
        area = polygon.area
        det = np.linalg.det(np.c_[face[:, :2], np.ones(3)])
        if abs(det) < 1e-14:
            degenerate += 1
            continue
        plane = np.linalg.solve(np.c_[face[:, :2], np.ones(3)], face[:, 2])
        lo = np.floor(face[:, :2].min(axis=0) / h).astype(int)
        hi = np.floor(face[:, :2].max(axis=0) / h).astype(int)
        total = 0.0
        for i in range(lo[0], hi[0] + 1):
            for j in range(lo[1], hi[1] + 1):
                p = polygon.intersection(box(i * h, j * h, (i + 1) * h, (j + 1) * h))
                if p.is_empty:
                    continue
                total += p.area
                if p.geom_type == 'Polygon':
                    xy = np.array(p.exterior.coords)[:, :2]
                elif p.geom_type == 'LineString':
                    xy = np.array(p.coords)[:, :2]
                elif p.geom_type == 'Point':
                    xy = np.array(p.coords)[:, :2]
                else:
                    raise ValueError('unexpected external clipping geometry ' + p.geom_type)
                roof = np.c_[xy, np.ones(len(xy))] @ plane
                clearance = height[i, j] * h - roof.max()
                min_clearance = min(min_clearance, float(clearance))
                checked += 1
        max_area_error = max(max_area_error, float(abs(total - area)))
    row = {'case': case, 'source_triangles': len(qtri), 'external_rectangle_intersections': checked, 'degenerate_projected_faces': degenerate, 'minimum_external_cover_clearance': min_clearance, 'max_projected_area_partition_error': max_area_error, 'gate': 'PASS' if degenerate == 0 and min_clearance >= -1e-05 and (max_area_error <= 1e-05) else 'FAIL_OR_UNKNOWN', 'method': 'GEOS/Shapely polygon rectangle intersection plus independent affine-plane reconstruction; original source facets, no source clipping code reused'}
    rows.append(row)
    print(json.dumps(row), flush=True)
write('raw/EXTERNAL_COVER.json', {'rows': rows, 'wall_seconds': time.perf_counter() - start, 'gate': 'PASS' if all((r['gate'] == 'PASS' for r in rows)) else 'FAIL_OR_UNKNOWN', 'independent_scientific_review': False})
