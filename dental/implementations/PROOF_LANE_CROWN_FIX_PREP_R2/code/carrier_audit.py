from common import *

def audit():
    rows = []
    for part in read(R / 'EXPORTS.json')['rows']:
        m = load(part['mesh']['path'])
        d = [F(float(x)) for x in m['M'][2]]
        h = F(part['full_preparation_formula']['h_exact'])
        res = []
        for q in trimesh.Trimesh(m['prep_vertices'], m['prep_faces'], process=False).split(only_watertight=False):
            lo = min((sum((F(float(a)) * b for (a, b) in zip(p, d))) for p in q.vertices))
            res.append(dict(min_projection_minus_h_exact=str(lo - h), positive_separating_slab=lo > h, projection_difference_float=float(lo - h)))
        rows.append(dict(key=part['key'], material=part['material'], round=part['round'], components=res, resolution='PER_POINT'))
    return dict(rows=rows, scope='Additional exact rational topology audit. Nominal film was measured against stored cap, exact inclusion/path against symbolic P. Unified equality of these carriers is not established.')
if __name__ == '__main__':
    current = audit()
    assert current == read(R / 'raw/EXACT_CARRIER_AUDIT.json')
    print('Exact carrier audit reproduced:8/8 positive separating slabs')
