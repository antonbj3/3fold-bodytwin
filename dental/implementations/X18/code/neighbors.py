"""Sampled exact roof-edge/neighbor triangle witnesses; full axial fit UNKNOWN."""
from geometry import *
from distance import closest

def nearest(points, tri):
    best = np.inf
    witness = None
    for (i, p) in enumerate(points):
        for start in range(0, len(tri), 2048):
            t = tri[start:start + 2048]
            (dd, q) = closest(np.broadcast_to(p, (len(t), 3)), t)
            k = int(np.argmin(dd))
            if dd[k] < best:
                best = float(dd[k])
                witness = dict(roof_sample_index=i, neighbor_triangle_index=start + k, roof_xyz_mm=p, neighbor_xyz_mm=q[k], midpoint_xyz_mm=(p + q[k]) / 2)
    return dict(minimum_gap_mm=best, witness=witness)

def run():
    allrows = []
    for rnd in ['R1', 'R2']:
        pfile = H / 'raw' / ('PREDICTIONS_' + rnd + '.json')
        if not pfile.exists():
            continue
        preds = json.loads(pfile.read_text())
        current = None
        for p in preds:
            if p['status'] != 'DESIGNED':
                continue
            if p['case'] != current:
                (data, _) = pair(p['case'])
                current = p['case']
            z = np.load(p['file'])
            row = dict(round=rnd, case=p['case'], fdi=p['fdi'], neighbors={})
            for fdi in [p['fdi'] - 1, p['fdi'] + 1]:
                tris = data['lower']['tri'][data['lower']['owner'] == fdi]
                if not len(tris):
                    row['neighbors'][str(fdi)] = dict(status='UNKNOWN_MISSING_NEIGHBOR')
                    continue
                target = data['lower']['tri'][data['lower']['owner'] == p['fdi']]
                (orig, faces) = query_height(target, z['xy'][z['border']], False)
                arms = {}
                for (arm, heights) in [('original', orig), ('practice', z['z_practice'][z['border']]), ('informed', z['z_informed'][z['border']])]:
                    points = np.c_[z['xy'][z['border']], heights]
                    points = points[np.isfinite(heights)]
                    arms[arm] = nearest(points, tris) if len(points) else dict(minimum_gap_mm=None, witness=None)
                    arms[arm]['within_0_1mm_proximity'] = arms[arm]['minimum_gap_mm'] is not None and arms[arm]['minimum_gap_mm'] <= 0.1
                row['neighbors'][str(fdi)] = dict(arms=arms, scope='Finite roof-edge samples only; contact on omitted axial crown sidewalls UNKNOWN')
            allrows.append(row)
    dump(H / 'raw/NEIGHBOR_RESULTS.json', dict(claim_type='capability', rows=allrows, external_referent=dict(kind='published_dataset', locator='Bits2Bites source member hashes in predictions manifests', compared_quantity='Same-site original and generated roof-edge distance witnesses to predicted adjacent teeth', refutes_us=True), axial_crown_contact='UNKNOWN: exported roof is not a whole clinically prepared crown'))
    print('Neighbor roof witnesses', len(allrows))
if __name__ == '__main__':
    run()
