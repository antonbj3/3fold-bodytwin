"""Adversarial check of digital surface discrepancy and conventional control output."""
from registration import *

def main():
    t = time.monotonic()
    w = np.load(DATA / 'r2_witnesses.npz')
    tri = stl(DATA / 'hao_demo/Demo/Demo_1/CBCT/reconstruction/out_smoothed.stl')
    cent = tri.mean(1)
    rows = []
    for jaw in ['lower', 'upper']:
        mask = cent[:, 2] < 39 if jaw == 'lower' else cent[:, 2] >= 39
        part = tri[mask]
        p = w[jaw + '_source_points']
        sel = np.isin(w[jaw + '_fdi'], [11, 21, 31, 41])
        T = w[jaw + '_control_T']
        x = transform(p[sel][::20], T)
        bad = x + np.array([10.0, 0.0, 0.0])
        (d, _, _) = nearest_surface(bad, part.reshape(-1, 3), np.arange(part.size // 3).reshape(-1, 3))
        rows.append({'jaw': jaw, 'check': 'conventional_ICP_control_pose_injection', 'injection': '10 mm x translation', 'p95_mm': float(np.percentile(d, 95)), 'gate_mm': 0.5, 'fault_rejected': bool(np.percentile(d, 95) > 0.5), 'resolution': 'PER_SURFACE_REGION'})
        x = w[jaw + '_transformed_points']
        q = w[jaw + '_closest_cbct_surface']
        dist = w[jaw + '_residual_mm']
        cell = w[jaw + '_surface_cell']
        ids = np.arange(0, len(p), max(1, len(p) // 100))
        face = part[cell[ids]]
        qq = q[ids]
        e = face[:, 1] - face[:, 0]
        f = face[:, 2] - face[:, 0]
        a = qq - face[:, 0]
        ee = np.einsum('ij,ij->i', e, e)
        ff = np.einsum('ij,ij->i', f, f)
        ef = np.einsum('ij,ij->i', e, f)
        ea = np.einsum('ij,ij->i', e, a)
        fa = np.einsum('ij,ij->i', f, a)
        det = ee * ff - ef * ef
        u = (ff * ea - ef * fa) / det
        v = (ee * fa - ef * ea) / det
        recon = face[:, 0] + u[:, None] * e + v[:, None] * f
        err = float(np.max(np.linalg.norm(recon - qq, axis=1)))
        within = bool(np.all((u >= -1e-07) & (v >= -1e-07) & (u + v <= 1 + 1e-07)))
        metricerr = float(np.max(np.abs(np.linalg.norm(x[ids] - qq, axis=1) - dist[ids])))
        badmetricerr = float(np.max(np.abs(np.linalg.norm(x[ids] - (qq + 10), axis=1) - dist[ids])))
        rows.append({'jaw': jaw, 'check': 'surface_witness_direct_triangle_and_distance', 'n': len(ids), 'triangle_reconstruction_error_mm': err, 'distance_identity_error_mm': metricerr, 'within_stated_triangle': within, 'pass': within and err <= 1e-08 and (metricerr <= 1e-09), 'injected_closest_point_wrong_by10_mm_rejected': badmetricerr > 1e-09, 'injected_metric_error_mm': badmetricerr, 'resolution': 'PER_POINT', 'scope': 'checks supplied witnesses are on indicated triangles and distances agree; not independent anatomy or a global nearest-distance proof'})
    dump(ROOT / 'raw/SURFACE_CONTROL_FAULTS.json', {'rows': rows, 'seconds': time.monotonic() - t, 'all_expected_digital_checks_pass': all((x.get('fault_rejected', x.get('pass', False)) for x in rows))})
    print('Surface control fault checks', rows)
if __name__ == '__main__':
    main()
