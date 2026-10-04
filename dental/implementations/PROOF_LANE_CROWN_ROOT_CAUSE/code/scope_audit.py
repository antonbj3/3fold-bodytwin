"""Descriptive scope audit; does not select candidates or change any gate."""
from joint_solid import *

def area(t):
    return float(np.linalg.norm(np.cross(t[:, 1] - t[:, 0], t[:, 2] - t[:, 0]), axis=1).sum() / 2)

def run():
    assert (R / 'PREREG_SCOPE2.json').exists()
    rows = []
    for rec in read(R / 'RESULTS_D.json')['rows']:
        if not rec.get('geometric_conjunction'):
            continue
        st = time.perf_counter()
        m = load(rec['mesh_path'])
        t = m['vertices'][m['faces']]
        roles = m['roles']
        key = rec['key']
        (an, aa) = (area(t[roles == 0]), area(t[roles == 2]))
        support = D / (key + '_scope_T.mesh')
        prep = D / (key + '_scope_Q.mesh')
        out = D / (key + '_scope_TminusQ.mesh')
        meshwrite(support, m['support_vertices'], m['support_faces'])
        meshwrite(prep, m['prep_vertices'], m['prep_faces'])
        p = subprocess.run([str(D / 'boolean'), str(support), str(prep), str(out)], capture_output=True, text=True, check=True)
        bo = json.loads(p.stdout)
        (v, f) = meshread(out)
        minus = trimesh.Trimesh(v, f, process=False)
        assert bo['success'] and bo['closed'] and (bo['exact_intersection_count'] == 0) and minus.is_watertight
        vt = rec['closure']['volume_mm3']
        vq = vt - minus.volume
        vk = rec['retained_digital_volume_mm3']
        assert 0 < vq < vk < vt
        row = dict(key=key, family=rec['family'], native_area_mm2=an, artificial_exterior_area_mm2=aa, artificial_exterior_fraction=aa / (aa + an), support_volume_mm3=vt, occupied_cavity_volume_mm3=vk, virtual_preparation_volume_mm3=float(vq), virtual_preparation_fraction=float(vq / vt), nominal_gap_volume_mm3=float(vk - vq), boolean=bo, seconds=time.perf_counter() - st, resolution='PER_SURFACE_REGION aggregated PER_TOOTH', scope='Artificial support T, never measured anatomical volume')
        rows.append(row)
        dump(R / 'RESULTS_SCOPE.json', dict(claim_type='capability', rows=rows, correction='C/D retained_digital_volume_mm3 denotes occupied cavity, including film; use virtual_preparation_volume_mm3 for Q intersect T.', external_referent=read(R / 'PREREG_SCOPE2.json')['external_referent']))
        print(key, 'artificial fraction', row['artificial_exterior_fraction'], 'Q/T', row['virtual_preparation_fraction'], flush=True)
    return rows
if __name__ == '__main__':
    run()
