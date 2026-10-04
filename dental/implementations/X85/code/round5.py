from common import *
from geometry import export_stl
import time, resource

def full_crown_force_query(crown, calibration=None):
    if calibration is None:
        return dict(status='UNKNOWN', reason='Same-specimen region-resolved actuation response and force observation absent', force_interval_N=None)
    if calibration.get('geometry_sha256') != crown['mesh_sha256'] or calibration.get('frame') != crown['frame']:
        raise ValueError('Full crown calibration geometry/frame mismatch')
    return dict(status='UNKNOWN', reason='Whole-crown actuation/basis/observation validation not yet implemented; no scalar-to-patch force transfer', force_interval_N=None)

def inherited_geometry():
    saved = sys.modules.get('common')
    paths = list(sys.path)
    try:
        vc = module('x85_v6_common_for_geometry', V6 / 'code/common.py')
        sys.modules['common'] = vc
        sys.path.insert(0, str(V6 / 'code'))
        return module('x85_v6_reuse_readonly', V6 / 'code/reuse.py')
    finally:
        if saved is None:
            sys.modules.pop('common', None)
        else:
            sys.modules['common'] = saved
        sys.path[:] = paths

def run():
    for (p, d) in read(ROOT / 'INPUT_LOCK_R5.json')['files'].items():
        if sha(p) != d['sha256']:
            raise ValueError('Full source input drift ' + p)
    start = time.perf_counter()
    pr = read(ROOT / 'PREREG_R5.json')
    FC = ROOT.parent / 'PROOF_LANE_FULL_CROWN_R2'
    FCDATA = DATA.parent / 'PROOF_LANE_FULL_CROWN_R2'
    refined = {(r['key'], r['participant']): r for r in read(FC / 'raw/WALL_REFINEMENT.json')['rows'] if r['round'] == 'R2'}
    rows = []
    pred = []
    cache = {}
    for r0 in pr['cohort']:
        r = dict(r0)
        path = FCDATA / 'R2' / r['participant'] / r['key'] / 'mesh.npz'
        with np.load(path, allow_pickle=False) as a:
            m = dict(a)
        v = m['vertices'].copy()
        f = m['faces']
        roles = m['face_roles']
        protected = np.unique(f[roles != 0])
        ext = np.unique(f[roles == 0])
        free = np.setdiff1d(ext, protected)
        with np.load(V4 / 'payload/whole_inputs' / r['key'] / 'preparation.npz', allow_pickle=False) as a:
            margin = float(a['margin_z'])
        (xmin, xmax) = (v[ext, 0].min(), v[ext, 0].max())
        b = (xmax - v[free, 0]) / (xmax - xmin)
        taper = np.clip((v[free, 2] - margin) / 0.6, 0, 1)
        delta = np.zeros(len(v))
        delta[free] = -0.02 * b * taper
        vnew = v.copy()
        vnew[:, 2] += delta
        dest = ROOT / 'exports' / ('FULL_' + r['key'])
        dest.mkdir(parents=True, exist_ok=True)
        for (name, vv) in [('unadjusted', v), ('regional20um_probe', vnew)]:
            export_stl(dest / (name + '.stl'), vv, f)
        artifact = DATA / ('FULL_' + r['key'] + '.npz')
        np.savez_compressed(artifact, vertices=vnew, faces=f, face_roles=roles, height_change_mm=delta)
        w = refined.get((r['key'], r['participant']))
        oldlower = None
        if w and w['mesh_sha256'] == sha(path):
            oldlower = w['refined']['lower_mm']
        displacement = float(np.max(abs(delta)))
        r.update(old_mesh_sha256=sha(path), artifact_path=str(artifact), artifact_sha256=sha(artifact), protected_vertex_identity_error=float(np.max(abs(vnew[protected] - v[protected]))), maximum_actual_vertex_displacement_um=displacement * 1000, all_point_displacement_bound_um=displacement * 1000, all_exterior_wall_lower_mm=None if oldlower is None else oldlower - displacement, wall_bound_kind='Previous same-hash all-exterior bound minus analytic barycentric motion bound; annulus excluded; rounding MISSING', force_query=full_crown_force_query(dict(mesh_sha256=sha(path), frame=r['key'])), milling='UNKNOWN_FULL_TOOL_HOLDER_FIXTURE', chairside_adjustment_saved_um=None)
        try:
            full_crown_force_query(dict(mesh_sha256=sha(path), frame=r['key']), dict(geometry_sha256='X82_OR_ROOF', frame='other'))
            r['wrong_force_binding_rejected'] = False
        except ValueError:
            r['wrong_force_binding_rejected'] = True
        rows.append(r)
        pred.append(dict(key=r['key'], artifact_sha256=r['artifact_sha256'], maximum_height_probe_mm=0.02, force_prediction='UNKNOWN', chairside_saved_adjustment_um=None))
        cache[r['key']] = (v, vnew, f, roles, margin)
    if not (ROOT / 'FROZEN_PREDICTIONS_R5.json').exists():
        freeze(ROOT / 'FROZEN_PREDICTIONS_R5.json', dict(predictions=pred, private_contact_queries_before_freeze=0, physical_measurement='NOT_RUN'))
    from review_integrity import assert_frozen_payload
    assert_frozen_payload(ROOT, 'FROZEN_PREDICTIONS_R5.json', 'predictions', clean(pred))
    q = v6_contact()
    reuse = inherited_geometry()
    VF = reuse.VF
    controls = []
    for r in rows:
        (v, vnew, f, roles, margin) = cache[r['key']]
        with np.load(V4 / 'payload/whole_private' / r['key'] / 'reference.npz', allow_pickle=False) as a:
            source = dict(a)
        tri = source['source_triangles']
        tri = tri[tri.mean(1)[:, 2] >= margin]
        (xy, ff) = reuse.grid(tri)
        ceiling = VF.height(source['antagonist_triangles'], xy, True)
        rg = ceiling - VF.height(tri, xy)
        oldg = ceiling - VF.height(v[f[roles == 0]], xy)
        newg = ceiling - VF.height(vnew[f[roles == 0]], xy)
        before = q.compare(xy, ff, oldg, rg)
        after = q.compare(xy, ff, newg, rg)
        change = q.compare(xy, ff, newg, oldg)
        parity = q.Q.contact_map(xy, ff, newg, oldg)
        error = abs(parity['contact_symdiff_mm2'] - change['symdiff_mm2']) if change['status'] == 'SCORED' else None
        r.update(contact_before=before, contact_after=after, contact_change=change, continuous_wall_gate='PASS' if r['all_exterior_wall_lower_mm'] is not None and r['all_exterior_wall_lower_mm'] >= 0.5 else 'UNKNOWN_OR_FAIL')
        controls.append(dict(key=r['key'], contact_polygon_parity_error_mm2=error, fault_plus3mm2_rejected=error is not None and abs(parity['contact_symdiff_mm2'] - (change['symdiff_mm2'] + 3)) > 1e-07, protected_vertices_identical=r['protected_vertex_identity_error'] == 0, force_binding_fault_rejected=r['wrong_force_binding_rejected']))
        print(r['key'], 'contact delta', change.get('symdiff_mm2'), 'wall', r['continuous_wall_gate'], flush=True)
    summary = dict(requested=len(rows), scored=len(rows), case_clusters=len({r['case_key'] for r in rows}), continuous_wall_passes=sum((r['continuous_wall_gate'] == 'PASS' for r in rows)), nonzero_contact_mask_changes=sum((r['contact_change'].get('symdiff_mm2', 0) > 1e-09 for r in rows)), force_transfer_rejections=sum((r['force_query']['status'] == 'UNKNOWN' for r in rows)), physical_validated_crowns=0, seconds=time.perf_counter() - start, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024)
    dump(ROOT / 'raw/R5_ROWS.json', rows)
    dump(ROOT / 'raw/R5_CONTROLS.json', controls)
    dump(ROOT / 'raw/R5_SUMMARY.json', summary)
    state('R5_DECIDED', str(summary), 'Seal one-command replay, external measurement transcription and10-crown lab protocol')
    print(json.dumps(summary))
if __name__ == '__main__':
    run()
