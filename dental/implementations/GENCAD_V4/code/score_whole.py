from common import *
sys.path.insert(0, str(DATA / 'participant_code'))
from legacy.geometry import height
from quality import contact_map
from score_panel import verified_npz
import trimesh
from skimage.measure import marching_cubes
from scipy.spatial import cKDTree

def samples(tri, n=1024):
    if not len(tri):
        return np.empty((0, 3))
    ar = np.linalg.norm(np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0]), axis=1) / 2
    if ar.sum() <= 0:
        raise ValueError('zero area')
    i = np.searchsorted(np.cumsum(ar), (np.arange(n) + 0.5) / n * ar.sum())
    return tri[i].mean(1)

def mesh(tri):
    return trimesh.Trimesh(tri.reshape(-1, 3), np.arange(len(tri) * 3).reshape(-1, 3), process=False)

def distances(m, p):
    if not len(p):
        return np.empty(0)
    return np.concatenate([trimesh.proximity.closest_point(m, pp)[1] for pp in np.array_split(p, max(1, int(np.ceil(len(p) / 128))))])

def topology(v, f):
    edges = np.sort(np.concatenate([f[:, [0, 1]], f[:, [1, 2]], f[:, [2, 0]]]), axis=1)
    ct = np.unique(edges, axis=0, return_counts=True)[1]
    m = trimesh.Trimesh(v, f, process=False)
    vol = float(m.volume)
    return dict(bad_edges=int(np.sum(ct != 2)), watertight=bool(m.is_watertight), winding_consistent=bool(m.is_winding_consistent), volume_mm3=vol, pass_gate=bool(np.all(ct == 2) and m.is_winding_consistent and (vol > 0)))

def run(track='R3'):
    from integrity import verify
    verify()
    predroot = DATA / ('prep_predictions' if track == 'R5' else 'whole_predictions')
    start = time.perf_counter()
    pr = read(ROOT / ('PREREG_R5.json' if track == 'R5' else 'PREREG_R3.json'))
    wf = read(ROOT / 'FROZEN_WHOLE_INPUTS.json')
    pf = read(ROOT / ('FROZEN_PREP_PREDICTIONS.json' if track == 'R5' else 'FROZEN_WHOLE_PREDICTIONS.json'))
    records = read(predroot / 'RECORDS.json')
    meta = {r['key']: r for r in wf['records']}
    rows = []
    control = []
    cached = {}
    for (idx, rec) in enumerate(records):
        row = dict(rec)
        key = rec['key']
        md = meta[key]
        row.update(split=md['split'], dataset=md['dataset'], resolution='PER_TOOTH')
        if rec['status'] != 'EXPORTED':
            rows.append(row)
            continue
        if key not in cached:
            a = verified_npz(DATA / 'whole_private' / key / 'reference.npz', wf['private_files'][key + '/reference.npz'])
            p = verified_npz(DATA / 'whole_inputs' / key / 'preparation.npz', wf['public_files'][key + '/preparation.npz'])
            tri = a['source_triangles']
            tri = tri[tri.mean(1)[:, 2] >= float(p['margin_z'])]
            if len(tri) == 0:
                raise ValueError('No measured source above margin')
            sm = mesh(tri)
            sp = samples(tri)
            s0 = distances(sm, sp)
            shift = distances(sm, sp + [0, 0, 2.0])
            control.append(dict(key=key, unshifted_p95_mm=float(np.quantile(s0, 0.95)), shifted_p95_mm=float(np.quantile(shift, 0.95)), source_shift_rejected=float(np.quantile(shift, 0.95)) > 0.35, original_accepted=float(np.max(s0)) < 1e-07))
            neighbors = {}
            for side in ['mesial', 'distal']:
                pp = samples(a[side + '_triangles'], 256)
                d = distances(sm, pp)
                if len(pp):
                    order = np.argsort(d, kind='stable')[:64]
                    neighbors[side] = (pp[order], d[order])
            h = float(p['step'])
            (pv, pf0, _, _) = marching_cubes(p['preparation'], 0.001, spacing=(h,) * 3, allow_degenerate=False)
            pv += p['origin']
            pm = trimesh.Trimesh(pv, pf0, process=False)
            v = tri.reshape(-1, 3)
            lo = v[:, :2].min(0)
            hi = v[:, :2].max(0)
            axes = [np.arange(lo[k], hi[k] + 0.125, 0.25) for k in range(2)]
            (xx, yy) = np.meshgrid(*axes, indexing='ij')
            xy = np.c_[xx.ravel(), yy.ravel()]
            ii = np.arange(len(xy)).reshape(xx.shape)
            ff = []
            for i in range(ii.shape[0] - 1):
                for j in range(ii.shape[1] - 1):
                    ff.extend([[ii[i, j], ii[i + 1, j], ii[i, j + 1]], [ii[i + 1, j], ii[i + 1, j + 1], ii[i, j + 1]]])
            ref = height(tri, xy)
            ceiling = height(a['antagonist_triangles'], xy, True)
            cached[key] = (a, p, sm, sp, neighbors, pm, xy, np.array(ff), ref, ceiling)
        (a, p, sm, sp, neighbors, pm, xy, ff, ref, ceiling) = cached[key]
        rel = rec['participant'] + '/' + key + '/mesh.npz'
        m = verified_npz(predroot / rel, pf['files'][rel])
        (v, f, roles) = (m['vertices'], m['faces'], m['face_roles'])
        tp = topology(v, f)
        ext = v[f[roles == 0]]
        if not len(ext):
            row.update(status='FAILED_NO_EXTERIOR')
            rows.append(row)
            continue
        ep = samples(ext)
        em = mesh(ext)
        forward = distances(sm, ep)
        reverse = distances(em, sp)
        err = float(max(np.quantile(forward, 0.95), np.quantile(reverse, 0.95)))
        proximal = {}
        raw = {}
        for (side, (pp, dd)) in neighbors.items():
            pred = distances(em, pp)
            proximal[side] = dict(mean_abs_gap_error_mm=float(np.mean(abs(pred - dd))), source_min_gap_mm=float(dd.min()), candidate_min_gap_mm=float(pred.min()), probes=len(pp), resolution='PER_SURFACE_REGION', scope='Unsigned mesh distance; overlap/contact force not identified')
            raw[side + '_probes'] = pp
            raw[side + '_source_gap'] = dd
            raw[side + '_candidate_gap'] = pred
        z = height(ext, xy)
        cm = contact_map(xy, ff, ceiling - z, ceiling - ref)
        fi = f[roles == 1]
        it = v[fi]
        ips = samples(it, 512)
        gap = distances(pm, ips)
        occ = ips[:, 2] >= np.quantile(ips[:, 2], 0.75)
        regional = {}
        for (name, mask0) in [('occlusal_height_quartile', occ), ('axial_and_margin_remainder', ~occ)]:
            gg = gap[mask0]
            regional[name] = dict(mean_mm=float(gg.mean()), p05_mm=float(np.quantile(gg, 0.05)), p95_mm=float(np.quantile(gg, 0.95)), n=len(gg), resolution='PER_SURFACE_REGION') if len(gg) else None
        trace = DATA / ('prep_scores' if track == 'R5' else 'whole_scores') / rec['participant']
        trace.mkdir(exist_ok=True, parents=True)
        np.savez_compressed(trace / (key + '.npz'), exterior_samples=ep, exterior_to_source_mm=forward, source_samples=sp, source_to_exterior_mm=reverse, intaglio_samples=ips, digital_gap_mm=gap, **raw)
        row.update(topology=tp, full_surface_p95_mm=err, anatomy_gate=err <= 0.35, proximal=proximal, contact=cm, digital_spacing=regional, digital_spacing_scope='Quarter-mm voxel-derived nominal CAD geometry, not calibrated physical cement film; regions geometrically defined and do not match published measurement locations', virtual_removed_volume_mm3=md['virtual_removed_volume_mm3'], cement_physical='UNKNOWN_PROTOCOL_AND_SPATIAL_FIELD', fracture_margin='UNKNOWN_SETUP_MISMATCH', complete_clinical_L1='UNKNOWN_MARGIN_INSERTION_COLLISION_MILLING', digital_complete_shell=tp['pass_gate'], sample_bound='UNKNOWN continuous coverage bound; p95 refers to area-stratified discrete probes')
        if idx < 8:
            control.append(dict(key=key, participant=rec['participant'], removed_face_rejected=not topology(v, f[:-1])['pass_gate']))
        rows.append(row)
        if idx % 16 == 0:
            print('whole scored', idx + 1, 'seconds', round(time.perf_counter() - start), flush=True)
    exports = [r for r in rows if r.get('full_surface_p95_mm') is not None]
    good = sum((r['digital_complete_shell'] for r in exports))
    controls_ok = all((all((v for (k, v) in r.items() if k.endswith('_rejected') or k.endswith('_accepted'))) for r in control))
    out = dict(round=track, claim_type='capability', decision='COMPLETE_VIRTUAL_SHELLS_SCORED' if good and controls_ok else 'WHOLE_SHELL_GATE_FAILED', requested=len(records), exported=len(exports), closed_shells=good, anatomy_pass=sum((r['anatomy_gate'] for r in exports)), rows=rows, controls=control, dropout_fraction=(len(records) - len(exports)) / len(records), external_referent=pr['external_referent'], seconds=time.perf_counter() - start, physical_cement_calibrated=False, fracture_margin_known=False, source_scope='Original registered full tooth with X11 inferred FDI. Virtual margin and source insertion envelope are closures; all generated axial faces included.')
    dump(ROOT / 'rounds' / f'{track}.json', out)
    state(track + '_DECIDED', out['decision'], 'Validate protocol-conditioned literature ports, integrity and copied replay')
    print(json.dumps({k: out[k] for k in ['decision', 'requested', 'exported', 'closed_shells', 'anatomy_pass']}))
if __name__ == '__main__':
    run('R5' if '--prep' in sys.argv else 'R3')
