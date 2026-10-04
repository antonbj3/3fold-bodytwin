"""Reference-independent geometric function. Reconstruction kept as a separate axis."""
from common import *
import importlib.util, collections
from scipy import ndimage
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from skimage.measure import marching_cubes
import trimesh
from contracts import reconstruction_tiers
spec = importlib.util.spec_from_file_location('v4_quality', V4 / 'code/quality.py')
quality = importlib.util.module_from_spec(spec)
spec.loader.exec_module(quality)
sys.path.insert(0, str(V4 / 'payload/participant_code'))
from legacy.geometry import height

def contact_components(xy, faces, gap, band=0.1):
    active = []
    areas = []
    edges = {}
    links = []
    for (fi, ids) in enumerate(faces):
        g = gap[ids]
        if not np.isfinite(g).all():
            continue
        ar = quality.area(quality.band_polygon(g, band), xy[ids])
        if ar <= 1e-12:
            continue
        idx = len(active)
        active.append(fi)
        areas.append(ar)
        for (a, b) in zip(ids, np.roll(ids, -1)):
            (ga, gb) = (gap[a], gap[b])
            if ga == gb:
                good = 0 <= ga <= band
            else:
                (t0, t1) = sorted([-ga / (gb - ga), (band - ga) / (gb - ga)])
                good = min(1, t1) > max(0, t0)
            if good:
                key = tuple(sorted((int(a), int(b))))
                if key in edges:
                    links.append((idx, edges[key]))
                else:
                    edges[key] = idx
    if not active:
        return dict(count=0, component_areas_mm2=[], contact_area_mm2=0.0)
    ij = np.array(links, int).reshape(-1, 2)
    graph = coo_matrix((np.ones(len(ij)), (ij[:, 0], ij[:, 1])), shape=(len(active), len(active))).tocsr()
    (n, lab) = connected_components(graph, directed=False)
    return dict(count=int(n), component_areas_mm2=np.bincount(lab, weights=areas).tolist(), contact_area_mm2=float(sum(areas)))

def functional_map(xy, faces, gap, band=0.1):
    finite = np.isfinite(gap)
    metric = quality.contact_map(xy, faces, gap, gap, band)
    if not metric['contact_support_mm2']:
        return dict(status='UNKNOWN_NO_OPPOSING_SUPPORT', contact_regions=None, contact_area_mm2=None, penetration_area_mm2=None, support_mm2=0.0)
    c = contact_components(xy, faces, gap, band)
    return dict(status='GEOMETRIC_ONLY', contact_regions=c['count'], contact_area_mm2=c['contact_area_mm2'], component_areas_mm2=c['component_areas_mm2'], penetration_area_mm2=metric['negative_gap_area_mm2'], min_gap_mm=float(gap[finite].min()), support_mm2=metric['contact_support_mm2'], resolution='PER_SURFACE_REGION')

def validate_function(r):
    if r['support_mm2'] == 0:
        return r['status'] == 'UNKNOWN_NO_OPPOSING_SUPPORT' and all((r[k] is None for k in ['contact_regions', 'contact_area_mm2', 'penetration_area_mm2']))
    return r['status'] == 'GEOMETRIC_ONLY' and all((np.isfinite(r[k]) and r[k] >= 0 for k in ['contact_regions', 'contact_area_mm2', 'penetration_area_mm2'])) and (r['contact_area_mm2'] <= r['support_mm2'] + 1e-08) and (r['penetration_area_mm2'] <= r['support_mm2'] + 1e-08)

def sample(tri, n=1024):
    if not len(tri):
        return np.empty((0, 3))
    a = np.linalg.norm(np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0]), axis=1) / 2
    ix = np.searchsorted(np.cumsum(a), (np.arange(n) + 0.5) / n * a.sum())
    return tri[ix].mean(1)

def make_mesh(tri):
    return trimesh.Trimesh(tri.reshape(-1, 3), np.arange(len(tri) * 3).reshape(-1, 3), process=False)

def distance(m, p):
    if not len(p):
        return np.array([])
    return np.concatenate([trimesh.proximity.closest_point(m, x)[1] for x in np.array_split(p, max(1, int(np.ceil(len(p) / 128))))])

def grid(tri):
    v = tri.reshape(-1, 3)
    a = [np.arange(v[:, k].min(), v[:, k].max() + 0.125, 0.25) for k in range(2)]
    (xx, yy) = np.meshgrid(*a, indexing='ij')
    ix = np.arange(xx.size).reshape(xx.shape)
    f = []
    for i in range(ix.shape[0] - 1):
        for j in range(ix.shape[1] - 1):
            f.extend([[ix[i, j], ix[i + 1, j], ix[i, j + 1]], [ix[i + 1, j], ix[i + 1, j + 1], ix[i, j + 1]]])
    return (np.c_[xx.ravel(), yy.ravel()], np.array(f))

def extract_external(rec):
    key = rec['key']
    p = npz(DATA / 'external_inputs' / f'{key}.npz')
    pred = npz(DATA / 'external_predictions' / f'{key}.npz')['tsdf']
    scale = float(p['scale'])
    coords = -1 + np.arange(64) * 2 / 64
    world = np.stack(np.meshgrid(coords, coords, coords, indexing='ij'), -1) / scale + p['center']
    crop = np.maximum(np.max(p['roi_low'] - world, axis=-1), np.max(world - p['roi_high'], axis=-1))
    exterior = np.maximum(pred / scale, crop)
    field = np.maximum(exterior, -p['cavity'])
    (v, f, _, _) = marching_cubes(field, 0, spacing=(2 / 64 / scale,) * 3, allow_degenerate=False)
    v += -1 / scale + p['center']
    pts = v[f].mean(1)
    ix = ((pts - p['center']) * scale + 1) * 32
    vals = np.stack([ndimage.map_coordinates(exterior, ix.T, order=1), ndimage.map_coordinates(-p['cavity'], ix.T, order=1)])
    roles = np.argmax(vals, axis=0).astype(np.int8)
    if trimesh.Trimesh(v, f, process=False).volume < 0:
        f = f[:, ::-1]
    dest = DATA / 'external_meshes' / key
    dest.mkdir(exist_ok=True, parents=True)
    np.savez_compressed(dest / 'mesh.npz', vertices=v, faces=f, face_roles=roles)
    trimesh.Trimesh(v @ p['source_R'].T + p['source_base'], f, process=False).export(dest / 'crown.stl')
    return (v, f, roles)

def score_mesh(rec, v, f, roles):
    from contracts import validate_external
    contract = dict(task_id=rec['key'], frame='site_' + rec['key'])
    if not validate_external(contract, dict(contract, units='mm', points=v)):
        raise ValueError('Invalid external native geometry contract')
    key = rec['key']
    a = npz(V4 / 'payload/whole_private' / key / 'reference.npz')
    prep = npz(V4 / 'payload/whole_inputs' / key / 'preparation.npz')
    tri = a['source_triangles']
    tri = tri[tri.mean(1)[:, 2] >= float(prep['margin_z'])]
    ext = v[f[roles == 0]]
    if not len(ext):
        return dict(key=key, status='NO_EXTERIOR', participant=rec['participant'])
    (xy, ff) = grid(tri)
    ceiling = height(a['antagonist_triangles'], xy, True)
    z = height(ext, xy)
    gap = ceiling - z
    function = {str(shift): functional_map(xy, ff, gap + shift) for shift in [-0.05, 0, 0.05]}
    em = make_mesh(ext)
    sm = make_mesh(tri)
    err = max(np.quantile(distance(sm, sample(ext)), 0.95), np.quantile(distance(em, sample(tri)), 0.95))
    proximal = {}
    for side in ['mesial', 'distal']:
        points = sample(a[side + '_triangles'], 256)
        if not len(points):
            proximal[side] = dict(status='UNKNOWN_NO_NEIGHBOUR')
            continue
        d = distance(em, points)
        proximal[side] = dict(unsigned_min_gap_mm=float(d.min()), unsigned_p05_gap_mm=float(np.quantile(d, 0.05)), resolution='PER_SURFACE_REGION', overlap_status='UNKNOWN_UNSIGNED_DISTANCE')
    m = trimesh.Trimesh(v, f, process=False)
    closed = bool(m.is_watertight and m.is_winding_consistent and (m.volume > 0))
    return dict(key=key, participant=rec['participant'], family=rec['family'], source_fdi=rec.get('source_fdi'), status='SCORED', resolution='PER_TOOTH', digital_closed_shell=closed, volume_mm3=float(m.volume), reconstruction_p95_mm=float(err), reconstruction_tiers=reconstruction_tiers(float(err)), function=function, proximal=proximal, virtual_removed_volume_mm3=rec.get('virtual_removed_volume_mm3'), source_virtual_volume_mm3=rec.get('virtual_source_volume_mm3'), substance_sparing='Fixed preparation across participants; no participant-specific sparing claim', continuous_distance_enclosure='MISSING; area-stratified triangle-centroid probes', contact_pose_uncertainty='Scenario +/-0.05mm; not a confidence interval')

def run():
    start = time.perf_counter()
    base = read(V4 / 'results.json')
    out = []
    meta = {r['key']: r for r in read(V4 / 'payload/whole_inputs/RECORDS.json')}
    cache = {}
    for (track, folder) in [('R3', 'whole_predictions'), ('R5', 'prep_predictions')]:
        for old in base['rounds'][track]['rows']:
            if old.get('full_surface_p95_mm') is None:
                out.append(dict(track=track, key=old['key'], participant=old['participant'], status=old['status'], reason=old.get('reason')))
                continue
            key = old['key']
            md = meta[key]
            if key not in cache:
                a = npz(V4 / 'payload/whole_private' / key / 'reference.npz')
                p = npz(V4 / 'payload/whole_inputs' / key / 'preparation.npz')
                tri = a['source_triangles']
                tri = tri[tri.mean(1)[:, 2] >= float(p['margin_z'])]
                (xy, ff) = grid(tri)
                ceiling = height(a['antagonist_triangles'], xy, True)
                cache[key] = (xy, ff, ceiling)
            (xy, ff, ceiling) = cache[key]
            m = npz(V4 / 'payload' / folder / old['participant'] / key / 'mesh.npz')
            ext = m['vertices'][m['faces'][m['face_roles'] == 0]]
            z = height(ext, xy)
            fs = {str(s): functional_map(xy, ff, ceiling - z + s) for s in [-0.05, 0, 0.05]}
            row = dict(track=track, key=key, participant=old['participant'], family=md['family'], source_fdi=md.get('source_fdi'), status='SCORED', resolution='PER_TOOTH', digital_closed_shell=old['digital_complete_shell'], volume_mm3=old['topology']['volume_mm3'], reconstruction_p95_mm=old['full_surface_p95_mm'], reconstruction_tiers=reconstruction_tiers(old['full_surface_p95_mm']), function=fs, proximal=old['proximal'], digital_spacing=old['digital_spacing'], virtual_removed_volume_mm3=md['virtual_removed_volume_mm3'], source_virtual_volume_mm3=md['virtual_source_volume_mm3'], substance_sparing='Fixed virtual preparation, not an individual outcome')
            out.append(row)
        print('functional completed', track, len(out), flush=True)
    external = []
    if (ROOT / 'FROZEN_PREDICTIONS.json').exists():
        frozen = read(ROOT / 'FROZEN_PREDICTIONS.json')
        for r in read(ROOT / 'EXTERNAL_INPUTS.json')['rows']:
            name = r['key'] + '.npz'
            assert sha(DATA / 'external_predictions' / name) == frozen['files'][name]['sha256']
            (v, f, roles) = extract_external(r)
            external.append(dict(track='ToothCraft', **score_mesh(dict(r, participant='ToothCraft_normal'), v, f, roles)))
        out.extend(external)
    dump(ROOT / 'raw/FUNCTIONAL_ROWS.json', out)
    result = dict(rows=len(out), external_rows=external, scored=sum((r['status'] == 'SCORED' for r in out)), contact_scope='Static geometric contact connected components; not loaded physical contact points', seconds=time.perf_counter() - start)
    dump(ROOT / 'raw/R4_FUNCTION.json', result)
    return out
if __name__ == '__main__':
    run()
