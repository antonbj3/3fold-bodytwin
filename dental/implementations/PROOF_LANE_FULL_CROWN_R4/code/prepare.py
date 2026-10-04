from geometry import *

def run():
    st = time.perf_counter()
    pr = read(ROOT / 'PREREG_A.json')
    rejected = []
    selected = []
    rows = []
    for pp in pr['selection']['ordered_paths']:
        p = Path(pp)
        lab = np.array(read(p.with_suffix('.json'))['labels'])
        fd = [11, 14, 16, 21, 24, 26]
        if any((np.sum(lab == x) < 100 for x in fd)):
            rejected.append(dict(path=pp, reason='required labels <100'))
            continue
        lower = list(p.parents[3].glob('data_part_*/lower/' + p.parent.name + '/*_lower.obj'))
        if len(lower) != 1:
            rejected.append(dict(path=pp, reason='lower absent/ambiguous'))
            continue
        (v, f) = read_obj(p)
        if any((np.sum(np.all(lab[f] == x, axis=1)) < 30 for x in fd)):
            rejected.append(dict(path=pp, reason='required triangles <30'))
            continue
        (lv, lf) = read_obj(lower[0])
        z = np.linalg.eigh(np.cov(v.T))[1][:, 0]
        if np.dot(z, v[lab > 0].mean(0) - v[lab == 0].mean(0)) < 0:
            z = -z
        x = np.array([1.0, 0.0, 0.0])
        x -= x @ z * z
        x /= np.linalg.norm(x)
        frame = np.array([x, np.cross(z, x), z])
        v = v @ frame.T
        lv = lv @ frame.T
        key = hashlib.sha256(p.parent.name.encode()).hexdigest()[:16]
        source = dict(case_key=key, upper_path=str(p), upper_sha256=sha(p), label_path=str(p.with_suffix('.json')), label_sha256=sha(p.with_suffix('.json')), lower_path=str(lower[0]), lower_sha256=sha(lower[0]), frame=frame, bite_registration='UNKNOWN; common raw coordinates only', arch_z_ranges_mm=[v[:, 2].min(), v[:, 2].max(), lv[:, 2].min(), lv[:, 2].max()])
        selected.append(source)
        for (fdi, family) in [(11, 'anterior'), (14, 'premolar'), (16, 'molar')]:
            rec = dict(case_key=key, key=key + '_' + family, family=family, fdi=fdi, status='PREPARED')
            try:
                (target, drop, components) = component(v[f[np.all(lab[f] == fdi, axis=1)]])
                (donor, dd, nc) = component(v[f[np.all(lab[f] == fdi + 10, axis=1)]])
                if max(drop, dd) > 0.05:
                    raise ValueError('component loss exceeds5%')
                margin = float(np.quantile(target.vertices[:, 2], 0.2))
                visible_hi = float(np.quantile(target.vertices[:, 2], 0.25))
                ts = cut(target, margin)
                loop = boundary(ts, margin)
                curve = ts.vertices[loop].copy()
                poly = Polygon(curve[:, :2])
                if not poly.is_valid:
                    raise ValueError('invalid marginal polygon')
                c = np.array(poly.centroid.coords[0])
                radius = poly.boundary.distance(Point(c)) - 0.8
                top = float(target.vertices[:, 2].max() - 1.5)
                if radius <= 0.7 or top - margin < 1:
                    raise ValueError('virtual die dimensions invalid')
                collar = target.triangles[target.triangles[:, :, 2].max(1) <= visible_hi]
                if len(collar) < 30:
                    raise ValueError('visible collar fewer30 triangles')
                sides = []
                for side_id in [21, 12] if fdi == 11 else [fdi - 1, fdi + 1]:
                    sides.append(v[f[np.all(lab[f] == side_id, axis=1)]])
                lo = target.bounds[0] - 3
                hi = target.bounds[1] + 3
                lt = lv[lf]
                keep = np.all(lt[:, :, :2].max(1) >= lo[:2], axis=1) & np.all(lt[:, :, :2].min(1) <= hi[:2], axis=1)
                origin = target.vertices.mean(0)
                origin[2] = margin
                public = dict(donor=donor.triangles - origin, collar=collar - origin, margin_curve=curve - origin, margin_z=np.array(0.0), cavity_center=c - origin[:2], cavity_radius=np.array(radius), cavity_top=np.array(top - margin), cavity_taper=np.array(np.tan(np.deg2rad(6.0))), mesial=sides[0] - origin, distal=sides[1] - origin, antagonist=lt[keep] - origin, origin_world=origin)
                for (folder, arr) in [('public', public), ('private', dict(target=target.triangles - origin, reference=ts.triangles - origin))]:
                    (DATA / folder).mkdir(exist_ok=True)
                    np.savez_compressed(DATA / folder / (rec['key'] + '.npz'), **arr)
                rec.update(target_component_drop=drop, donor_component_drop=dd, target_components=components, donor_components=nc, public_path=DATA / 'public' / (rec['key'] + '.npz'), private_path=DATA / 'private' / (rec['key'] + '.npz'), source_faces=len(target.faces), donor_faces=len(donor.faces), visible_collar_faces=len(collar), antagonist_faces=int(keep.sum()), reference_faces=len(ts.faces))
            except Exception as e:
                rec.update(status='PREPARATION_REJECTED', reason=repr(e))
            rows.append(rec)
            print(rec['key'], rec['status'], rec.get('reason', ''), flush=True)
        if len(selected) == 6:
            break
    assert len(selected) == 6
    lock = dict(sources=selected, records=rows, eligibility_rejections=rejected, seconds=time.perf_counter() - st, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024)
    for r in rows:
        if r['status'] == 'PREPARED':
            r.update(public_sha256=sha(r['public_path']), private_sha256=sha(r['private_path']))
    freeze(ROOT / 'FROZEN_TEST_SET.json', lock)
    state('DATA_PREPARED', dict(requested=18, prepared=sum((r['status'] == 'PREPARED' for r in rows))), 'Generate collar-only and local-boundary variants before held-surface scoring')
if __name__ == '__main__':
    run()
