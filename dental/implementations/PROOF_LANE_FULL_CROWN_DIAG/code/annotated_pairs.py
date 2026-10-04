from common import *
from local_pairs import component, fit, sample, metrics, distance
import resource

def read_obj(p):
    vv = []
    ff = []
    for line in Path(p).open():
        if line.startswith('v '):
            vv.append([float(x) for x in line.split()[1:4]])
        elif line.startswith('f '):
            q = [int(x.split('/')[0]) - 1 for x in line.split()[1:]]
            if min(q) < 0:
                raise ValueError('negative OBJ indices')
            for j in range(1, len(q) - 1):
                ff.append([q[0], q[j], q[j + 1]])
    return (np.array(vv), np.array(ff, int))

def run():
    pr = read(ROOT / 'PREREG_ANNOTATED.json')
    rows = []
    screen = []
    sources = []
    start = time.perf_counter()
    cases = 0
    for pp in pr['selection']['paths_in_order']:
        p = Path(pp)
        labels = np.array(read(p.with_suffix('.json'))['labels'])
        fd = pr['selection']['teeth']
        if any((np.sum(labels == x) < 100 for x in fd)):
            screen.append(dict(path=pp, reason='Missing100 vertices of required tooth'))
            continue
        (v, f) = read_obj(p)
        assert len(v) == len(labels)
        if any((np.sum(np.all(labels[f] == x, axis=1)) < 30 for x in fd)):
            screen.append(dict(path=pp, reason='Missing30 wholly labelled triangles'))
            continue
        z = np.linalg.eigh(np.cov(v.T))[1][:, 0]
        if np.any(labels == 0) and np.dot(z, v[labels > 0].mean(0) - v[labels == 0].mean(0)) < 0:
            z = -z
        x = np.array([1.0, 0, 0])
        x -= np.dot(x, z) * z
        x /= np.linalg.norm(x)
        y = np.cross(z, x)
        v = v @ np.array([x, y, z]).T
        case = hashlib.sha256(p.parent.name.encode()).hexdigest()[:16]
        cases += 1
        sources.append(dict(path=pp, obj_sha256=sha(p), labels_path=str(p.with_suffix('.json')), labels_sha256=sha(p.with_suffix('.json')), case_key=case))
        for (fdi, fam) in [(11, 'anterior_crown'), (14, 'premolar_crown'), (16, 'molar_crown')]:
            (t, dt, nt) = component(v[f[np.all(labels[f] == fdi, axis=1)]])
            (s, ds, ns) = component(v[f[np.all(labels[f] == fdi + 10, axis=1)]])
            s[:, :, 0] *= -1
            (R, u, starts) = fit(sample(s, 1536), sample(t, 1536))
            s = s @ R.T + u
            mm = {str(n): metrics(s, t, n) for n in [2048, 8192]}
            selfmax = float(distance(t, sample(t, 2048)).max())
            valid = ds <= 0.05 and dt <= 0.05
            row = dict(case_key=case, key=case + '_' + fam, family=fam, source_fdi=fdi, resolution='PER_TOOTH', metrics=mm, registration_start_rms_mm=starts, homolog_dropped_area_fraction=ds, target_dropped_area_fraction=dt, components=[ns, nt], valid_for_envelope=valid, rejection_reason=None if valid else 'component area loss >5%', self_max_mm=selfmax, identity_control_pass=selfmax <= 1e-10, probe_sensitivity_p95_mm=abs(mm['2048']['p95_mm'] - mm['8192']['p95_mm']), clinical_acceptance='UNKNOWN')
            rows.append(row)
            np.savez_compressed(DATA / (row['key'] + '_annotated.npz'), homolog=s, target=t)
            print(row['key'], mm['8192'], valid, flush=True)
            dump(ROOT / 'raw/ANNOTATED_PAIRS.json', dict(rows=rows, sources=sources, screening_rejections=screen, seconds=time.perf_counter() - start, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024))
            state('ANNOTATED_MEASURING', str(len(rows)) + '/18 pairs', 'Freeze independent descriptive thresholds before R3 scoring')
        if cases == 6:
            break
    assert cases == 6, 'Insufficient eligible arches'
if __name__ == '__main__':
    run()
