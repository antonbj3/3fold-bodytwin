"""Inject faults into actual geometry, then invoke the actual C/D scorer."""
from joint_solid import *
import score_ray as sc
import score_joint as sj

def native_identity(m, source):
    canonical = lambda tt: sorted((tuple(sorted((tuple(x) for x in t))) for t in tt))
    return canonical(m['vertices'][m['faces'][m['roles'] == 0]]) == canonical(source)

def main():
    latest = read(R / 'raw/DEMO_LATEST.json')
    work = Path(latest['work_directory'])
    sc.R = work
    sc.D = work / 'data'
    sj.D = work / 'data'
    rows = read(work / 'RESULTS_D.json')['rows']
    controls = []
    for row in rows:
        m = load(row['mesh_path'])
        source = m['native_triangles']
        good = native_identity(m, source)
        bad = {k: v.copy() for (k, v) in m.items()}
        idx = bad['faces'][bad['roles'] == 0][0, 0]
        bad['vertices'][idx, 0] += 0.1
        moved_rejected = not native_identity(bad, source)
        bad = {k: v.copy() for (k, v) in m.items()}
        bad['prep_vertices'][:, 2] += 0.2
        path = work / 'data' / (row['key'] + '_fault_gap.npz')
        np.savez_compressed(path, **bad)
        rec = {k: v for (k, v) in row.items() if k in read(R / 'FROZEN_PREDICTIONS_D.json')['rows'][0]}
        rec.update(mesh_path=str(path), mesh_sha256=sha(path))
        score = sc.evaluate_c(rec)
        gap_rejected = not score['gates']['gap']
        bad = {k: v.copy() for (k, v) in m.items()}
        bad['faces'] = bad['faces'][1:]
        bad['roles'] = bad['roles'][1:]
        path = work / 'data' / (row['key'] + '_fault_open.npz')
        np.savez_compressed(path, **bad)
        rec.update(mesh_path=str(path), mesh_sha256=sha(path))
        score2 = sc.evaluate_c(rec)
        record = dict(key=row['key'], native_identity_independent=good, moved_native_facet_rejected=moved_rejected, shifted_prep_rejected=gap_rejected, shifted_prep_max_gap_error_mm=score['film']['max_residual_mm'], deleted_facet_rejected=not score2['gates']['closed'])
        assert all((record[k] for k in ['native_identity_independent', 'moved_native_facet_rejected', 'shifted_prep_rejected', 'deleted_facet_rejected']))
        controls.append(record)
    dump(R / 'raw/GEOMETRY_FAULT_CONTROLS.json', dict(all_pass=True, rows=controls, scope='Independent native-triangle identity plus injections through actual film and closed-mesh scorer; no physical validation'))
    print('All actual-geometry fault controls rejected the injected faults', flush=True)
if __name__ == '__main__':
    main()
