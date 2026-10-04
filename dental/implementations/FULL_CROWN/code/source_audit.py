from dental_release.paths import expand as _release_expand
from fc_common import *
import trimesh
sys.path.insert(0, _release_expand('@DENTAL_WORK_ROOT@/PROOF_LANE_GENCAD_V5/deps'))
from prescan_local_holes import main_component, cap
F = functional()

def run():
    rows = []
    for rec in read(DATA / 'public/RECORDS.json'):
        key = rec['key']
        p = npz(V4 / 'payload/whole_inputs' / key / 'preparation.npz')
        a = npz(V4 / 'payload/whole_private' / key / 'reference.npz')
        tri = a['source_triangles']
        m = trimesh.Trimesh(tri.reshape(-1, 3), np.arange(len(tri) * 3).reshape(-1, 3), process=True)
        dup = int(len(m.faces) - m.unique_faces().sum())
        deg = int(len(m.faces) - m.nondegenerate_faces().sum())
        ed = np.sort(m.edges, axis=1)
        (_, ct) = np.unique(ed, axis=0, return_counts=True)
        parts = m.split(only_watertight=False)
        areas = [q.area for q in parts]
        row = dict(key=key, family=rec['family'], resolution='PER_TOOTH', faces=len(tri), duplicate_faces=dup, degenerate_faces=deg, boundary_edges=int((ct == 1).sum()), nonmanifold_edges=int((ct > 2).sum()), face_components=len(parts), largest_component_area_fraction=float(max(areas) / sum(areas)) if areas else None, reference_is_closed=bool(m.is_watertight), source_fdi_status='INFERRED_X11_NOT_INDEPENDENTLY_ANNOTATED')
        try:
            (main, _, _) = main_component(tri)
            (cm, ci) = cap(main)
            row['closure_loops'] = ci['loops']
            ctri = cm.triangles[len(main.faces):]
            ref = tri[tri.mean(1)[:, 2] >= p['margin_z']]
            for (tag, participant) in [('R6', 'prescan_signed_solid_die_shell'), ('R7', 'prescan_local_holes_die_shell')]:
                fp = DATA / (tag + '_predictions') / participant / key / 'mesh.npz'
                if not fp.exists():
                    continue
                pred = npz(fp)
                pts = F.sample(pred['vertices'][pred['faces'][pred['face_roles'] == 0]])
                ds = F.distance(F.make_mesh(ref), pts)
                dc = F.distance(F.make_mesh(ctri), pts) if len(ctri) else np.full(len(pts), np.inf)
                tail = ds >= np.quantile(ds, 0.95)
                row[tag + '_generated_tail_closer_to_local_closure_fraction'] = float(np.mean(dc[tail] < ds[tail]))
                row[tag + '_all_generated_closer_to_local_closure_fraction'] = float(np.mean(dc < ds))
        except Exception as e:
            row['closure_audit_error'] = str(e)
        rows.append(row)
    dump(ROOT / 'raw/SOURCE_GEOMETRY_AUDIT.json', dict(rows=rows, claim_type='capability', interpretation='Open/nonmanifold inferred label regions are not closed anatomical crown truth. This audit does not prove scanner error or clinical invalidity; closure-facet attribution is a geometric diagnostic.', external_referent=REFERENT))
if __name__ == '__main__':
    run()
