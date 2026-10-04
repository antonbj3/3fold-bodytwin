from common import *
import score
from loop_crown import clip

def assess(row):
    r = score.assess(row)
    m = load(row['mesh_path'])
    (v, f, roles) = (m['vertices'], m['faces'], m['roles'])
    M = m['M']
    b = float(m['base'])
    rec = next((x for x in inputs() if x['key'] == row['key']))
    S = m['target']
    (sv, sf) = triangles_mesh(S)
    (lowv, lowf) = clip(sv, sf, b, M, above=False)
    visible_low = lowv[lowf]
    source_preserved = row['source_closure'].get('method') != 'winding_isosurface'
    assembly = np.r_[v[f[roles == 0]], visible_low]
    ashape = shape_score(assembly, S)
    rim = v[m['margin_vertex_ids']]
    sd = Distance(sv, sf)
    d = sd.query(rim)[0]
    plane = abs((rim @ M.T)[:, 2] - b)
    shoulder = m['shoulder_gap_mm']
    if 'prep_cap_vertices' in m:
        sid = m['shoulder_vertex_ids']
        shoulder = Distance(m['prep_cap_vertices'], m['prep_cap_faces']).query(v[sid])[0]
    r['original_gate_record'] = dict(r['gates'])
    r['gates']['margin'] = 'PASS' if d.max() <= 0.025 and plane.max() <= 0.025 else 'FAIL'
    r['margin'].update(actual_constructed_rim_max_source_distance_mm=float(d.max()), actual_rim_plane_residual_mm=float(plane.max()), source_curve_vertices=len(rim), full_curve_enclosure='Every piecewise linear source-section segment is retained; intersection float error not an empirical scanner bound.')
    if row['round'] in ['D', 'E']:
        r['shoulder']['rounded_radius_mm'] = None
        r['shoulder']['reason'] = 'Adaptive roof and shoulder rounding radius is not certified; physical material compliance UNKNOWN'
    r['assembly_shape'] = ashape
    r['assembly_shape']['source_preserving_native_lower'] = source_preserved
    r['shoulder_nominal_film'] = dict(min_mm=float(shoulder.min()), max_mm=float(shoulder.max()), field='PL interpolation between0um outer rim and50um inner shoulder seam', pass_design=bool(shoulder.min() >= -1e-08 and shoulder.max() <= 0.05 + 1e-08), resolution='PER_SURFACE_REGION', physical_status='PHENOMENOLOGICAL; no measured seated-film field')
    r['assembly_geometric_gates'] = dict(native_subset=r['gates']['native_containment'], restored_form='PASS' if ashape['p95_mm'] <= read(R / 'PREREG_B.json')['metrics']['shape_p95_mm'][row['family']] else 'FAIL', native_lower_model='PASS' if source_preserved else 'UNKNOWN_REGULARIZED', margin=r['gates']['margin'], wall=r['gates']['wall'], active_gap=r['gates']['nominal_gap'], shoulder_gap='PASS' if r['shoulder_nominal_film']['pass_design'] else 'FAIL', isolated_insertion=r['gates']['insertion'], mesh=r['gates']['mesh'], STL=r['gates']['STL'])
    r['restored_assembly_geometric_conjunction'] = all((x == 'PASS' for x in r['assembly_geometric_gates'].values()))
    r['complete_digital_crown'] = all((r['gates'][k] == 'PASS' for k in ['native_containment', 'form', 'margin', 'wall', 'nominal_gap', 'insertion', 'mesh', 'material_preparation', 'milling', 'STL']))
    r['complete_chain'] = all((x == 'PASS' for x in r['gates'].values()))
    r['full_preparation_subset'] = dict(status='PASS_SET_IDENTITY' if r['gates']['native_containment'] == 'PASS' else 'UNKNOWN', formula='P=(T intersect {z<=b}) union Q', certificate=row['containment'].get('path', row['containment'].get('components')), scope='Exact logical CSG set, not an unverified float Boolean or STL of P. Q certified on stored coordinates; T closed native completion.')
    path = D / (Path(row['mesh_path']).stem + '_assembly_surface.npz')
    np.savez_compressed(path, triangles=assembly)
    r['assembly_surface'] = dict(path=str(path), sha256=sha(path), scope='Visible outer source+prosthesis surface; not closed export')
    dump(R / 'exports' / (Path(row['mesh_path']).stem + '_CERTIFICATE.json'), r)
    return r

def run(tag='C_PILOT'):
    rows = []
    for row in read(R / f'FROZEN_PREDICTIONS_{tag}.json')['rows']:
        tick = time.perf_counter()
        if row['status'] == 'GENERATED':
            try:
                row = assess(row)
            except Exception as e:
                row = dict(**row, score_error=repr(e), complete_chain=False, complete_digital_crown=False, restored_assembly_geometric_conjunction=False)
                __import__('traceback').print_exc()
        else:
            row = dict(**row, complete_chain=False, complete_digital_crown=False, restored_assembly_geometric_conjunction=False)
        row['score_seconds'] = time.perf_counter() - tick
        rows.append(row)
        dump(R / f'RESULTS_{tag}.json', dict(rows=rows, claim_type='capability', external_referent=read(R / 'PREREG_A.json')['external_referent']))
        print(row['key'], row['material'], row.get('assembly_geometric_gates'), row.get('score_error'), flush=True)
if __name__ == '__main__':
    run(sys.argv[1] if len(sys.argv) > 1 else 'C_PILOT')
