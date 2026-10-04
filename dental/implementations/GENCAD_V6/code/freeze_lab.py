from common import *
import shutil
rows = read(ROOT / 'raw/R2_ROWS.json')
loads = {r['uid']: r for r in read(ROOT / 'raw/R3_ROWS.json')}
key = '89d83f17897eb151_molar_crown'
pair = []
for (label, name) in [('A', 'field_generator'), ('B', 'x42_shape')]:
    r = next((r for r in rows if r['key'] == key and r['participant'] == name and (r['track'] == 'V5B_R3')))
    mesh = Path(r['mesh_path'])
    dest = ROOT / 'exports' / label
    dest.mkdir(exist_ok=True, parents=True)
    for src in [mesh, mesh.with_name('crown.stl')]:
        shutil.copy2(src, dest / src.name)
    pair.append(dict(label=label, uid=r['uid'], mesh=str(dest / 'mesh.npz'), mesh_sha256=sha(dest / 'mesh.npz'), stl=str(dest / 'crown.stl'), stl_sha256=sha(dest / 'crown.stl'), geometry_accepted_by_parent=r['parent_geometry_pass'], local_contact_prediction=r['contact'], gap_scenario=r['pattern_enclosure'], per_tooth_load_sets=loads[r['uid']]['classes'], physical_prediction='UNKNOWN', source_mesh=str(mesh), source_mesh_sha256=sha(mesh), field_path=r['field_path'], field_sha256=r['field_sha256']))
freeze(ROOT / 'FROZEN_MEASUREMENT_PREDICTIONS.json', dict(claim_type='capability', key=key, pair=pair, preparation_path=str(V4 / 'payload/whole_inputs' / key / 'preparation.npz'), preparation_sha256=sha(V4 / 'payload/whole_inputs' / key / 'preparation.npz'), reference_path=str(V4 / 'payload/whole_private' / key / 'reference.npz'), reference_sha256=sha(V4 / 'payload/whole_private' / key / 'reference.npz'), measurement='NOT_RUN', coordinate_contract='mesh.npz uses site local mm; STL is predecessor supplied physical-world export; preparation source_R/source_base define local-to-world', protocol='LAB_PROTOCOL.md', scope='Digital predictions frozen before future matched observation; no absolute physical force prediction'))
(ROOT / 'exports/README.md').write_text('Two existing complete crown proposals on one registered lower first-molar site. A: field_generator; B: x42_shape. mesh.npz preserves indexed coordinates/roles; crown.stl is an exact copy of the parent export. Dimensions are mm. Source transforms and predictions are pinned in FROZEN_MEASUREMENT_PREDICTIONS.json. Neither an accepted clinical design nor a measured fabricated specimen. See LAB_PROTOCOL.md.\n')
