from common import ROOT, dump, sha, stamp, check_frozen, state
from lxml import etree
if (ROOT / 'PREREG_R4.json').exists():
    raise SystemExit('R4 already frozen.')
old = check_frozen('PREREG_R3.json')
text = ' '.join(etree.parse(str(ROOT / 'raw/bond.xml')).getroot().itertext())
assert '4.91' in text and '6.2' in text and ('0.089' in text)
p = dict(old)
p.update(id='X74-restorative-R4', frozen_utc=stamp(), parent='PREREG_R3.json; rounds/R3.json; HANDOFF_R3.md', changed_operation='Only complete source geometries are generated. Verified flat optical tiles and shear substrate disc are exportable; size-5 gelatin capsule and missing crown dome are acquired/scanned, never invented.', source_geometry={'optical_tile_mm': [8.0, 8.0, 1.0], 'substrate_tile_mm': [10.0, 10.0, 2.0], 'optical_locator': 'PMC10703855 methods', 'bond_disc_diameter_mm': 14.0, 'bond_disc_thickness_mm': 2.5, 'capsule_nominal_diameter_mm': 4.91, 'capsule_nominal_height_mm': 6.2, 'capsule_wall_thickness_mm': 0.089, 'capsule_shape': 'UNKNOWN: buy the specified size-5 gelatin capsule and scan; nominal dimensions alone do not specify the dome.', 'cover_thickness_mm': 1.0, 'cover_internal_space_mm': 0.3, 'cover_shape': 'UNKNOWN until actual capsule/cast is scanned; no cover STL generated', 'bond_locator': 'PMC10829558 Specimen preparation'}, prediction_freeze='FROZEN_LAB_PREDICTIONS.json carries full tested protocols and source values, immutable before any future lab observation.')
dump(ROOT / 'PREREG_R4.json', p)
(ROOT / 'PREREG_R4.json.sha256').write_text(sha(ROOT / 'PREREG_R4.json') + '\n')
dump(ROOT / 'DECOMPOSITION_R4.json', {'parent_sha256': sha(ROOT / 'DECOMPOSITION_R2_R3.json'), 'revised_leaf': {'status': 'UNKNOWN', 'name': 'capsule/cover dome surface', 'stop': 'Nominal height, diameter and wall thickness do not determine the dome surface or bonded internal area; acquire and scan the actual item.'}, 'constructible_leaves': {'flat optical tiles': 'EXTERNALLY_MEASURED source dimensions', 'flat shear substrate disc': 'EXTERNALLY_MEASURED source dimensions', 'tangent rigid frame': 'DERIVED_UNDER_ASSUMPTIONS: orthonormal frame preserves dimensions; STL rounding tested, not fabricated accuracy'}, 'source_failure': 'R3 wrong nominal capsule values refuted by primary methods. All failed files and criteria remain frozen.'})
state(status='R4_FROZEN', latest_gate='R3 source-dimension failure preserved', next_operation='Export verified tile/disc geometries and typed future lab scorer; no unmeasured dome geometry')
print('R4 frozen', sha(ROOT / 'PREREG_R4.json'))
