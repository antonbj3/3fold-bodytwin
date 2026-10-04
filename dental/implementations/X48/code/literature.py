"""Source tables remain product-, measurement- and window-specific."""
from common import *
import xml.etree.ElementTree as ET, collections, csv

def run():
    material = read(DENTAL / 'notes/aligner_intervals/ALIGNER_INTERVALS.json')
    x = ET.parse(ROOT / 'inputs/bilello_2022.xml').getroot()
    doi = x.find(".//article-id[@pub-id-type='doi']").text
    table = x.find(".//table-wrap[@id='Tab3']")
    parsed = []
    missing = []
    for row in table.findall('.//tr')[1:]:
        cells = [''.join(c.itertext()) for c in row]
        for (quantity, v) in zip(['axial_rotation', 'intrusion', 'buccolingual_tip'], cells[1:]):
            r = dict(tooth_group=cells[0], motion=quantity, locator=f'https://doi.org/{doi}#Tab3', source_file='inputs/bilello_2022.xml', source_sha256=sha(ROOT / 'inputs/bilello_2022.xml'), source_table='Table3', source='Bilello 2022', product='Invisalign', patients=10, time_window='end of treatment INCLUDING refinements (mean25 refinement trays per arch)', resolution='POPULATION', metric='reported mean accuracy percentage; source remaps >100% values', comparable_to_X48=False, reason='X48 is weekly Naturaligner crown pose relative to moving incisors; source is end-treatment component accuracy in occlusal planes, with refinement and auxiliaries')
            if v == '(**)':
                missing.append({**r, 'status': 'SOURCE_INSUFFICIENT_SAMPLE'})
            else:
                parsed.append({**r, 'value_pct': float(v), 'status': 'SOURCE_REPORTED_MEASUREMENT'})
    abstract = [dict(source='Haouili2020', locator='https://pubmed.ncbi.nlm.nih.gov/32620479/', doi='10.1016/j.ajodo.2019.12.015', source_location='primary abstract Results', patients=38, product='Invisalign Full/Teen', motion=k, value_pct=v, resolution='POPULATION', time_window='end treatment', comparable_to_X48=False, reason='different material, treatment duration, absolute reference and component metric') for (k, v) in [('all_movements', 50.0), ('buccolingual_tip', 56.0), ('axial_rotation', 46.0)]]
    abstract.extend([dict(source='Sachdev2021', locator='https://pubmed.ncbi.nlm.nih.gov/34625386/', doi='10.1016/j.ejwf.2021.08.003', source_location='primary abstract Results', patients=30, product='in-house clear aligners', motion=k, value_pct=v, resolution='POPULATION', time_window='end treatment', comparable_to_X48=False, reason='anterior-only cohort, fixed posterior/rugae reference, different product/window/metric') for (k, v) in [('all_movements', 56.18), ('mesiodistal_translation', 72.33), ('intrusion', 43.28)]])
    write(ROOT / 'raw/LITERATURE_MOTION_REFERENTS.json', dict(supplied_bank=dict(file=str(DENTAL / 'notes/aligner_intervals/ALIGNER_INTERVALS.json'), records=len(material), usable_movement_records=0, rejected_movement_records=len(material), rejection_reason='wrong quantity: thickness, stress retention, modulus and Tg; no achieved/planned motion entries', rejected_fraction=1.0, quantities=dict(collections.Counter((r['quantity'] for r in material)))), table_rows=parsed, source_unavailable_cells=missing, primary_abstract_rows=abstract, movement_claims_not_available=dict(extrusion='no numeric primary extrusion value extracted from these sources', bodily='root displacement unobserved; no bodily claim from IOS crowns'), external_referent=dict(kind='independent_measurement', locator=f'https://doi.org/{doi}#Tab3', compared_quantity='reported planned/achieved movement accuracy by tooth group and motion; explicit mismatch to X48 relative weekly components', refutes_us=True), direct_validation_eligible_rows=0))
    motion = read(ROOT / 'raw/R2_RELATIVE_POSE.json')
    exports = []
    for r in motion:
        for (k, a) in r['signed_components']['achieved'].items():
            p = r['signed_components']['planned'][k]
            motion_type = 'intrusion' if p >= 0 else 'extrusion'
            if k != 'intrusion_extrusion':
                motion_type = k
            exports.append(dict(patient=r['patient'], jaw=r['jaw'], fdi=r['fdi'], anchor_fdi=r['anchor_fdi'], week=r['week'], interval=r['interval'], motion=motion_type, component=k, achieved=a, planned=p, unit='deg' if k in ['buccolingual_tip', 'mesiodistal_tip', 'axial_rotation'] else 'mm', achieved_over_planned=a / p if abs(p) > 1e-12 else None, resolution='PER_TOOTH', quantity_definition='relative-to-named-moving-incisor crown component', status='DIAGNOSTIC_ONLY; physical error bound UNKNOWN', rigorous_signed_component_enclosure='MISSING', absolute_motion='UNKNOWN', body_translation='UNKNOWN_ROOT_NOT_OBSERVED'))
    with (ROOT / 'raw/movement_types.csv').open('w') as o:
        w = csv.DictWriter(o, fieldnames=list(exports[0]))
        w.writeheader()
        w.writerows(exports)
    (rows, ix) = pose_index()
    floor = []
    for r in rows:
        if r['source'] != 'Sirona' or r['week'] == 0:
            continue
        floor.append(dict(patient=r['patient'], jaw=r['jaw'], fdi=r['fdi'], week=r['week'], resolution='PER_TOOTH', fit_surface_plane_RMS_mm=r['surface_error']['trimmed_plane_rms_mm'], fit_surface_point_RMS_mm=r['surface_error']['trimmed_point_rms_mm'], sampling_spread_mm=r['sampling_spread_mm'], point_plane_centroid_delta_mm=r['point_control_centroid_delta_mm'], point_plane_rotation_delta_deg=r['point_control_rotation_delta_deg'], reference_choice_centroid_spread_mm=r['reference_spread_mm'], reference_choice_rotation_spread_deg=r['reference_rotation_spread_deg'], patient_repeat_scan_floor='UNKNOWN', source_label='published transferred crown segmentation'))
    with (ROOT / 'raw/registration_floor.csv').open('w') as o:
        w = csv.DictWriter(o, fieldnames=list(floor[0]))
        w.writeheader()
        w.writerows(floor)
    print('Primary table extracted', len(parsed), 'numeric cells', len(missing), 'source-unavailable; diagnostic movement rows', len(exports))
if __name__ == '__main__':
    run()
