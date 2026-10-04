"""Small, source-located local inputs. No dataset download or tree copy."""
from dental_release.paths import expand as _release_expand
import hashlib, json, time
from pathlib import Path
from lxml import etree
ROOT = Path(__file__).resolve().parents[1]
DENT = ROOT.parents[1]
FIELD = Path(_release_expand('@DENTAL_EXTERNAL_ROOT@/engines/3fold-motion-engine/_private/romi_collab/build/SOL_FALT_ALIGNERGREN_20261001'))
SOURCES = {'park': Path(_release_expand('@DENTAL_CORPUS_ROOT@/europepmc/fulltext/PMC10322848.xml')), 'cho': FIELD / 'raw/Cho_fulltext.xml', 'arch': DENT / 'results/LANE_X19_ALIGNER_FORCE/inputs/ARCH_GEOMETRY.json', 'kaur': DENT / 'results/LANE_X19_ALIGNER_FORCE/inputs/SENSOR_KAUR_2021.json', 'intervals': DENT / 'notes/aligner_intervals/ALIGNER_INTERVALS.json', 'continuous_gap': DENT / 'results/LANE_NEXT_P_OCCLUSION_VALIDATION/code/continuous_gap.py', 'field_port': FIELD / 'PORT.json'}

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def write(p, data):
    (ROOT / p).write_text(json.dumps(data, indent=2, ensure_ascii=False, allow_nan=False) + '\n')

def main():
    start = time.perf_counter()
    import re
    tree = etree.parse(str(SOURCES['park']))
    records = []
    for (tab, quantity) in [('Tab3', 'finished_thickness'), ('Tab4', 'passive_gap')]:
        for (i, row) in enumerate(tree.xpath(f'//table-wrap[@id="{tab}"]//tr'), 1):
            cells = [''.join(c.itertext()).strip() for c in row]
            if not cells or cells[0] not in ['Anterior teeth', 'Posterior teeth', 'Palatogingival', 'Palatal', 'Incisal/occlusal', 'Buccal', 'Buccogingival']:
                continue
            for (j, arm) in enumerate(['TS', 'TM', 'PA', 'PC'], 1):
                vals = re.findall('[0-9]+(?:\\.[0-9]+)?', cells[j])
                records.append(dict(quantity=quantity, region=cells[0], arm=arm, original_text=cells[j], median_mm=float(vals[0]) / 1000, q1_mm=float(vals[1]) / 1000, q3_mm=float(vals[2]) / 1000, locator=f'''{SOURCES['park']}#(//table-wrap[@id="{tab}"]//tr)[{i}]/*[{j + 1}]''', resolution_level='PER_SURFACE_REGION', time_scale='HANDOVER', source_doi='10.1038/s41598-023-36851-5', joint_tooth_surface_map=False))
    write('inputs/PARK_REGIONS.json', records)
    arch = json.loads(SOURCES['arch'].read_text())
    selected = [a for a in arch['arches'] if a['dataset'] == 'Teeth3DS' and a['jaw'] == 'upper']
    write('inputs/ARCHES.json', selected)
    write('inputs/KAUR_REFERENT.json', json.loads(SOURCES['kaur'].read_text()))
    cho = etree.parse(str(SOURCES['cho']))
    rows = []
    for (tab, fdi) in [('T4', 12), ('T5', 11)]:
        for (i, tr) in enumerate(cho.xpath(f'//table-wrap[@id="{tab}"]//tr'), 1):
            cells = [''.join(c.itertext()).strip() for c in tr]
            if len(cells) < 2 or cells[1] != 'Amount':
                continue
            component = cells[0][:2]
            factor = 0.001 if component.startswith('F') else 1.0
            numbers = [c for c in cells[2:] if c]
            if len(numbers) != 15:
                raise ValueError((tab, i, cells, numbers))
            for (k, material) in enumerate(['Zendura', 'Trioclear', 'Graphy']):
                v = numbers[k * 5:k * 5 + 5]
                for (j, delta) in enumerate([0.3, 0.6]):
                    rows.append(dict(fdi=fdi, material=material, activation_mm=delta, component=component, mean=float(v[j * 2].replace(',', '').replace('–', '-')) * factor, sd=float(v[j * 2 + 1].replace(',', '').replace('–', '-')) * factor, unit='N' if component.startswith('F') else 'Nmm', original_mean=v[j * 2], original_unit='mN' if component.startswith('F') else 'mN*M (m interpreted as metre; corroborated by Dahlberg N mm text)', locator=f'''{SOURCES['cho']}#(//table-wrap[@id="{tab}"]//tr)[{i}]''', source_doi='10.4041/kjod25.003', resolution_level='PER_TOOTH', time_scale='SIMULTANEOUS', moment_origin='root apex (source calls it estimated centre of resistance)', baseline='sensor zeroed after seating separate passive baseline aligner', acquisition='last 10 s of 2 min; 37 C; six specimens per material/activation group'))
    write('inputs/CHO_SENSOR.json', rows)
    intervals = json.loads(SOURCES['intervals'].read_text())
    relaxation = []
    for i in [101, 122, 131, 137, 138, 139]:
        record = dict(intervals[i], source_index=i, resolution_level='PHENOMENOLOGICAL', time_scale='HANDOVER')
        record['can_scale_shell_modulus_conditionally'] = i in [101, 122]
        record['replacement_measurement'] = 'Matched formed aligner, same seating/activation, same temperature and elapsed time, six-axis force trace'
        relaxation.append(record)
    write('inputs/RELAXATION.json', relaxation)
    manifest = [dict(key=k, path=str(p), sha256=sha(p), bytes=p.stat().st_size) for (k, p) in SOURCES.items()]
    write('SOURCE_MANIFEST.json', manifest)
    write('raw/PREPARATION.json', dict(wall_s=time.perf_counter() - start, source_files=len(manifest), park_candidate_rows=len(records), park_primary_kept=sum((r['arm'] == 'TS' for r in records)), park_rejected=sum((r['arm'] != 'TS' for r in records)), park_rejection_reason='Different material/process arms; not fused into Duran TS', arches_candidates=len(arch['arches']), arches_kept=len(selected), arches_rejected=len(arch['arches']) - len(selected), arch_rejection_reasons='Predicted FDI geometry or lower arch; source-labelled upper arches only', tooth_type_surface_cells_requested=4 * 5, tooth_type_surface_cells_measured=0, joint_table_unavailable_fraction=1.0, cho_rows=len(rows), cho_duplicate_tables_not_counted=['T6', 'T7'], kaur_force_rows=3, kaur_moment_rows_rejected=3, kaur_moment_rejection_reason='Moment origin/table unavailable in primary abstract', own_mesh_copies=0, license={'Park': 'CC BY 4.0', 'Cho': 'CC BY-NC 4.0', 'Kaur': 'Copyright; only transcribed factual mean/SD', 'arch': 'Private local derived Teeth3DS landmarks; source licence must be checked before redistribution'}))
    print('Local inputs prepared:', len(records), 'Park regional rows;', len(rows), 'Cho sensor rows;', len(selected), 'upper arches')
if __name__ == '__main__':
    main()
