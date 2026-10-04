"""One-time R3 refinement. Preserve scientific values and previous draft bytes."""
import datetime, hashlib, json
from pathlib import Path
R = Path(__file__).resolve().parent

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def write(p, x):
    p.write_text(json.dumps(x, ensure_ascii=False, indent=2) + '\n')
previous = json.loads((R / 'raw/R3_PREVIOUS_VERIFIED_OBSERVATIONS.json').read_text())
assert sha(R / 'VERIFIED_OBSERVATIONS.json') == sha(R / 'raw/R3_PREVIOUS_VERIFIED_OBSERVATIONS.json')
diagnostic = {'P10_Donnermeyer', 'P11_Jang', 'P12_Allihaibi', 'P13_Piecha', 'P14_Raghav', 'P15_Hazard', 'P16_Ma', 'P20_Mayo'}
count_fields = {'P01_Chew': 'teeth', 'P03_Gurler': 'teeth', 'P04_Ramani': 'participants', 'P05_Hoang': 'studies', 'P06_Leong': 'reviews', 'P07_Sulaiman': 'teeth', 'P08_Li': 'teeth', 'P09_Linas': 'teeth', 'P10_Donnermeyer': 'studies/publications', 'P11_Jang': 'patients', 'P12_Allihaibi': 'teeth/roots', 'P13_Piecha': 'reader-case pairs', 'P14_Raghav': 'patients', 'P15_Hazard': 'patients/teeth', 'P16_Ma': 'images', 'P17_Li_review': 'studies/articles', 'P18_dogs': 'cases/dogs/teeth', 'P19_Zanini': 'publications', 'P20_Mayo': 'cases'}
dimensionless = {'OR_FP_PP', 'pain_bleeding_rho', 'age_outcome_V', 'necrosis_OR', 'accuracy', 'sensitivity', 'specificity', 'internal_accuracy', 'external_accuracy', 'ICC_lateral', 'ICC_AP'}
percentage = {'PP_full', 'PP_deep', 'ITT_full', 'ITT_deep', 'success_percent', 'effective_percent', 'uncertain_percent', 'ineffective_percent', 'PA_AP', 'CBCT_AP', 'clinician_sensitivity', 'AI_sensitivity', 'conventional_accuracy', 'digital_accuracy', 'ultrasound_accuracy', 'intact_percent', 'dehiscence_percent', 'fenestration_percent'}
for row in previous:
    id = row['id']
    if id in diagnostic:
        row['time_scale'] = 'SIMULTANEOUS'
    if id in {'C06_replica_mass', 'C11_postcore'}:
        row['observation_resolution'] = row['edge_resolution'] = 'PER_TOOTH'
    if id in {'C12_KATANA', 'C13_Ivoclar'}:
        row['statistic_resolution'] = 'PER_POINT'
    row['value_metadata'] = {}
    for key in row['values']:
        unit = row['unit']
        resolution = row['observation_resolution']
        time = row['time_scale']
        if id.startswith('P'):
            unit = count_fields.get(id, row['unit'])
            if key in percentage:
                unit = '%'
            if key in dimensionless:
                unit = 'dimensionless'
            if id == 'P04_Ramani' and key == 'success':
                unit = '%'
            if id == 'P05_Hoang' and key != 'studies':
                unit = 'dimensionless'
            if id == 'P10_Donnermeyer':
                unit = 'studies' if key == 'studies' else 'publications'
            if id == 'P12_Allihaibi' and key in {'teeth', 'roots', 'PARL_teeth'}:
                unit = 'roots' if key == 'roots' else 'teeth'
            if id == 'P15_Hazard' and key not in dimensionless:
                unit = 'patients' if key == 'patients' else 'teeth'
            if id == 'P17_Li_review':
                unit = 'studies' if key == 'studies' else 'articles'
            if id == 'P18_dogs':
                unit = 'dogs' if key == 'dogs' else 'teeth' if key in {'teeth', 'success'} else 'cases'
            if id == 'P07_Sulaiman' and key in {'pain_bleeding_rho', 'age_outcome_V'}:
                time = 'SIMULTANEOUS' if key == 'pain_bleeding_rho' else 'HANDOVER'
        if id == 'F05_maxshare' and 'teeth' in key:
            unit = 'teeth'
            resolution = 'PER_ARCH'
        if id == 'F06_regions_conflict':
            resolution = 'PER_SURFACE_REGION'
        if id == 'F07_reference':
            unit = 'dimensionless' if key.startswith('ICC') else '%'
        if id == 'F08_Hattori':
            resolution = 'PER_ARCH' if key.startswith('arch') else 'PER_POINT'
        if id in {'G01_Varga', 'G02_Wu'}:
            unit = 'degree' if 'angle' in key else 'mm'
        if id == 'M01_RMS' and key == 'whole':
            resolution = 'PER_TOOTH'
        row['value_metadata'][key] = {'unit': unit, 'observation_resolution': resolution, 'statistic_resolution': row['statistic_resolution'], 'time_scale': time, 'available_raw_individual_observations': False}
    if id.startswith('P'):
        row['quantity'] = row['protocol']
current = json.loads((R / 'VERIFIED_OBSERVATIONS.json').read_text())
assert len(previous) == len(current)
assert all((a['values'] == b['values'] for (a, b) in zip(previous, current)))
write(R / 'VERIFIED_OBSERVATIONS.json', previous)
files = ['VERIFIED_OBSERVATIONS.json', 'SOURCE_AUDIT.json', 'PREREG_R1.json', 'PREREG_R2.json', 'PREREG_R3.json']
write(R / 'CURATED_SOURCE_LOCK.json', {'timestamp_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'sha256': {name: sha(R / name) for name in files}, 'review_state': 'PENDING_INDEPENDENT_REVIEW', 'previous_lock': 'raw/R3_PREVIOUS_CURATED_SOURCE_LOCK.json'})
draft = R.parent.parent / 'notes/expansion/src_verified_20261004.jsonl'
receipt = json.loads((R / 'raw/DRAFT_RECEIPT.json').read_text())
assert sha(draft) == receipt['sha256'] == sha(R / 'raw/R3_PREVIOUS_src_verified_20261004.jsonl')
draft.write_text(''.join((json.dumps(row, ensure_ascii=False, sort_keys=True) + '\n' for row in previous)))
write(R / 'raw/R3_REFINEMENT.json', {'unchanged_value_bundles': len(previous), 'field_contract_count': sum((len(x['values']) for x in previous)), 'corrected_same_visit_time_bundles': sorted(diagnostic), 'previous_draft_sha256': receipt['sha256'], 'refined_draft_sha256': sha(draft), 'scientific_gate_changed': False})
print('Refined', len(previous), 'bundles without changing values.')
