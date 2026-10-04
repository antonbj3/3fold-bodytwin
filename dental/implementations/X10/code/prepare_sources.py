import csv, json, re, time
from common import *
start = time.perf_counter()
manifest = []
for pmc in ['PMC10729975', 'PMC10265779', 'PMC10048864', 'PMC10923227']:
    (path, ts, license_text) = tables(pmc)
    write('raw/' + pmc + '_tables.json', ts)
    manifest.append({'path': str(path), 'sha256': sha(path), 'bytes': path.stat().st_size, 'license_statement': license_text, 'derived_file': 'raw/' + pmc + '_tables.json'})
global_values = {'Primescan': (17.3, 4.9), 'Trios 4': (20.8, 6.2), 'Medit i500': (25.2, 7.3), 'CS3600': (26.9, 15.9), 'Trios 3': (27.7, 6.8), 'Omnicam 4.6': (57.5, 3.2)}
write('raw/global_scanners.json', {'doi': '10.3390/dj9070075', 'locator': 'Table 2', 'verification_url': 'https://pmc.ncbi.nlm.nih.gov/articles/PMC8303663/', 'local_seed_path': str(PROJECT / 'notes/chains_brainstorm/04_manufacturing_materials.md'), 'estimand': 'full arch best-fit surface deviation summary, not a local bound or signed sigma', 'values': {k: {'mean_um': a, 'between_scan_sd_um': b} for (k, (a, b)) in global_values.items()}})
ts = load('raw/PMC10729975_tables.json')

def matrix(table_number):
    rows = ts[table_number - 1]['rows']
    out = {}
    for r in rows:
        if r and r[0] in ['Omni', 'Prime', 'Trios 3', 'Trios 4']:
            assert len(r) == 8, r
            out[r[0]] = [float(x) for x in r[1:]]
    assert len(out) == 4
    return out
write('raw/implant_medians.json', {'doi': '10.1371/journal.pone.0295790', 'n_per_cell': 10, 'replicates_are_not_independent_patient_trials': True, 'position_3d_um': matrix(5), 'angular_deg': matrix(8), 'cross_arch_span_um': matrix(11), 'table_locators': {'position_3d_um': 'Table 5', 'angular_deg': 'Table 8', 'cross_arch_span_um': 'Table 11'}, 'estimand': 'group median geometric error against D2000 reference; no tail, direction or paired samples supplied'})
(d, h) = clinical_primary_record('25693497')
write('raw/clinical_primary_record.json', d)
manifest.append({'locator': 'local cadcam.jsonl exact line PMID25693497', 'sha256_line': h, 'doi': d['doi'], 'derived_file': 'raw/clinical_primary_record.json', 'license': 'primary abstract metadata; numerical facts only used; publisher full text not redistributed'})
assert '149' in d['abstractText'] and '112' in d['abstractText']
clinical = [{'scanner': 'Omnicam 4.6', 'actual_name': 'CEREC AC Omnicam', 'value_um': 149.0, 'summary': 'median', 'iqr_um': [114.0, 218.0], 'doi': d['doi'], 'locator': 'primary abstract Results', 'version_equivalence': 'UNKNOWN', 'n_teeth': 49, 'n_patients': 24}, {'scanner': 'Trios 3', 'actual_name': 'Heraeus Cara TRIOS', 'value_um': 112.0, 'summary': 'median', 'iqr_um': [94.0, 149.0], 'doi': d['doi'], 'locator': 'primary abstract Results', 'version_equivalence': 'UNKNOWN', 'n_teeth': 49, 'n_patients': 24}]
review = load('raw/PMC10923227_tables.json')[0]
target_row = [r for r in review['rows'] if r and 'Lee et al., 2020' in ' '.join(r)]
assert len(target_row) == 1
numbers = ' '.join(target_row[0])
assert '49.1' in numbers and '56.5' in numbers
write('raw/clinical_2020_transcription.json', {'primary_doi': '10.3390/jcm9124035', 'primary_pmc': 'PMC7764839', 'local_table_source': 'PMC10923227 Table 1', 'row': target_row[0], 'status': 'SECONDARY_TABLE_TRANSCRIPTION; primary fulltext absent locally; primary abstract verifies design, not numbers'})
for (scanner, v, sd) in [('Medit i500', 49.1, 7.7), ('CS3600', 56.5, 12.7)]:
    clinical.append({'scanner': scanner, 'actual_name': scanner, 'value_um': v, 'sd_um': sd, 'summary': 'mean', 'doi': '10.3390/jcm9124035', 'locator': 'local PMC10923227 Table 1 Lee2020 row', 'version_equivalence': 'UNKNOWN', 'n_patients': 20, 'source_level': 'secondary transcription of primary measurement'})
write('raw/clinical_gaps.json', clinical)
write('SOURCE_MANIFEST.json', {'sources': manifest, 'preparation_seconds': time.perf_counter() - start, 'data_policy': 'Existing local corpus only for external measurement tables. Public web used for literature verification; no raw dataset download.', 'no_personal_data_extracted': True})
print('Prepared local tables and clinical numeric summaries.')
