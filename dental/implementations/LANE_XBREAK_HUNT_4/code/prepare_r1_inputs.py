import json, hashlib, datetime
from pathlib import Path
p = json.loads(Path('PREREG_R1.json').read_text())
p['round'] = 'R1v2'
p['frozen_at_utc'] = datetime.datetime.now(datetime.timezone.utc).isoformat()
p['supersedes'] = 'R1 input contract invalid BEFORE numerical exploration: Table1 has class totals, not test marginals; Figure1 carries marginals and missing tests.'
p['input'] = 'PMC4884138 Table1 reference totals + Figure1 exact marginal counts. Preserve cold missing=6 and percussion missing=1 as third categories; 18 atoms. Table2 was visually read during source acquisition and is NOT blinded; not consumed by R1.'
p['strongest_equally_informed_control'] = 'Independent 18-atom Charnes-Cooper LP including missing test categories. Same information as analytic Frechet calculation.'
p['metrics']['facit'] = 'Source Figure1 and later Table2 cross-check, same cohort; no independent validation.'
Path('PREREG_R1v2.json').write_text(json.dumps(p, indent=2) + '\n')
Path('PREREG_R1v2.sha256').write_text(hashlib.sha256(Path('PREREG_R1v2.json').read_bytes()).hexdigest() + '\n')
d = json.loads(Path('DECOMPOSITION_R1.json').read_text())
d['round'] = 'R1v2'
d['representation'] = '18 (D,C,P) cells, C/P categories 0,1,missing. Known reference totals and two ternary marginal tables.'
d['leaves'].append({'leaf': 'missing_index_tests', 'status': 'EXTERNALLY_MEASURED', 'basis': 'Figure1 cold response n702 and percussion n707 versus reference N708. Missing=6 and1.', 'stop_argument': 'No missing-at-random assumption; unknown joint allocation retained.'})
Path('DECOMPOSITION_R1v2.json').write_text(json.dumps(d, indent=2) + '\n')
s = {'source': {'pmid': '27118600', 'doi': '10.1016/j.joen.2016.03.016', 'pmcid': 'PMC4884138', 'figure': 'Figure1', 'table_reference_totals': 'Table1', 'source_sha256': hashlib.sha256(Path('sources/pbrn.pdf').read_bytes()).hexdigest(), 'extraction': 'Manual transcription of original PDF page938 Figure1; HTML image independently available. Review pending.'}, 'N': 708, 'labels': {'D0': 'Bleeding observed in chamber', 'D1': 'No bleeding observed in chamber', 'C1': 'No response to cold', 'C0': 'Response to cold', 'C2': 'Cold response missing', 'P1': 'Pain on percussion', 'P0': 'No pain on percussion', 'P2': 'Percussion missing'}, 'class_totals': [359, 349], 'cold_margins': [[284, 73, 2], [38, 307, 4]], 'percussion_margins': [[148, 210, 1], [99, 250, 0]], 'selection': 'All scheduled for initial orthograde RCT; one tooth per patient. 62 treating dentists; source-cluster independence UNKNOWN.'}
Path('raw/R1_SOURCE_COUNTS.json').write_text(json.dumps(s, indent=2) + '\n')
