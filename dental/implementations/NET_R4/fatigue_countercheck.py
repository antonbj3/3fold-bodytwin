import xml.etree.ElementTree as ET
from common import *

def run():
    root = ET.parse(source_path('fatigue_counterexample'))
    doi = root.find('.//article-meta/article-id[@pub-id-type="doi"]').text
    rows = []
    group = None
    for table in root.findall('.//table-wrap'):
        if not table.get('id', '').endswith('t003'):
            continue
        for tr in table.findall('.//tr'):
            cells = [' '.join(''.join(x.itertext()).split()) for x in tr.findall('td')]
            if len(cells) == 1 and 'Implant' in cells[0]:
                group = cells[0]
            if len(cells) == 4 and '*' in cells[1]:
                rows.append({'system': group, 'loading_level_percent': int(cells[0]), 'peak_N': int(cells[1].split()[0]), 'cycles': [int(x.replace(',', '').strip()) for x in cells[2].split(';')], 'table': table.get('id'), 'resolution_level': 'PER_TOOTH', 'specimen_ids': 'not published; table positions retained'})

    def universal40(rs):
        return all((x['loading_level_percent'] == 40 for x in rs))
    choi = [x for x in read('SOURCE_OBSERVATIONS.json')['records'][2]['reported']['fatigue_rows'] if not any(x['events'])]
    mutated = [dict(x, loading_level_percent=50) for x in choi]
    out = {'source': {'kind': 'independent_measurement', 'locator': 'https://doi.org/' + doi + '#Table3', 'local_file': entry('fatigue_counterexample'), 'compared_quantity': 'Tested load fraction with three5e6-cycle runouts per system', 'refutes_us': True}, 'rows': rows, 'Choi_local_40percent_matches': universal40(choi), 'external_universal_40percent_rejected': not universal40(rows), 'injected_50percent_rejected': not universal40(mutated), 'test_count': 1, 'scope': 'Rejects a universal exact40% rule; does not infer a continuous population endurance limit or determine the target-system fatigue model.', 'all_pass': len(rows) == 3 and (not universal40(rows)) and (not universal40(mutated))}
    write('FATIGUE_COUNTERCHECK.json', out)
    return out
if __name__ == '__main__':
    print(run()['all_pass'])
