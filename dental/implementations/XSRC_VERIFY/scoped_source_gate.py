"""Independent Table 2 gate for sampled C04; compares labelled cells, not tokens."""
import copy, json, re, sys, xml.etree.ElementTree as E
from decimal import Decimal
from pathlib import Path

def reference(source_root):
    root = E.parse(source_root / 'sources/PMC12696814.xml').getroot()
    tab = next((t for t in root.findall('.//table-wrap') if t.get('id') == 't2'))
    out = {}
    suffix = {'o': 'occlusal', 'cp': 'cusp', 'ax': 'axial', 'ch': 'chamfer'}
    for tr in tab.findall('.//tr')[1:]:
        cells = [' '.join(''.join(c.itertext()).split()) for c in tr]
        label = cells[0]
        material = next((m for m in ['PIC', 'ZR', 'LD'] if label.startswith(m)))
        suf = label[len(material):]
        numeric = [re.findall('\\d+(?:\\.\\d+)?', x)[:2] for x in cells if '±' in x]
        out[material + '_' + suffix[suf]] = numeric[0]
        if len(numeric) == 2:
            out['marginal_' + material] = numeric[1]
    return out

def check(row, source_root):
    expected = reference(source_root)
    return row['id'] == 'C04_regions' and row['unit'] == 'um' and (row['table_locator'] == 'Table2 t2') and (set(row['values']) == set(expected)) and all((list(map(Decimal, row['values'][k])) == list(map(Decimal, v)) for (k, v) in expected.items()))
if __name__ == '__main__':
    source_root = Path(sys.argv[1])
    row = next((r for r in json.loads((source_root / 'VERIFIED_OBSERVATIONS.json').read_text()) if r['id'] == 'C04_regions'))
    mutant = copy.deepcopy(row)
    (mutant['values']['PIC_occlusal'], mutant['values']['PIC_axial']) = (mutant['values']['PIC_axial'], mutant['values']['PIC_occlusal'])
    result = {'valid_row_accepted': check(row, source_root), 'known_numbers_wrong_region_rejected': not check(mutant, source_root), 'conclusion_changed': False, 'scope': 'C04 Table2 only; remaining manual source semantics separately reviewed.'}
    assert result['valid_row_accepted'] and result['known_numbers_wrong_region_rejected']
    print(json.dumps(result))
