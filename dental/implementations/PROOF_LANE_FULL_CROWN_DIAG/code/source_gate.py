"""Bind the two advertised source controls to extracted table cells."""
import copy
import re
import xml.etree.ElementTree as ET
from pathlib import Path
from bs4 import BeautifulSoup

def extracted_values(source_dir):
    root = ET.parse(Path(source_dir) / 'technician_repeat.xml')
    table = root.find('.//table-wrap[@id="tbl1"]')
    rows = [[' '.join(c.itertext()).strip() for c in r] for r in table.findall('.//tr')]
    row = next((r for r in rows if r[0] == 'Digital design'))
    means = [float(re.search('\\d+(?:\\.\\d+)?', c).group()) for c in row[1:3]]
    soup = BeautifulSoup((Path(source_dir) / 'madcrowner.html').read_text(), 'html.parser')
    tables = [soup.find('table', id='S4.T4.6.1')]
    assert len(tables) == 1
    assert tables[0] is not None and 'Fidelity' in tables[0].get_text()
    rows = [[c.get_text(' ', strip=True) for c in r.find_all(['td', 'th'])] for r in tables[0].find_all('tr')]
    row = next((r for r in rows if r and r[0] == 'MADCrowner'))
    values = [float(row[i]) for i in (3, 6, 9, 12)]
    return {'technician_within': means, 'MADCrowner_mesh': values}

def source_matches(records, source_dir):
    expected = extracted_values(source_dir)
    rows = {r['id']: r for r in records}
    return rows['technician_within']['value'] == expected['technician_within'] and rows['technician_within']['unit'] == 'um' and (rows['MADCrowner_mesh']['value'] == expected['MADCrowner_mesh']) and (rows['MADCrowner_mesh']['unit'] == ['mm2', 'mm2', 'mm', 'dimensionless'])

def corrected_controls(records, source_dir):
    assert source_matches(records, source_dir)
    wrong_rms = copy.deepcopy(records)
    next((r for r in wrong_rms if r['id'] == 'technician_within'))['value'][0] = 450
    wrong_hd = copy.deepcopy(records)
    next((r for r in wrong_hd if r['id'] == 'MADCrowner_mesh'))['value'][2] = 10.046
    return {'source_cells_match_report': True, 'wrong_within_RMS_450um_rejected': not source_matches(wrong_rms, source_dir), 'wrong_MADCrowner_HD_10mm_rejected': not source_matches(wrong_hd, source_dir)}
