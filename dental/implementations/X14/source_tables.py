"""Extract primary JATS tables, preserving rows, units and source lineage."""
from pathlib import Path
import hashlib
import re
import xml.etree.ElementTree as ET
ROOT = Path(__file__).resolve().parent

def expanded_rows(table):
    pending = {}
    output = []
    for tr in table.findall('./tbody/tr'):
        row = {}
        for (col, (remaining, value)) in list(pending.items()):
            row[col] = value
            if remaining == 1:
                del pending[col]
            else:
                pending[col] = (remaining - 1, value)
        col = 0
        for td in tr.findall('td'):
            while col in row:
                col += 1
            value = ' '.join(''.join(td.itertext()).split())
            (rs, cs) = (int(td.get('rowspan', 1)), int(td.get('colspan', 1)))
            for j in range(cs):
                row[col + j] = value
                if rs > 1:
                    pending[col + j] = (rs - 1, value)
            col += cs
        assert sorted(row) == list(range(len(row)))
        output.append([row[j] for j in range(len(row))])
    return output

def hirano_angles(table_number):
    path = ROOT / 'sources/Hirano2025.xml'
    table_id = f'materials-18-04234-t00{table_number}'
    tree = ET.parse(path)
    rows = expanded_rows(tree.find(f".//table-wrap[@id='{table_id}']/table"))
    maker = None
    result = []
    for (i, row) in enumerate(rows):
        if row[0]:
            maker = row[0]
        (pre, post) = [float(x.replace('−', '-')) for x in row[2:4]]
        result.append(dict(material=('n' if table_number == 3 else 'c') + '-MCL-' + maker, manufacturer=maker, area=row[1].removeprefix('Area '), pre_deg=pre, post_deg=post, delta_deg=post - pre, n_specimens=7, statistic='tabulated aggregate, statistic unnamed', sd_deg=None, source_table_id=table_id, source_row=i + 1, locator=f'doi:10.3390/ma18184234 Table{table_number} row{i + 1}', source_sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    return result

def shrinkage_cells():
    path = ROOT / 'sources/Shrinkage2025.xml'
    tree = ET.parse(path)
    table_id = 'materials-18-03217-t002'
    rows = expanded_rows(tree.find(f".//table-wrap[@id='{table_id}']/table"))
    result = []
    for (i, row) in enumerate(rows):
        assert len(row) == 7, (i, row)
        for (method, cell) in zip(['micrometer', 'light_microscopy', 'surface_scan'], row[4:]):
            numbers = re.findall('[-−]?\\d+(?:\\.\\d+)?', cell)
            assert len(numbers) in (1, 2), (i, cell)
            mean = float(numbers[0].replace('−', '-'))
            sd = float(numbers[1]) if len(numbers) == 2 else None
            result.append(dict(material=row[0], position=row[1], horizontal=row[2], axis=row[3], method=method, mean_shrinkage_pct=mean, sd_shrinkage_pct=sd, n_per_cell=None, text=cell, source_table_id=table_id, source_row=i + 1, locator=f'doi:10.3390/ma18143217 Table2 row{i + 1}, {method}', source_sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    assert len(result) == 459
    return result
