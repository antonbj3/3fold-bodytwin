from dental_release.paths import expand as _release_expand
from common import *
import xml.etree.ElementTree as ET
from decimal import Decimal
SOURCE = Path(_release_expand('@DENTAL_CORPUS_ROOT@/europepmc/fulltext/PMC10958416.xml'))

def source_table():
    table = ET.parse(SOURCE).getroot().findall('.//table-wrap')[0]
    grid = []
    pending = {}
    for tr in table.findall('.//tr'):
        row = {}
        nxt = {}
        for (k, (left, txt)) in pending.items():
            row[k] = txt
            if left > 1:
                nxt[k] = (left - 1, txt)
        col = 0
        for cell in tr:
            while col in row:
                col += 1
            txt = ' '.join(''.join(cell.itertext()).split())
            rs = int(cell.get('rowspan', '1'))
            cs = int(cell.get('colspan', '1'))
            for k in range(col, col + cs):
                row[k] = txt
                if rs > 1:
                    nxt[k] = (rs - 1, txt)
            col += cs
        grid.append([row.get(k, '') for k in range(20)])
        pending = nxt
    return grid

def run():
    if not (ROOT / 'EXTERNAL_SOURCE_LOCK.json').exists():
        freeze(ROOT / 'EXTERNAL_SOURCE_LOCK.json', dict(locator='https://doi.org/10.1016/j.heliyon.2024.e28130 Table1', path=SOURCE, sha256=sha(SOURCE), verification='Retrospective transcription, table inspected before this lock; not a prospective prediction', compared_quantity='Repeated same-tooth contact count and relative T-Scan signal, not absolute axial force'))
    if sha(SOURCE) != read(ROOT / 'EXTERNAL_SOURCE_LOCK.json')['sha256']:
        raise ValueError('External source drift')
    g = source_table()
    teeth = list(map(int, g[1][4:]))
    rows = [dict(day=x[0], scale=x[1], method=x[2], time=x[3], values=dict(zip(teeth, x[4:]))) for x in g[2:]]
    kept = []
    for time in ['9 a.m.', '16 p.m.']:
        cr = next((r for r in rows if r['day'] == 'Initial' and r['method'] == 'Silicone (<20 μm)' and (r['time'] == time)))
        fr = next((r for r in rows if r['day'] == 'Initial' and r['method'] == 'T-Scan' and (r['time'] == time)))
        kept.append(dict(tooth=16, day='Initial', time=time, count=int(cr['values'][16]), TScan_relative_signal_percent=float(fr['values'][16]), count_locator='Table1 Silicone(<20um) initial ' + time + ' tooth16', signal_locator='Table1 T-Scan initial ' + time + ' tooth16'))
    identity = kept[0]['count'] - kept[1]['count']
    difference = float(abs(Decimal(str(kept[0]['TScan_relative_signal_percent'])) - Decimal(str(kept[1]['TScan_relative_signal_percent']))))
    result = dict(external_referent=dict(kind='independent_measurement', locator='https://doi.org/10.1016/j.heliyon.2024.e28130 Table1', local_locator=SOURCE, sha256=sha(SOURCE), compared_quantity='Same-tooth repeated silicone contact counts and relative T-Scan force signal', refutes_us=True), selected_cells=kept, contact_count_identity_error=identity, relative_force_difference_percentage_points=difference, resolution='PER_TOOTH', absolute_region_force_admitted=0, missing_force_calibration='UNKNOWN: relative signals/counts are not absolute regional reactions', fault_plus1_signal_rejected=kept[0]['TScan_relative_signal_percent'] + 1 != float(next((r for r in rows if r['day'] == 'Initial' and r['method'] == 'T-Scan' and (r['time'] == '9 a.m.')))['values'][16]), attrition=dict(source_cells=len(rows) * 16, selected_paired_contact_and_signal_cells=4, not_selected_cells=len(rows) * 16 - 4, admitted_absolute_force_cells=0, reason='Bounded explicit counterexample; unselected source cells are not missing or new independent validation'))
    dump(ROOT / 'raw/EXTERNAL_MEASUREMENT.json', result)
    return result
if __name__ == '__main__':
    print(json.dumps(clean(run())))
