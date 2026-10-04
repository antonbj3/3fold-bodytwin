"""Exact decimal rounding enclosures for paired group temperature means."""
from dental_release.paths import expand as _release_expand
from pathlib import Path
from decimal import Decimal
from lxml import etree
import json, re
P = Path(__file__).resolve().parent
S = Path(__import__('os').environ.get('XFACIT_SOURCE_DIR', _release_expand('@DENTAL_WORK_ROOT@/XFACIT-harvest/sources')))
root = etree.parse(str(S / 'PMC12306102.xml'))

def first(cell):
    text = ' '.join(''.join(cell.itertext()).split())
    return re.findall('(?<![\\w.])[-+]?\\d+(?:\\.\\d+)?', text)[0]

def balance(initial, final, delta):
    vals = [Decimal(x) for x in [initial, final, delta]]
    halfs = [Decimal(5).scaleb(v.as_tuple().exponent - 1) for v in vals]
    residual = vals[1] - vals[0] - vals[2]
    radius = sum(halfs)
    return {'initial': initial, 'final': final, 'reported_delta': delta, 'residual_K': str(residual), 'rounding_residual_interval_K': [str(residual - radius), str(residual + radius)], 'consistent_with_same_paired_observations': abs(residual) <= radius, 'resolution_level': 'POPULATION', 'assumption': 'All three means use the same paired sample; nearest printed rounding.'}
t = root.xpath('.//table-wrap[@id="Tab2"]')[0]
rows = t.xpath('.//tr')
checks = []
for (j, group) in enumerate(['A', 'B', 'C', 'D']):
    col = 3 if j == 0 else 2
    selected = [(base + j, col) for base in [5, 9, 13]]
    printed = [first(list(rows[r])[c]).replace(' ', '') for (r, c) in selected]
    printed = []
    for (r, c) in selected:
        text = ''.join(''.join(list(rows[r])[c].itertext()).split())
        printed.append(re.findall('[-+]?\\d+(?:\\.\\d+)?', text)[0])
    check = balance(*printed)
    check['group'] = group
    check['selector'] = {'table_id': 'Tab2', 'raw_rows_cols': selected}
    checks.append(check)
assert len(checks) == 4
control_good = balance('20.00', '23.00', '3.00')
control_bad = balance('20.00', '23.00', '999.00')
out = {'source': 'https://doi.org/10.1186/s12903-025-06593-z', 'locator': 'Table2[Tab2], paired group means', 'checks': checks, 'quarantine': any((not c['consistent_with_same_paired_observations'] for c in checks)), 'fault_control': {'consistent_triple_accepted': control_good['consistent_with_same_paired_observations'], 'injected_999_delta_rejected': not control_bad['consistent_with_same_paired_observations']}, 'scope': 'Source-table consistency only; no inference about apparatus bias or physical heat model.', 'physical_rigorous_enclosure': 'UNKNOWN; enclosure covers printing-rounding arithmetic only.'}
(P / 'SOURCE_BALANCE.json').write_text(json.dumps(out, indent=2) + '\n')
assert out['quarantine'] and all(out['fault_control'].values())
print('paired source balance:', checks)
