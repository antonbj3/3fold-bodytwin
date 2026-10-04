from geometry import *
from compare_pressure import compare
import csv
r = json.loads((H / 'rounds/R3.json').read_text())['rows'][0]
f = json.loads((H / 'FROZEN_PREDICTIONS_R3.json').read_text())['files'][0]
checks = []
for name in ['wrong_hash', 'missing_channel', 'negative_force']:
    p = H / 'raw/failed_checks' / ('FAULT_PRESSURE_' + name + '.csv')
    rows = []
    for k in range(r['patch_count']):
        if name == 'missing_channel' and k == r['patch_count'] - 1:
            continue
        rows.append(dict(case=r['case'], fdi=r['fdi'], patch_id=k, force_N=-1 if name == 'negative_force' and k == 0 else 100 / r['patch_count'], basis_sha256='0' * 64 if name == 'wrong_hash' else f['sha256']))
    with p.open('w') as fh:
        w = csv.DictWriter(fh, fieldnames=['case', 'fdi', 'patch_id', 'force_N', 'basis_sha256'])
        w.writeheader()
        w.writerows(rows)
    rejected = False
    reason = ''
    try:
        compare(p)
    except (AssertionError, ValueError) as e:
        rejected = True
        reason = str(e)
    checks.append(dict(name=name, rejected=rejected, reason=reason, input_status='INJECTED_INVALID_SOFTWARE_FIXTURE, NOT a measurement'))
dump(H / 'raw/PRESSURE_PORT_VERIFICATION.json', dict(checks=checks, physical_measurements=0, pass_=all((c['rejected'] for c in checks))))
assert all((c['rejected'] for c in checks))
