"""Read one recorded per-crown result; does not imply physical calibration."""
from common import *
import argparse
p = argparse.ArgumentParser()
p.add_argument('--uid')
args = p.parse_args()
rows = read(ROOT / 'raw/R2_ROWS.json')
loads = {r['uid']: r for r in read(ROOT / 'raw/R3_ROWS.json')}
if args.uid is None:
    print('\n'.join((r['uid'] for r in rows)))
else:
    found = next((r for r in rows if r['uid'] == args.uid), None)
    if found is None:
        raise SystemExit('Unknown uid; run without --uid to list')
    l = loads[args.uid]
    print(json.dumps(dict(contact=found, load_sets=l.get('classes'), robust_load_sets=l.get('robust_classes'), physical_force='UNKNOWN'), indent=2))
