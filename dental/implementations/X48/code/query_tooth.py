"""One tooth/week query with moving-reference and uncertainty fields retained."""
from common import *
import argparse
p = argparse.ArgumentParser()
p.add_argument('--patient', choices=['A', 'B', '3485', '6457'], required=True)
p.add_argument('--tooth', type=int, required=True)
p.add_argument('--week', type=int, required=True)
p.add_argument('--interval', choices=['weekly', 'cumulative'], default='weekly')
a = p.parse_args()
patient = {'A': '3485', 'B': '6457'}.get(a.patient, a.patient)
x = [r for r in read(ROOT / 'raw/R4_TYPE_ENCLOSURES.json') if r['patient'] == patient and r['fdi'] == a.tooth and (r['week'] == a.week) and (r['interval'] == a.interval)]
if not x:
    p.error('No matching measured tooth/time; nothing imputed')
print(json.dumps(dict(patient=a.patient, fdi=a.tooth, week=a.week, interval=a.interval, claim_type='information_link', absolute_anatomical_motion='UNKNOWN', body_translation='UNKNOWN_ROOT_NOT_OBSERVED', relative_type_observations=x), indent=2))
