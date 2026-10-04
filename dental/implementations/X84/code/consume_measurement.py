"""Run: python3 code/consume_measurement.py measurement.json [--output path]."""
from pathlib import Path
import argparse, json
from patient_geometry import P, sha, write
from load_port import bind, fixed_load_decision
p = argparse.ArgumentParser()
p.add_argument('measurement')
p.add_argument('--output', default=str(P / 'exports/NEW_MEASURED_LOAD_CASE.json'))
a = p.parse_args()
q = json.loads(Path(a.measurement).read_text())
inventory = sorted((v['fdi'] for v in json.loads((P / 'raw/SOURCE_MANIFEST.json').read_text())['stls']))
out = bind(q, 'LiuHao2023Demo1', sha(P / 'raw/SOURCE_MANIFEST.json'), expected_source_fdi=inventory)
if out.get('forces'):
    row = next((r for r in out['forces'] if r['fdi'] == 36))
    out['tooth36_fixed100N_query'] = fixed_load_decision(row['force_interval_N'])
else:
    out['tooth36_fixed100N_query'] = fixed_load_decision(None)
write(Path(a.output), out)
print(json.dumps(out, indent=2))
