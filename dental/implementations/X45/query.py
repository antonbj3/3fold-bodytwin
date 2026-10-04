"""Query final bounded chain ports after ./run_all.sh."""
import argparse, json
from pathlib import Path
import ports
ROOT = Path(__file__).resolve().parent
ap = argparse.ArgumentParser()
ap.add_argument('chain', choices=['healing', 'wear', 'coupon'])
ap.add_argument('--group', default='OD')
ap.add_argument('--day', type=float, default=21)
ap.add_argument('--material', default='Admira fusion flow')
ap.add_argument('--cycles', type=int, default=40000)
ap.add_argument('--grade', default='XT', choices=['T', 'ST', 'XT', 'P'])
ap.add_argument('--state', default='aged', choices=['aged', 'unaged'])
ap.add_argument('--probability', type=float, default=0.05)
a = ap.parse_args()
raw = json.loads((ROOT / 'raw/PRIMARY_MEASUREMENTS.json').read_text())
try:
    if a.chain == 'healing':
        knots = [p for p in raw['ISQ'][a.group] if p['day'] in [0, 14, 28, 90]]
        result = ports.interpolate_enclosed([p['day'] for p in knots], [p['mean_ISQ'] for p in knots], a.day)
        result.update(source='doi:10.1111/cid.13140, Table2', group=a.group, day=a.day, physical_BIC='UNKNOWN')
    elif a.chain == 'wear':
        knots = [p for p in raw['WEAR'][a.material] if p['cycles'] in [5000, 20000, 120000]]
        result = ports.wear_envelope([p['cycles'] for p in knots], [p['mean_mg'] for p in knots], a.cycles)
        result.update(source='doi:10.18502/fid.v20i10.12609, Table3', material=a.material, cycles=a.cycles, volume='UNKNOWN: density missing')
    else:
        state = json.loads((ROOT / 'MATERIAL_STATE_REGISTRY.json').read_text())[a.grade][a.state]
        result = ports.material_query(state, a.probability)
        result.update(source='doi:10.1055/s-0042-1755630, Table4', grade=a.grade, state=a.state, probability=a.probability)
    print(json.dumps(result, indent=2, ensure_ascii=False))
except (ports.PortError, KeyError) as e:
    print(json.dumps({'outcome': 'UNKNOWN_OR_REFUSED', 'reason': str(e)}, indent=2))
    raise SystemExit(2)
