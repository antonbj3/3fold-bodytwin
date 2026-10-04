"""One-command independent-measurement input; never substitutes simulated data."""
import argparse, json
from pathlib import Path
from calibration_port import fit

def main():
    a = argparse.ArgumentParser()
    a.add_argument('measurement_json')
    a.add_argument('--length-scale-mm', type=float, required=True)
    a.add_argument('--residual-tolerance-N', type=float, default=0.001)
    args = a.parse_args()
    rows = json.loads(Path(args.measurement_json).read_text())
    if not rows or any((r.get('source_kind') != 'independent_measurement' for r in rows)):
        print(json.dumps(dict(status='REQUIRES_INDEPENDENT_MEASUREMENT', detail='Do not send simulator fixture to this laboratory entry point')))
        return
    print(json.dumps(fit(rows, args.length_scale_mm, args.residual_tolerance_N), indent=2))
if __name__ == '__main__':
    main()
