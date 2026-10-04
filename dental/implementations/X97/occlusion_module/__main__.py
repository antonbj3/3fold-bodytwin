"""JSON CLI with optional relative NPZ array references."""
import argparse, json
from pathlib import Path
import numpy as np
from .api import analyze

def load_packet(path):
    path = Path(path)

    def resolve(x):
        if isinstance(x, dict) and set(x) == {'npz', 'key'}:
            p = (path.parent / x['npz']).resolve()
            with np.load(p, allow_pickle=False) as a:
                return a[x['key']].copy()
        if isinstance(x, dict):
            return {k: resolve(v) for (k, v) in x.items()}
        if isinstance(x, list):
            return [resolve(v) for v in x]
        return x
    return resolve(json.loads(path.read_text()))

def main():
    p = argparse.ArgumentParser(description='Registered bite + research crown -> occlusion report (mm, N).')
    p.add_argument('input', help='JSON packet, arrays inline or {npz,key} relative references')
    p.add_argument('--out', required=True)
    a = p.parse_args()
    packet = load_packet(a.input)
    out = analyze(packet['bite'], packet['crown'], packet.get('request'))
    dest = Path(a.out)
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(json.dumps(out, indent=2, allow_nan=False) + '\n')
    print(json.dumps({'status': out['status'], 'out': str(dest), 'physical_validation': out['physical_validation']}))
if __name__ == '__main__':
    main()
