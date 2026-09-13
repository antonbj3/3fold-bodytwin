"""Bound absorption termination before the fixed detector window."""
import json
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
PACKET=1<<30


def measure():
    paths=np.load(ROOT/'reports/photon_detector_paths_l4.npz')['paths_0']
    properties=np.array([[.005,1,.01,1.37],[.01,2,.8,1.4]],dtype=np.float32).astype(float)
    depth=paths@properties[:,0];time=paths@properties[:,3]/299.792458
    maximum_depth=5.*299.792458*np.max(properties[:,0]/properties[:,3])
    threshold=np.log(2.*PACKET)
    return {'detected_paths':len(paths),'window_ns':5.,'conservative_depth_bound':float(maximum_depth),
        'rounding_zero_depth':float(threshold),'minimum_continuous_weight_before_window':float(PACKET*np.exp(-maximum_depth)),
        'saved_max_depth':float(depth.max()),'saved_max_time_ns':float(time.max())}


def main():
    legs=[measure(),measure()]
    gates={'tables_repeat_exact':legs[0]==legs[1],
        'no_absorption_termination_before_window':all(r['conservative_depth_bound']<r['rounding_zero_depth'] for r in legs),
        'saved_paths_within_bound':all(r['saved_max_depth']<=r['conservative_depth_bound'] and r['saved_max_time_ns']<=5. for r in legs)}
    report={'legs':legs,'gates':gates,'scope':'path reuse for unchanged scattering/refraction/source/detector and5ns window only'}
    (ROOT/'reports/photon_path_reuse.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report));return 0 if all(gates.values()) else 1


if __name__=='__main__':raise SystemExit(main())
