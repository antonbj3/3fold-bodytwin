"""Measure static versus normalized modulation observables on frozen paths."""
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
from bodytwin.geometry.optics.nirs_forward_v1 import signal,BLOOD_FRACTION,sao2_hill


def measure():
    paths=np.load(ROOT/'reports/photon_detector_paths_l4.npz')['paths_0']
    rows=[];hashes=[]
    for po2 in (40.,60.,80.,100.):
        sat=float(sao2_hill(po2))
        base,w=signal(paths,sat)
        low,wl=signal(paths,sat,blood_fraction=BLOOD_FRACTION*.99)
        high,wh=signal(paths,sat,blood_fraction=BLOOD_FRACTION*1.01)
        hashes.append(hashlib.sha256(np.stack((w,wl,wh)).tobytes()).hexdigest())
        original=None
        for gain in (.5,1.,2.):
            g=np.array([gain,1.])
            a=np.array(low['intensities'])*g;b=np.array(high['intensities'])*g
            modulation=(a-b)/((a+b)/2)
            intensity=np.array(base['intensities'])*g
            row=dict(po2=po2,sao2=sat,gain760=gain,static_ratio=float(intensity[0]/intensity[1]),
                     expected_static_ratio=base['ratio_760_850']*gain,
                     normalized_modulation_ratio=float(modulation[0]/modulation[1]),
                     modulation=modulation.tolist())
            rows.append(row)
    return rows,hashes


if __name__=='__main__':
    a,ha=measure();b,hb=measure()
    gates=dict(full_repeat=a==b and ha==hb,
               static_scales=all(r['static_ratio']==r['expected_static_ratio'] for r in a),
               modulation_invariant=all(len({r['normalized_modulation_ratio'] for r in a if r['po2']==p})==1 for p in (40.,60.,80.,100.)),
               matching_calibration_anchor=False)
    report=dict(rows=a,weight_hashes=ha,gates=gates,accepted_matching_anchors=0,
                clinical_calibration_certified=False,scope='Synthetic gain/observable counterexample; no validated pulse physiology.')
    (ROOT/'reports/nirs_observable_contract.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report));raise SystemExit(0 if all(gates.values()) else 2)
