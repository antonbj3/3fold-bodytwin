"""Measure existing chain seams before designing a data-defined runner."""
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
from bodytwin.chains.nirs_chain_v1 import leg,propagate,NOMINAL
from bodytwin.geometry.optics.nirs_forward_v1 import sao2_hill,signal,BLOOD_FRACTION


def main():
    rows=[]
    for name,path,report in [('box','photon_detector_paths_l4.npz','nirs_chain_v1.json'),('curved','photon_curved_detector_paths_l4.npz','nirs_curved_chain_v1.json')]:
        paths=np.load(ROOT/'reports'/path)['paths_0'];legs=[leg(paths),leg(paths)]
        po2,hb,fraction=NOMINAL;sat=float(sao2_hill(po2));out,_=signal(paths,sat,hb,BLOOD_FRACTION*fraction)
        manual=np.asarray([[po2,hb,fraction,sat,*np.array(out['mua_per_mm']).ravel(),*out['intensities'],out['ratio_760_850']]])
        archived=json.loads((ROOT/'reports'/report).read_text())['legs'][0]
        rows.append(dict(name=name,path_shape=list(paths.shape),path_hash=hashlib.sha256(paths.tobytes()).hexdigest(),legs=legs,seam_shape=list(manual.shape),seam_exact=manual.tobytes()==propagate(paths,NOMINAL[None,:]).tobytes(),archived_exact=legs[0]==archived))
    gates=dict(repeat=all(r['legs'][0]==r['legs'][1] for r in rows),physical_nulls=all(r['legs'][0]['finite_physical'] and r['legs'][0]['zero_band_exact'] and r['legs'][0]['zero_blood_ratio']==1 for r in rows),seam_exact=all(r['seam_exact'] for r in rows))
    result=dict(rows=rows,gates=gates,scope='Synthetic fixed-path conditional chains; no clinical calibration.')
    (ROOT/'reports/nirs_spec_mechanism.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result));return 0 if all(gates.values()) else 2

if __name__=='__main__':raise SystemExit(main())
