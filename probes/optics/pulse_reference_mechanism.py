"""Inventory pinned public MC data and check a conventional inverse expression."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from scipy.io import loadmat
PINS={'demo_mc_data.mat':'49ae7d8e41c5e96464575a7d5a578bc08e03933c7d11815e8a4e2fb56b6242c6',
      'ext_coef.mat':'5abb6f2906a45ec1dc77f6e20bd8893847197293564bc230d3dadcaf58b9920a',
      'self_calibrated_algorithm_demo.m':'0c30c1822408b972fa7e7c0fe2e5fbf3224e4f8881d325fc747439aeef9dd63b'}


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('source_dir',type=Path);p.add_argument('--output',required=True,type=Path);a=p.parse_args()
    for name,h in PINS.items():
        if hashlib.sha256((a.source_dir/name).read_bytes()).hexdigest()!=h:raise ValueError('Source hash mismatch')
    data=loadmat(a.source_dir/'demo_mc_data.mat');coeff=loadmat(a.source_dir/'ext_coef.mat')
    intensity=data['I_r_all'];od=np.log(intensity[:,:,:,0]/intensity[:,:,:,1])
    wavelengths=coeff['wavelength'].ravel();idx=[int(np.flatnonzero(wavelengths==w)[0]) for w in (760,840)]
    e=np.column_stack((coeff['ext_hbo'].ravel()[idx],coeff['ext_hb'].ravel()[idx]))*2.3026/1000
    legs=[];arrays={}
    for leg in (0,1):
        rows=[]
        for sat in (.2,.4,.6,.8,1.):
            mixture=e[:,0]*sat+e[:,1]*(1-sat);ratio=mixture[0]/mixture[1]
            denominator=(e[0,0]-e[0,1])-ratio*(e[1,0]-e[1,1])
            demo=(e[1,0]*ratio-e[0,1])/denominator
            derived=(e[1,1]*ratio-e[0,1])/denominator
            rows.append(dict(saturation=sat,ratio=float(ratio),demo_prediction=float(demo),derived_prediction=float(derived),demo_error=float(abs(demo-sat)),derived_error=float(abs(derived-sat))))
        legs.append(rows);arrays[f'pulse_od_{leg}']=od.copy()
    gates=dict(full_repeat=legs[0]==legs[1] and arrays['pulse_od_0'].tobytes()==arrays['pulse_od_1'].tobytes(),
               positive_input=bool(np.isfinite(intensity).all() and np.min(intensity)>0 and np.isfinite(od).all() and np.min(od)>0),
               derived_consistency=all(r['derived_error']<=1e-12 for r in legs[0]),demo_consistency=all(r['demo_error']<=1e-12 for r in legs[0]))
    result=dict(legs=legs,gates=gates,pins=PINS,intensity_shape=list(intensity.shape),intensity_count=int(intensity.size),
                wavelengths=data['wavelengths_all'].ravel().tolist(),saturations=data['oxy_sim_all'].ravel().tolist(),distances=data['r_all'].ravel().tolist(),
                pulse_od_min=float(od.min()),pulse_od_max=float(od.max()),pulse_od_hash=hashlib.sha256(od.tobytes()).hexdigest(),
                source='https://doi.org/10.1184/R1/24530353.v1',scope='Public simulated data inventory and literal conventional-expression control; not clinical calibration or a test of the separate self-calibrated estimator.')
    a.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result));return 0 if all(gates.values()) else 2


if __name__=='__main__':raise SystemExit(main())
