"""Conditional synthetic oxygen/Hb/tissue-to-detector uncertainty chain.

Input bands are declared engineering assumptions. Geometry, extinction and
finite-photon library uncertainty are excluded; clinical calibration is absent.
"""
import hashlib
import itertools
import json
from pathlib import Path
import sys
import numpy as np
ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT/'src'))
from bodytwin.geometry.optics.nirs_forward_v1 import signal,sao2_hill,BLOOD_FRACTION
NOMINAL=np.array([60.,15.,1.]);BANDS=np.array([2.,.75,.1]);SEED=20260912;DRAWS=256


def propagate(paths,inputs):
    records=[];cache={}
    for po2,hb,fraction in inputs:
        sat=float(sao2_hill(po2));blend=BLOOD_FRACTION*fraction
        key=(sat,float(hb),np.asarray(blend).tobytes());out=cache.get(key)
        if out is None:
            out,_=signal(paths,sat,hb,blend);cache[key]=out
        records.append([po2,hb,fraction,sat,*np.array(out['mua_per_mm']).ravel(),*out['intensities'],out['ratio_760_850']])
    return np.asarray(records)


def leg(paths):
    random=np.random.default_rng(SEED).uniform(-1.,1.,(DRAWS,3))
    joint=propagate(paths,NOMINAL+BANDS*random);single=[];hashes={}
    for axis in range(3):
        offsets=np.zeros((DRAWS,3));offsets[:,axis]=BANDS[axis]*random[:,axis]
        values=propagate(paths,NOMINAL+offsets);single.append(float(np.ptp(values[:,-1])))
        hashes[f'isolated_{axis}']=hashlib.sha256(values.tobytes()).hexdigest()
    nominal=propagate(paths,NOMINAL[None,:]);zero=propagate(paths,NOMINAL+0.*random)
    grid=propagate(paths,np.array([NOMINAL+BANDS*np.array(x) for x in itertools.product((-1,0,1),repeat=3)]))
    grid_single=[]
    for axis in range(3):
        values=[]
        for side in (-1,0,1):
            x=NOMINAL.copy();x[axis]+=side*BANDS[axis];values.append(x)
        grid_single.append(float(np.ptp(propagate(paths,np.array(values))[:,-1])))
    blood_null=signal(paths,float(sao2_hill(60.)),15.,np.zeros(2))[0]['ratio_760_850']
    oxygen_low=signal(paths,float(sao2_hill(58.)))[0]['ratio_760_850']
    oxygen_high=signal(paths,float(sao2_hill(62.)))[0]['ratio_760_850']
    hashes.update(joint=hashlib.sha256(joint.tobytes()).hexdigest(),grid=hashlib.sha256(grid.tobytes()).hexdigest(),
                  zero=hashlib.sha256(zero.tobytes()).hexdigest())
    return {'hashes':hashes,'nominal_ratio':float(nominal[0,-1]),'joint_span':float(np.ptp(joint[:,-1])),
        'isolated_spans':single,'sum_isolated_spans':sum(single),'sample_interval_95':np.quantile(joint[:,-1],[.025,.975]).tolist(),
        'grid_span':float(np.ptp(grid[:,-1])),'grid_isolated_span_sum':sum(grid_single),
        'random_within_grid':bool(joint[:,-1].min()>=grid[:,-1].min() and joint[:,-1].max()<=grid[:,-1].max()),
        'finite_physical':bool(np.isfinite(joint).all() and np.all(joint[:,3]>0) and np.all(joint[:,3]<1) and np.all(joint[:,8:10]>0) and np.all(joint[:,8:10]<1)),
        'zero_band_exact':bool(np.all(zero==nominal)),'zero_blood_ratio':float(blood_null),
        'oxygen_direction':bool(oxygen_low<nominal[0,-1]<oxygen_high)}


def main():
    paths=np.load(ROOT/'reports/photon_detector_paths_l4.npz')['paths_0'];legs=[leg(paths),leg(paths)]
    gates={'full_chain_repeat_exact':legs[0]==legs[1],'finite_physical':all(r['finite_physical'] for r in legs),
        'zero_band':all(r['zero_band_exact'] for r in legs),'zero_blood':all(r['zero_blood_ratio']==1. for r in legs),
        'oxygen_direction':all(r['oxygen_direction'] for r in legs),
        'sample_composition':all(r['joint_span']<=r['sum_isolated_spans'] for r in legs),
        'grid_composition':all(r['grid_span']<=r['grid_isolated_span_sum'] for r in legs),
        'samples_inside_measured_grid':all(r['random_within_grid'] for r in legs)}
    report={'legs':legs,'gates':gates,'seed':SEED,'draws':DRAWS,'nominal_inputs':NOMINAL.tolist(),'uniform_half_bands':BANDS.tolist(),
        'clinical_calibration_certified':False,'scope':'conditional synthetic three-input propagation; excludes geometry/extinction/finite-photon uncertainty and clinical calibration'}
    (ROOT/'reports/nirs_chain_v1.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report));return 0 if all(gates.values()) else 1


if __name__=='__main__':raise SystemExit(main())
