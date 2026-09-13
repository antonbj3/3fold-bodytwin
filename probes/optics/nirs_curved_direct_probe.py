"""Compare fixed-path wavelength prediction with independent full transport calls."""
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
import warp as wp
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
from bodytwin.geometry.layered_cylinder_v1 import layered_cylinder
from bodytwin.geometry.optics.tissue_photon_paths_v1 import simulate,PACKET
from bodytwin.geometry.optics.nirs_forward_v1 import absorption,sao2_hill,LAUNCHED


def main():
    wp.init();m=layered_cylinder();paths=np.load(ROOT/'reports/photon_curved_detector_paths_l4.npz')['paths_0'];cases=[];ratios=[]
    for po2 in (40.,60.,80.,100.):
        mua,_=absorption(float(sao2_hill(po2)));intensities=[]
        for wavelength,row in zip((760,850),mua):
            prop=np.array([[0,0,1,1],[row[0],1,.01,1.37],[row[1],2,.8,1.4]],dtype=np.float32)
            predicted=np.floor(PACKET*np.exp(-(paths@prop[1:,0].astype(float)))+.5).astype(np.int64)
            predicted_integer=int(predicted.sum());legs=[]
            for _ in range(2):
                out=simulate(m['vertices'],m['faces'],prop,m['front'],m['back'],source=[32,32,2.0001],direction=[0,0,1],
                    initial_region=1,origin=[0,0,0],shape=[64,64,124],reflect=True)
                escaped=out['terminal'][:,1]>0
                selected=escaped&(out['exits'][:,6]<=5.)&(np.sum((out['exits'][:,:3]-[42,32,2])**2,axis=1)<=4.)
                integer=int(out['terminal'][selected,1].sum());sums=out['terminal'][:,:3].sum(axis=0,dtype=np.int64)
                legs.append({'hashes':{k:hashlib.sha256(a.tobytes()).hexdigest() for k,a in out.items()},
                    'detected_photons':int(selected.sum()),'detected_integer':integer,'prediction_integer':predicted_integer,
                    'prediction_error_integer':abs(integer-predicted_integer),'absorbed_integer':int(sums[0]),'escaped_integer':int(sums[1]),
                    'residual_integer':int(sums[2]),'leaks':int(out['terminal'][:,3].sum()),'event_caps':int(out['counters'][2]),
                    'finite_nonnegative':bool(all(np.isfinite(out[k]).all() and np.min(out[k])>=0 for k in ('absorption','terminal','paths')))})
            cases.append({'po2':po2,'wavelength_nm':wavelength,'stored_properties':prop.astype(float).tolist(),'legs':legs})
            intensities.append(legs[0]['detected_integer']/(LAUNCHED*PACKET))
            print(json.dumps({'po2':po2,'wavelength':wavelength,'prediction_error_integer':legs[0]['prediction_error_integer'],
                              'repeat_exact':legs[0]==legs[1]}),flush=True)
        ratios.append({'po2':po2,'intensities':intensities,'ratio_760_850':intensities[0]/intensities[1]})
    gates={'all_full_arrays_repeat_exact':all(c['legs'][0]==c['legs'][1] for c in cases),
        'exact_energy':all(r['absorbed_integer']+r['escaped_integer']==LAUNCHED*PACKET for c in cases for r in c['legs']),
        'zero_failed_packets':all(r['residual_integer']==r['leaks']==r['event_caps']==0 for c in cases for r in c['legs']),
        'finite_nonnegative':all(r['finite_nonnegative'] for c in cases for r in c['legs']),
        'path_prediction_integer_bound':all(r['prediction_error_integer']<=3590 for c in cases for r in c['legs']),
        'ratio_increases':all(a['ratio_760_850']<b['ratio_760_850'] for a,b in zip(ratios,ratios[1:]))}
    report={'cases':cases,'ratios':ratios,'gates':gates,'integer_bound':3590,'scope':'synthetic fixed-scattering wavelength transport; clinical slope calibration remains unverified'}
    (ROOT/'reports/nirs_curved_direct_transport.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report));return 0 if all(gates.values()) else 1


if __name__=='__main__':raise SystemExit(main())
