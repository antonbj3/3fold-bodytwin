"""Synthetic two-wavelength detector forward model from certified MC paths.

Scattering, refraction, geometry and detector remain fixed across wavelengths.
This experimental model has no clinical pulse-oximetry slope certificate.
"""
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
ROOT=Path(__file__).resolve().parents[4];sys.path.insert(0,str(ROOT/'src'))
from bodytwin.cells.haematology.blood_oxygen_transport import sao2_hill,HB_G_DL

EXTINCTION=np.array([[586.,1548.52],[1058.,691.32]])
BLOOD_FRACTION=np.array([.02,.04])
BACKGROUND=np.array([.002,.004])
LAUNCHED=1000000


def absorption(saturation,hb_g_dl=HB_G_DL,blood_fraction=BLOOD_FRACTION):
    fraction=np.asarray(blood_fraction,dtype=float)
    if not np.isfinite([saturation,hb_g_dl]).all() or not 0<=saturation<=1 or hb_g_dl<=0:
        raise ValueError('Finite saturation in[0,1] and positive hemoglobin required')
    if fraction.shape!=(2,) or not np.isfinite(fraction).all() or np.any((fraction<0)|(fraction>1)):
        raise ValueError('Two blood fractions in[0,1] required')
    # g/dL->g/L and cm^-1->mm^-1 factors cancel; extinction is decadic per mol/L.
    scale=np.log(10.)*hb_g_dl/64500.
    blood=scale*(EXTINCTION[:,0]*saturation+EXTINCTION[:,1]*(1-saturation))
    derivative=scale*(EXTINCTION[:,0]-EXTINCTION[:,1])
    return BACKGROUND[None,:]+blood[:,None]*fraction[None,:],derivative[:,None]*fraction[None,:]


def signal(paths,saturation,hb_g_dl=HB_G_DL,blood_fraction=BLOOD_FRACTION):
    paths=np.asarray(paths,dtype=float)
    if paths.ndim!=2 or paths.shape[1]!=2 or len(paths)==0 or not np.isfinite(paths).all() or np.any(paths<0):
        raise ValueError('Finite nonnegative two-region detector paths required')
    mua,dmua=absorption(saturation,hb_g_dl,blood_fraction)
    weights=np.exp(-(paths@mua.T));total=weights.sum(axis=0);intensity=total/LAUNCHED
    mean_paths=(weights.T@paths)/total[:,None]
    dlog=-np.sum(mean_paths*dmua,axis=1)
    return {'sao2':float(saturation),'hb_g_dl':float(hb_g_dl),'mua_per_mm':mua.tolist(),
        'intensities':intensity.tolist(),'ratio_760_850':float(intensity[0]/intensity[1]),
        'd_log_ratio_d_sao2':float(dlog[0]-dlog[1])},weights


def main():
    paths=np.load(ROOT/'reports/photon_detector_paths_l4.npz')['paths_0'];legs=[];array_hashes=[]
    for _ in range(2):
        rows=[];arrays=[]
        for po2 in (40.,60.,80.,100.):
            sat=float(sao2_hill(po2));row,w=signal(paths,sat);step=1e-5
            plus=signal(paths,sat+step)[0]['ratio_760_850'];minus=signal(paths,sat-step)[0]['ratio_760_850']
            finite_difference=(np.log(plus)-np.log(minus))/(2*step)
            row.update(po2=po2,derivative_fd=float(finite_difference),
                derivative_relative_error=float(abs(finite_difference-row['d_log_ratio_d_sao2'])/abs(row['d_log_ratio_d_sao2'])))
            rows.append(row);arrays.append(w)
        legs.append(rows);array_hashes.append(hashlib.sha256(np.array(arrays).tobytes()).hexdigest())
    gates={'full_weights_and_results_repeat_exact':legs[0]==legs[1] and array_hashes[0]==array_hashes[1],
        'physical_intensities':all(np.isfinite(r['intensities']).all() and 0<min(r['intensities'])<=max(r['intensities'])<1 for r in legs[0]),
        'ratio_increases':all(a['ratio_760_850']<b['ratio_760_850'] for a,b in zip(legs[0],legs[0][1:])),
        'positive_analytic_slope':all(r['d_log_ratio_d_sao2']>0 for r in legs[0]),
        'derivative_identity':all(r['derivative_relative_error']<=1e-7 for r in legs[0])}
    report={'legs':legs,'weight_array_hashes':array_hashes,'engineering_gates':gates,'clinical_calibration_certified':False,
        'accepted_matching_clinical_calibration_anchors':0,'launched_photons':LAUNCHED,'detected_paths':len(paths),
        'blood_fractions':BLOOD_FRACTION.tolist(),'background_mua':BACKGROUND.tolist(),
        'extinction_source':'https://omlc.org/spectra/hemoglobin/summary.html',
        'oxygen_cell_sha256':hashlib.sha256((ROOT/'src/bodytwin/cells/haematology/blood_oxygen_transport.py').read_bytes()).hexdigest(),
        'scope':'synthetic fixed optical geometry and path library; direct transport comparison and clinical calibration are separate gates'}
    (ROOT/'reports/nirs_forward_v1.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report));return 0 if all(gates.values()) else 1


if __name__=='__main__':raise SystemExit(main())
