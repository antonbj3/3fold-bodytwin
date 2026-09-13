"""Measure the existing oxygen cell to published extinction-coefficient seam."""
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'src'))
from bodytwin.cells.haematology.blood_oxygen_transport import sao2_hill,HB_G_DL

EXTINCTION=np.array([[586.,1548.52],[1058.,691.32]])


def measure():
    # Extinction columns are oxygenated/deoxygenated, in cm^-1 per mol/L.
    # g/dL -> g/L -> mol/L; cm^-1 -> mm^-1; decadic -> natural attenuation.
    concentration=HB_G_DL*10./64500.
    rows=[]
    for po2 in (40.,60.,80.,100.):
        sat=float(sao2_hill(po2))
        mua=np.log(10.)*concentration*(EXTINCTION[:,0]*sat+EXTINCTION[:,1]*(1.-sat))/10.
        derivative=np.log(10.)*concentration*(EXTINCTION[:,0]-EXTINCTION[:,1])/10.
        rows.append({'po2_mmhg':po2,'sao2':sat,'mua_760_per_mm':float(mua[0]),'mua_850_per_mm':float(mua[1]),
                     'd_mua_760_d_sao2':float(derivative[0]),'d_mua_850_d_sao2':float(derivative[1]),
                     'equal_path_log_intensity_ratio_slope_per_mm':float(derivative[1]-derivative[0])})
    return rows


def main():
    legs=[measure(),measure()]
    gates={'two_tables_identical':legs[0]==legs[1],
           'physical_ranges':all(0<r['sao2']<1 and r['mua_760_per_mm']>0 and r['mua_850_per_mm']>0 for r in legs[0]),
           'spectral_direction':all(r['d_mua_760_d_sao2']<0<r['d_mua_850_d_sao2'] and r['equal_path_log_intensity_ratio_slope_per_mm']>0 for r in legs[0])}
    source=ROOT/'src/bodytwin/cells/haematology/blood_oxygen_transport.py'
    report={'legs':legs,'gates':gates,'hb_g_dl':HB_G_DL,'extinction_cm_inverse_per_molar':EXTINCTION.tolist(),
            'oxygen_cell_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
            'source':'https://omlc.org/spectra/hemoglobin/summary.html',
            'scope':'whole-blood absorption arithmetic; equal-path analytic diagnostic, not tissue MC or pulse calibration'}
    (ROOT/'reports/hemoglobin_mechanism.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report),flush=True)
    return 0 if all(gates.values()) else 1


if __name__=='__main__':raise SystemExit(main())
