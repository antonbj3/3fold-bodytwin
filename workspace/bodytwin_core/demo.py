"""Run the sourced 1–2-person the collaborator box-lift integration from the workspace root."""
import argparse
import json
import time
import numpy as np
from .config import paths
from . import geometry, solver, whatif, population, nullmodels

def timed(label, call, timings):
    start=time.perf_counter()
    value=call()
    timings[label]=round(time.perf_counter()-start,3)
    return value

def run(subjects=('z001',), config=None):
    if len(subjects) not in (1,2):raise ValueError('Choose one or two subjects')
    p=paths(config)
    n2=json.loads((p['n2b_results']/'n2b_results.json').read_text())
    band_source=n2['bands']['ens_s9p5_p2_lm_M1000']['hip_r']
    null_source=n2['nullmodels']
    for subject in subjects:
        tag='vsd_'+subject
        if tag not in ('vsd_z001','vsd_z009'):
            raise ValueError('the collaborator solver data supports z001 and z009 in this demo')
        timings={}
        geo=timed('instantiate_s',lambda:geometry.instantiate(subject,config=config),timings)
        force,receipts=timed('solve_141_frames_s',lambda:solver.solve_lift(tag,config),timings)
        curve=timed('hip_force_s',lambda:solver.hip_curve(tag,force,config),timings)
        t=int(np.argmax(curve));peak=float(curve[t])
        curves={param:whatif.reported_curve(tag,param,config) for param in ('CCD','AV')}
        for param,row in curves.items():
            if abs(float(row['full_peak_N']['0'])-peak)>0.05:
                raise RuntimeError(f'{tag} {param}: frozen what-if baseline differs from solved peak')
        whatif_peaks={param:{k:row['full_peak_N'][k] for k in ('-10','0','10') if k in row['full_peak_N']}
                      for param,row in curves.items()}
        derivatives={}
        for param in ('CCD','AV'):
            derivatives[param]=timed(f'{param}_implicit_derivative_s',
                                     lambda param=param:whatif.reported_derivative(tag,param,t,config),timings)
        band=timed('population_band_s',lambda:{
            'source':'N2b A363, 1000-individual sigma=9.5 mm attachment ensemble; reference-relative transfer',
            'lower_N':peak*(1+band_source['peak_p2_5_pct']/100),
            'upper_N':peak*(1+band_source['peak_p97_5_pct']/100),
            'lower_pct':band_source['peak_p2_5_pct'],
            'upper_pct':band_source['peak_p97_5_pct'],
            'N50_kappa':population.hip_band(peak,config)['kappa'],
            'N50_note':'individual N50 band needs posterior-draw SD; none supplied for this the collaborator individual'},timings)
        null=timed('nullmodels_s',lambda:{
            'B24_lifting_peak_N':nullmodels.b24_lifting(null_source['BW_N'],config),
            'N1_John_GRF_proxy_peak_N':{side:nullmodels.n1_knee(value,2.2) for side,value in null_source['GRF_peak_N'].items()},
            'N1g_external_knee_gait_RMSE_BW':0.378006,
            'N1g_John_lift_N':'unavailable: no activity-group coefficient calibrated on this lift',
            'scope_note':'B24 spine facets and N1/N1g knee gait do not measure this the collaborator hip-force truth'},timings)
        out={'subject':subject,'geometry_accepted':bool(geo.accepted),
             'geometry_registry_ids':len(geo.instance.reg.entities),'frames':len(force),
             'solver_kkt_pass_1e10':sum(r['kkt']<=1e-10 for r in receipts),
             'peak_frame':t,'hip_peak_N':peak,'whatif_peak_N_by_degree':whatif_peaks,'derivatives':derivatives,
             'whatif_full_peak_fd_N_per_deg':{param:whatif.reported_curve(tag,param,config)['central_fd_peak_N_per_deg'] for param in ('CCD','AV')},
             'derivative_note':'Archived implicit values are at the listed stable frame; full-peak finite differences are separate. Active-set changes prevent an implicit peak value where reported.',
             'band':band,
             'nullmodels':null,'timings_s':timings}
        print(json.dumps(out,indent=2,default=float),flush=True)

def main():
    a=argparse.ArgumentParser();a.add_argument('--subjects',nargs='+',default=['z001'],choices=['z001','z009']);a.add_argument('--config');a.add_argument('--determined',action='store_true')
    args=a.parse_args()
    if args.determined:
        from .determined_demo import run as determined_run
        determined_run()
    else:run(tuple(args.subjects),args.config)
if __name__=='__main__':main()
