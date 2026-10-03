#!/usr/bin/env python3
"""Only dermal new area changes; all layers contribute reversible release."""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[key]='2'
import time
import json
import resource
import numpy as np
from layered_skin_r2 import Layers, write_json, LANE

def trajectory(nx, coupling=2e7):
    model=Layers(nx,nx//2,0,coupling)
    rows=[]; U=None
    for a in (0,.008,.012,.016,.020):
        U,g,row=model.solve([.020,a,0],initial=U)
        with (LANE/f'r2/dermis_n{nx}_k{int(coupling)}_a{int(a*1000)}.npz').open('xb') as f:
            np.savez_compressed(f,U=U,gap=g,x=model.x)
        if rows:
            prev=rows[-1]; da=a-prev['lengths_m'][1]
            release=(prev['energy_J']-row['energy_J'])/da
            R=5e-6; q=2*R/.005
            prework=15000*(2*q-q*q)
            row['all_layer_release_force_N']=release
            row['release_equivalent_J_m2']=release/.0015
            row['Gamma0_dermis_only_J_m2']=150
            row['conditional_surface_minus_release_force_N']=.0015*150-release
            row['conditional_plus_bridge_force_N']=.0015*(150+prework)-release
            row['empirical_missing_other_process_and_friction_N']=None
            row['accumulated_dermis_surface_work_J']=.0015*150*a
            row['energy_balance_closure_error_J']=abs(row['conditional_surface_minus_release_force_N']*da-(row['energy_J']-prev['energy_J'])-.0015*150*da)
        rows.append(row)
    return rows

def main():
    start=time.process_time(); timings={}; results={}
    for nx in (64,96):
        t=time.process_time(); results[str(nx)]=trajectory(nx); timings[f'trajectory_n{nx}']=time.process_time()-t
    t=time.process_time(); results['decoupled64']=trajectory(64,0); timings['decoupled64']=time.process_time()-t
    t=time.process_time(); model=Layers(128,64,0)
    U,g,epi128=model.solve([.020,0,0]); timings['epi128']=time.process_time()-t
    with (LANE/'r2/epidermis_only_n128.npz').open('xb') as f:
        np.savez_compressed(f,U=U,gap=g,x=model.x)
    # An unstrained incision has no prestress-release/opening in this instrument.
    _,_,no_pre=Layers(32,16,0).solve([.020]*3,stretch=1.)
    mesh=[]
    for a,b in zip(results['64'],results['96']):
        mesh.append({'a_m':b['lengths_m'][1],
                     'gap_difference_relative':[abs(x-y)/max(abs(y),1e-20) for x,y in zip(a['max_gap_m'],b['max_gap_m'])],
                     'release_difference_J_m2':abs(a['release_equivalent_J_m2']-b['release_equivalent_J_m2']) if 'release_equivalent_J_m2' in a else None})
    epi96=results['96'][0]['max_gap_m'][0]
    difference=abs(epi128['max_gap_m'][0]-epi96)/epi128['max_gap_m'][0]
    result={'scope':'SYNTHETIC_FINITE_STRAIN_LAYER_ENERGY_ACQUISITION',
            'coupled_dermis_trajectories':results,'epidermis_only_n128':epi128,
            'no_prestretch_check':no_pre,'dermis_64_to96_refinement':mesh,
            'epidermis_96_to128_relative_gap_change':difference,
            'epidermis_refinement_gate':'PASS_DISCRETE_REFINEMENT' if difference<=.02 else 'FAIL',
            'max_energy_identity_error_J':max(row.get('energy_balance_closure_error_J',0) for rows in results.values() for row in rows),
            'release_remaining_synthetic_precision_budget_J_m2':14.875,
            'release_error_bound_scope':'Two-mesh difference is a sensitivity estimate, not a continuum bound. Empirical release error unknown.',
            'same_information_FE_control':'TIE by identical equations; no measured runtime advantage',
            'CPU_by_step_s':timings,'CPU_total_s':time.process_time()-start,
            'peak_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    write_json(LANE/'r2/layer_extension_results.json',result)
    print(json.dumps({'CPU_s':result['CPU_total_s'],'peak_RSS_KiB':result['peak_RSS_KiB'],'epi_refinement':difference,'release_refinement':mesh}))

if __name__=='__main__':main()
