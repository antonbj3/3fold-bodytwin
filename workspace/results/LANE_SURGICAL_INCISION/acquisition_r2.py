#!/usr/bin/env python3
"""Energy-balance identifiability and blade-arrival bridge-work instrument.
All generated observations are synthetic and explicitly named as such.
"""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[key]='2'
import json
import time
from pathlib import Path
import numpy as np
from scipy.optimize import linprog
from layered_skin_r2 import write_json, LANE

def interval(X,y,noise):
    c=np.array([1.,0,0,0])
    A=np.vstack((X,-X)); b=np.r_[y+noise,-y+noise]
    low=linprog(c,A_ub=A,b_ub=b,bounds=[(0,None)]*4,method='highs')
    high=linprog(-c,A_ub=A,b_ub=b,bounds=[(0,None)]*4,method='highs')
    if not low.success or not high.success:
        raise RuntimeError((low.message,high.message))
    return {'Gamma0_min_J_m2':float(low.fun),'Gamma0_max_J_m2':float(-high.fun),
            'halfwidth_J_m2':float((-high.fun-low.fun)/2),'rank':int(np.linalg.matrix_rank(X)),
            'min_witness':low.x.tolist(),'max_witness':high.x.tolist(),'LP_solves':2,
            'gate':'PASS_SYNTHETIC' if (-high.fun-low.fun)/2<=30+1e-8 else 'FAIL'}

def main():
    start=time.process_time()
    p=json.loads((LANE/'PREREG_R2_ACQUISITION.json').read_text())
    targets=json.loads((LANE/'SOURCE_TARGETS_R2.json').read_text())
    tear=targets['primary_world_targets']['PISSARENKO_2020']
    world={'review_state':'PENDING_INDEPENDENT_REVIEW','unfitted_Gtear_candidate_J_m2':15000,
           'mode_III_relative_error':abs(15000-20600)/20600,
           'mode_I_relative_error':abs(15000-30380)/30380,
           'world_tear_gate':'FAIL','skin_intrinsic_cut_intercept_gate':'UNKNOWN',
           'paired_skin_cut_tear_ratio':None,'radius_transition_world_gate':'UNKNOWN',
           'strongest_per_observable_calibrated_control':{'mode_III_fit_J_m2':20600,'mode_I_fit_J_m2':30380,'training_error':0.,'heldout_error':None},
           'force_thresholds_N':{'smooth_grade3_lower_bound':1.9*9.80665,'serrated_grade2_or3':.7*9.80665},
           'force_threshold_ratio_gate':'UNKNOWN_NONMATCHED_ENDPOINTS',
           'table_reconstruction':[]}
    for mode,orient,h,F,J in tear['rows']:
        reconstructed=2*F/(h*1e-3)
        world['table_reconstruction'].append({'mode':mode,'orientation':orient,'h_mm':h,'F_N':F,
              'reported_J_J_m2':J*1000,'derived_2F_over_h_J_m2':reconstructed,'relative_difference':abs(reconstructed-J*1000)/(J*1000)})
    world['max_table_reconstruction_relative_difference']=max(r['relative_difference'] for r in world['table_reconstruction'])
    write_json(LANE/'r2/world_crossprediction.json',world)

    phi,tau,l,rf=.3,1e4,.005,5e-6
    Gp=phi*tau*l*l/rf; Tp=2*phi*tau*l/rf
    R=np.r_[0.,rf,2*rf,np.geomspace(1e-7,l/2,100)]
    delta=2*R
    q=np.minimum(delta/l,1)
    before=Gp*(2*q-q*q)
    cut=150+before
    # Analytic bridge-work integral checked separately by numerical quadrature.
    integration=[]
    for d in (2*rf,4*rf,.001,l):
        grid=np.linspace(0,d,10001)
        numeric=np.trapezoid(Tp*np.maximum(1-grid/l,0),grid)
        expected=Gp*(2*min(d/l,1)-min(d/l,1)**2)
        integration.append(abs(numeric-expected)/max(expected,1e-12))
    bridge={'scope':'DERIVED_UNDER_MONOTONE_ALIGNED_BRIDGE_AND_SYNTHETIC_DELTA_EQUALS_2R_CLOSURE',
            'Gp_J_m2':Gp,'matrix_Gamma0_J_m2':150,'tear_total_if_matrix_work_added_J_m2':Gp+150,
            'peak_traction_Pa':Tp,'fiber_peak_stress_Pa':Tp/phi,'strength_limit_Pa':1e8,
            'R_equal_fiber_radius':{'R_m':rf,'cut_work_J_m2':float(cut[1]),'tear_to_cut_ratio':float((Gp+150)/cut[1])},
            'R_equal_fiber_diameter':{'R_m':2*rf,'cut_work_J_m2':float(cut[2]),'tear_to_cut_ratio':float((Gp+150)/cut[2])},
            'R_for_half_pullout_work_m':float(l/2*(1-np.sqrt(.5))),
            'half_work_R_over_fiber_diameter':float(l/2*(1-np.sqrt(.5))/(2*rf)),
            'quadrature_max_relative_error':max(integration),
            'interpretation':'Fiber diameter alone does not set the transition in this explicit bridge construction. Anchorage length and actual opening at blade arrival enter. This is not a fitted or measured skin transition.',
            'intrinsic_intercept_warning':'R=0 removes radius-dependent work only if pre-cleavage process work also vanishes. A constant process term survives the radius-zero limit.'}
    write_json(LANE/'r2/bridge_work.json',bridge)
    with (LANE/'r2/bridge_work_raw.npz').open('xb') as f:
        np.savez(f,radius_m=R,delta_m=delta,precleavage_work_J_m2=before,total_cut_work_J_m2=cut)

    radii=np.array(p['radii_m']); z=radii/20e-6
    q=2*radii/l; shape=2*q-q*q
    X=np.column_stack((np.ones(5),np.ones(5),z,shape))
    theta=np.array([150.,900.,60.,Gp])
    y=X@theta
    null=np.array([1.,-1.,0,0])
    base=interval(X,y,np.full(5,10.))
    # Same-path recut isolates contact, but not first-pass constant process work.
    Xcontact=np.column_stack((np.zeros(5),np.zeros(5),z,np.zeros(5)))
    ycontact=Xcontact@theta
    X2=np.vstack((X,Xcontact)); y2=np.r_[y,ycontact]; e2=np.r_[np.full(5,10.),np.full(5,5.)]
    contact=interval(X2,y2,e2)
    candidates={}
    lp=4
    for worknoise in (20.,5.):
        rows=[]
        for i,r in enumerate(radii):
            workrow=np.array([0.,1.,0.,shape[i]])
            X3=np.vstack((X2,workrow)); y3=np.r_[y2,workrow@theta]; e3=np.r_[e2,worknoise]
            record=interval(X3,y3,e3); record['radius_m']=float(r); record['work_error_bound_J_m2']=worknoise
            rows.append(record); lp+=2
        candidates[str(int(worknoise))]=rows
    best20=min(candidates['20'],key=lambda r:r['halfwidth_J_m2'])
    best5=min(candidates['5'],key=lambda r:r['halfwidth_J_m2'])
    result={'scope':'SYNTHETIC_ACQUISITION_IDENTIFIABILITY_NOT_EMPIRICAL_DATA',
            'parameter_columns':['Gamma0','constant_process_work','contact_amplitude','bridge_work_amplitude'],
            'radius_series_only':base,'radius_series_plus_contact':contact,
            'radius_null_direction':null.tolist(),'null_direction_max_effect':float(np.abs(X@null).max()),
            'independent_bridge_work_20J':best20,'refined_independent_bridge_work_5J':best5,
            'all_candidate_designs':candidates,'LP_solves':lp,
            'same_information_control':'Identical positive LP inverse: TIE, cost/gain 1x by construction; no separate speedup claimed',
            'required_measurement':'Integrate independently observed bridge traction against measured local opening before cleavage, on same skin/sample and loading path; bound viscoelastic loss and prestress-release mismatch. Recut alone is insufficient.',
            'worst_case_energy_error_budget_J_m2':{'cutting':10,'contact_at_best_radius':float(5*radii.min()/20e-6),'bridge_work':5,'additional_release_or_model_error_remaining':None},
            'total_CPU_s':time.process_time()-start,
            'empirical_feasibility_of_acquisition_precision':'UNKNOWN'}
    write_json(LANE/'r2/acquisition_results.json',result)
    with (LANE/'r2/acquisition_synthetic_raw.npz').open('xb') as f:
        np.savez(f,radius_m=radii,X_cut=X,y_cut_SYNTHETIC=y,X_contact=Xcontact,y_contact_SYNTHETIC=ycontact,theta_SYNTHETIC=theta)
    print(json.dumps({'world_tear_errors':[world['mode_III_relative_error'],world['mode_I_relative_error']],
                      'radius_interval':base,'contact_interval':contact,'work20':best20,'work5':best5}))

if __name__=='__main__':main()
