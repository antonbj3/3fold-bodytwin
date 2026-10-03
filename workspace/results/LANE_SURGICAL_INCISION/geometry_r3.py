#!/usr/bin/env python3
"""Frozen published skin geometry tests. No fabricated cutting observations.

Normal pressure-width and additive nose-pressure are constitutive instruments,
not a fitted chemical cleavage model. All outputs refuse existing destinations.
"""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[key] = '2'
import argparse
import hashlib
import json
import math
import resource
import time
from pathlib import Path
import numpy as np
from scipy.optimize import brentq

LANE = Path(__file__).resolve().parent
LBF = 4.4482216152605

def write_json(path, value):
    with Path(path).open('x') as f:
        json.dump(value, f, indent=2, allow_nan=False)
        f.write('\n')

def halfwidth(depth, angle, radius):
    theta = math.radians(angle/2)
    switch = radius*(1-math.sin(theta))
    if radius > 0 and depth < switch:
        return math.sqrt(max(0, 2*radius*depth-depth*depth))
    return depth*math.tan(theta) + radius*math.cos(theta)/(1+math.sin(theta))

def interpolated_force(depth_pct, curve):
    d = np.array(curve['depth_pct'], dtype=float)
    f = np.array(curve['force_lbf'], dtype=float)
    # Full-thickness values do not establish a unique depth/force relation.
    valid = d < 100
    d, f = d[valid], f[valid]
    if not d[0] <= depth_pct <= d[-1]:
        return None
    return float(np.interp(depth_pct, d, f))

def inverse_on_support(target_force, fun, lo, hi):
    if target_force < fun(lo):
        return {'interval_pct':[0,lo], 'point_pct':None, 'scope':'LEFT_CENSORED_REFERENCE_SUPPORT'}
    if target_force > fun(hi):
        return {'interval_pct':[hi,100], 'point_pct':None, 'scope':'RIGHT_CENSORED_REFERENCE_SUPPORT'}
    d = brentq(lambda d: fun(d)-target_force, lo, hi, xtol=1e-10)
    return {'interval_pct':[d,d], 'point_pct':d, 'scope':'INTERPOLATED_ON_REFERENCE_SUPPORT'}

def distance_to_interval(y, band):
    return max(band[0]-y, y-band[1], 0)

def run(out):
    start, wall = time.process_time(), time.perf_counter()
    data = json.loads((LANE/'SOURCE_TARGETS_R3_v2.json').read_text())
    nbs = data['NBS_1973']
    curves, edges = nbs['curves_specimen18'], nbs['edges']
    thickness = nbs['specimen']['average_skin_thickness_mm']*1e-3
    a = curves['A']
    force_rows=[]
    for name in ('B','E','F','G'):
        curve, edge = curves[name], edges[name]
        for force, depth, endpoint in zip(curve['force_lbf'],curve['depth_pct'],curve['endpoint']):
            if endpoint != 'true_cut' or not 23 <= depth <= 96:
                continue
            def prediction(d):
                x = halfwidth(d*.01*thickness,edge['angle_deg'],edge['nominal_radius_m'])
                xa = halfwidth(d*.01*thickness,60,0)
                return interpolated_force(d,a)*x/xa
            pred = prediction(depth)
            eps = nbs['rounding_allowance_pp']
            band = [prediction(max(23,depth-eps)),prediction(min(96,depth+eps))]
            force_rows.append({'edge':name,'true_cut':True,'observed_depth_pct':depth,
                'observed_normal_force_N':force*LBF,'predicted_normal_force_N':pred*LBF,
                'prediction_tabular_rounding_band_N':[x*LBF for x in band],
                'relative_error':abs(pred-force)/force,
                'minimum_relative_error_with_depth_rounding':distance_to_interval(force,band)/force})
    # New radius, same angle, same specimen: B -> K uses B as single calibration geometry.
    b = curves['B']
    def rounded_contact_force(d):
        width_ratio=halfwidth(d*.01*thickness,90,50.8e-6)/halfwidth(d*.01*thickness,90,0)
        return interpolated_force(d,b)*width_ratio
    contact_k = inverse_on_support(20,rounded_contact_force,7,87)
    k_observed, k_eps = 40., 2.
    contact_k.update(observed_depth_pct=k_observed,observed_reading_band_pct=[38,42],
                     normal_force_N=20*LBF,
                     minimum_relative_depth_error=distance_to_interval(k_observed,contact_k['interval_pct'])/k_observed,
                     minimum_relative_error_with_reading=distance_to_interval(42,contact_k['interval_pct'])/42,
                     gate='FAIL')
    first={'scope':'REAL_PUBLISHED_NORMAL_FORCE_AND_CUT_DEPTH; CONSTITUTIVE_COMMON_PRESSURE_TEST',
           'calibration_geometry':'A60degree nominalR0; B90degree nominalR0 for isolated B->K radius test',
           'heldout_force_rows':force_rows,
           'max_relative_force_error':max(x['relative_error'] for x in force_rows),
           'rounded_same_angle_K':contact_k,
           'strongest_same_information_control':'Identical width law/pressure pullback; TIE. Full FE not executed.',
           'gate':'FAIL'}
    write_json(out/'pressure_width_results.json',first)

    # Executed next construction: new blunt tip pressure port. Independent acquisition
    # of one F datum plus B reference; heterogeneous F-radius is only a scenario.
    radius_f, radius_k = 25.4e-6, 50.8e-6
    calibration_depth, calibration_force = 65.,12.
    sharp_force = interpolated_force(calibration_depth,b)
    tip_offset = calibration_force-sharp_force
    coefficient = tip_offset*LBF/radius_f
    def offset_b_force(d):
        return interpolated_force(d,b)+2*tip_offset
    pred_k = inverse_on_support(20,offset_b_force,7,87)
    pred_k.update(observed_depth_pct=40,observed_reading_band_pct=[38,42],
                  relative_error=abs(pred_k['point_pct']-40)/40,
                  min_relative_error_with_reading=distance_to_interval(42,pred_k['interval_pct'])/42,
                  gate='FAIL')
    heldout_f=[]
    for force, observed in [(16,80),(20,71)]:
        pred=inverse_on_support(force,lambda d:interpolated_force(d,b)+tip_offset,7,87)
        pred.update(normal_force_lbf=force,observed_depth_pct=observed,
                    min_relative_depth_error=distance_to_interval(observed,pred['interval_pct'])/observed)
        heldout_f.append(pred)
    second={'scope':'PUBLISHED_FORCE_DEPTH_WITH_UNVERIFIED_EFFECTIVE_F_EDGE_RADIUS_SCENARIO',
        'changed_operation':'additive normal tip-pressure offset proportional to R rather than common flank pressure',
        'training':'B fullcurve plus F12lbf/65% anchor; F effectiveR25.4um scenario, NOT a uniform measured radius',
        'normal_offset_lbf_at_F':tip_offset,'normal_pressure_coefficient_N_m':coefficient,
        'K_heldout':pred_k,'F_same_geometry_heldout_loads':heldout_f,
        'gate':'FAIL','strongest_control':'same additive contact law and same datum: TIE',
        'why_not_empirical_radius_validation':'F local curvature varies; conditional pointprediction not validated across radius despite matched skin. Missing current tip radius/contact mechanism cannot be fitted from geometry label.'}
    write_json(out/'tip_pressure_results.json',second)

    needle=data['SHERGOLD_2005']
    d1,d2=needle['shank_diameter_m'];p1,p2=needle['pressure_MPa'];e1,e2=needle['digitization_bound_MPa']
    pred_p=p1*d1/d2
    extreme_errors=[abs(v-o)/o for v in [(p1-e1)*d1/d2,(p1+e1)*d1/d2] for o in [p2-e2,p2+e2]]
    needleres={'scope':'APPROXIMATE_PRIMARY_FIGURE_DIGITIZATION; NEW_SHANK_GEOMETRY_NOT_EDGE_RADIUS',
        'training_pressure_MPa':p1,'heldout_pressure_MPa':p2,'predicted_heldout_pressure_MPa':pred_p,
        'training_force_N':p1*1e6*math.pi*d1*d1/4,
        'heldout_observed_force_N':p2*1e6*math.pi*d2*d2/4,
        'heldout_predicted_force_N':pred_p*1e6*math.pi*d2*d2/4,
        'central_relative_error':abs(pred_p-p2)/p2,
        'extreme_digitization_relative_error':max(extreme_errors),
        'central_gate':'PASS','robust_digitization_gate':'FAIL',
        'effective_fracture_geometry_product_J_m2':p1*1e6*math.pi*d1/4,
        'chemical_Gamma0_J_m2':None,'skin_scalpel_edge_radius_prediction':'UNKNOWN',
        'control':'same F=cD scaling: TIE; constant a/R and neglected bulk contribution are unverified'}
    write_json(out/'needle_geometry_results.json',needleres)

    # Power identifies which force must be acquired; never infer tangential work
    # from a force component orthogonal to the prescribed cutting velocity.
    depth=.40*thickness
    gamma_candidates=[150.,2500.]
    power={'scope':'DERIVED_UNDER_QUASISTATIC_STEADY_DEPTH_AND_PROJECTED_AREA_dA=d*dx',
        'observed_normal_force_N':20*LBF,'observed_cut_depth_m':depth,
        'tangential_velocity_m_s':.0254,'steady_normal_velocity_m_s':0,
        'normal_force_power_W':0,'unmeasured_tangential_force_N':None,
        'work_relation':'Ft dx + Fn dz = dPsi + (Gamma0+B_before)dA + dD_contact',
        'nonidentifiability_witnesses':[{'effective_work_J_m2':g,
            'surface_work_force_N':g*depth,'surface_power_W':g*depth*.0254,
            'is_measurement':False,'normal_observables_changed':False} for g in gamma_candidates],
        'invalid_normal_force_divided_by_depth_J_m2':20*LBF/depth,
        'energy_observation_jacobian_rank_without_Ft':0,
        'effective_skin_sharp_cut_intercept_J_m2':None,'chemical_Gamma0_J_m2':None,
        'required_observations':['synchronous tangential force and cut-normal displacement at blade','incremental true crack area, separate crush endpoint','actual edge curvature along the contacting line and bevel','recoverable prestress/bulk work and first-pass/contact work','independent pre-cleavage full interface work if chemical Gamma0 required'],
        'strongest_control':'same virtual-work identity: TIE'}
    write_json(out/'power_port_results.json',power)

    peer=json.loads((LANE/'r3/BINDINGS_PORTS_R3_SNAPSHOT.json').read_text())['process_volume_R3']
    coupling={'source':'r3/BINDINGS_PORTS_R3_SNAPSHOT.json','peer_gate':peer['gate'],
        'peer_h_pz_m':peer['h_pz_m'],'peer_micro_cell_h_m':peer['cell_hypothesis']['h_m'],
        'NBS_K_radius_m':50.8e-6,'radius_over_micro_cell_h':50.8e-6/peer['cell_hypothesis']['h_m'],
        'crossprediction':'UNKNOWN: different species, missing Ft/work and measured h_pz; do not map normal groove width/depth to V/A',
        'sharp_suppression_counterexample_reused':peer['same_radius_contact_counterexample'],
        'mechanism_constraint':'delta_b/contact/wedge forcing survives rho->0; R2 delta_b=2R and Gamma150 scenario remain unvalidated',
        'gain_vs_uniform_FE':None,'strongest_same_information':'TIE'}
    write_json(out/'peer_coupling.json',coupling)
    summary={'operation':'Published skin sharpness/geometry acquisition -> heldout common-pressure prediction -> additive nose-pressure acquisition -> power-conjugate work port',
        'gate':'FAIL','skin_chemical_Gamma0':'UNKNOWN','uniform_fine_FE_comparison':'NOT_EXECUTED',
        'strongest_same_information_control':'TIE','prediction_results':{
            'width_max_relative_force_error':first['max_relative_force_error'],
            'K_width_min_relative_depth_error':contact_k['minimum_relative_depth_error'],
            'K_tip_predicted_depth_pct':pred_k['point_pct'],
            'K_tip_relative_depth_error':pred_k['relative_error'],
            'needle_central_relative_error':needleres['central_relative_error'],
            'needle_digitization_max_relative_error':needleres['extreme_digitization_relative_error']},
        'cost':{'instrumented_CPU_s':time.process_time()-start,'instrumented_wall_s':time.perf_counter()-wall,
            'peak_RSS_MiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,
            'new_global_FE_solves':0,'new_LP_solves':0,'actual_calibration_geometries_step1':1,
            'actual_calibration_geometries_step2':2,'literature_acquisition_research_IO_s':None},
        'source_sha256':hashlib.sha256((LANE/'SOURCE_TARGETS_R3_v2.json').read_bytes()).hexdigest(),
        'review_state':'PENDING_INDEPENDENT_REVIEW','scientific_admission':False}
    write_json(out/'summary.json',summary)
    print(json.dumps(summary,indent=2))

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--out',type=Path,default=LANE/'r3/geometry')
    args=parser.parse_args();args.out.mkdir(parents=True,exist_ok=False);run(args.out)

if __name__=='__main__':main()
