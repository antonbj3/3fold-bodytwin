#!/usr/bin/env python3
"""A measured-force/unknown-local-radius acquisition ray, not a new skin fit."""
import itertools
import json
import time
from pathlib import Path
import numpy as np
from geometry_r3 import LANE, LBF, halfwidth, interpolated_force, inverse_on_support, write_json

def main():
    t=time.process_time()
    d=json.loads((LANE/'SOURCE_TARGETS_R3_v2.json').read_text())['NBS_1973']
    curves,edges=d['curves_specimen18'],d['edges'];b=curves['B'];a=curves['A']
    rounding=[]
    for row in json.loads((LANE/'r3/geometry/pressure_width_results.json').read_text())['heldout_force_rows']:
        vals=[]
        for signs in itertools.product((-1,1),repeat=4):
            ca={'depth_pct':[x+.5*s for x,s in zip(a['depth_pct'][:4],signs)],'force_lbf':a['force_lbf'][:4]}
            for delta in (-.5,.5):
                dep=row['observed_depth_pct']+delta;fa=interpolated_force(dep,ca)
                if fa is None:continue
                ed=edges[row['edge']];h=d['specimen']['average_skin_thickness_mm']*1e-3
                vals.append(fa*halfwidth(dep*.01*h,ed['angle_deg'],ed['nominal_radius_m'])/halfwidth(dep*.01*h,60,0)*LBF)
        low,high=min(vals),max(vals);obs=row['observed_normal_force_N']
        rounding.append({'edge':row['edge'],'depth_pct':row['observed_depth_pct'],
             'shared_reference_and_target_rounding_band_N':[low,high],
             'minimum_relative_error':max(low-obs,obs-high,0)/obs})
    offset=12-interpolated_force(65,b);rk=50.8
    required= rk*offset/(20-interpolated_force(40,b))
    good_r=[rk*offset/(20-interpolated_force(depth,b)) for depth in (32,48)]
    reading_r=[rk*offset/(20-interpolated_force(depth,b)) for depth in (38,42)]
    derivative=22*rk*offset/required**2 # dDepth_pct / dR_F_um, B2..4lbf branch
    scenarios=[]
    for center in (25.4,required):
        for eps in (.1,.5,1.,2.):
            band=[]
            for radius in (center-eps,center+eps):
                pred=inverse_on_support(20,lambda dep:interpolated_force(dep,b)+rk/radius*offset,7,87)
                band.extend(pred['interval_pct'])
            scenarios.append({'center_R_F_um':center,'radius_error_um':eps,
                'predicted_K_depth_interval_pct':[min(band),max(band)],
                'center_fitted_from_K':center==required,
                'status':'DIAGNOSTIC_ONLY' if center==required else 'CONDITIONAL_SENSITIVITY'})
    conditional_corner=[]
    for signs in itertools.product((-1,1),repeat=4):
        cb={'depth_pct':[x+.5*s for x,s in zip(b['depth_pct'][:4],signs)],'force_lbf':b['force_lbf'][:4]}
        for depth in (64.5,65.5):
            o=12-interpolated_force(depth,cb)
            pred=inverse_on_support(20,lambda dep:interpolated_force(dep,cb)+2*o,cb['depth_pct'][0],cb['depth_pct'][-1])
            conditional_corner.extend(pred['interval_pct'])
    result={'scope':'EMPIRICAL_REFERENCE_FORCE_CURVE_WITH_UNKNOWN_CONTACT_CURVATURE; NO_REFIT_PREDICTION',
        'shared_rounding_enclosures':rounding,
        'maximum_min_relative_force_error_with_joint_rounding':max(x['minimum_relative_error'] for x in rounding),
        'K_nominal_radius_ratio2_depth_band_pct':[min(conditional_corner),max(conditional_corner)],
        'DIAGNOSTIC_required_F_radius_um':required,
        'DIAGNOSTIC_F_radius_for_K_20pct_depth_budget_um':good_r,
        'DIAGNOSTIC_F_radius_for_K_digitization_band_um':reading_r,
        'depth_sensitivity_pp_per_um_at_diagnostic_radius':derivative,
        'necessary_radius_error_um_if_other_errors_zero':8/derivative,
        'unknown_radius_scenarios':scenarios,
        'contact_gauge':'Only k*R_F is acquired from F. k->c*k,R_F->R_F/c preservesF. Another radiusquery needslocalmetrology.',
        'normal_power_Jacobian_rank_Gamma':0,
        'empirical_radius_law_gate':'UNKNOWN; actual localR_Funmeasured and can explainK without contradictingapprox25.4um description',
        'nominal_ratio2_conditional_gate':'FAIL_PRESERVED',
        'control':'Same identified contactproduct and inverse geometry ray: TIE; acquisition precision not a measuredachievable guarantee',
        'cost_CPU_s':time.process_time()-t,
        'review_state':'PENDING_INDEPENDENT_REVIEW','scientific_admission':False}
    write_json(LANE/'r3/metrology_results.json',result)
    print(json.dumps({k:v for k,v in result.items() if k not in ('shared_rounding_enclosures','unknown_radius_scenarios')},indent=2))

if __name__=='__main__':main()
