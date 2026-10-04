#!/usr/bin/env python3
"""Scientific checks of fixed geometry and newly acquired observable ports."""
import itertools
import hashlib
import json
import math
import time
import numpy as np
from geometry_r3 import LANE, halfwidth, interpolated_force, write_json

def main():
    start = time.process_time()
    data = json.loads((LANE/'SOURCE_TARGETS_R3_v2.json').read_text())['NBS_1973']
    b = data['curves_specimen18']['B']
    geom = json.loads((LANE/'r3/geometry/tip_pressure_results.json').read_text())
    met = json.loads((LANE/'r3/metrology_results.json').read_text())
    continuity=[]
    for angle in (15,60,90,105):
        theta=math.radians(angle/2);radius=50.8e-6
        z=radius*(1-math.sin(theta))
        circle=math.sqrt(2*radius*z-z*z)
        flank=z*math.tan(theta)+radius*math.cos(theta)/(1+math.sin(theta))
        dc=(radius-z)/circle
        continuity.append({'angle_deg':angle,'width_discontinuity_m':abs(circle-flank),
                           'slope_discontinuity':abs(dc-math.tan(theta))})
    # Independent analytical inversion on B's 2..4 lbf branch, no root solver.
    offset=12-(2+(65-22)*2/(66-22))
    analytick=22+(20-2*offset-2)*22
    savedk=geom['K_heldout']['point_pct']
    rf=met['DIAGNOSTIC_required_F_radius_um']
    def inverse_radius(r):return 22+(20-50.8*offset/r-2)*22
    eps=1e-4
    derivative=(inverse_radius(rf+eps)-inverse_radius(rf-eps))/(2*eps)
    expected_derivative=met['depth_sensitivity_pp_per_um_at_diagnostic_radius']
    # F constant offset is a measurable product; its value is independent of assumed RF.
    offset_ranges=[]
    for load,depth in ((12,65),(16,80),(20,71)):
        vals=[]
        for signs in itertools.product((-1,1),repeat=4):
            cb={'depth_pct':[x+.5*s for x,s in zip(b['depth_pct'][:4],signs)],'force_lbf':b['force_lbf'][:4]}
            for delta in (-.5,.5):
                vals.append(load-interpolated_force(depth+delta,cb))
        offset_ranges.append({'load_lbf':load,'depth_pct':depth,
            'identified_offset_lbf':load-interpolated_force(depth,b),
            'shared_reference_and_depth_rounding_interval_lbf':[min(vals),max(vals)]})
    commonlo=max(r['shared_reference_and_depth_rounding_interval_lbf'][0] for r in offset_ranges)
    commonhi=min(r['shared_reference_and_depth_rounding_interval_lbf'][1] for r in offset_ranges)
    # Interior sampling checks corner band; sampling is not an enclosure proof.
    rng=np.random.default_rng(6103);inside=[]
    for _ in range(1000):
        cb={'depth_pct':(np.array(b['depth_pct'][:4])+rng.uniform(-.5,.5,4)).tolist(),'force_lbf':b['force_lbf'][:4]}
        depth=65+rng.uniform(-.5,.5);o=12-interpolated_force(depth,cb)
        ftarget=20-2*o
        inside.append(float(np.interp(ftarget,cb['force_lbf'],cb['depth_pct'])))
    band=met['K_nominal_radius_ratio2_depth_band_pct']
    pdf=LANE/'literature/r3/nbs_skin_cutting.pdf'
    manifest=json.loads((pdf.with_suffix('.acquisition.json')).read_text())
    summary=json.loads((LANE/'r3/geometry/summary.json').read_text())
    checks={
        'circle_flank_continuity':max(x['width_discontinuity_m'] for x in continuity)<1e-18,
        'circle_flank_slope_continuity':max(x['slope_discontinuity'] for x in continuity)<1e-12,
        'independent_K_analytic_inverse':abs(analytick-savedk)<1e-8,
        'finite_difference_radius_derivative':abs(derivative-expected_derivative)<1e-6,
        'sampled_shared_rounding_inside_corner_band':min(inside)>=band[0]-1e-8 and max(inside)<=band[1]+1e-8,
        'NBS_primary_PDF_hash_matches':hashlib.sha256(pdf.read_bytes()).hexdigest()==manifest['sha256'],
        'source_freeze_hash_matches':hashlib.sha256((LANE/'SOURCE_TARGETS_R3_v2.json').read_bytes()).hexdigest()==summary['source_sha256'],
        'F_nonmonotone_data_preserved':data['curves_specimen18']['F']['depth_pct']==[14,32,65,80,71],
        'K_grooves_not_fracture_observations':data['curves_specimen18']['K']['endpoint']==['groove','groove','true_cut']
    }
    result={'checks':checks,'all_checks_pass':all(checks.values()),'geometry_continuity':continuity,
       'analytic_K_depth_pct':analytick,'finite_difference_radius_sensitivity_pp_per_um':derivative,
       'K_shared_rounding_interior_sample_range_pct':[min(inside),max(inside)],
       'F_constant_offset_identification':offset_ranges,
       'constant_F_offset_intersection_lbf':[commonlo,commonhi],
       'constant_F_offset_gate':'FAIL_WITH_RESOLUTION_ALLOWANCES' if commonlo>commonhi else 'UNKNOWN',
       'failure_meaning':'No constant offset can explain these reported F loads/depths, irrespective of RF scenario. Force calibration, specimen variability and biological repeatability are unbounded; no statistical rejection.',
       'source_visual_check':'Table B E/F/G specimen18 transcribed from PDF page71; Fig17 K specimen18 at20lb open circle40%, lower12/16lb filled grooves. Fig12 F4lb/14% open true-cut. Reviewed rendered primary pages.',
       'cost_CPU_s':time.process_time()-start,
       'numerical_limits':'Binary64 and deterministic sampling; corner sensitivity bands are resolution scenarios, not verified outward-rounding or biological confidence intervals.',
       'review_state':'PENDING_INDEPENDENT_REVIEW','scientific_admission':False}
    write_json(LANE/'VERIFICATION_R3.json',result)
    print(json.dumps(result,indent=2))
    assert result['all_checks_pass']

if __name__=='__main__':main()
