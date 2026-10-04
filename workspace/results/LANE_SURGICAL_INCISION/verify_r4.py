"""Scientific bookkeeping/scope checks, not independent empirical admission."""
from pathlib import Path
import json,hashlib,itertools,math,datetime
P=Path(__file__).resolve().parent
checks=[]
def read(n):return json.loads((P/n).read_text())
def check(name,ok,detail=None):checks.append({'name':name,'passed':bool(ok),'detail':detail})
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
D=read('SOURCE_TARGETS_R4_v2.json');R=read('r4/summary.json');F=read('r4/fracture_transfer.json');T=read('r4/total_force_transfer.json');H=read('r4/human_and_crossstudy.json');O=read('r4/piglet_geometry_history.json');W=read('r4/work_area_requirements.json');P0=read('PORTS_R4_FINAL.json')
check('Prediction uses frozen corrected source',sha(P/'SOURCE_TARGETS_R4_v2.json')==R['source_sha256']==read('FREEZE_RECEIPT_R4_v2.json')['sha256'])
check('Original frozen source preserved',sha(P/'SOURCE_TARGETS_R4.json')==read('FREEZE_RECEIPT_R4.json')['source_targets_sha256'])
check('Prereg unchanged',sha(P/'PREREG_R4.json')==read('FREEZE_RECEIPT_R4.json')['PREREG_sha256'])
check('Series prereg unchanged',sha(P/'PREREG_R4_SERIES.json')==read('FREEZE_RECEIPT_R4.json')['PREREG_series_sha256'])
for meta in read('SOURCE_MANIFEST_R4.json')['acquisitions']:
 n=Path(meta['metadata']).name.replace('.acquisition.json','.response');q=P/'literature/r4'/n
 if q.exists():check('Acquired response integrity '+n,sha(q)==meta.get('sha256'))
check('Two independent skin series beyond initial needle figure available',R['independent_primary_studies_numeric']>=2 and (P/'literature/r4/barnett2016.response').read_bytes().startswith(b'%PDF'))
B=D['Barnett2016'];J0=B['first_minus_repeat_effective_J_means_J_m2']['16'][1]
check('No heldout J refit',all(r['predicted_effective_J_J_m2']==J0 for r in F['heldout_tools']))
check('All new-tool central gates fail20percent',all(r['J_relative_error']>.2 and r['central_gate']=='FAIL' for r in F['heldout_tools']))
check('21G near-threshold failure remains uncertainty-qualified',F['heldout_tools'][1]['same_crack_port_reading_error']['best_relative_error']<.2<F['heldout_tools'][1]['same_crack_port_reading_error']['worst_relative_error'])
check('Reported bands give empty common J intersection',F['intersection_endpoints_J_m2'][0]>F['intersection_endpoints_J_m2'][1] and not F['common_J_feasible'])
check('Same-tool speed central+bounded-reading errors below20percent',all(r['relative_error']<=.2 and r['reading_error']['worst_relative_error']<=.2 for r in F['heldout_speed_same_tool']))
check('Strong rich control heldout27 passes20percent with declared training difference',all(r['central_gate']=='PASS' for r in T['published_richer_control_heldout27']) and not T['published_model_REIMPLEMENTED'])
check('Human R3 11percent and 38percent limits retained',math.isclose(H['reused_human_heldout'][0]['relative_error'],.11428571428571425,abs_tol=1e-12) and H['reused_human_heldout'][0]['reading_error']['worst_relative_error']>.38)
check('Owen medians and range correction are primary literal values',D['Owen2022']['experiments']['2']['21_new_late']['range_N']==[.4,.58] and D['Owen2022']['experiments']['4']['21_reuse100']['median_N']==1.28)
check('Owen missing diameter/Gamma not silently measured',D['Owen2022']['diameter_m'] is None and O['effective_J_J_m2'] is None)
check('Same-hole repeat and new-site tool reuse distinguished','New site every puncture' in O['same_cohort_history_adverse_test'][0]['interpretation'])
check('Best oracle constant still fails20percent, no holdout re-fit',W['single_common_J_oracle_diagnostic']['minimum_possible_max_relative_error']>.45)
# Exact ratio error under all bounded numerator/area/precleavage corner combinations.
r=W['prospective_measurement_budget']['input'];cap=r['effective_work_upper_for_design_J_m2'];eps=[r['external_work_bound_over_A_J_m2'],r['recoverable_work_bound_over_A_J_m2'],r['contact_work_bound_over_A_J_m2']];ra=r['area_relative_bound'];eb=r['precleavage_work_bound_J_m2'];bound=W['prospective_measurement_budget']['Gamma0_total_bound_J_m2'];errors=[]
for signs in itertools.product([-1,1],repeat=5):
 num=cap+sum(s*e for s,e in zip(signs,eps));area=1+signs[3]*ra;deltaB=signs[4]*eb
 errors.append(abs(num/area-cap-deltaB))
check('Prospective exact quotient bound covers all bounded corners',max(errors)<=bound+1e-10,{'max_corner':max(errors),'bound':bound,'measurement_achieved':False})
check('Prospective budget is unachieved conditional target',bound<30 and W['prospective_measurement_budget']['achieved_measurement_precision'] is None)
check('PORTS.json and final snapshot exactly identical',read('PORTS.json')==P0)
check('No empirical cleavage/cut substitute',P0['Gamma0_intrinsic_J_m2'] is None and P0['Gamma_cut_empirical_J_m2'] is None and P0['Gamma_cut_scenario']['status']=='SYNTHETIC_NOT_CI')
check('No per-tool friction fabricated',P0['friction_and_spreading']['per_tool_force_N'] is None)
check('Every biological damage/perfusion width unknown',all(x.get('cell_viability_zone_halfwidth_m') is None and x.get('perfusion_loss_zone_halfwidth_m') is None for x in P0['damage_zone_width_by_tool']))
check('Trace does not become biological injury',all(not x.get('transfer_to_RESPONSE_biological_damage_allowed',False) for x in P0['damage_zone_width_by_tool']))
check('No vascular/hemostasis/strength port inferred',all(P0[k] is None for k in ['vessel_map','bleeding_flow_m3_s','hemostasis_time_s','healing_strength_law']))
check('No empirical layer Gamma injected',all(v is None for v in P0['layer_specific_Gamma_cut_J_m2'].values()))
check('No new FE or false cost advantage',R['cost']['new_global_FE_solves']==0 and R['cost']['new_LP_solves']==0 and R['gain']==1 and not R['cost']['timing_gain_claim'])
check('Unknown acquisition cost not zero',R['cost']['acquisition_research_IO_total_s'] is None)
check('No admission and independent review pending',not P0['scientific_admission'] and P0['review_state']=='PENDING_INDEPENDENT_REVIEW')
for name in ['SOURCE_REVIEW_R4.md','MEASUREMENT_SPEC_R4.md','r4/needle_series_overview.pdf','r4/needle_series_overview.png','CHECKPOINT_R3.json','r3/geometry/needle_geometry_results.json','SOURCE_TARGETS_R3_v2.json']:
 check('Required/existing artifact preserved '+name,(P/name).exists())
result={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'checks':checks,'count':len(checks),'passed':sum(c['passed'] for c in checks),'gate':'PASS' if all(c['passed'] for c in checks) else 'FAIL','scope':'Bookkeeping, preserved physical port definitions and conditional bounded-error calculation; NOT independent scientific review/empirical biological validation','review_state':'PENDING_INDEPENDENT_REVIEW','scientific_admission':False}
with (P/'VERIFICATION_R4.json').open('x') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps({k:result[k] for k in ['count','passed','gate']}));
if result['gate']!='PASS':
 print(json.dumps([c for c in checks if not c['passed']],indent=2));raise SystemExit(1)
