"""Append/finalize once, preserving prior checkpoint metadata."""
from pathlib import Path
import json,datetime,hashlib,shutil
P=Path(__file__).resolve().parent
UTC=datetime.datetime.now(datetime.timezone.utc).isoformat()
def read(n):return json.loads((P/n).read_text())
def write(n,v):
 with (P/n).open('x') as f:json.dump(v,f,indent=2,allow_nan=False);f.write('\n')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
S=read('r4/summary.json');V=read('VERIFICATION_R4.json');W=read('r4/work_area_requirements.json')
if V['gate']!='PASS':raise RuntimeError('verification required')
for n in ['CHECKPOINT_R4_FINAL.json','night_rounds/r4.json','WORK_STATUS_R3_BEFORE_R4.json','NEXT_ROUND_R3_BEFORE_R4.md','ARTIFACT_MANIFEST_R4_FINAL.json','FINAL_HANDOFF_CHECK_R4.json']:
 if (P/n).exists():raise FileExistsError(P/n)
section=(P/'ROUND4_SECTION.md').read_text();old=(P/'RESULTS.md').read_text()
if "## Natt 1/10, round 4" in old:raise RuntimeError('section already exists')
shutil.copyfile(P/'WORK_STATUS.json',P/'WORK_STATUS_R3_BEFORE_R4.json');shutil.copyfile(P/'NEXT_ROUND.md',P/'NEXT_ROUND_R3_BEFORE_R4.md')
files=['RESULTS.md','WORK_STATUS.json','NEXT_ROUND.md','CHECKPOINT_R4_FINAL.json','PORTS.json','PORTS_R4_FINAL.json','MEASUREMENT_SPEC_R4.md','PREREG_R4.json','PREREG_R4_SERIES.json','PREREG_R4_WORK_AREA.json','SOURCE_TARGETS_R4.json','SOURCE_TARGETS_R4_v2.json','SOURCE_REVIEW_R4.md','SOURCE_MANIFEST_R4.json','FREEZE_RECEIPT_R4.json','FREEZE_RECEIPT_R4_v2.json','PEER_CONTEXT_R4.json','DECOMPOSITION_R4.json','ATTEMPTS_R4.json','COST_ACCOUNTING_R4.json','VERIFICATION_R4.json','GRAPH_FEEDBACK_R4.json','COMMANDS_R4.md','NEXT_ROUND_R4_FINAL.md','r4/summary.json','r4/fracture_transfer.json','r4/total_force_transfer.json','r4/human_and_crossstudy.json','r4/piglet_geometry_history.json','r4/paired_work_observation.json','r4/work_area_requirements.json','r4/needle_series_overview.pdf','r4/needle_series_overview.png','needle_series_r4.py','work_area_r4.py','freeze_series_r4.py','acquire_r4.py','handoff_r4.py','verify_r4.py','plot_r4.py','finalize_r4.py','literature/r4/']
nextop='Matched Ft/Fn/x/z, true3D crack area, measured local edge radius, independent release/contact/repeat drift and terminal precleavage work plus biological/vascular damage acquisition; MEASUREMENT_SPEC_R4.md'
obstacle='Diameter/bevel/contact co-vary; tool-dependent effective paired-work J, unknown true3D area/repeat nuisance drift, Gamma0+B null and unknown viability/perfusion widths'
cp={'round':4,'status':'LANE_SCOPE_CLOSED_DATA_ACQUISITION_REQUIRED','round_complete':True,'parent_goal_status':'OPEN','innovation_gate':'FAIL','strongest_equally_informed_control_gate':'TIE','review_state':'PENDING_INDEPENDENT_REVIEW','scientific_admission':False,'main_results':files,'new_capability':'Three public independent skin needle series frozen; calibrated single-gauge tool/rate/crossstudy tests, conditional paired-work identification, area/nuisance budget and consumer handoff','new_gate':'Same-tool effectiveJ speed PASS; new-tool common-effectiveJ FAIL; empirical cleavage/cut/damage UNKNOWN','portable_operator':'needle_series_r4.py --out <NEW lane-local output>','final_ports':'PORTS_R4_FINAL.json','sources':'SOURCE_TARGETS_R4_v2.json','verification':'VERIFICATION_R4.json','costs':'COST_ACCOUNTING_R4.json','next':'NEXT_ROUND.md','next_operation':nextop,'preserve':['r1','r2','r3','SOURCE_TARGETS_R4.json','SOURCE_TARGETS_R4_v2.json','r4','WORK_STATUS_R3_BEFORE_R4.json','NEXT_ROUND_R3_BEFORE_R4.md'],'updated_utc':UTC}
write('CHECKPOINT_R4_FINAL.json',cp)
write('night_rounds/r4.json',{'round':4,'operation':S['operation']+' -> exact nuisance/area measurement budget -> final RESPONSE ports','control':'Same-information conventional paired-work/area and geometry scaling:TIE. Richer published all-four-gauge Barnett control heldout27 passes; no uniform-fine FE executed.','outcome':{**S['metrics'],'common_J_oracle_minimax_relative_error':W['single_common_J_oracle_diagnostic']['minimum_possible_max_relative_error'],'prospective_Gamma0_bound_J_m2':W['prospective_measurement_budget']['Gamma0_total_bound_J_m2'],'prospective_measurement_executed':False,'Gamma_pierce_effective':'tool-specific approximate primary paired-work assay values; not intrinsicGamma','Gamma_cut_empirical':None,'biological_damage_width':None,'verification':'42/42 PASS bookkeeping/model-scope','innovation':'FAIL','same_information_control':'TIE','lane_scope_closed':True,'parent_goal_status':'OPEN'},'gate':'FAIL','obstacle':obstacle,'next':nextop,'files':files,'review_state':'PENDING_INDEPENDENT_REVIEW','scientific_admission':False,'updated_utc':UTC})
with (P/'RESULTS.md').open('a') as f:f.write('\n'+section)
(P/'NEXT_ROUND.md').write_text((P/'NEXT_ROUND_R4_FINAL.md').read_text())
status={'lane':'LANE_SURGICAL_INCISION','round':4,'status':'LANE_SCOPE_CLOSED_DATA_ACQUISITION_REQUIRED','round_complete':True,'lane_scope_complete':True,'parent_goal_status':'OPEN','review_state':'PENDING_INDEPENDENT_REVIEW','scientific_admission':False,'innovation_gate':'FAIL','strongest_equally_informed_control_gate':'TIE','same_tool_effective_J_speed_gate':'PASS_FOR_REPORTED_ASSAY_READINGS','new_tool_common_effective_J_gate':'FAIL','empirical_skin_Gamma0_Gamma_cut_and_biological_damage':'UNKNOWN','completed_capability':cp['new_capability'],'independent_numeric_primary_studies':3,'last_checkpoint':'CHECKPOINT_R4_FINAL.json','next_round':'NEXT_ROUND.md','next_operation':nextop,'blockers':[obstacle,'Per-tool phase3 friction and synchronous first/repeat force-depth curves not acquired; source printed regressions unreproducible as stated','True skin scalpel work and terminalprocess/active-area acquisition absent; biological/vascular closures unknown'],'defaults':'PREREG_R4_SERIES.json; SOURCE_TARGETS_R4_v2.json; null empirical unknowns; syntheticcut150..380 opt-in only','results':'RESULTS.md','final_ports':'PORTS.json','measurement_spec':'MEASUREMENT_SPEC_R4.md','costs':'COST_ACCOUNTING_R4.json','verification':'VERIFICATION_R4.json','updated_utc':UTC}
(P/'WORK_STATUS.json').write_text(json.dumps(status,indent=2,allow_nan=False)+'\n')
missing=[n for n in files if not (P/n).exists()]
checks={'all_required_files_exist':not missing,'no_missing_files':missing,'one_R4_section_only':(P/'RESULTS.md').read_text().count("## Natt 1/10, round 4")==1,'RESULTS_old_content_prefix_preserved':(P/'RESULTS.md').read_text().startswith(old),'NEXT_matches_final':(P/'NEXT_ROUND.md').read_text()==(P/'NEXT_ROUND_R4_FINAL.md').read_text(),'status_round_matches_night':read('WORK_STATUS.json')['round']==read('night_rounds/r4.json')['round']==4,'previous_status_checkpoint_round3':read('WORK_STATUS_R3_BEFORE_R4.json')['round']==3,'PORTS_snapshots_equal':(P/'PORTS.json').read_bytes()==(P/'PORTS_R4_FINAL.json').read_bytes(),'biological_claims_not_admitted':read('WORK_STATUS.json')['scientific_admission'] is False}
checks['gate']='PASS' if all(v for k,v in checks.items() if k not in ['no_missing_files']) else 'FAIL'
write('FINAL_HANDOFF_CHECK_R4.json',checks)
artifacts=[]
for f in sorted(P.rglob('*')):
 if not f.is_file():continue
 rel=str(f.relative_to(P))
 if rel.startswith('r4/mpl_config/'):continue
 if rel.startswith(('r4/','literature/r4/')) or 'R4' in f.name or f.name.endswith('_r4.py') or rel in ['WORK_STATUS.json','NEXT_ROUND.md','RESULTS.md','PORTS.json','night_rounds/r4.json']:
  artifacts.append({'path':rel,'bytes':f.stat().st_size,'sha256':sha(f)})
write('ARTIFACT_MANIFEST_R4_FINAL.json',{'round':4,'utc':UTC,'entries':artifacts,'review_state':'PENDING_INDEPENDENT_REVIEW','scientific_admission':False,'manifest_self_excluded':True,'matplotlib_cache_not_scientific_artifacts':True})
print(json.dumps({'round':4,'status':status['status'],'gate':S['gate'],'control':'TIE','handoff_checks':checks['gate'],'artifact_count':len(artifacts)},indent=2))
if checks['gate']!='PASS':raise SystemExit(1)
