"""Retrospective computational holdouts on frozen PRIMARY skin needle series.
No FE, fabricated force curves, or per-heldout calibration. Exclusive outputs.
"""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[k]='2'
import json,math,time,resource,hashlib,argparse
from pathlib import Path
P=Path(__file__).resolve().parent

def write(p,v):
 with p.open('x') as f:json.dump(v,f,indent=2,allow_nan=False);f.write('\n')
def error(pred,obs):return abs(pred-obs)/abs(obs)
def envelope(pred,obs):
 return {'best_relative_error':0 if max(pred[0],obs[0])<=min(pred[1],obs[1]) else min(error(x,y) for x in pred for y in obs),'worst_relative_error':max(error(x,y) for x in pred for y in obs),'kind':'rectangular bounded reading; NOT probability/CI'}
def y_to_band(fig,pair,top_value,zero,top,bound):
 return [(zero-max(pair))*top_value/(zero-top)-bound,(zero-min(pair))*top_value/(zero-top)+bound]

def run(out):
 cpu=time.process_time();wall=time.perf_counter();data=json.loads((P/'SOURCE_TARGETS_R4_v2.json').read_text());b=data['Barnett2016'];s=data['Shergold2005_reused'];o=data['Owen2022'];i=1
 D=b['diameter_mm'];J=b['first_minus_repeat_effective_J_means_J_m2'];a=b['crack_length_means_mm'];F=b['total_force_means_N'];J0=J['16'][i];F0=F['16'][i];epsJ=b['Fig9']['reading_bound_J_m2'];epsF=b['Fig15']['reading_bound_N']
 rows=[]
 for g in ['18','21','25']:
  obsJ=J[g][i];predPt=J0*a[g]*1e-3;obsPt=obsJ*a[g]*1e-3
  rows.append({'gauge':int(g),'diameter_mm':D[g],'speed_mm_s':20,'heldout_effective_J_J_m2':obsJ,'predicted_effective_J_J_m2':J0,'J_relative_error':error(J0,obsJ),'predicted_propagation_force_N':predPt,'reference_propagation_force_J_times_a_N':obsPt,'force_reference_scope':'constructed from independently reported steady-work J and unloaded all-speed crack mean, not synchronous peak force','same_crack_port_reading_error':envelope([J0-epsJ,J0+epsJ],[obsJ-epsJ,obsJ+epsJ]),'central_gate':'PASS' if error(J0,obsJ)<=.2 else 'FAIL','extra_work_required_vs_single_J_J_m2':obsJ-J0})
 speed=[]
 for idx,v in enumerate(b['speeds_mm_s']):
  if idx==i:continue
  speed.append({'gauge':16,'heldout_speed_mm_s':v,'observed_J_J_m2':J['16'][idx],'predicted_J_J_m2':J0,'relative_error':error(J0,J['16'][idx]),'reading_error':envelope([J0-epsJ,J0+epsJ],[J['16'][idx]-epsJ,J['16'][idx]+epsJ])})
 bands={g:y_to_band(b['Fig9'],pair,4000,608,15,epsJ) for g,pair in b['Fig9']['at20_bar_y_px'].items()}
 common=[max(x[0] for x in bands.values()),min(x[1] for x in bands.values())]
 fracture={'training':'Gauge16 at20mm/s ONLY','training_J_J_m2':J0,'heldout_tools':rows,'heldout_speed_same_tool':speed,'common_J_feasible':common[0]<=common[1],'intersection_endpoints_J_m2':common,'per_tool_plotted_bar_plus_reading_bands_J_m2':bands,'error_bar_statistic':'UNSPECIFIED in Fig9 caption; descriptive plotting bands, not CI','control':'Conventional first-minus-repeat work/area, same inputs -> identical values:TIE','gate':'FAIL','intrinsic_Gamma0_J_m2':None,'chemical_identifiability':'Effective J includes radius-independent irreversible work. Gamma0+B_before null remains.'}
 write(out/'fracture_transfer.json',fracture)
 force_rows=[]
 for g in ['18','21','25','27']:
  obs=F[g][i] if g!='27' else (441-b['Fig16']['mean_y_px'][i])*2/426
  ratio=D[g]/D['16'];models={name:F0*ratio**power for name,power in [('constant',0),('linear_D',1),('quadratic_D',2)]}
  force_rows.append({'gauge':int(g),'speed_mm_s':20,'observed_peak_N':obs,'models':{name:{'prediction_N':v,'relative_error':error(v,obs),'central_gate':'PASS' if error(v,obs)<=.2 else 'FAIL'} for name,v in models.items()},'linear_reading_error':envelope([(F0-epsF)*ratio,(F0+epsF)*ratio],[obs-epsF,obs+epsF])})
 pub=[]
 for v,y,ym in zip(b['speeds_mm_s'],b['Fig16']['mean_y_px'],b['Fig16']['model_line_y_px']):
  obs=(441-y)*2/426;pred=(441-ym)*2/426
  pub.append({'gauge':27,'speed_mm_s':v,'observed_N':obs,'published_model_digitized_N':pred,'relative_error':error(pred,obs),'absolute_error_N':abs(pred-obs),'reading_error':envelope([pred-.03,pred+.03],[obs-.03,obs+.03]),'central_gate':'PASS' if error(pred,obs)<=.2 else 'FAIL'})
 write(out/'total_force_transfer.json',{'training':'One gauge16 peak at20mm/s; no other force anchors','heldout':force_rows,'published_richer_control_heldout27':pub,'published_control_information':'All four gauge J/friction/crack/tension datasets AND contact-factor fit to force per gauge; cannot label fit over gauges16..25 our heldout success','same_information_control':'identical one-anchor power law/doublepass algebra:TIE','published_model_REIMPLEMENTED':False,'published_model_read_from_figure':True})
 c=s['pressure_MPa'][0]*1e6*math.pi*s['shank_diameter_m'][0]/4
 human=[]
 for idx in [1]:
  d=s['shank_diameter_m'][idx];pred=c*d;obs=s['pressure_MPa'][idx]*1e6*math.pi*d*d/4
  r=d/s['shank_diameter_m'][0];p0=s['pressure_MPa'][0];e0=s['digitization_bound_MPa'][0];pi=s['pressure_MPa'][idx];ei=s['digitization_bound_MPa'][idx]
  pband=[(p0-e0)*1e6*math.pi*s['shank_diameter_m'][0]**2/4*r,(p0+e0)*1e6*math.pi*s['shank_diameter_m'][0]**2/4*r];oband=[(pi-ei)*1e6*math.pi*d*d/4,(pi+ei)*1e6*math.pi*d*d/4]
  human.append({'diameter_m':d,'predicted_friction_corrected_force_N':pred,'observed_friction_corrected_force_N':obs,'relative_error':error(pred,obs),'reading_error':envelope(pband,oband),'central_gate':'PASS','robust_reading_gate':'FAIL'})
 cross=[]
 for g in ['16','18','21','25']:
  pred=c*D[g]*1e-3;obs=J[g][i]*a[g]*1e-3
  cross.append({'gauge':int(g),'predicted_force_N':pred,'observed_effective_J_times_a_N':obs,'relative_error':error(pred,obs),'gate':'PASS' if error(pred,obs)<=.2 else 'FAIL','scope':'Cross-study HUMAN friction-corrected includes spreading vs PORCINE first-minus-repeat removes spreading. Test of deliberately strong port-equivalence hypothesis; NOT species-invariant posterior.'})
 write(out/'human_and_crossstudy.json',{'reused_human_heldout':human,'human_F_over_D_N_m':c,'human_c_is_Gamma_times_a_over_D_plus_bulk_over_D':True,'human_Gamma_J_m2':None,'cross_study_unmodified_c':cross,'conclusion':'R3~1.8kJ/m2 is effective force/diameter geometry product, not independently measured Gamma_pierce; crossstudy failures do not by themselves refute a species-specific physical model'})
 # Owen: independent first-use anchor transferred without a cohort refit.
 orows=[];ratio=o['nominal_diameter_scenario_mm']['18']/o['nominal_diameter_scenario_mm']['21'];O0=o['experiments']['2']['21_new']['median_N']
 for exp,key in [('2','18_new'),('3','21_new'),('4','21_new'),('4','18_new')]:
  obs=o['experiments'][exp][key];r=ratio if key.startswith('18') else 1
  pred=O0*r;orows.append({'experiment':int(exp),'tool':key,'observed_median_N':obs['median_N'],'predicted_median_N':pred,'relative_error':error(pred,obs['median_N']),'central_gate':'PASS' if error(pred,obs['median_N'])<=.2 else 'FAIL','observed_range_N':obs['range_N'],'diameter_scope':'nominal scenario borrowed from Barnett, true Owen diameters unknown' if key.startswith('18') else 'same21G nominal diameter; tool/cohort differences unmodelled'})
 reuse=[]
 for key in ['21_reuse12','21_reuse36','21_reuse100','21_blunt']:
  new=o['experiments']['4']['21_new'];obs=o['experiments']['4'][key];pred=new['median_N']
  reuse.append({'tool_history':key,'diameter_changed':False,'new_median_N':pred,'history_median_N':obs['median_N'],'history_to_new_ratio':obs['median_N']/pred,'no_history_pred_relative_error':error(pred,obs['median_N']),'central_gate':'PASS' if error(pred,obs['median_N'])<=.2 else 'FAIL','interpretation':'Tool wear/residue changes. New site every puncture; not same-hole repeat assay; radius not measured.'})
 write(out/'piglet_geometry_history.json',{'source':'Owen2022','single_calibration':'experiment2 first-use21G0.45N; independent-cohortheldout no refit','heldout':orows,'same_cohort_history_adverse_test':reuse,'effective_J_J_m2':None,'why':'No displacement work, first-minus-same-hole-repeat force or crack area; cannot separate friction/deformation/Gamma','n_conflicts_preserved':{'exp2_18G_methods':5,'exp2_18G_results_caption':10,'exp4_21blunt_methods_caption':6,'exp4_21blunt_results':7},'range_not_CI':True})
 # Changed observation: local paired-work cancels only equal nuisance histories.
 diagnostics=[]
 for g in ['16','18','21','25']:
  Pt=J[g][i]*a[g]*1e-3;total=F[g][i];diagnostics.append({'gauge':int(g),'peak_total_N':total,'steady_J_times_unloaded_mean_a_N':Pt,'naive_peak_remaining_N':total-Pt,'cannot_claim_friction_force':True,'reason':'Peak and steady interval not synchronous; unloaded a from different aggregate trials; no covariance. Negative center at25 is not a proved energy violation.'})
 eq=[]
 for g in ['16','18','21','25']:
  d=D[g]*1e-3;v=.02;j=1.305e7*d+1.52e4*v+3.389e10*d*d-2.222e7*d*v-1.003e13*d**3+8.006e9*d*d*v
  friction=.1493+365.7*d+.4463*math.log(v)
  eq.append({'gauge':int(g),'printed_eq12_SI_J_J_m2':j,'fig9_at20_J_J_m2':J[g][i],'printed_eq12_relative_discrepancy':error(j,J[g][i]),'printed_eq13_SI_friction_N':friction,'printed_eq13_v_mm_s_friction_N':.1493+365.7*d+.4463*math.log(20)})
 paired={'operation':'Replace diameter/peak proxy with integral(P_first-P_repeat)dl divided by matched new crack area; no synthetic solve','effective_identified':'Published Barnett first-minus-repeat effective work, under equal nuisance histories','relation':'W_delta=Gamma_pierce_eff*DeltaA + DeltaPsi_difference + DeltaD_contact_difference + DeltaD_other_difference','intrinsic_nonidentifiability':'Gamma_pierce_eff=Gamma0+B_before. Independent pre-cleavage terminal work required. Work difference contains nuisance drift if second pass changes conditioning/contact.','peak_vs_steady_diagnostics':diagnostics,'friction_port':{'Barnett_studywide_model_fraction':.21,'per_tool_phase3_measured_force_N':None,'per_tool_spreading_force_N':None,'reason':'3D figure and unreproducible printed regressions do not supply frozen matched per-case friction/spreading curves'},'printed_source_regression_reproduction':eq,'source_regression_gate':'FAIL_AS_PRINTED_SI; no silent sign/coefficient/unit correction','measured_biological_damage_width_m':None,'strongest_control':'Same-information paired work/area conventional fracture assay:identical:TIE'}
 write(out/'paired_work_observation.json',paired)
 summary={'round':4,'operation':'Primary needle series plus paired-work/area observation and geometry/history adverse tests','gate':'FAIL','same_information_control_gate':'TIE','gain':1.0,'native_common_Gamma_crossprediction':'FAIL for four-tool effective-work closure; intrinsicGamma UNKNOWN','independent_primary_studies_numeric':3,'metrics':{'human_heldout_central_error':human[0]['relative_error'],'human_worst_reading_error':human[0]['reading_error']['worst_relative_error'],'Barnett_new_tool_J_errors':[r['J_relative_error'] for r in rows],'Barnett_same_tool_speed_max_J_error':max(r['relative_error'] for r in speed),'Barnett_single_J_intersection_empty':not fracture['common_J_feasible'],'published_gauge27_max_relative_error':max(r['relative_error'] for r in pub),'Owen21to18_exp2_error':orows[0]['relative_error'],'Owen_cohort_exp4_same21G_error':orows[2]['relative_error'],'Owen100_to_new_ratio':reuse[2]['history_to_new_ratio']},'cost':{'numeric_CPU_s':time.process_time()-cpu,'numeric_wall_s':time.perf_counter()-wall,'peak_RSS_MiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,'new_global_FE_solves':0,'new_LP_solves':0,'fallbacks':0,'acquisition_research_IO_total_s':None,'physical_acquisition_executed':False,'timing_gain_claim':False},'source_sha256':hashlib.sha256((P/'SOURCE_TARGETS_R4_v2.json').read_bytes()).hexdigest(),'review_state':'PENDING_INDEPENDENT_REVIEW','scientific_admission':False}
 write(out/'summary.json',summary);print(json.dumps(summary,indent=2))
if __name__=='__main__':
 q=argparse.ArgumentParser();q.add_argument('--out',type=Path,default=P/'r4');args=q.parse_args();args.out.mkdir(parents=True,exist_ok=False);run(args.out)
