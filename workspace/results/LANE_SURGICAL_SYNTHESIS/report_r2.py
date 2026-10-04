"""Build R2 review artifacts after computation and tests; exclusive outputs."""
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:
    os.environ[k]='2'
import datetime, hashlib, json, re
from pathlib import Path
import numpy as np
from surgical_chain.chain import write
ROOT=Path(__file__).resolve().parent
def read(p):return json.loads((ROOT/p).read_text())
def save_text(p,text):
    with (ROOT/p).open('x') as f:f.write(text)

def main():
    a=read('r2/attribution_v1/ATTRIBUTION.json');v=read('r2/measurement_v1/MEASUREMENT_VALUE.json')
    nonlinear=read('r2/nonlinear_information_v1/VALIDATION.json')
    refinement=read('r2/measurement_v1/REFINEMENT.json')
    log=(ROOT/'TESTS_R2_v1.log').read_text();match=re.search(r'Ran (\d+) tests',log)
    assert match and log.rstrip().endswith('OK');count=int(match[1])
    rows={r['case']:r for r in a['rows']};rank=v['variants']['n256_noise1.0_floorFalse']
    cost=v['costs'];now=datetime.datetime.now(datetime.timezone.utc).isoformat()
    source_integrity=[]
    for row in read('CODE_MANIFEST_R1.json')['code_files']:
        h=hashlib.sha256((ROOT/row['file']).read_bytes()).hexdigest()
        source_integrity.append(dict(row,unchanged=h==row['sha256']))
    assert all(r['unchanged'] for r in source_integrity)
    verification={'review_state':'PENDING_INDEPENDENT_REVIEW','scientific_admission':False,'passed':count,'count':count,
       'original_R1_tests':51,'new_R2_tests':count-51,'gate':'PASS','log':'TESTS_R2_v1.log',
       'attribution_gate':a['gate'],'numerical_refinement_gate':refinement['gate'],
       'R1_code_integrity':source_integrity,'scope':'Mathematical/source/diagnostic regressions; no native biological validation'}
    write(ROOT/'VERIFICATION_R2.json',verification)
    head="# Misattribution for the surgical chain, R2\n\nPENDING_INDEPENDENT_REVIEW. The same nine conserved targets are used; RMSE is counted on the seven lines with `independent_day=true`, with replicas remaining. ΔRMSE = 43,703696 − new RMSE; positive number is improvement. No parameter has been customized in R2.\n\nThe isotropic R1 observer is exactly `S=75 A m`, where `A=trace(U+I+M)` and `m=trace(M)/A`. R3 uses `S=75 C_R3 x_HP`. R3's HP curve is external input and C_R3 a model state that inherits force calibration from RESPONSE R1. The D/H/A reservoir in synthesis R1 has no power consumer.\n\n| Ersatt port | RMSE pp | ΔRMSE pp |\n|---|---:|---:|\n"
    names={'baseline':'Ingen: ursprunglig kedja','factor_A1_m0':"Amount only A → R3 C",'factor_A0_m1':'Enbart mognad m → R3 HP-normalisering',
      'factor_A1_m1':"Both A and m → R3",'current_FV_C_only':"Amount A → the chain",'DHA_HP_only':'m → isolerad D/H/A-HP med R3-skala',
      'bridge_reference_identity':'Bridge → R3 referens1','angle_reference_identity':'Riktning → isotrop R3 referens',
      'early_state_hazard_R3':'Tidig state/hazard → R3','upstream_all_R3':"All six upstream R3 ports"}
    for key in ['gap','biological_width','perfusion_width','face_conductance','k_deposit','k_mature_legacy']:
        names['upstream_'+key]=key+' → R3'
    for name,label in names.items():
        row=rows[name];head+=f"| {label} | {row['RMSE_pp']:.6f} | {row['delta_RMSE_pp']:+.6f} |\n"
    s=a['two_factor_Shapley_RMSE_reduction_pp']
    head+=f"""\nThe port of quantity carries {s['abundance']:.6f} pp av den totalt {sum(s.values()):.6f} pp. major improvement according to two factor-Shapley; maturity carries {s['maturity']:.6f} pp. Interaktionen i RMSE is {a['factorial_interaction_RMSE_pp']:.6f} pp. Shapley is here descriptive error attribute, no biological causal effect. Absolute port changes cannot be equated with independent acquired information.\n\nThe exact obstacle is in the coupling of the fixtures. FV C is formed with `k_deposit F gate (1-C)` and lost with `.008 macrophage C`. U/I/M forms simultaneously with `birth_scale F gate (1-C_FV)` and further loses `turnover_UI (U+I)`. Internal maturation fluxes cancel each other out, but the total U/I/M pool is never forced to FV C. When FV C fills, the U/I/M's birthgate closes even if its own inventory is empty.\n\nDag90 is the chain's FV C=0,640170, U/I/M sum=0,313012 and R3 C=0,998320. Just importing R3's k_deposit≈5/day to the chain drives C to≈1 and the U/I/M sum till0,009881, which degrades RMSE till59,827729 pp. All upstream R3 ports simultaneously reproduce R3's FV C at the saved times, but the remaining U/I/M observer provides endast0,620360% strength dag90. This distinguishes a misconnected reference set from an early O2 error.\n\nR3 calibration ancestry: `surgical_chain/vendor/response_r1.py:main` uses `least_squares` vs older styrka3% dag7 and20% dag21 and get k_deposit≈5/day. R3 did not do **new** strengthfit; its quantity port is nevertheless not independent force-free information. Synthesis's k_deposit=0,1 is another synthetic default choice. The statement "same chemistry input" thus does not describe this comparison.\n\nThe port change is a run diagnostic reconstruction, not a changed biological model. Native traction/fracture law remains null, both proxy branches retain previous FAIL and strongest equally informed port/FV-/cohort control is TIE. Data: [ATTRIBUTION.json](r2/attribution_v1/ATTRIBUTION.json), [raw curves](r2/attribution_v1/ATTRIBUTION_CURVES.npz), [inventory diagnosis](r2/attribution_v1/INVENTORY_DIAGNOSIS.json).\n"""
    save_text('ATTRIBUTION_R2.md',head)
    spec="# Quantitative measurement specification, R2\n\nPENDING_INDEPENDENT_REVIEW. This is a computed **planning analysis of unchanged synthetic closures**. Native strength uncertainty, full hemostasis time and edge necrosis are UNKNOWN/null. The table states prediction errors for model proxies; it is no biological posterior or achieved laboratory precision. No acquisition has been ordered.\n\n256 frozen Sobol scenarios with seed6122 are propagated through the entire chain. [Design and precision budgets](MEASUREMENT_DESIGN_R2.json) state all intervals and costs. Biological and perfusion widths share a synthetic donor driver; other drivers are independent **only as a planning assumption**. Measurement error has separate random error and shared calibration error. Optimal linear prediction uses `V_after=V_y−C_yz(C_zz+R)^−1 C_zy`; interpretation as posterior variance requires a Gaussian approximation to these moments. We also show an estimator trial trained on128 scenarios and evaluated on128 new ones.\n\nR1's six items remain. M1/M2/M3 take already represented physical port observations. M5 assumes in a conditional information scenario that chemistry/direction/connected crossings and mechanics can identify the synthetic Q observer, bridge and a constant reference scale. That assay mapping is missing in the world. M4 lacks linkage from chemistry/collagen moles to U/I/M and M6 from work to biological/perfusion damage; therefore their zero values mean only zero **implemented** downstream information value. The value when the unknown laws are identified is UNKNOWN.\n\nCost means chosen relative work units for preparation, instruments and analysis, with0,5–2× sensitivity. Actual hours/prices are null. The score is the equally weighted mean reduction of variance in strength28/90, min(hemostasis time,120) and damage proxy width, divided by cost. With a constant censored outcome its sampling variance is counted as0 without giving the unknown late time variance0.\n\n| Planning rank | R1 item | Cost | SD strength28 before→after pp | SD strength90 before→after pp | SD min(T,120) before→after min | SD damage proxy before→after µm | Mean variance reduction/cost |\n|---:|---|---:|---:|---:|---:|---:|---:|\n"
    for i,row in enumerate(rank,1):
        prior=row['prior_SD'];post=row['conditional_linear_SD']
        fmt=lambda j:f'{prior[j]:.3f}→{post[j]:.3f}'
        spec+=f"| {i} | {row['id']}: {row['title']} | {row['workload_units']:.0f} | {fmt(0)} | {fmt(1)} | {fmt(2)} | {fmt(4)} | {row['value_per_cost']:.5f} |\n"
    spec+=f"\n{v['censored_samples']}/256 scenarios do not reach haemostasis before120mine and are censored; the proportion of closed is{v['closure_by120_probability_under_assumed_prior']:.4f}. Samplings-SD for the limited time is{v['prior_SD'][2]:.3f}min. This does not identify later haemostasis: the requested unlimited time SD and information value remains null. r/Q/plugg on a longer measured clock before a new time priority.120My sightings can't tell different late events.\n\nThe damage proxy is the Q036 translation's width where `1−exp(−hazard)>0,20`; the hypoxia proxy uses10Torr. These declared thresholds are no necrosis assay. Native edge necrosis needs separate viability with a typed cell, time and local registration. For all six items its posterior-SD and SD reduction are null.\n\n## Precision per delobservation\n\nσ states a chosen future random-error target; “shared” is a calibration error shared across a package. These are requirements for the calculation, not verified assay accuracy. All subitems from R1 that lack an observation operator are in the design's `unrepresented` and are commented on above; they have UNKNOWN native information value.\n\n| Paket | Observation | σ / gemensamt | ΔSD styrka28 pp | ΔSD styrka90 pp | ΔSD skadeproxy µm |\n|---|---|---:|---:|---:|---:|\n"
    for key,obs in v['individual_observations'].items():
        for row in obs:
            delta=np.array(row['prior_SD'])-row['conditional_linear_SD']
            spec+=f"| {key} | {row['observation']} | {row['sigma']:.6g} / {row['common_sigma']:.6g} | {delta[0]:.4f} | {delta[1]:.4f} | {delta[4]:.4f} |\n"
    spec+="\nSI units and the observations' exact meaning: gap/widths/path in m, pressure in Pa (driving) or Torr (pressure20/100/300/face), μ in Pa·s, r/density-factor and coverage/seal/bridge/reference/trace_M dimensionless, initial_flow m³/s, face_flux mlO2/(m² min), oxygen_storage mlO2/(m³ Torr), rates1/min or1/dygn according to the port table, D/H mol/molreferenskollagen, strength21 pp. M5's mechanical reference is an explicitly synthetic port acquisition at21dygn, not a late strength-target fit.\n\n## Sensitivity, dependencies and useful prioritization\n\n"
    for label,rr in v['variants'].items():
        spec+=f"- {label}: "+' > '.join(x['id'] for x in rr)+'.\n'
    spec+='\n'
    for row in rank:
        q=row['heldout_train128_test128_MSE_reduction_fraction']
        spec+=f"{row['id']}: nytt128-scenarioprov ger varians/MSE-minskning28/90/skadeproxy {100*q[0]:.1f}/{100*q[1]:.1f}/{100*q[4]:.1f}%.\n\n"
    spec+="Precision0,5×/2× and128/256 compare the same frozen prior. Separate cost intervals can change rank: for each pair dominance is counted only when the candidate's lower value/cost exceeds the control's upper one. The additional floors20pp and50µm are open sensitivity assumptions for an unresolved observation law, no measured error bounds. They can neither be identified nor reduced by the modeled port measurements.\n\nThe concrete next acquisition shall combine an amount inventory on a fixed reference with the early measurement ports that carry the planning table. Native strength still requires M4's absolute chemistry/reference **and** M5's connected mechanical observations with the same cohort. A low model score for M4 is an indicator of a missing consumer, not a reason to forgo chemistry. M6 is needed for new tools when the transfer law is filled. The parent goal's empirical prioritization is therefore still conditional.\n\nAll raw scenarios, subobservations, covariance/noise and costs are in [MEASUREMENT_VALUE.json](r2/measurement_v1/MEASUREMENT_VALUE.json) and sample000…255.json. Large field data are on `external_mount`. The strongest control receives the same acquisition, moments and optimal linear estimator: **TIE**.\n"
    spec+="\n## Executed nonlinear information trial after linear loss\n\nThe preserved128→128 trial worsened hypoxia prediction MSE by68,16% forM2 and32,09% forM1. The moment table's sampling-based contraction therefore cannot alone guarantee predictive benefit for this threshold observable. We changed the observation calculation to a positive Gaussian kernel with training-based, measurement-noised leave-one-out over five frozen bandwidths. Physics, biology and measurement budgets are unchanged.64 new Sobol scenarios (seed6124), each with32 independent future measurement-noise draws, evaluate prediction risk. Training uses only the earlier256 scenarios. The initial CV seed version is preserved; training6126 and future measurement noise6125 are separated in V2 before evaluation.\n\nHere RMS is **prediction error on new scenarios**, which differs from the earlier Gaussian moment-SD. No native posterior intervals follow.\n\n| Package | New64: RMS strength28 pp | RMS strength90 pp | RMS damage proxy µm | RMS hypoxia proxy µm | Hypoxia gate against prior | Strength/damage against linear control |\n|---|---:|---:|---:|---:|---|---|\n"
    for q in nonlinear['rows']:
        err=q['kernel_RMS_prediction_error']
        spec+=f"| {q['id']} | {err[0]:.3f} | {err[1]:.3f} | {err[4]:.3f} | {err[5]:.3f} | {q['kernel_vs_prior_hypoxia_gate']} | {q['strength_injury_vs_linear_gate']} |\n"
    spec+=f"\nFryst pivotgrind: **{nonlinear['gate']}**. The new observed value/cost order is "+' > '.join(q['id'] for q in nonlinear['rank'])+"All original linear errors, the initial training, V2-coefficienter, nya64-fields and measuring noises/residuals remain in [VALIDATION.json](r2/nonlinear_information_v1/VALIDATION.json) and its raw/model-filer. Strongest conventional kernel control has the same input and operation: TIEThis is a test of an information operator on synthetic data, no new biological ability.\n"
    save_text('MEASUREMENT_SPEC_R2.md',spec)
    write(ROOT/'DECOMPOSITION_R2.json',{'parent_goals':['GRAPH_INVERSE','MULTIPHYSICS'],
       'decision':'Locate strength-loss port and next acquisition, without fitting new biology',
       'leaves':[{'name':'Isotropic projection/factorial attribution','status':'DERIVED_UNDER_ASSUMPTIONS','stopping':'Existing tensors isotropic; no microscopic force law needed to diagnose algebra','test':'factorization residual and exact R3 reconstruction'},
          {'name':'FV and UIM birth/removal ledgers','status':'CONSTITUTIVE_CLOSURE','stopping':'Read existing source and cancellation of U/I/M transitions; absolute tissue mass unknown','test':'k_deposit substitution saturates FV C and starves observer'},
          {'name':'R3 HP/rat-strength curves','status':'EXTERNALLY_MEASURED','stopping':'Separate source assays and inherited modelfit, no cross-cohort common posterior','test':'Preserved seven independent-day row scores'},
          {'name':'Noisy experimental design','status':'DERIVED_UNDER_ASSUMPTIONS','stopping':'Finite synthetic prior and Gaussian moment planning; native posterior unspecified','test':'Scalar analytic update, shared-calibration regression and heldout estimator'},
          {'name':'Native chemical/traction/viability observation maps and labcost','status':'UNKNOWN','stopping':'No matched measurements or identified conversion; retain null','test':'Acquire fixed reference and registered assay map; do not infer from model SD'}],
       'source_context':'Read peer nightR1 NEXT_ROUND and receipts; matched conventional kernels are TIE, no certified coupledcostgain transferable',
       'review_state':'PENDING_INDEPENDENT_REVIEW'})
    summary={'round':2,'updated_utc':now,'review_state':'PENDING_INDEPENDENT_REVIEW','scientific_admission':False,
        'attribution':{'baseline_RMSE_pp':a['baseline_RMSE_pp'],'abundance_only_RMSE_pp':rows['factor_A1_m0']['RMSE_pp'],
           'maturity_only_RMSE_pp':rows['factor_A0_m1']['RMSE_pp'],'both_RMSE_pp':rows['factor_A1_m1']['RMSE_pp'],
           'shapley_pp':s,'gate':a['gate'],'root_obstacle':'Independent FV and UIM reference inventories; FV saturation gates observerbirth; externalHP not coupled chemistry'},
        'measurement_planning':{'sample_count':256,'rank':[x['id'] for x in rank],'prior_SD':v['prior_SD'],
           'censored_samples':v['censored_samples'],'native_prediction_uncertainty':v['native_outcomes'],
           'scope':'Synthetic prior and linear-MSE design; assay maps, fullhemostasis, native necrosis and laboratorycost UNKNOWN'},
        'tests':{'passed':count,'count':count},'numerical_refinement_gate':refinement['gate'],
        'integrated_strength_proxy_gate':'FAIL','integrated_strength_proxy_RMSE_pp':43.703695914540624,
        'strongest_control_gate':'TIE','empirical_joint_gate':'UNKNOWN','new_biology_added':False,
        'R1_global_context_import':'Verified four coordinator receipts under tasks/graph_runs, pending independent review',
        'solver_costs':cost,'attribution_wall_s':a['wall_s'],'attribution_CPU_s':a['CPU_s']}
    summary['nonlinear_information_pivot']=nonlinear
    summary['all_nonlinear_scenario_solves']=320
    write(ROOT/'r2/ROUND_SUMMARY.json',summary)
    write(ROOT/'CHECKPOINT_R2.json',{'round':2,'review_state':'PENDING_INDEPENDENT_REVIEW','scientific_admission':False,
       'physical_state_checkpoint':'r1/nominal_complete_v4/checkpoint21.npz','attribution':'r2/attribution_v1/ATTRIBUTION.json',
       'measurement_moments':'r2/measurement_v1/MEASUREMENT_VALUE.json','design':'MEASUREMENT_DESIGN_R2.json',
       'sampler_state':'MEASUREMENT_PRIOR_R2.npz','unchanged_R1_source_code':True,'integrated_proxy_gate':'FAIL',
       'native_uncertainty':'UNKNOWN','next':'Single reference inventory and shared formation/removal ports; separate native assay map prerequisite'})
    save_text('COMMANDS_R2.md',"# Reproduktion R2\n\nThreads2 is set in the scripts. Choose **new** output/file names for reruns; existing destination files are rejected. PREREG/design must be verified before a new run. Large field data is on games-240, not shared_data.\n\n```bash\npython results/LANE_SURGICAL_SYNTHESIS/diagnose_r2.py --out results/LANE_SURGICAL_SYNTHESIS/r2/attribution_REPLAY\npython results/LANE_SURGICAL_SYNTHESIS/tests_r2.py --out results/LANE_SURGICAL_SYNTHESIS/r2/tests_REPLAY\n```\n\nEnsemble: `plan_measurements_r2.py batch --start 0 --stop 32 --out <new lane directory> --store <new games-240 lane directory>` runs original sample IDs to exclusive new destinations. Reuse the frozen R2 design and prior. Run batch intervals0:32,32:96,96:160,160:224,224:256 and then summary and refinement with the same new --out/--store. Do not mix priors or results between rounds.\n\n`report_r2.py` requires a complete ensemble, refinement and OK log and writes exclusive report files. OriginalR1 code is verified against CODE_MANIFEST_R1.json. Result binding uses the lane's local graph, kindreview; canonicalimport is performed by the coordinator.\n\nNonlinear information test: nonlinear_information_r2.py prior → fit → build --start 0 --stop 32 → build --start 32 --stop 64 → evaluate. V2 uses independent training/measurement noise. Its fixed exclusive destinations require a new lane-local replay copy with new OUT/STORE for reruns; existing data must not be overwritten.\n")
    print(json.dumps({'tests':count,'attribution_gate':a['gate'],'planning_rank':[x['id'] for x in rank],
        'nonlinear_pivot_gate':nonlinear['gate'],'new64_rank':[q['id'] for q in nonlinear['rank']]},ensure_ascii=False))

if __name__=='__main__':main()
