#!/usr/bin/env python3
"""Append round-three research handoff, preserving prior raw results."""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[key]='2'
import datetime
import hashlib
import json
import time
from pathlib import Path

LANE=Path(__file__).resolve().parent
os.environ['MPLCONFIGDIR']=str(LANE/'r3/mplconfig')

def read(path):return json.loads((LANE/path).read_text())
def write(path,obj):
    with (LANE/path).open('x') as f:
        json.dump(obj,f,indent=2,allow_nan=False);f.write('\n')

def plot(data, geometry, metrology):
    start=time.process_time()
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    from geometry_r3 import interpolated_force
    nbs=data['NBS_1973'];curves=nbs['curves_specimen18']
    fig,axes=plt.subplots(1,3,figsize=(13,4.6),layout='constrained')
    for name,color in [('B','#3276aa'),('F','#df8a28'),('K','#8248a8')]:
        c=curves[name]
        for force,depth,endpoint in zip(c['force_lbf'],c['depth_pct'],c['endpoint']):
            axes[0].plot(force,depth,'o',color=color,markerfacecolor='white' if endpoint=='true_cut' else color)
        axes[0].plot([],[],'o',color=color,markerfacecolor='white',label=name)
    axes[0].set(xlabel='Applied normal force (lbf)',ylabel='Residual injury depth (% thickness)',
                title='Primary skin specimen 18',ylim=(0,105),xlim=(0,21))
    axes[0].legend(title='Edge');axes[0].text(.04,.05,'Open: true cut; filled: groove\nNo measured tangential force',transform=axes[0].transAxes,fontsize=9)
    axes[1].errorbar([0],[40],yerr=[2],fmt='o',color='black',label='K observation / reading allowance')
    axes[1].errorbar([1],[64],yerr=[[1.5],[1.5]],fmt='s',color='#df8a28',label='Offset law, RF=25.4 µm scenario')
    axes[1].vlines(2,87,100,color='#3276aa',linewidth=4,label='Width law: support censored')
    axes[1].set(xticks=[0,1,2],xticklabels=['Measured K','Offset','Width'],ylabel='Depth (%)',ylim=(0,105),title='Heldout K at 20 lbf')
    axes[1].legend(fontsize=8,loc='lower right')
    b=curves['B'];offset=geometry['normal_offset_lbf_at_F'];r=np.linspace(22,27,300)
    # Analytic inverse through observed B support; no extrapolation.
    target=20-50.8*offset/r
    dep=np.interp(target,b['force_lbf'][:4],b['depth_pct'][:4],left=np.nan,right=np.nan)
    axes[2].plot(r,dep,color='#df8a28')
    axes[2].axhspan(38,42,color='black',alpha=.12,label='K reading allowance')
    axes[2].axvline(25.4,color='#3276aa',ls='--',label='Unverified RF scenario')
    axes[2].axvline(metrology['DIAGNOSTIC_required_F_radius_um'],color='#8248a8',ls=':',label='Inverse RF: diagnostic only')
    axes[2].set(xlabel='Assumed actual RF (µm)',ylabel='K depth predicted by offset law (%)',title='Missing radius metrology')
    axes[2].legend(fontsize=8,loc='lower right')
    fig.suptitle('NBSIR 73-262: published normal-force / injury-depth data; PENDING_INDEPENDENT_REVIEW',fontsize=11)
    for ax in axes:ax.grid(alpha=.16)
    for ext in ('pdf','png'):
        out=LANE/f'r3/round3_overview.{ext}'
        if out.exists():raise FileExistsError(out)
        fig.savefig(out,dpi=160)
    plt.close(fig)
    return time.process_time()-start

def main():
    now=datetime.datetime.now(datetime.timezone.utc).isoformat()
    data=read('SOURCE_TARGETS_R3_v2.json');geom=read('r3/geometry/summary.json')
    tip=read('r3/geometry/tip_pressure_results.json');met=read('r3/metrology_results.json')
    ver=read('VERIFICATION_R3.json')
    assert ver['all_checks_pass']
    for name in ('WORK_STATUS.json','NEXT_ROUND.md'):
        assert (LANE/'r3'/('prior_'+name)).exists()
    marker="## Night 1/10, Round 3"
    assert marker not in (LANE/'RESULTS.md').read_text()
    figure_cpu=plot(data,tip,met)
    instrumented_cpu=geom['cost']['instrumented_CPU_s']+met['cost_CPU_s']+ver['cost_CPU_s']
    acquisition=[]
    for p in sorted((LANE/'literature/r3').glob('*.acquisition.json')):
        acquisition.append({'manifest':str(p.relative_to(LANE)),**json.loads(p.read_text())})
    write('SOURCE_MANIFEST_R3.json',{'acquisitions':acquisition,
        'primary_data_contract':'SOURCE_TARGETS_R3_v2.json','review_state':'PENDING_INDEPENDENT_REVIEW',
        'scientific_admission':False,'failed_responses_preserved':True})
    write('COST_ACCOUNTING_R3.json',{
        'scope':'New R3 calculations only; R1/R2 are sunk work, separately preserved and not omitted from any claimed end-to-end gain.',
        'numerical_instrumented_CPU_lower_bound_s':instrumented_cpu,
        'geometry_CPU_s':geom['cost']['instrumented_CPU_s'],
        'metrology_CPU_s':met['cost_CPU_s'],'verification_CPU_s':ver['cost_CPU_s'],
        'figure_setup_and_render_CPU_s':figure_cpu,
        'combined_instrumented_CPU_lower_bound_s':instrumented_cpu+figure_cpu,
        'source_acquisition_research_IO_and_finalization_s':None,
        'import_setup_s':None,'full_end_to_end_cost_s':None,
        'reported_geometry_peak_RSS_MiB':geom['cost']['peak_RSS_MiB'],
        'peak_RSS_scope':'Geometry process measurement only; not complete round resource telemetry.',
        'literature_r3_stored_bytes':sum(p.stat().st_size for p in (LANE/'literature/r3').iterdir() if p.is_file()),
        'HTTP_acquisition_manifests':len(acquisition),
        'valid_primary_PDF_acquisitions':sum(bool(a.get('valid_pdf')) for a in acquisition),
        'new_global_FE_solves':0,'new_LP_solves':0,'new_history_integrations':0,
        'uniform_fine_FE_comparison':'NOT_EXECUTED','strongest_same_information_control':'TIE',
        'warm_query_speedup':None,'order_of_magnitude_gain_gate':'FAIL',
        'reason':'The observable acquisition is decisive, but no coupled nonlinear state/history query or end-to-end advantage was demonstrated.',
        'review_state':'PENDING_INDEPENDENT_REVIEW','scientific_admission':False})
    write('DECOMPOSITION_R3.json',{
        'root':'Predict changed blade geometry, gap/damage and bleeding without full global state/history recomputation',
        'leaves':[
            {'operation':'Acquire same-skin published geometry series','information':'NBS specimen18, nominal edge radius/angle, Fn and residual injury depth','status':'ACQUIRED_WITH_LIMITS'},
            {'operation':'Transfer mean normal contact pressure using blade width','information':'A reference curve, nominal geometry, assumed depth/contact mapping','status':'FAIL','relative_error':geom['prediction_results']['width_max_relative_force_error']},
            {'operation':'Acquire constant nose-pressure product and transfer radius','information':'B reference+F12lb anchor','status':'CONDITIONAL_FAIL_EMPIRICAL_RADIUS_UNKNOWN','K_depth_prediction_pct':64},
            {'operation':'Identify actual curvature separately from contact coefficient','information':'Only k*RF acquired; local RF absent','status':'UNKNOWN','required_RF_precision_um_conditional':met['necessary_radius_error_um_if_other_errors_zero']},
            {'operation':'Test same-F offset with changed load','information':'F12/16/20lb and B shared curve','status':ver['constant_F_offset_gate']},
            {'operation':'Acquire work conjugate to cutting motion','information':'Fn known but vn=0; Ft absent','status':'UNKNOWN','direct_energy_observation_rank_Gamma':0},
            {'operation':'Separate intrinsic cleavage from pre-cut process','information':'Effective work includes Gamma0+B0 plus release/contact','status':'UNKNOWN_EXACT_R2_NULL_REMAINS'},
            {'operation':'Map advancing fracture work to microstructure','information':'BINDINGS r3 process zone null; cell contour/partial loop are scenarios','status':'UNKNOWN'},
            {'operation':'Local nonlinear history update vs strongest FE','information':'R1/R2 preserved, no new adequate acquisition','status':'NOT_TESTED_R3'},
            {'operation':'Cell damage/vessel bleeding','information':'No common measured ports','status':'UNKNOWN'}],
        'control':'Same information and same contact/work laws; algebraic TIE. Uniform fine FE not run.',
        'review_state':'PENDING_INDEPENDENT_REVIEW','scientific_admission':False})
    write('ATTEMPTS_R3.json',{
        'sequence':[
            {'prereg':'PREREG_R3.json','operation':'Acquire primary edge sharpness/work data','outcome':'Real human skin normal-force/depth series acquired; work series absent','gate':'UNKNOWN_FOR_WORK'},
            {'prereg':'PREREG_R3_GEOMETRY.json','operation':'A->B/E/F/G common contact-width pressure and B->K radius','gate':'FAIL','max_force_error':geom['prediction_results']['width_max_relative_force_error']},
            {'prereg':'PREREG_R3_GEOMETRY.json','operation':'Execute new additive nose-pressure acquisition after first failure','gate':'FAIL_CONDITIONAL','K_error':geom['prediction_results']['K_tip_relative_depth_error']},
            {'prereg':'PREREG_R3_METROLOGY.json','operation':'Replace assumed radius by joint contact product and resolve acquisition precision','gate':'UNKNOWN_EMPIRICAL','diagnostic_RF_um':met['DIAGNOSTIC_required_F_radius_um'],'heldout_reuse':'Diagnostic inverse only; no refit pass'},
            {'prereg':'PREREG_R3_VERIFICATION.json','operation':'Execute same-F load falsifier; validate analytical inverse/derivative/data integrity','gate':'FAIL_CONSTANT_OFFSET','verification':'PASS'},
            {'prereg':'PREREG_R3_GEOMETRY.json','operation':'Expose missing force-work conjugate and consume BINDINGS r3','gate':'UNKNOWN_WORK_PORT','process_volume':'UNMEASURED'}],
        'auxiliary_needle_central_gate':'PASS_CONDITIONAL','auxiliary_digitization_gate':'FAIL',
        'strongest_control':'TIE','innovation_gate':'FAIL','parent_goal':'OPEN',
        'review_state':'PENDING_INDEPENDENT_REVIEW','scientific_admission':False})
    write('NEXT_WORK_PORT_R3.json',{
        'status':'UNACQUIRED_TEMPLATE_NOT_DATA','primary_source':None,'species':None,'specimen_id':None,
        'training_geometry':{'actual_contact_radius_m':None,'radius_uncertainty_m':None,'bevel_and_contact_map':None},
        'heldout_geometry':{'actual_contact_radius_m':None,'radius_uncertainty_m':None,'bevel_and_contact_map':None},
        'observations':{'synchronous_Ft_N':None,'synchronous_Fn_N':None,'blade_tangent_displacement_m':None,
            'blade_normal_displacement_m':None,'incremental_true_crack_area_m2':None,
            'recoverable_energy_change_J':None,'contact_dissipation_J':None,
            'independent_precleavage_process_work_J':None,'joint_uncertainty':None},
        'energy_convention':'Integral(Ft dx+Fn dz)=DeltaPsi+(Gamma0+B_before)*DeltaA+D_contact; crack area counts one projected interface, not arbitrarily doubled faces.',
        'missing_port_rule':'Return UNKNOWN, not a fitted Gamma, if work/area and separation ports absent.',
        'review_state':'PENDING_INDEPENDENT_REVIEW','scientific_admission':False})
    write('GRAPH_FEEDBACK_R3.json',{
        'target_id':'BT-CTX-SURG-INCISION','kind':'review','lane':'LANE_SURGICAL_INCISION','round':3,
        'review_state':'PENDING_INDEPENDENT_REVIEW','scientific_admission':False,'negative_result':True,
        'feedback':'Actual primary human-skin sharpness series acquired. Contact-width transfer fails 91.12%; added nose-pressure scenario fails 60%. True local radius remains unmeasured; same-F constant offset fails resolution scenarios. Published normal force is not conjugate to cutting motion; Gamma/work cannot be inferred. BINDINGS h_pz null and same-radius forcing counterexample consumed. No new coupled query gain; strongest same-information TIE.',
        'result_path':str(LANE/'RESULTS.md'),
        'binding_status':'PREPARED_LOCALLY_COORDINATOR_REVIEW_REQUIRED_NO_GRAPH_MUTATION'})
    next_text="# Next round — continue from R3\n\nEverything PENDING_INDEPENDENT_REVIEW. The parent target is OPEN; innovation gate FAIL, strongest equally informed control TIE. No new FE/LP is needed before a real work port exists.\n\n## Read first\n\n`CHECKPOINT_R3.json`, the latest section in `RESULTS.md`, `SOURCE_REVIEW_R3.md`, `SOURCE_TARGETS_R3_v2.json`, `r3/metrology_results.json`, `VERIFICATION_R3.json`, `NEXT_WORK_PORT_R3.json` and frozen `r3/BINDINGS_PORTS_R3_SNAPSHOT.json`. Read current steering and other lanes' latest ports; their empirical h_pz and skin intercept were null in R3. R1/R2's executable layer/history operators and all misses remain.\n\n## Next key operation\n\nChange the observation operation from normal load/damage depth to **the work's two force components and measured contact/crack kinematics**. The published NBS series is actual skin cutting with varied geometry, but Fn has no effect in direct work balance when vn=0. Acquire an original that has simultaneous Ft, Fn and displacement, real new crack area and actual local radius. Prioritise the Doran original through a legitimate author/institution host, later instrumented skin-cutting studies or available data files; never use secondarily cited forces as primary curves. DOI/full text remaining in `SOURCE_MANIFEST_R3.json` and saved errors are starting points, not evidence.\n\n1. Freeze the source's numerical data and their joint uncertainty before a new prediction. Fill `NEXT_WORK_PORT_R3.json` in a **new** file. Leave missing fields null. Matched material, body site/prestretch/direction, first cut versus later passage and speed must be explicitly stated. Distinguish a real cut from a crushed track. Residual damage depth is not automatically loaded penetration or new crack area.\n2. Integrate `Wext=∫Ft dx+∫Fn dz`. Keep a clear interval for `Geff=(Wext−ΔPsi−Dcontact)/ΔA`. Use a shared reference/material uncertainty for all geometries, not a separately fitted Γ for each observable. A full equally informed cohesive contact model is the control; count its acquisitions, solutions and fallback.\n3. If chemical Γ0 is requested, independent full first-pass pre-cleavage work is needed: `Γ0=Geff−B_before`. R2's Γ0+B0 null direction does not disappear with more radii without a mechanism that identifies B0. Synthetic 150 or the scissors value 2500 J/m² must not become measured Γ0.\n4. Freeze one training geometry and another radius held out from all parameter fits. Predict measured work/force, then a new gap through the preserved layer operator if the ports match. Run the actual control only when the data makes this meaningful. No history-free/10× gain without full cost and physical error budget.\n\n## If only NBS data remains available\n\nPreserve R3's negative results; do not make a new fit and call it prediction. A->other geometries misses by up to91,12%; B->K width model gives at least87% depth against40±2%. B+F12lb gives K64% under the RF25,4µm scenario, with shared rounding sensitivity62,5–65,5%. Actual RF is indeterminate. The23,787µm inverted from K is diagnostic, not a new test hit. Necessary local RF error is approximately0,503µm for an8pp depth budget **if** this closure and RK were exact and other errors zero; no achieved measurement precision is claimed.\n\nFurthermore, no constant F offset can match12/16/20lb with the declared resolution scenarios:8,000–8,091 /9,143–9,524 /14,857–15,238lb. This is a failed closure test that is not solved by changing RF. Moving forward requires variation/contact kinematics along the actual cutting sites and synchronous force–motion; a load-dependent fitted offset does not make the model predictive.\n\n## BINDINGS-kopplingen\n\nRequest deltaV/DeltaA for an actually advancing fracture front with full first-pass energy, anchor detachment/fibre rupture/delamination and the wedge's imposed contact work. Microcell contour length is not measured process-zone depth, partial hysteresis is not full fracture energy. The same radius can give different work. Do not transfer a porcine microcell scenario to human autopsy skin as a shared posterior.\n\n## Preservation and running\n\nAll R3 results are written exclusively. `geometry_r3.py --out <ny path i lane>` can be reproduced to a new directory; metrology_r3.py and verify_r3.py refuse existing default files. Reproduction requires a new output destination in a copy, no overwriting. Run with OMP/OPENBLAS/MKL/NUMEXPR_NUM_THREADS=2. `finalize_r3.py` must not be rerun: it rejects existing artifacts/sections.\n\nChange only the lane's own results or an allowed shared-data directory. No subagents, graph mutations, cloud or product code. Local GRAPH_FEEDBACK_R3.json awaits the coordinator's review.\n"
    section="\n## Night 1/10, Round 3\n\n**Capability before:** R2 had synthetic linked inventory equilibrium and a conditional purchase that separated Γ0/B0. Actual skin cutting with varied edge radius and a measured working port was missing. The parentage, cell damage and bleeding were open.\n\nDesired ability beyond current → predict a new blade radius cutting work and cut response out of common material ports, without full history calculation.\nConflicting concurrent requirements → real skin series, correct conjugate force, joint contact/process uncertainty, and a geometry not used for parameter fit.\nPrecise obstacle → published cutting data can measure normal load and residual damage depth but still lack cutting work; nominal radius does not identify actual contact curvature.\nChanged operation → acquisitiona primary skin series, try contact width, switch to separate nose work/pressure port and execute metrology inverse plus work two component balance.\n\n**New operation:** An actual published human skin series was found in [Sorrells & Berger, NBSIR73-262](https://www.govinfo.gov/content/pkg/GOVPUB-C13-7440d7745c063aa551900802884c8d51/pdf/GOVPUB-C13-7440d7745c063aa551900802884c8d51.pdf). It provides normal load and residual damage depth for different angles/radii on the same autopsy specimen. First, A's pressure/contact width transfer curve was frozen to B/E/F/G; B->K isolates the same angle with new nominal radius. After FAIL a new gate was executed: B's curve plus F12lb/65% anchors identify an additive normal force product kRF. When this prediction also failed, the uncertain nominal radius was replaced with a jointly identified product inverse and a precision requirement for the next purchase. The K-inverse is used diagnostically only; the original missed prediction remains. Five PREREG files and one source patch are preserved. This is analysis with omitted geometries in already published data, not a new blind or prospective study.\n\n**Reasonable Default Choices:** Sample18: Female Abdominal Abduction Skin69years,1,2mm,8days after death, subcutaneous tissue removed, wood mandrel,1inch/s. Piecewise linear interpolation only within measured support; full-thickness depth is censored. True-cut is separated from groove. ±0,5percentages table rounding and ±2 for K-reading are selected resolution scenarios, not repeatability ranges. Force uncertainty and pre-extension are missing. Contact width treats residual damage depth as loaded penetration under an explicit, untested closure; also depth of cut may include crushing. Q:s25,4µm is a scenario from the report's approximate radius in some places, not a measured homogeneous blade. 20% gate, two CPU threads, no new FE/LP or story solutions.\n\n**Strongest control:** Equally informed phenomenological contact/cohesive model gets same reference curve and extra F-anchor. The same contact/inverse/work algebra matches: **TIE**. No uniform fine FE was run in R3, and no measured errors or cost gains against such FE are claimed. R1/R2's full synthetic equilibrium and costs are preserved.\n\n**Faktiskt utfall:**\n\n| Trial | Actual numbers | Gate and scope |\n|---|---|---|\n| Primary skin series | B90°/nominaltR0 provides full transection at20lb; K90°/50,8µm gives40±2% at the same load/sample | Real geometry effect; cutting work **UNKNOWN** |\n| A->B/E/F/G common contact pressure model | Largest normal-force error91,12%; at least90,91% after common rounding scenarios | **FAIL** the construction, not general contact theory |\n| B->K, contact width only | Predicted depth≥87% on reference support against40%; central minimum relative error117,5% | **FAIL**, extrapolation avoided |\n| New additive nose port B+F12lb | kRF=8,045lb; K64% against40%,60% error. Shared rounding62,5–65,5% | **FAIL** under RF25,4µm/RK50,8µm-scenario |\n| Actual RF unidentified | K requires diagnosticallyRF23,787µm;32–48% depth corresponds toRF23,294–24,302µm | Empirical radius law **UNKNOWN**; no postfitPASS |\n| Required radius precision |15,891percent units/µm locally;8pp budget requires approximately0,503µm if other errors were zero | Conditional purchase requirement, no precision achieved |\n| F's additional loads, independent of the assumed RF | Offset interval12/16/20lb:8,000–8,091 /9,143–9,524 /14,857–15,238lb; empty common intersection | **FAIL** constant offset with selected resolution scenarios; statistical variance unknown |\n| [Shergold2005:s](https://www.fleckmech.org/papers/190.pdf)0,3->0,6mm needle, F=cD | Predicted3,9MPa against roughly read3,5; centrally11,43% error, extreme reading error38,33% | Central **PASS_CONDITIONAL**, robust **FAIL**; no scalpel radius/Γ0 |\n| Conjugate working port | K:Fn88,964N,d0,48mm,vt0,0254m/s; desktop tv0 givesFn·vn=0 | Ft, efficient cutting work and chemicalΓ0 **UNKNOWN** |\n\nAt the K depth gives effective work candidates150/2500J/m² surface termar0,072/1,2N in the tangent direction. These are explicit accounting evidence of a missing observation, no measured Γ or full resolved contact fields. `Fn/d≈185342J/m²` is not conjugate cutting work. Direct energy information from the normal force on the stationary path has rank0; an independent contact model may add information but is not identified here. Even measured tangential work leaves the Γ0+B0 zero direction if pre-cleavage/contact/release work is missing.\n\n**Cross linking and source corrections:** BINDINGS PORTS_R3 was consumed and hash snapshotted. Measured h_pz and the skin's shear intercept is null. SameR1µm but different wedge load gives their synthetic150 respectively379,275J/m²; small radius thus does not guarantee low process work. Cell contour14,701µm and partial loop is not a measured human process zone. The Comley2010 original gives trouser tear, not the suggested wedge-cut series; Comley2011 penetrates adipose, McCarthy measures elastomers, and Davis needle tip area. They are not mixed for skin calibration. Pereira/Doran full text downloads failed; secondary quoted figures are not used. See SOURCE_REVIEW_R3.md and SOURCE_MANIFEST_R3.json. F4lb was corrected before calculation to true-cut in v2; original freeze and all HTTP errors remain.\n\n**Verification and cost:** Independent analytical inverse gives K64%; circle/flange tangent and numerical radius derivative are correct, PDF/data hashes are correct, source non-monotonic F data and groove symbols are preserved. Common reference rounds and1000deterministic internal samples were checked; this is sensitivity analysis in binary64, no rounded biological error margin. New numerical calculations including verification have instrumented CPU lower bound **0,04318s**,0FE,0LP. Figure setup/rendering is counted separately in COST_ACCOUNTING_R3.json; source purchase, import, research time and full final cost is UNKNOWN, not zero. No hot linked issue or order of magnitude gain has been demonstrated. [Figure with primary data and missed predictions](r3/round3_overview.pdf).\n\n**Remaining Obstacles:** Correct force–motion couple, actual new crack area, current local edge curvature, joint contact/preload energy, and first pass irreversible interface work are missing. A radius series in another observable space does not solve this. Cell damage, bleeding, hemostasis and global history-free nonlinear query remain UNKNOWN/not demonstrated.\n\n**Next design change:** Acquire synchronously Ft/Fn/displacement and actual contact/crack kinematics, integrate the two-component work and separate ΔPsi,Dcontact and B_before with common error budget. Use NEXT_WORK_PORT_R3.json as blank observation list. Freeze one geometry for calibration and predict another with measured local radius; give no Γ-fit from normal load or groove depth. BINDINGS next gate must carry full first-pass energy for an advancing rupture front. Innovation gate **FAIL**, strongest control **TIE**, empirical skin cut work/Γ0 **UNKNOWN**. All PENDING_INDEPENDENT_REVIEW, no scientific admission.\n"
    with (LANE/'RESULTS.md').open('a') as f:f.write(section)
    (LANE/'NEXT_ROUND.md').write_text(next_text)
    status=read('WORK_STATUS.json')
    status.update(round=3,status='ROUND_COMPLETE',parent_goal_status='OPEN',
        review_state='PENDING_INDEPENDENT_REVIEW',scientific_admission=False,
        innovation_gate='FAIL',strongest_equally_informed_control_gate='TIE',
        empirical_skin_validation='UNKNOWN_FOR_WORK_AND_Gamma0',
        completed_capability='Primary published same-skin geometry/normal-force/injury-depth series acquired; heldout contact transfer, additional nose-pressure acquisition, diagnostic radius metrology and power-port obstruction executed',
        last_checkpoint='CHECKPOINT_R3.json',next_round='NEXT_ROUND.md',
        next_operation='Synchronous two-component blade work, actual contact curvature/true crack area, and independent release/contact/pre-cleavage acquisition',
        blockers=['Published normal force not conjugate to stationary-depth cutting motion; Ft absent',
            'Actual local edge radius and contact/loaded-depth mapping unmeasured',
            'Constant nose offset fails reported same-F load data under resolution scenarios',
            'Gamma0+B0 null remains without independent first-pass process work',
            'BINDINGS h_pz and empirical skin work null; biological/vascular closures unmeasured'],
        defaults='SOURCE_TARGETS_R3_v2.json',results='RESULTS.md',costs='COST_ACCOUNTING_R3.json',
        verification='VERIFICATION_R3.json',updated_utc=now)
    (LANE/'WORK_STATUS.json').write_text(json.dumps(status,indent=2)+'\n')
    checkpoint={'round':3,'status':'ROUND_COMPLETE_PARENT_GOAL_OPEN','gate':'FAIL','control_gate':'TIE',
        'empirical_intrinsic_crossprediction':'UNKNOWN','review_state':'PENDING_INDEPENDENT_REVIEW','scientific_admission':False,
        'main_results':['r3/geometry/summary.json','r3/geometry/power_port_results.json','r3/metrology_results.json','VERIFICATION_R3.json'],
        'portable_operators':['geometry_r3.py','metrology_r3.py','verify_r3.py'],
        'data':'SOURCE_TARGETS_R3_v2.json','preserve':['r1','r2','SOURCE_TARGETS_R3.json','SOURCE_TARGETS_R3_v2.json','r3'],
        'verification':'VERIFICATION_R3.json','costs':'COST_ACCOUNTING_R3.json',
        'next':'NEXT_ROUND.md','observations_template':'NEXT_WORK_PORT_R3.json','updated_utc':now}
    write('CHECKPOINT_R3.json',checkpoint)
    files=['RESULTS.md','WORK_STATUS.json','NEXT_ROUND.md','CHECKPOINT_R3.json','SOURCE_REVIEW_R3.md',
        'SOURCE_TARGETS_R3.json','SOURCE_TARGETS_R3_v2.json','SOURCE_MANIFEST_R3.json',
        'PREREG_R3.json','PREREG_R3_GEOMETRY.json','PREREG_R3_DATA_CORRECTION.json','PREREG_R3_METROLOGY.json','PREREG_R3_VERIFICATION.json',
        'geometry_r3.py','metrology_r3.py','verify_r3.py','finalize_r3.py','r3/geometry/summary.json',
        'r3/geometry/pressure_width_results.json','r3/geometry/tip_pressure_results.json','r3/geometry/needle_geometry_results.json',
        'r3/geometry/power_port_results.json','r3/geometry/peer_coupling.json','r3/metrology_results.json',
        'r3/BINDINGS_PORTS_R3_SNAPSHOT.json','VERIFICATION_R3.json','COST_ACCOUNTING_R3.json',
        'DECOMPOSITION_R3.json','ATTEMPTS_R3.json','NEXT_WORK_PORT_R3.json','GRAPH_FEEDBACK_R3.json',
        'r3/round3_overview.pdf','r3/round3_overview.png']
    write('night_rounds/r3.json',{
        'round':3,'operation':'Acquire primary human-skin edge geometry series -> heldout contact-width prediction -> additional nose-pressure acquisition -> diagnostic actual-radius metrology and conjugate-work port',
        'control':'Same-information phenomenological contact/cohesive and exact inverse/work algebra; TIE. Full uniform fine FE not executed.',
        'outcome':{
            'same_skin_series':'NBSIR73-262 specimen18, normal load / residual injury depth, no tangential work',
            **geom['prediction_results'],
            'width_max_min_force_error_joint_rounding':met['maximum_min_relative_force_error_with_joint_rounding'],
            'K_offset_joint_rounding_depth_band_pct':met['K_nominal_radius_ratio2_depth_band_pct'],
            'actual_radius_law':'UNKNOWN; conditional25.4um prediction FAIL retained',
            'diagnostic_RF_um':met['DIAGNOSTIC_required_F_radius_um'],
            'necessary_RF_precision_um_conditional':met['necessary_radius_error_um_if_other_errors_zero'],
            'constant_F_offset_gate':ver['constant_F_offset_gate'],
            'skin_effective_work_and_chemical_Gamma0':'UNKNOWN','BINDINGS_h_pz_m':None,
            'strongest_same_information_control':'TIE','new_FE_solves':0,'new_LP_solves':0,
            'numerical_instrumented_CPU_lower_bound_s':instrumented_cpu},
        'gate':'FAIL','obstacle':'Measured normal force is orthogonal to stationary cutting motion; true crack-area/work/contact/radius/pre-cleavage ports missing. Constant contact offset also fails reported loads; no local coupled-history gain demonstrated.',
        'next':'Acquire Ft/Fn and blade motion, actual curvature and advancing crack area, independent recoverable/contact/precleavage work; then freeze training geometry and predict heldout geometry with shared uncertainty.',
        'files':files+['ARTIFACT_MANIFEST_R3.json'],'review_state':'PENDING_INDEPENDENT_REVIEW',
        'scientific_admission':False,'updated_utc':now})
    immutable=[p for p in files if p not in ('RESULTS.md','WORK_STATUS.json','NEXT_ROUND.md')]
    immutable.append('night_rounds/r3.json')
    write('ARTIFACT_MANIFEST_R3.json',{
        'immutable_artifacts':[{'path':p,'sha256':hashlib.sha256((LANE/p).read_bytes()).hexdigest(),
            'bytes':(LANE/p).stat().st_size} for p in immutable],
        'mutable_handoff_files':['RESULTS.md','WORK_STATUS.json','NEXT_ROUND.md'],
        'primary_PDF_hashes':[{'manifest':a['manifest'],'sha256':a.get('sha256')} for a in acquisition if a.get('valid_pdf')],
        'review_state':'PENDING_INDEPENDENT_REVIEW','scientific_admission':False,'utc':now})
    print(json.dumps({'round':3,'status':'ROUND_COMPLETE','gate':'FAIL','control':'TIE',
                      'numerical_CPU_lower_bound_s':instrumented_cpu,'figure_CPU_s':figure_cpu,
                      'artifacts':len(immutable)},indent=2))

if __name__=='__main__':main()
