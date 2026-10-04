#!/usr/bin/env python3
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[key]='2'
import json
import time
import hashlib
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from layered_skin_r2 import LANE,write_json,Layers,verify

def read(path):return json.loads((LANE/path).read_text())

def main():
    t0=time.process_time(); now=datetime.now(timezone.utc).isoformat()
    w=read('r2/world_crossprediction.json'); b=read('r2/bridge_work.json'); acq=read('r2/acquisition_results.json')
    ext=read('r2/layer_extension_results.json'); targeted=read('r2/targeted_refinement_n128.json')
    layers=read('r2/layers_v2_n96/summary.json'); coarse=read('r2/layers_v2_n64/summary.json')
    energy=read('r2/energy_balance_check_corrected.json')
    verifications={str(n):read(f'r2/layers_v2_n{n}/verification.json') for n in (32,64,96)}
    verifications['angle45']=verify(Layers(32,16,45))
    fullref=[]
    for r,s in zip(coarse['records'],layers['records']):
        if r['kind']=='full_depth' and r['lengths_m'][1]>0:
            fullref.append(max(abs(a-b)/max(abs(b),1e-20) for a,b in zip(r['max_gap_m'],s['max_gap_m'])))
    verification={'review_state':'PENDING_INDEPENDENT_REVIEW','constitutive_directional_checks':verifications,
                  'max_global_equilibrium_residual_relative':max(layers['max_residual_relative'],ext['epidermis_only_n128']['residual_relative'],targeted['residual_relative']),
                  'minimum_final_full_model_detF':min(r['minJ'] for r in layers['records']),
                  'full_depth_gap_64_to96_max_relative_change':max(fullref),
                  'epidermis_gap_96_to128_relative_change':ext['epidermis_96_to128_relative_gap_change'],
                  'dermis8_gap_96_to128_relative_change':targeted['96_to128_relative_gap_change'],
                  'corrected_energy_identity_error_J':energy['max_absolute_error_J'],
                  'initial_energy_check_FAIL_preserved':True,
                  'skin_accuracy_gate':'UNKNOWN','continuum_enclosure':'UNKNOWN',
                  'interpretation':'Residuals/finite differences/refinement support numerical instrument only. Energy identity is bookkeeping, not independent empirical validation.'}
    write_json(LANE/'VERIFICATION_R2.json',verification)
    runnames=['layers_n32','layers_n64','layers_n96','layers_v2_n32','layers_v2_n64','layers_v2_n96']
    cpus={name:read(f'r2/{name}/summary.json')['CPU_s'] for name in runnames}
    cpus.update(layer_extension=ext['CPU_total_s'],targeted_refinement=targeted['CPU_total_s'],acquisition=acq['total_CPU_s'])
    cost={'instrumented_CPU_s':cpus,'instrumented_CPU_lower_bound_s':sum(cpus.values()),
          'preserved_initial_construction_cost_included':True,'equilibrium_query_count':150,
          'LP_solve_count':acq['LP_solves'],'peak_instrumented_RSS_KiB':ext['peak_RSS_KiB'],
          'threads_per_process':2,'individual_command_under_10min_2GiB':True,
          'global_sparse_equilibrium_solves_required':True,'local_state_update_gain':'NOT_DEMONSTRATED',
          'best_same_information_control_gain':1.0,
          'source_acquisition_IO_research_time_s':None,'process_setup_not_instrumented_s':None,
          'unknown_cost_policy':'unknown is not zero; no complete acquisition-cost advantage claimed'}
    write_json(LANE/'COST_ACCOUNTING_R2.json',cost)
    attempts=[
      {'operation':'same-skin cut/tear crossprediction from nominated micro/fiber ports','parents':['BINDINGS/MECHANICAL_RESULTS_R1','PISSARENKO_2020'],
       'outcome':'15kJ/m2 tear candidate FAIL 27.18%/50.63%; intrinsic skin intercept UNKNOWN. 150 is a synthetic pullout ceiling, not cleavage.',
       'obstacle':'No measured common cutting/intercept/microstructure parameter set','next':'integrate bridge work up to blade arrival','evidence':'r2/world_crossprediction.json'},
      {'operation':'pre-cleavage bridge-work truncation','outcome':'Conditional R=10um gives269.76J/m2 vs15150J/m2 total tear; half-pullout radius732um=73.2 fiber diameters.',
       'control':'ordinary same-information cohesive bridge law TIE','obstacle':'delta_b(R,friction,prestretch) and recruitment are unmeasured','next':'layer-release and independent work acquisition','evidence':'r2/bridge_work.json'},
      {'operation':'three coupled finite-strain vector layers','outcome':'Resolved full/partial-depth gap, anisotropy and interlayer prestress release; numerical scenario only',
       'preserved_failures':'initial oblique-family symmetry; first force quadrature length; mis-signed energy reporting check;64->96partial gap>2%',
       'next_executed':'mirror fiber families;physical crack lengths;correct energy check;targeted128 refinement',
       'control':'identical conventional FE TIE','evidence':'r2/layer_extension_results.json'},
      {'operation':'independent pre-cleavage work acquisition after radius/contact nullspace','outcome':'Γ interval0..1060.52J/m2 after recut;work20J gives119.875..180.125 FAIL;work5J gives134.875..165.125 PASS_SYNTHETIC',
       'control':'same positive LP inverse TIE','obstacle':'real microscopic work precision, missing layer closures and empirical release uncertainty',
       'next':'same-sample traction/opening work plus contact and prestretch observations','evidence':'r2/acquisition_results.json'}]
    write_json(LANE/'ATTEMPTS_R2.json',{'round':2,'review_state':'PENDING_INDEPENDENT_REVIEW','attempts':attempts})

    fig,axes=plt.subplots(2,2,figsize=(11,8),layout='constrained')
    ax=axes[0,0]
    ax.errorbar([0,1],[20.60,30.38],yerr=[2.15,4.90],fmt='o',label='Porcine skin, Table 1: mean ± SD')
    ax.axhline(15,color='tab:red',linestyle='--',label='Unfitted fiber scenario')
    ax.set(xticks=[0,1],xticklabels=['Mode III','Mode I'],ylabel='Tear work (kJ/m²)',title='World test: frozen tear candidate FAIL')
    ax.legend(fontsize=8)
    raw=np.load(LANE/'r2/bridge_work_raw.npz'); idx=np.argsort(raw['radius_m'])
    axes[0,1].semilogx(raw['radius_m'][idx][1:]*1e6,raw['total_cut_work_J_m2'][idx][1:]/1000)
    axes[0,1].axvline(10,color='gray',ls=':');axes[0,1].axhline(15.15,color='gray',ls='--')
    axes[0,1].set(xlabel='Blade radius (µm)',ylabel='Conditional cut work (kJ/m²)',title='Synthetic bridge closure: δ before blade = 2R')
    ax=axes[1,0]
    for key,label in [('96','Bonded layers'),('decoupled64','No layer transfer')]:
        rows=ext['coupled_dermis_trajectories'][key]
        ax.plot([r['lengths_m'][1]*1000 for r in rows],[r['max_gap_m'][1]*1000 for r in rows],'o-',label=label)
    ax.set(xlabel='Dermal cut length (mm)',ylabel='Dermal gap (mm)',title='Synthetic: cut epidermis, intact subcutis')
    ax.legend(fontsize=8)
    ax=axes[1,1]
    cases=[acq['radius_series_only'],acq['radius_series_plus_contact'],acq['independent_bridge_work_20J'],acq['refined_independent_bridge_work_5J']]
    for i,c in enumerate(cases):
        ax.plot([c['Gamma0_min_J_m2'],c['Gamma0_max_J_m2']],[i,i],lw=5,solid_capstyle='round')
    ax.axvline(150,color='gray',ls='--')
    ax.set(yticks=range(4),yticklabels=['Radius series','+ recut/contact','+ bridge work ±20','+ bridge work ±5'],xlabel='Feasible Γ₀ (J/m²)',title='Synthetic acquisition; positive constraints')
    fig.suptitle('INCISION R2 — empirical gate FAIL; model/acquisition results are scenarios',fontsize=12)
    for suffix in ('png','pdf'):
        path=LANE/f'r2/round2_overview.{suffix}'
        if path.exists():raise FileExistsError(path)
        fig.savefig(path,dpi=180)
    plt.close(fig)

    results=f'''
## Night 1/10, round 2

**Capability before:** R1 gave a synthetic scalar incision operator and numerical capsules; strongest equally informed control TIE, skin crossprediction UNKNOWN. Round guidance took priority over R1’s Green-column proposal: no further ROM/certification hunt.

Desired capability beyond the current state → cutting force, tear toughness and layer gap from the same binding/fibre parameters and uncertainty.
Simultaneous conflicting requirements → separate cleavage, process work, blade contact and prestress energy while the layers exchange load.
Precise obstacle → a radius-series intercept identifies the sum Γ0+B0, and the proposed microlevel150J/m² is not chemical cleavage energy.
Changed operation → integrate the fibre bridges’ work until blade arrival, couple to the layers’ released energy and test independent acquisition of this work.

**New operation:** Four frozen PREREG files. An alternative to the binary cutting/tearing hypothesis was executed: `T(δ)=2φτ(l−δ)_+/rf`, `W_before=Gp(2q−q²)`, `q=min(δ_b/l,1)`. A three-layer model with two displacement components, compressible neo-Hookean matrix, two recruited ±θ fibre families and interlayer slip actually solves equilibrium, gap and total stored energy. During dermal incision extension the epidermis is already cut and subcutis intact; all layers contribute to `F=h_dermis(Γ0+W_process)+dU_all/da+F_contact`. A positive interval inverse was then executed with radius series, separate contact and independent pre-cleavage work, including improved acquisition after the first missed precision.

**Reasonable defaults:** All layer stiffnesses and geometries are scenarios: h=0,1/1,5/3mm, E=100/200/20kPa, ν=0,3, total fibre coefficient20/100/0kPa, slip stiffness20MPa/m and0 as ablation. Domain64×32mm, prescribed biaxial stretch1,21 from a separate porcine reference, and1/1,10 as interventions. Mesh32×16,64×32,96×48; targeted128×64 refinement. Fibre-fracture ports are φ0,3,τ10kPa,l5mm,rf5µm,σf100MPa. `δ_b=2R` is an explicit synthetic closure, no measured skin-contact law. Unknown empirical layer energies, cell/vessel ports and contact work remain null.

**Strongest control:** Phenomenological cohesive model with separately calibrated cutting/tearing energy per observable, identical layer equilibrium and the same acquisition. Its fit to the tear means is exact on training data; no held-out result is claimed. The same bridge-law/FE/positive LP gets all parameters and matches: **TIE,1×**. The microparametrisation has not predicted both fracture modes against the world. No full global state/history-free capability or new warm cost gain has been demonstrated.

**Actual outcome:**

| Test | Numbers | Gate/scope |
|---|---|---|
| Unchanged fibre prediction15kJ/m² against skin | Mode III20,60±2,15; Mode I30,38±4,90kJ/m²; error27,18%/50,63% | **FAIL**, frozen20% requirement |
| Skin Γ0 and cut/tear ratio | No identified scalpel intercept/radius series in the same skin; no invented intercept | **UNKNOWN** |
| Origin of150J/m² | BINDINGS nominal fibril pullout12;150 is the `φσf l/2` ceiling under synthetic parameters | Not chemical binding prediction |
| Blade-arrival bridge | R=5/10µm gives209,94/269,76J/m²; half-pullout work only atR732µm,73,2 fibre diameters | Derived under scenario; radius transition in skin **UNKNOWN** |
| Layer gap at20mm,θ0 | Only epidermis cut:0,458mm; epi+dermis cut:dermis3,718mm; all cut:dermis6,150mm | Synthetic computed capability |
| Layer transfer | Epi+dermis,64 mesh: dermis gap3,685mm with slip coupling,5,873mm without | The mechanism affects the answer; no world match |
| Fibre orientation | Full incision20mm,96 mesh: dermis gap6,150/8,240/8,220mm at0/45/90° | Synthetic anisotropy |
| Dermal extension16→20mm | Released energy gives0,21639N; Γ0=150 gives surface term0,225N and remaining0,00861N before unknown process/contact work | Energy balance under fixed boundary displacement |
| Radius series and separate contact | Γ0 possible0–1060,52J/m²; Jacobian rank3/4 and exact null direction(1,−1,0,0) | **FAIL** identification |
| New independent work acquisition±20J/m² | Γ0=119,875–180,125; half-width30,125J/m² | **FAIL**, frozen30 budget |
| Work acquisition refined to±5J/m² | Γ0=134,875–165,125; half-width15,125J/m² | **PASS_SYNTHETIC**; realisable skin precision **UNKNOWN** |

The world anchor is [Pissarenko2020, original table1](https://meyersgroup.ucsd.edu/papers/journals/Meyers%20477.pdf), the same piglet for the tear modes. Microscopy shows bridging, delamination and fibre fracture; tearing is not solely pullout of intact fibres. The table’s2F/h is reconstructed within0,498% from rounded forces/thicknesses. This is reuse of already available observations, no new independent measurement. [Barnett2016](https://pure.psu.edu/en/publications/fracture-mechanics-model-of-needle-cutting-tissue/) separates porcine needle trials’ fracture/deformation/friction (means61/18/21%), but gives no scalpel-radius intercept for Γ0. [PMID32162009](https://pubmed.ncbi.nlm.nih.gov/32162009/) states knife thresholds equivalent to18,63N and6,86N respectively with different injury grades; these give neither a comparable energy ratio nor Γ0. [Smith2013](https://stars.library.ucf.edu/honorstheses1990-2015/1543/) verifies12µm scalpel tip in the abstract, but numerical PDF data could not be acquired. [Chu2019](https://www.nature.com/articles/s41598-019-52054-3) measures force in muscle and skin-graft morphology separately; muscle force gains are not transferred to skin.

**Numerical confidence and preserved misses:** Equilibrium residual≤{verification['max_global_equilibrium_residual_relative']:.3g}; directional gradient/Hessian checks≤{max(max(v['energy_gradient_directional_relative_error'],v['gradient_hessian_directional_relative_error']) for v in verifications.values()):.3g}. Full-depth gap change64→96 at most{100*max(fullref):.3f}%. Epidermis-only64→96 missed the2% requirement with2,676%;96→128 gives0,963%. Dermis8mm marginally missed2% at64→96 with2,081%;96→128 gives0,962%. No continuous or biological error bound is claimed. The first oblique fibre family and first force-length convention are preserved; the final model has symmetric families and physical tip separation. The energy check’s first report used the wrong sign on ΔU and gave0,00309J; the check is corrected without changed solutions/forces, now≤2,17×10⁻¹⁹J. It is a bookkeeping identity, not independent physical validation. All raw initial files remain.

**Remaining obstacles:** Molecular cleavage energy, blade contact’s `δ_b(R)` and pre-cleavage work are unidentified. Γ0 can change while B0 changes oppositely without any effect in radius series or post-cut gap. Zero radius does not guarantee zero radius-independent irreversible work. In the preregistered synthetic inverse budget,14,875J/m² remains for release/model error; the mesh difference in dermis16→20mm is1,040J/m², but this is a sensitivity estimate and no empirical error certificate. Layer-FE is a compressible2D membrane stack with assumed slip, prescribed incision and no 3D blade contact/fracture selection. Bulk stiffness and microscopic bridge work have no jointly measured parameter posterior. Cell damage, bleeding and hemostasis remain UNKNOWN.

**Cost:** All six preserved/final layer sweeps,15 dermal cases, two targeted128 cases, no-prestretch and24 LP solutions count. Instrumented CPU lower bound is **{cost['instrumented_CPU_lower_bound_s']:.2f}s**; largest instrumented RSS is **{ext['peak_RSS_KiB']/1024:.1f}MiB**. Source acquisition, I/O and researcher time are UNKNOWN and do not count as zero. Queries use full sparse Newton; no local-update gain. Strongest same-information control is1×. See COST_ACCOUNTING_R2.json and [the four-panel figure](r2/round2_overview.pdf).

**Next construction change:** Acquire same-specimen pre-cleavage traction/opening work and contact/release, together with cut/tear curves and microstructure. Compute `Γ0=F/h+G_release−G_contact−W_before` with a common error budget. Start by testing whether the extra±5J/m² work port is measurable; at nominal6MPa traction,5J/m² corresponds to opening uncertainty≤0,83µm if force error were zero. Executable priorities and starting commands are in NEXT_ROUND.md. The empirical innovation outcome is **FAIL**, strongest algorithm control **TIE**, skin’s joint crossprediction **UNKNOWN**. All results PENDING_INDEPENDENT_REVIEW; no scientific admission.
'''
    # Append only: retain the complete previous section verbatim.
    with (LANE/'RESULTS.md').open('a') as f:f.write(results)
    write_json(LANE/'CHECKPOINT_R2.json',{'round':2,'status':'ROUND_COMPLETE_PARENT_GOAL_OPEN','review_state':'PENDING_INDEPENDENT_REVIEW',
        'gate':'FAIL','control_gate':'TIE','empirical_intrinsic_crossprediction':'UNKNOWN',
        'portable_operators':['layered_skin_r2.py','layer_extension_r2.py','acquisition_r2.py'],
        'main_results':['r2/world_crossprediction.json','r2/layer_extension_results.json','r2/acquisition_results.json'],
        'preserve':['r1','r2/layers_n32','r2/layers_n64','r2/layers_n96','r2/layered_skin_r2_initial_snapshot.py','r2/layer_extension_r2_initial_snapshot.py'],
        'next':'NEXT_ROUND.md','updated_utc':now})
    # These two designated status/handoff files are intentionally updated. Their
    # exact R1 versions were archived before any mutation.
    status={'lane':'LANE_SURGICAL_INCISION','round':2,'status':'ROUND_COMPLETE','parent_goal_status':'OPEN',
        'review_state':'PENDING_INDEPENDENT_REVIEW','scientific_admission':False,
        'innovation_gate':'FAIL','strongest_equally_informed_control_gate':'TIE','empirical_skin_validation':'UNKNOWN',
        'completed_capability':'Synthetic coupled finite-strain layered cut/gap/release plus positive acquisition inverse; preserves world failure and micro-port correction',
        'last_checkpoint':'CHECKPOINT_R2.json','next_round':'NEXT_ROUND.md',
        'next_operation':'same-specimen pre-cleavage traction-opening work plus contact/prestress release to break Gamma0/B0 ambiguity',
        'blockers':['150 is not measured chemical cleavage','no same-skin cut-radius and tear data','delta_b and pre-cut process work unknown','layer constitutive/3D contact and biological ports unmeasured'],
        'defaults':'SOURCE_TARGETS_R2.json','results':'RESULTS.md','costs':'COST_ACCOUNTING_R2.json','updated_utc':now}
    with (LANE/'WORK_STATUS.json').open('w') as f:json.dump(status,f,indent=2);f.write('\n')
    nexttext='''# Next load-bearing operation after round 2

PENDING_CONSTRUCTION / PENDING_INDEPENDENT_REVIEW. Parents GRAPH_INVERSE, MULTIPHYSICS, BT-CTX-SURG-INCISION and BINDINGS R1; no graphs/sources have been changed. Empirical gate FAIL; best same-information control TIE; skin’s Γ0 and joint crossprediction UNKNOWN. New guidance always takes precedence over this handoff.

Capability → cutting/tearing force and layer gap from a common measured binding/fibre state.
Conflict → radius series is easier to measure than separate irreversible works, but the intercept contains both.
Obstacle → Γ0/B0 has an exact null direction. 150J/m² is short-fibril pullout’s synthetic upper ceiling; nominal12, not chemical cleavage. 15kJ/m² misses both skin tear modes.
Changed operation → same-specimen independent work on fibre bridges until blade arrival and contact/prestress energy; then an unfitted held-out blade prediction.

## Start directly

1. Read CHECKPOINT_R2.json, SOURCE_INTERPRETATION_CORRECTIONS_R2.md, DECOMPOSITION_R2.json and r2/acquisition_results.json. Preserve all R1/R2 files. Do not continue Green/residual certification as the main track.
2. Primary acquisition: published **skin**, same specimen/cohort and orientation, cut/tear work or complete force curves plus thickness, radius/angle, local opening and traction or another independent process-work observation. Measure subsequent reinsertion for contact but do not assume it subtracts constant pre-cut damage or restores prestress. Data from muscle/gel, scissor/needle/scalpel and other species are kept separate. Smith2013 PDF had access errors; no numerical data were acquired. The Barnett file with .pdf suffix is an HTML error response, not evidence.
3. Freeze new world ordinates before model/fit. Report the exact missing observation if matched open data are not found. Sources are in SOURCE_TARGETS_R2.json; Pissarenko Table1 is already verified and is not a new held-out test.
4. Execute the work port’s measurability test: `Gamma0=F_cut/h + G_release - G_contact - W_before`. R2 synthetic precision: cutoff±10J/m², contact coefficient±5, bridge-work±5 gaveΓ134.875..165.125 atR1µm. With bridge-work±20 it gave119.875..180.125 and FAIL. ForΓ150 with20% budget,14.875J/m² remains for ALL release/closure/geometry errors. Nominal T6MPa requires opening uncertainty≤0.83µm for5J/m² if traction were exact; actual precision UNKNOWN. Use deterministic intervals or measured common covariance, no independent marginals as an invented posterior.
5. Mechanical consumer: `layered_skin_r2.py` has Layers(nx,ny,angle,coupling), energy(U,hessian=True), solve([a_epi,a_dermis,a_subcutis],stretch=1.21). It is full sparse Newton in synthetic compressible2D layers with±theta fibres and reversible slip. Prescribed incision and additiveΓ mean that post-cut gap alone does not identifyΓ. No warm global-state-free capability. Existing final results r2/layers_v2_n96 and r2/layer_extension_results.json may be reused; importing the module writes nothing.
6. Prioritise actual blade arrival `delta_b(R,wedge,friction,prestretch)` and recruitment/damage history over a larger mesh. R2’s explicit `delta_b=2R` gives half pullout work atR732µm,73 fibre diameters: fibre diameter alone does not derive the transition. Replace this closure with an independently measured/local contact law and test held-out sharpness without a new Γ fit. If it does not exist, parametrised intervals with UNKNOWN world gate.
7. Matched control gets the same microscopy/work/contact data and a mode-separated cohesive law with full work curve. Measure prediction error and count of **independent** empirical acquisitions, not just training fit. Algebraic identity is TIE. Actual binding capability requires exceeding even this strong control contract without extra free fitting.

## Numerical raw data and limitations

Derma16→20mm over intact subcutis releases0.21639N in all layers; Γ150*1.5mm=0.225N leaves0.00861N before unknown process/contact. 64→96 release difference1.040J/m² is a sensitivity, no rigorous continuous or empirical error bound. Epi-only96→128 gap change0.963%,dermis8mm0.962%; stop routine mesh increases. Final verification is in VERIFICATION_R2.json. Initial raw files, oblique-family/length convention and mis-signed energy report are preserved; corrected_energy_check is bookkeeping, no biological validation.

Reproduce new data in a NEW output directory: `OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 NUMEXPR_NUM_THREADS=2 python results/LANE_SURGICAL_INCISION/layered_skin_r2.py --nx 64 --out results/LANE_SURGICAL_INCISION/r3_layers_probe`. Do not run acquisition/finalize over R2; they write exclusive R2 files. For a new acquisition round, copy the module/PREREG to a new suffix and change the destination before running. COMMANDS_R2.md has the original commands.

Biological cell, vessel, bleeding and hemostasis ports are UNKNOWN. No clinical admission. Only this lane/shared-data is writable; no subagents/cloud/queues/source-graph changes. Graph feedback GRAPH_FEEDBACK_R2.json is local material for the coordinator’s ordinary review binding.
'''
    with (LANE/'NEXT_ROUND.md').open('w') as f:f.write(nexttext)
    roundreport={'round':2,'operation':'skin cut/tear crossprediction -> pre-cleavage bridge-work -> coupled vector layered prestress release -> independent microscopic-work acquisition',
        'control':'per-observable calibrated cohesive energies; same-information bridge/FE/positive-LP inverse',
        'outcome':{'skin_tear_relative_errors':[w['mode_III_relative_error'],w['mode_I_relative_error']],
            'skin_intrinsic_intercept':'UNKNOWN','paired_empirical_crossprediction':'UNKNOWN','micro150_origin':'synthetic short-fibril sliding ceiling, nominal12',
            'radius_plus_contact_Gamma0_interval_J_m2':[0,acq['radius_series_plus_contact']['Gamma0_max_J_m2']],
            'refined_synthetic_Gamma0_interval_J_m2':[134.875,165.125],'work20_gate':'FAIL','work5_gate':'PASS_SYNTHETIC',
            'synthetic_layered_gap_and_release':'EXECUTED','same_information_control':'TIE',
            'instrumented_CPU_lower_bound_s':cost['instrumented_CPU_lower_bound_s']},
        'gate':'FAIL','obstacle':'Micro cleavage/work/contact/recruitment and common skin ports unmeasured; intrinsic/process constant null direction; radius alone lacks blade-arrival map',
        'next':'independent same-specimen traction-opening pre-cleavage work and contact/prestress-release acquisition; heldout blade-sharpness prediction',
        'files':['RESULTS.md','WORK_STATUS.json','NEXT_ROUND.md','CHECKPOINT_R2.json','SOURCE_TARGETS_R2.json','DECOMPOSITION_R2.json',
                 'PREREG_R2.json','PREREG_R2_ACQUISITION.json','PREREG_R2_LAYER_EXTENSION.json','SOURCE_INTERPRETATION_CORRECTIONS_R2.md',
                 'r2/world_crossprediction.json','r2/bridge_work.json','r2/acquisition_results.json','r2/layers_v2_n96/summary.json',
                 'r2/layer_extension_results.json','r2/energy_balance_check_corrected.json','r2/targeted_refinement_n128.json',
                 'VERIFICATION_R2.json','COST_ACCOUNTING_R2.json','ATTEMPTS_R2.json','r2/round2_overview.pdf','GRAPH_FEEDBACK_R2.json'],
        'review_state':'PENDING_INDEPENDENT_REVIEW','scientific_admission':False,'updated_utc':now}
    write_json(LANE/'night_rounds/r2.json',roundreport)
    write_json(LANE/'GRAPH_FEEDBACK_R2.json',{'target_id':'BT-CTX-SURG-INCISION','kind':'review','lane':'LANE_SURGICAL_INCISION','round':2,
        'review_state':'PENDING_INDEPENDENT_REVIEW','scientific_admission':False,'negative_result':True,
        'feedback':'Unfitted 15kJ skin tear candidate fails both exact modes; nominated150 micro work is a synthetic ceiling. New blade-arrival, layered release and acquisition constructions executed; real skin cleavage unknown, strongest same-information TIE.',
        'result_path':str(LANE/'RESULTS.md'),'binding_status':'PREPARED_LOCALLY_COORDINATOR_REVIEW_REQUIRED_NO_GRAPH_MUTATION'})
    print(json.dumps({'cost_CPU_lower_bound_s':cost['instrumented_CPU_lower_bound_s'],'full_gap_refinement_max':max(fullref),'finalize_CPU_s':time.process_time()-t0}))

if __name__=='__main__':main()
