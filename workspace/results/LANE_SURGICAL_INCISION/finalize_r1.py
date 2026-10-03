#!/usr/bin/env python3
"""Lane-local final artifact synthesis from executed, preserved raw experiments."""
import ast
import hashlib
import json
import platform
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import scipy
from incision_r1 import LANE,ROOT,write_json
from compiled_query_r1 import evaluate


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def read(p):
    return json.loads((LANE/p).read_text())


def main():
    stamp=datetime.now(timezone.utc).isoformat()
    linear=read("r1/physical_patch/summary.json")
    nonlinear=read("r1/nonlinear/summary.json")
    compiler=read("r1/residual_compiler/summary.json")
    centered=read("r1/centered_residual/summary.json")
    hx=read("r1/HX_CHAIN_BOUND_REFINEMENT.json")
    checks={}
    for pre in ("PREREG_R1","PREREG_R1_EXTENSION","PREREG_R1_RESIDUAL_COMPILER","PREREG_R1_CENTERED_RESIDUAL"):
        checks[pre+"_hash"]=sha(LANE/(pre+".json"))==(LANE/(pre+".sha256")).read_text().split()[0]
    checks["exact_linear_all_meshes"]=all(r["methods"]["condensed"]["comparison"]["gate"]=="PASS" for r in linear["mechanical_sizes"])
    checks["exact_wake_capsule_all_meshes"]=all(r["methods"]["capsule"]["comparison"]["gate"]=="PASS" for r in linear["mechanical_sizes"])
    checks["frozen_wake_failure_preserved"]=all(r["methods"]["frozen_wake"]["comparison"]["gate"]=="FAIL" for r in linear["mechanical_sizes"])
    checks["nonlinear_no_unbudgeted_accept"]=all(not p["certificate_accepted"] and p["fallback"] is not None and p["returned_gap_max_error_m"]<=1e-7 for r in nonlinear["rows"] for p in r["patches"])
    checks["nonlinear_bounds_contain_errors"]=all(p["bound_contains_observed_error"] for r in nonlinear["rows"] for p in r["patches"])
    checks["centered_all_bounds_contain_errors"]=centered["all_bounds_contain_observed_errors"]
    checks["centered_polynomial_identity_abs"]=centered["all_polynomial_identities_max_abs_N"]<1e-10
    checks["centered_gram_numeric_allowance_contains_direct"]=all(r["compiled_norm_contains_direct_with_allowance"] for r in centered["cases"])
    checks["centered_no_false_accept"]=all(not r["accepted"] or r["actual_gap_error_m"]<=1e-7 for r in centered["cases"])
    checks["centered_all_rejections_have_costed_fallback"]=all(r["accepted"] or (r["fallback"] is not None and r["fallback_cpu_s"]>0 and r["returned_gap_error_m"]<=1e-7) for r in centered["cases"])
    with np.load(LANE/"r1/centered_residual/operator_tier2.npz") as f:
        op={k:f[k] for k in f.files}
    q9=evaluate(op,9)
    q8=evaluate(op,8)
    checks["portable_support9_accept"]=q9["status"]=="CONDITIONAL_DISCRETE_ACCEPT"
    checks["portable_new_support8_reject"]=q8["status"]=="REQUIRES_FULL_FALLBACK" and q8["gap_m"] is None
    checks["portable_query_no_global_solve_or_reconstruction"]=all(q["global_solves"]==q["global_state_reconstructions"]==0 for q in (q9,q8))
    checks["physical_observation_strip_not_censored"]=all(not r["methods"]["condensed"]["rows"][-1]["zone_censored_at_patch_edge"] for r in linear["mechanical_sizes"])
    checks["HX_reported_normalized_states_bounded"]=all(0<=s["reported_projected_state"][name]<=1 for r in hx["cases"] for s in r["snapshots"] for name in ("q","o","A","Gly","M","D","E","F"))
    for p in LANE.glob("*.py"):
        ast.parse(p.read_text(),filename=str(p))
    checks["all_lane_python_ast_parses"]=True
    if not all(checks.values()):
        raise RuntimeError({k:v for k,v in checks.items() if not v})
    write_json(LANE/"VERIFICATION_R1.json",{"checks":checks,"all_passed":True,"meaning":"Implementation/numerical discriminants passed; innovation and empirical gates did not become PASS","updated_utc":stamp})
    stages=[]
    for path in ("r1/first/summary.json","r1/physical_patch/summary.json","r1/nonlinear/summary.json","r1/residual_compiler/summary.json","r1/centered_residual/summary.json"):
        s=read(path)
        stages.append({"stage":path,"measured_instrumented_CPU_s":s["total_cpu_s"],"includes":"Executed acquisition/queries/validation/failures/fallbacks within script timer; excludes imports, process launch and some final serialization"})
    stages.append({"stage":"r1/HX_CHAIN_BOUND_REFINEMENT.json","measured_instrumented_CPU_s":hx["cpu_s"],"includes":"refined Q036 RHS solves"})
    linear_cost=[]
    for r in linear["mechanical_sizes"]:
        linear_cost.append({"nx":r["nx"],"global_dofs":r["global_free_dofs"],"port_dofs":r["seam_dofs"],"operator_bytes":r["operator_bytes"],"setup_CPU_s":r["build_cpu_s"]+r["schur_and_observation_acquisition_cpu_s"],"full22step_CPU_s":r["methods"]["full"]["cpu_s"],"condensed22step_CPU_s":r["methods"]["condensed"]["cpu_s"],"capsule22step_CPU_s":r["methods"]["capsule"]["cpu_s"],"gain_vs_repeated_full":r["methods"]["condensed"]["speedup_vs_repeated_full_query_only"],"same_information_best_control_gain":1.0})
    cost={"stages":stages,"instrumented_experiment_CPU_lower_bound_s":sum(s["measured_instrumented_CPU_s"] for s in stages),"linear":linear_cost,"nonlinear_patch_total_fallback_CPU_s":sum(p["fallback_cpu_s"] for r in nonlinear["rows"] for p in r["patches"]),"nonlinear_patch_total_setup_CPU_s":sum(p["setup_cpu_s"] for r in nonlinear["rows"] for p in r["patches"]),"nonlinear_patch_total_cert_CPU_s":sum(p["certificate_cpu_s"] for r in nonlinear["rows"] for p in r["patches"]),"centered_setup_CPU_s":centered["initial_setup_cpu_s"]+centered["compile"]["total_acquisition_compile_cpu_s"],"centered_validation_CPU_s":sum(r["validation_cpu_s"] for r in centered["cases"]),"centered_fallback_CPU_s":sum(r["fallback_cpu_s"] for r in centered["cases"]),"centered_query_benchmark_repetitions_per_case":25,"centered_operator_bytes":centered["compile"]["operator_bytes"],"centered_global_Gram_acquisition_RHS_solves":centered["compile"]["gram_global_rhs_solves"],"best_control_break_even":None,"best_control_break_even_reason":"Identical eligible same-information conventional operator: no positive cost difference; no claimed finite break-even","complete_cold_cost_advantage":"NOT_SHOWN","unmeasured_costs":"Process launch/imports/some serialization, literature retrieval, human-equivalent code construction, all broader physical model calibration: UNKNOWN; never zero","max_observed_mechanical_run_RSS_MiB":linear["max_rss_MiB"],"RSS_scope":"Measured by mechanical run; other process peaks not instrumented. Their largest compiler arrays were bounded by known dimensions, no heavy job used.","review":"PENDING_INDEPENDENT_REVIEW"}
    write_json(LANE/"COST_ACCOUNTING_R1.json",cost)
    references={"frozen_comparison_prereg":"PREREG_R1.json","gel_pdf":{"path":"literature/goda_2024.pdf","sha256":sha(LANE/"literature/goda_2024.pdf"),"url":"https://arxiv.org/abs/2412.02021","status":"PRIMARY_PUBLIC_PDF_ACQUIRED"},"post_computation_context_only":[{"url":"https://pmc.ncbi.nlm.nih.gov/articles/PMC4520804/","measurement":"Porcine in vivo/ex vivo control-patch prestretch, Table1: average stretches1.21/1.21, average area stretch1.44; single-animal prototype, excision-dependent stress-free reference","status":"EXTERNALLY_MEASURED_WITH_SCOPE","not_fitted_or_used_as_gate":True},{"url":"https://pubmed.ncbi.nlm.nih.gov/21605167/","measurement":"Human ex vivo standardized linear excisional wound tool uses2mm blade separation, wound gap measured histologically; not spontaneous scalpel-gap data","status":"PRIMARY_PUBLIC_MEASUREMENT_DESIGN","not_fitted_or_used_as_gate":True}],"unavailable_data":"No matched public skin radius-force-depth/gap/cell-zone/vessel data acquired; not replaced by generated points"}
    write_json(LANE/"WORLD_REFERENCES_R1.json",references)
    parents=[ROOT/"results/BT-HX-Q033/model.py",ROOT/"results/BT-HX-Q033/PREREG.md",ROOT/"results/BT-HX-Q036/model.py",ROOT/"results/BT-HX-Q049/model.py",ROOT/"tasks/build_night/COMMON.md",ROOT/"tasks/build_night/steer/LANE_SURGICAL_INCISION.md",ROOT/"tasks/lanes/LANE_SURGICAL_INCISION.md"]
    private=Path("source_repository/scripts/msk")
    parents += [private/name for name in ("skin_pulp_mechanics.py","collagen_triple_helix_thermal_stability.py","tendon_collagen_hierarchical_mechanics.py","myofascial_transmission.py","coagulation_hemostasis.py","platelet_hemostasis.py","wound_healing_cascade.py","acute_phase_inflammation.py","bone_fracture_toughness_lefm.py","skin_barrier_tewl.py")]
    write_json(LANE/"SOURCE_MANIFEST_R1.json",{"sources":[{"path":str(p),"sha256":sha(p),"access":"read-only"} for p in parents],"software":{"python":platform.python_version(),"numpy":np.__version__,"scipy":scipy.__version__},"threads":{"OMP":2,"OPENBLAS":2,"MKL":2,"NUMEXPR":2},"forbidden_internal_data_read":False,"updated_utc":stamp})
    attempts=[
        {"operation":"exact fine exterior Schur with cohesive max-history","evidence":"r1/physical_patch/summary.json","outcome":"Numerical PASS/TIE","obstacle":"linear exterior and fixed seam; standard equally informed substructuring matches","next":"nonlinear exterior"},
        {"operation":"freeze completed wake, then execute exact wake capsule","evidence":"r1/first/summary.json","outcome":"freeze FAIL, capsule same-model PASS/TIE","obstacle":"wake response grows and nonlinear remaining seam not constant size","next":"physical nonlinear patch with charged rebuild"},
        {"operation":"same energy tear-to-cut then cleavage/failure-work ports","evidence":"r1/physical_patch/world_transfer.json","outcome":"single energy FAIL62.5%, separate-assay slope PASS10.03%","obstacle":"published gel relation; skin strength and FPZ length absent","next":"measured shared skin ports"},
        {"operation":"distributed convex nonlinear exterior with tip and wake patch","evidence":"r1/nonlinear/summary.json","outcome":"12/12 cheap patches FAIL; charged full fallback","obstacle":"nonlocal nonlinear deformation plus global residual certification","next":"polynomial residual Gram compiler"},
        {"operation":"cubic residual Gram with 5/17 global modes","evidence":"r1/residual_compiler/summary.json","outcome":"no query global solve;0/8 new supports accepted","obstacle":"new topology subspace error and artificial top-boundary cancellation","next":"affine physical prestress center"},
        {"operation":"centered residual capsule, topology/history discriminant","evidence":"r1/centered_residual/summary.json","outcome":"4/4 acquired supports accepted;0/4 new supports;0/8 changed histories;all bounds contain errors","obstacle":"new tip cusp and branch/history sensitivity; no strong-control gain","next":"local Green support enrichment with shared historical residual sector"}
    ]
    write_json(LANE/"ATTEMPTS_R1.json",{"parent_ids":["BT-HX-Q033","BT-HX-Q036","BT-HX-Q049","BT-CTX-SURG-INCISION","GRAPH_INVERSE","MULTIPHYSICS"],"attempts":attempts,"review":"PENDING_INDEPENDENT_REVIEW"})
    checkpoint={"round":1,"status":"ROUND_COMPLETE_PARENT_RESEARCH_OPEN","review":"PENDING_INDEPENDENT_REVIEW","gate":"FAIL","strongest_matched_control":"TIE","original_lane_was_new":True,"current_operation":"affine prestress residual Gram with frozen cohesive-history support","current_operator":"r1/centered_residual/operator_tier2.npz","portable_consumer":"compiled_query_r1.py","actual_fallback":"centered_residual_r1.py","mechanical_chain":"incision_r1.py / r1/physical_patch/summary.json","physiology_refinement":"r1/HX_CHAIN_BOUND_REFINEMENT.json","next":"NEXT_ROUND.md","preserve":["r1/first/","r1/physical_patch/","r1/nonlinear/","r1/residual_compiler/","r1/centered_residual/"],"updated_utc":stamp}
    write_json(LANE/"CHECKPOINT_R1.json",checkpoint)
    status={"lane":"LANE_SURGICAL_INCISION","round":1,"status":"ROUND_COMPLETE","parent_goal_status":"OPEN","review_state":"PENDING_INDEPENDENT_REVIEW","scientific_admission":False,"innovation_gate":"FAIL","strongest_equally_informed_control_gate":"TIE","empirical_skin_validation":"UNKNOWN","completed_capability":"Runnable synthetic incision->gap+mechanical strain/cohesive zone->Q033 exudate and Q036 normalized initial-state/history; frozen-history nonlinear residual capsule with actual fallback","last_checkpoint":"CHECKPOINT_R1.json","next_round":"NEXT_ROUND.md","next_operation":"shared deformation basis plus acquired local new-tip Green columns and signed history residual closure","blockers":["new topology outside17D subspace","cohesive branch/history enclosure absent","vector layered biological model and measured shared ports absent"],"defaults":"PREREG_R1.json / DECOMPOSITION_R1.json","results":"RESULTS.md","costs":"COST_ACCOUNTING_R1.json","updated_utc":stamp}
    (LANE/"WORK_STATUS.json").write_text(json.dumps(status,indent=2)+"\n")
    feedback={"target_id":"BT-CTX-SURG-INCISION","lane":"LANE_SURGICAL_INCISION","kind":"review","result_file":"results/LANE_SURGICAL_INCISION/RESULTS.md","result_sha256":sha(LANE/"RESULTS.md"),"review_state":"PENDING_INDEPENDENT_REVIEW","scientific_admission":False,"negative_result":True,"outcome":"Broad no-global-recompute/large-strong-control-gain gateFAIL; discrete conditional suboperators and known-assay two-port slope test useful, strongest mathematical controlsTIE","measured_quantity":"synthetic gap errors, cohesive/strain proxies, CPU costs; published gel regression transfer","units":"m,1,N,CPU-s,J/m2,J/m3,Pa","uncertainty":"unmeasured biological closures; floating allowance not independently verified; no continuum or evolving-cohesive branch certificate","population_regime":"synthetic scalar P1 plus frozen-history convex nonlinear bulk; PAAm published assays are separate","binding_status":"READY_FOR_COORDINATOR_BINDING_NOT_DISPATCHED","scope_reason":"COMMON restricts all writes to lane/shared-data; ./graph dispatch and feedback would write tasks/graph_runs outside allowed scope"}
    write_json(LANE/"GRAPH_FEEDBACK.json",feedback)
    round_record={"round":1,"operation":"exact incision Schur/history -> completed-wake elimination -> nonlinear exterior patch -> cubic residual Gram -> affine prestress recentering and topology/history discriminant","control":"same-law uniform fine P1 FE plus eligible identical same-information substructuring / polynomial residual Gram certifier; vector layered skin FE not yet implemented","outcome":{"numerical_linear":"PASS in same discrete model","frozen_wake_max_gap_error_m":max(r["methods"]["frozen_wake"]["comparison"]["max_gap_error_m"] for r in linear["mechanical_sizes"]),"nonlinear_patch_accepts":0,"nonlinear_patch_cases":12,"centered_accepted":centered["accepted_counts"],"centered_cases":16,"all_centered_bounds_contain_errors":True,"single_toughness_transfer_relative_error":.625,"separate_failure_work_slope_relative_error":3108/31000,"empirical_skin":"UNKNOWN","strongest_control":"TIE","global_compiled_query_solves":0,"instrumented_experiment_CPU_lower_bound_s":cost["instrumented_experiment_CPU_lower_bound_s"]},"gate":"FAIL","obstacle":"new tip topology and changed cohesive history outside global trial subspace; capsule arithmetic/branch/continuum certification and empirical layered skin ports incomplete; strongest equal-information methods match","next":"local new-tip Green support enrichment with shared historical nonlinear residual sector and full charged control; see NEXT_ROUND.md","files":["RESULTS.md","WORK_STATUS.json","NEXT_ROUND.md","CHECKPOINT_R1.json","DECOMPOSITION_R1.json","ATTEMPTS_R1.json","COST_ACCOUNTING_R1.json","VERIFICATION_R1.json","SOURCE_REVIEW_R1.md","SOURCE_MANIFEST_R1.json","WORLD_REFERENCES_R1.json","GRAPH_FEEDBACK.json","COMMANDS_R1.md","compiled_query_r1.py","r1/physical_patch/summary.json","r1/nonlinear/summary.json","r1/residual_compiler/summary.json","r1/centered_residual/summary.json","r1/centered_residual/operator_tier2.npz","r1/HX_CHAIN_BOUND_REFINEMENT.json"],"review":"PENDING_INDEPENDENT_REVIEW","updated_utc":stamp}
    write_json(LANE/"night_rounds/r1.json",round_record)
    artifacts=[]
    for p in sorted(LANE.rglob("*")):
        if p.is_file() and p.name not in ("ARTIFACT_MANIFEST_R1.json","codex_r1.log") and "__pycache__" not in str(p):
            artifacts.append({"path":str(p.relative_to(LANE)),"sha256":sha(p),"bytes":p.stat().st_size})
    write_json(LANE/"ARTIFACT_MANIFEST_R1.json",{"artifacts":artifacts,"mutable_metadata":["WORK_STATUS.json","RESULTS.md","NEXT_ROUND.md"],"updated_utc":stamp})
    print(json.dumps({"verification":checks,"gate":"FAIL","strongest_control":"TIE","artifacts":len(artifacts),"measured_experiment_CPU_lower_bound_s":cost["instrumented_experiment_CPU_lower_bound_s"]}))


if __name__=="__main__":
    main()
