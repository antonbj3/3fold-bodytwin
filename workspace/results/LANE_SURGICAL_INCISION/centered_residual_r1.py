#!/usr/bin/env python3
"""Affine prestress compiler and topology/history discriminating experiment."""
import hashlib
import json
import time
import numpy as np
from incision_r1 import LANE,Membrane,write_json
from nonlinear_exterior_r1 import Nonlinear
from residual_compiler_r1 import compile_operator,query,monomials


def main():
    out=LANE/"r1/centered_residual"
    out.mkdir(exist_ok=False)
    expected=(LANE/"PREREG_R1_CENTERED_RESIDUAL.sha256").read_text().split()[0]
    if hashlib.sha256((LANE/"PREREG_R1_CENTERED_RESIDUAL.json").read_bytes()).hexdigest()!=expected:
        raise RuntimeError("centered prereg drift")
    all0=time.process_time()
    acq0=time.process_time()
    m=Membrane(64,32)
    nl=Nonlinear(m)
    initial=time.process_time()-acq0
    op,P,C,compileinfo=compile_operator(nl,[0]+list(range(1,32,2)),out,2,centered=True)
    with np.load(out/"operator_tier2.npz") as f:
        portable={k:f[k] for k in f.files}
    rows=[]
    original_kappa=nl.kappa.copy()
    cases=[(length,1.,"unseen_topology") for length in (8,16,24,32)]
    cases += [(length,mult,"matched_topology" if mult==1 else "new_history_matched_topology") for length in (1,9,17,25) for mult in (.95,1.,1.05)]
    for length,mult,kind in cases:
        cut=(m.x>=.008)&(m.x<.008+length*1e-3-1e-12)
        supplied_kappa=mult*original_kappa
        times=[]
        for _ in range(25):
            gap,z,info=query(portable,cut,supplied_kappa)
            times.append(info["cpu_s"])
        info.update(cut_length_mm=length,kappa_multiplier=mult,kind=kind,warm_median_cpu_s=float(np.median(times)))
        validate0=time.process_time()
        U=nl.uref+P @ z
        v,_=monomials(z,op["combos"])
        seam=op["k0"]*np.exp(-supplied_kappa/op["delta0"])*(~cut)*U[:m.m]
        rpoly=C @ v;rpoly[:m.m]+=seam
        rdirect,_=nl.bulk(U);rdirect[:m.m]+=seam
        # Residual absolute identity is meaningful also at near-zero residuals.
        info["polynomial_residual_max_abs_N"]=float(np.max(np.abs(rpoly-rdirect)))
        direct2=float(rdirect @ nl.lu0.solve(rdirect))
        info["direct_residual_dual_norm2"]=direct2
        info["compiled_norm_contains_direct_with_allowance"]=bool(direct2<=max(0,info["norm2"])+info["arithmetic_allowance_norm2"]+1e-20)
        nl.kappa=supplied_kappa
        truth,full=nl.full(cut,nl.uref)
        error=float(np.max(np.abs(gap-2*truth[:m.m])))
        info["actual_gap_error_m"]=error
        info["bound_contains_observed_error"]=bool(error<=info["gap_bound_m"]+1e-12)
        info["arithmetic_only_gap_floor_m"]=float(op["compliance"]*np.sqrt(info["arithmetic_allowance_norm2"]))
        info["full_control"]=full
        info["validation_cpu_s"]=time.process_time()-validate0
        if info["accepted"]:
            info["fallback"]=None;info["fallback_cpu_s"]=0.;info["returned_gap_error_m"]=error
        else:
            f0=time.process_time()
            fu,finfo=nl.full(cut,U)
            info["fallback"]=finfo
            info["fallback_cpu_s"]=time.process_time()-f0
            info["returned_gap_error_m"]=float(np.max(np.abs(2*fu[:m.m]-2*truth[:m.m])))
        info["charged_warm_cpu_s"]=info["warm_median_cpu_s"]+info["fallback_cpu_s"]
        info["warm_gain_vs_full_including_fallback"]=full["cpu_s"]/info["charged_warm_cpu_s"]
        rows.append(info)
        np.savez(out/f"validation_l{length}_k{mult}.npz",gap=gap,reference_gap=2*truth[:m.m],z=z)
    nl.kappa=original_kappa
    summary={"prereg_sha256":expected,"initial_setup_cpu_s":initial,"compile":compileinfo,"cases":rows,
        "total_cpu_s":time.process_time()-all0,"accepted_counts":{kind:sum(r["accepted"] for r in rows if r["kind"]==kind) for kind in set(r["kind"] for r in rows)},
        "all_bounds_contain_observed_errors":all(r["bound_contains_observed_error"] for r in rows),
        "all_polynomial_identities_max_abs_N":max(r["polynomial_residual_max_abs_N"] for r in rows),
        "strongest_control":"TIE with equally informed affine polynomial residual Gram certifier",
        "evolving_history_and_empirical_skin":"UNKNOWN","review":"PENDING_INDEPENDENT_REVIEW"}
    write_json(out/"summary.json",summary)
    print(json.dumps({"output":str(out),"accepted":summary["accepted_counts"],"all_bounds_contain_errors":summary["all_bounds_contain_observed_errors"],"max_polynomial_identity_abs_N":summary["all_polynomial_identities_max_abs_N"],"total_cpu_s":summary["total_cpu_s"]}),flush=True)


if __name__=="__main__":
    main()
