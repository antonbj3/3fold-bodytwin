#!/usr/bin/env python3
"""Exact polynomial residual algebra; discrete convex conditional certificate.

Compiler is intrusive and conventional; it does not assert methodological
novelty. The new capability tested here is a nonlinear exterior residual query
without full-field reconstruction, under new cohesive seam supports.
"""
import os
for key in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","NUMEXPR_NUM_THREADS"):
    os.environ[key]="2"
import hashlib
import itertools
import time
import numpy as np
from incision_r1 import LANE,Membrane,write_json
from nonlinear_exterior_r1 import Nonlinear


def monomials(z,combos):
    a=np.r_[1.,z]
    vals=a[combos]
    v=np.prod(vals,axis=1)
    deriv=np.zeros((len(combos),len(z)))
    for col in range(3):
        remaining=[i for i in range(3) if i!=col]
        term=vals[:,remaining[0]]*vals[:,remaining[1]]
        valid=combos[:,col]>0
        np.add.at(deriv,(np.flatnonzero(valid),combos[valid,col]-1),term[valid])
    return v,deriv


def query(op,cut,kappa):
    """Only compiled small matrices, seam data and reduced coordinates."""
    t0=time.process_time()
    PB,RC=op["PB"],op["RC"]
    diag=op["k0"]*np.exp(-kappa/op["delta0"])*(~cut)
    K=PB.T @ (diag[:,None]*PB)
    z=op["zref"].copy()
    for it in range(50):
        v,dv=monomials(z,op["combos"])
        r=RC @ v+K @ z
        rel=np.linalg.norm(r)/op["reduced_scale"]
        if rel<1e-10:
            break
        step=np.linalg.solve(RC @ dv+K,-r)
        alpha=1.0
        for _ in range(20):
            nv,_=monomials(z+alpha*step,op["combos"])
            nr=RC @ nv+K @ (z+alpha*step)
            if np.linalg.norm(nr)<np.linalg.norm(r):
                z=z+alpha*step
                break
            alpha*=.5
        else:
            raise RuntimeError("reduced nonlinear line search failed")
    else:
        raise RuntimeError("reduced nonlinear Newton failed")
    # Supported seam/bulk cross term is part of the joint certificate.
    s=diag*(PB @ z)
    bulk=float(v @ op["H"] @ v)
    cross=float(2*s @ op["D"] @ v)
    seam=float(s @ op["G"] @ s)
    norm2=bulk+cross+seam
    absolute=float(np.abs(v) @ op["H_abs"] @ np.abs(v)+2*np.abs(s) @ op["D_abs"] @ np.abs(v)+np.abs(s) @ op["G_abs"] @ np.abs(s))
    allowance=128*np.finfo(float).eps*absolute
    bound=float(op["compliance"]*np.sqrt(max(0,norm2)+allowance))
    no_cross_bound=float(op["compliance"]*np.sqrt(max(0,bulk+seam)+allowance))
    return 2*PB @ z,z,{"cpu_s":time.process_time()-t0,"gap_bound_m":bound,"norm2":norm2,"arithmetic_allowance_norm2":allowance,
        "residual_bulk_norm2":bulk,"signed_cross_norm2":cross,"seam_norm2":seam,
        "bound_if_cross_term_omitted_m":no_cross_bound,
        "accepted":bool(bound<=1e-7),"newton_iterations":it+1,"reduced_residual_relative":float(rel),
        "global_state_reconstructions":0,"global_solves":0,"coefficient_rebuilds":0}


def compile_operator(nl,lengths,out,tier):
    m=nl.m
    t0=time.process_time()
    acquisition=[]
    states=[]
    for length in lengths:
        cut=(m.x>=.008)&(m.x<.008+length*1e-3-1e-12)&(length>0)
        U,info=nl.full(cut,nl.uref)
        acquisition.append({"cut_length_mm":length,**info})
        states.append(U)
    basis0=time.process_time()
    P,sing,_=np.linalg.svd(np.stack(states,axis=1),full_matrices=False)
    d=int(np.sum(sing>sing[0]*1e-11))
    P=P[:,:d]
    combos=np.array(list(itertools.combinations_with_replacement(range(d+1),3)),dtype=int)
    # Gradient expansion including the known Dirichlet boundary in coordinate0.
    allP=np.vstack((P,np.zeros((m.nx+1,d))))
    E=np.einsum("eik,eia->eka",nl.grad,allP[nl.tris])
    fullconst=np.r_[np.zeros(m.n),nl.top]
    E0=np.einsum("eik,ei->ek",nl.grad,fullconst[nl.tris])
    E=np.concatenate((E0[:,:,None],E),axis=2)
    C=np.zeros((m.n,len(combos)))
    for col,triple in enumerate(combos):
        # Sum ordered tensor entries belonging to one commutative monomial.
        local=np.zeros((len(nl.tris),3))
        for i,j,k in set(itertools.permutations(tuple(triple))):
            term=np.sum(E[:,:,i]*E[:,:,j],axis=1)[:,None]*E[:,:,k]
            local+=nl.beta*nl.coeff[:,None]*np.einsum("eik,ek->ei",nl.grad,term)
        vals=np.zeros(len(m.coords))
        np.add.at(vals,nl.tris.ravel(),local.ravel())
        C[:,col]=vals[:m.n]
    cmap={tuple(t):i for i,t in enumerate(combos)}
    C[:,cmap[(0,0,0)]]-=nl.strain*m.f_unit
    AP=m.A @ P
    for i in range(d):
        C[:,cmap[(0,0,i+1)]]+=AP[:,i]
    coeff_cpu=time.process_time()-basis0
    gram0=time.process_time()
    invC=nl.lu0.solve(C)
    B=np.zeros((m.n,m.m));B[:m.m]=np.eye(m.m)
    invB=nl.lu0.solve(B)
    H=C.T @ invC
    H=(H+H.T)/2
    D=invC[:m.m]
    G=invB[:m.m]
    G=(G+G.T)/2
    RC=P.T @ C
    op={"PB":P[:m.m].copy(),"RC":RC,"H":H,"D":D,"G":G,
        "H_abs":np.abs(H),"D_abs":np.abs(D),"G_abs":np.abs(G),
        "combos":combos,"zref":P.T @ nl.uref,"delta0":np.array(m.mat.delta0),
        "k0":2*m.mat.K*m.mat.h*m.weights,"compliance":np.array(2*np.sqrt(np.max(np.diag(G)))),
        "reduced_scale":np.array(max(np.linalg.norm(P.T @ (nl.strain*m.f_unit)),1e-20)),
        "x":m.x,"kappa":nl.kappa}
    gram_cpu=time.process_time()-gram0
    np.savez(out/f"operator_tier{tier}.npz",**op)
    return op,P,C,{"tier":tier,"training":acquisition,"basis_dof":d,"residual_monomials":len(combos),"coefficient_compile_cpu_s":coeff_cpu,
        "gram_cpu_s":gram_cpu,"gram_global_rhs_solves":len(combos)+m.m,
        "total_acquisition_compile_cpu_s":time.process_time()-t0,"operator_bytes":sum(v.nbytes for v in op.values()),
        "smallest_retained_basis_singular_value":float(sing[d-1])}


def main():
    out=LANE/"r1/residual_compiler"
    out.mkdir(exist_ok=False)
    expected=(LANE/"PREREG_R1_RESIDUAL_COMPILER.sha256").read_text().split()[0]
    if hashlib.sha256((LANE/"PREREG_R1_RESIDUAL_COMPILER.json").read_bytes()).hexdigest()!=expected:
        raise RuntimeError("compiler prereg drift")
    all0=time.process_time()
    setup0=time.process_time()
    m=Membrane(64,32)
    nl=Nonlinear(m)
    initial_setup=time.process_time()-setup0
    tiers=[[0,4,12,20,28],[0]+list(range(1,32,2))]
    reports=[]
    for tier,lengths in enumerate(tiers,1):
        op,P,C,report=compile_operator(nl,lengths,out,tier)
        # Reopen saved artifact as the only data supplied to warm query.
        with np.load(out/f"operator_tier{tier}.npz") as f:
            portable={k:f[k] for k in f.files}
        rows=[]
        for length in (8,16,24,32):
            cut=(m.x>=.008)&(m.x<.008+length*1e-3-1e-12)
            times=[]
            for _ in range(25):
                gap,z,info=query(portable,cut,nl.kappa)
                times.append(info["cpu_s"])
            info["warm_median_cpu_s"]=float(np.median(times))
            info["cut_length_mm"]=length
            validation0=time.process_time()
            U=P @ z
            v,_=monomials(z,op["combos"])
            seam=op["k0"]*np.exp(-nl.kappa/op["delta0"])*(~cut)*U[:m.m]
            rpoly=C @ v
            rpoly[:m.m]+=seam
            rdirect,_=nl.bulk(U)
            rdirect[:m.m]+=seam
            info["polynomial_residual_max_absolute_error_N"]=float(np.max(np.abs(rpoly-rdirect)))
            info["polynomial_residual_relative_error"]=float(np.linalg.norm(rpoly-rdirect)/max(np.linalg.norm(rdirect),1e-20))
            exactdual=float(rdirect @ nl.lu0.solve(rdirect))
            info["direct_residual_dual_norm2"]=exactdual
            info["gram_contains_direct_norm2_with_allowance"]=bool(exactdual<=max(0,info["norm2"])+info["arithmetic_allowance_norm2"]+1e-20)
            truth,fullinfo=nl.full(cut,nl.uref)
            error=float(np.max(np.abs(gap-2*truth[:m.m])))
            info["actual_gap_error_m"]=error
            info["bound_contains_observed_error"]=bool(error<=info["gap_bound_m"]+1e-12)
            info["full_control"]=fullinfo
            info["validation_cpu_s"]=time.process_time()-validation0
            if not info["accepted"]:
                fall0=time.process_time()
                fallback,finfo=nl.full(cut,U)
                info["fallback_cpu_s"]=time.process_time()-fall0
                info["fallback"]=finfo
                info["returned_gap_error_m"]=float(np.max(np.abs(2*fallback[:m.m]-2*truth[:m.m])))
            else:
                info["fallback_cpu_s"]=0.
                info["fallback"]=None
                info["returned_gap_error_m"]=error
            info["charged_warm_query_including_fallback_cpu_s"]=info["warm_median_cpu_s"]+info["fallback_cpu_s"]
            info["speedup_vs_full_query_including_fallback"]=fullinfo["cpu_s"]/info["charged_warm_query_including_fallback_cpu_s"]
            rows.append(info)
            np.savez(out/f"validation_tier{tier}_l{length}.npz",candidate_gap=gap,reference_gap=2*truth[:m.m],z=z)
        report["heldouts"]=rows
        report["accepted_count"]=sum(r["accepted"] for r in rows)
        report["worst_gap_bound_m"]=max(r["gap_bound_m"] for r in rows)
        report["strongest_conventional_control"]="TIE: identical intrusive polynomial residual Gram / reduced-basis certifier, same acquisitions"
        write_json(out/f"tier{tier}.json",report)
        reports.append(report)
        print(__import__("json").dumps({"tier":tier,"basis":report["basis_dof"],"monomials":report["residual_monomials"],"acquisition_cpu_s":report["total_acquisition_compile_cpu_s"],"heldouts":[{k:r[k] for k in ("cut_length_mm","gap_bound_m","actual_gap_error_m","accepted","warm_median_cpu_s","polynomial_residual_relative_error")} for r in rows]}),flush=True)
        if report["accepted_count"]==4:
            break
    summary={"prereg_sha256":expected,"initial_mesh_and_nonlinear_prestress_setup_cpu_s":initial_setup,"tiers":reports,
        "total_cpu_s":time.process_time()-all0,"whole_evolving_history_certificate":"UNKNOWN",
        "conditional_discrete_certificate":"exact-arithmetic convex frozen-history residual theorem; binary64 assembly/evaluation allowance is engineering, pending independent numerical review",
        "review":"PENDING_INDEPENDENT_REVIEW","gate":"TIE",
        "next":"Enclose Gram arithmetic and propagate cohesive history branch tubes; acquire source-linked cleavage/strength ports and a vector layered nonlinear skin model"}
    write_json(out/"summary.json",summary)
    print(__import__("json").dumps({"out":str(out),"total_cpu_s":summary["total_cpu_s"]}),flush=True)


if __name__=="__main__":
    main()
