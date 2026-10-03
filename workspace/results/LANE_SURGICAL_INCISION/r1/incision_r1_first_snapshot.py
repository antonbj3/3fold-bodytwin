#!/usr/bin/env python3
"""Synthetic scalar P1 opening + irreversible cohesive interface experiment.

All writes stay in this lane. Parent HX modules are imported read-only. No
biological validation is implied. See PREREG_R1.json for numerical contracts.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
import os
import resource
import sys
import time
from dataclasses import dataclass
from pathlib import Path

for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[key] = "2"

import numpy as np
from scipy.linalg import cho_factor, cho_solve
from scipy.sparse import coo_matrix, diags
from scipy.sparse.linalg import splu

LANE = Path(__file__).resolve().parent
ROOT = LANE.parent.parent


def write_json(path, value):
    with Path(path).open("x") as f:
        json.dump(value, f, indent=2, allow_nan=False)
        f.write("\n")


def parent_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


@dataclass
class Material:
    E: float = 200000.0
    h: float = 0.002
    G: float = 200.0
    strength: float = 100000.0
    ax: float = 1.7

    @property
    def delta0(self):
        return self.G / (math.e * self.strength)

    @property
    def K(self):
        return self.G / self.delta0**2


class Membrane:
    """Upper half domain, scalar u; full opening delta=2u on y=0.

    Energy is twice the upper-half bulk energy plus cohesive interface energy.
    Scalar elliptic proxy deliberately omits vector incompressibility and shear.
    """

    def __init__(self, nx, ny, mat=None):
        t0, w0 = time.process_time(), time.perf_counter()
        self.mat = Material() if mat is None else mat
        self.nx, self.ny = nx, ny
        self.L, self.H = 0.064, 0.016
        self.dx, self.dy = self.L/nx, self.H/ny
        self.x = np.linspace(0, self.L, nx+1)
        self.coords = np.array([(x,y) for y in np.linspace(0,self.H,ny+1) for x in self.x])
        nall = len(self.coords)
        self.n = (nx+1)*ny
        self.m = nx+1
        self.b = np.arange(self.m)
        self.i = np.arange(self.m, self.n)
        self.weights = np.full(self.m, self.dx)
        self.weights[[0,-1]] *= 0.5
        triangles = []
        for j in range(ny):
            for k in range(nx):
                a = j*(nx+1)+k
                triangles.extend(((a,a+1,a+nx+2),(a,a+nx+2,a+nx+1)))
        triangles = np.array(triangles)
        xyz = self.coords[triangles]
        e1, e2 = xyz[:,1]-xyz[:,0], xyz[:,2]-xyz[:,0]
        area = (e1[:,0]*e2[:,1]-e1[:,1]*e2[:,0])/2
        grad = np.stack((np.stack((xyz[:,1,1]-xyz[:,2,1], xyz[:,2,0]-xyz[:,1,0]),axis=1),
                         np.stack((xyz[:,2,1]-xyz[:,0,1], xyz[:,0,0]-xyz[:,2,0]),axis=1),
                         np.stack((xyz[:,0,1]-xyz[:,1,1], xyz[:,1,0]-xyz[:,0,0]),axis=1)),axis=1)/(2*area[:,None,None])
        mid = xyz.mean(axis=1)
        modulus = self.mat.E*(1+0.22*np.sin(2*np.pi*mid[:,0]/self.L)*np.cos(np.pi*mid[:,1]/self.H))
        loc = self.mat.h*modulus[:,None,None]*area[:,None,None]*(self.mat.ax*grad[:,:,0,None]*grad[:,None,:,0]+grad[:,:,1,None]*grad[:,None,:,1])
        rr = np.repeat(triangles,3,axis=1).ravel()
        cc = np.tile(triangles,(1,3)).ravel()
        all_A = coo_matrix((loc.ravel(),(rr,cc)),shape=(nall,nall)).tocsc()
        self.A = all_A[:self.n,:self.n].tocsc()
        # Top displacement = strain * half-height; side boundaries natural.
        self.f_unit = -np.asarray(all_A[:self.n,self.n:] @ np.full(nx+1,self.H)).ravel()
        self.top_constant_unit = float(np.full(nx+1,self.H) @ (all_A[self.n:,self.n:] @ np.full(nx+1,self.H)))
        self.build_cpu = time.process_time()-t0
        self.build_wall = time.perf_counter()-w0
        t0, w0 = time.process_time(), time.perf_counter()
        Aii = self.A[self.i,:][:,self.i].tocsc()
        Aib = self.A[self.i,:][:,self.b].toarray()
        self.lu = splu(Aii)
        # Stored acquisition uses global fine solves ONCE and must be charged.
        R = -self.lu.solve(Aib)
        c = self.lu.solve(self.f_unit[self.i])
        self.S = self.A[self.b,:][:,self.b].toarray() + Aib.T @ R
        self.S = (self.S+self.S.T)/2
        self.fb = self.f_unit[self.b]-Aib.T @ c
        # Scalar energies without reconstructing a full state during queries.
        self.energy_constant_unit = self.top_constant_unit-float(self.f_unit[self.i] @ c)
        # Finite observation contract: seam gap and four near-seam rows.
        nrows = min(5,ny)
        self.strip_R = np.vstack((np.eye(self.m),R[:(nrows-1)*self.m]))
        self.strip_c = np.r_[np.zeros(self.m),c[:(nrows-1)*self.m]]
        self.nrows = nrows
        self.schur_cpu = time.process_time()-t0
        self.schur_wall = time.perf_counter()-w0
        self.acquisition_bytes = self.S.nbytes+self.fb.nbytes+self.strip_R.nbytes+self.strip_c.nbytes
        # R is not retained. Candidate queries cannot reconstruct exterior state.
        del R, c
        self.triangles, self.grad, self.area = triangles,grad,area

    def stiffness(self, kappa, cut):
        # delta=2u, so upper-half boundary force derivative is 2 K h dx.
        return 2*self.mat.K*self.mat.h*self.weights*np.exp(-kappa/self.mat.delta0)*(~cut)

    def equilibrium(self, kappa, cut, strain, mode="condensed", old_u=None):
        stiff = self.stiffness(kappa,cut)
        if mode == "full":
            diag = np.zeros(self.n)
            diag[:self.m] = stiff
            U = splu(self.A+diags(diag,format="csc")).solve(strain*self.f_unit)
            return U[:self.m], U, 0.0
        if mode in ("capsule","frozen_wake"):
            cut_x = self.x[cut]
            tip = cut_x.max() if len(cut_x) else -1.0
            wake = np.flatnonzero(cut & (self.x < tip-3*self.dx))
        else:
            wake = np.array([],dtype=int)
        active = np.setdiff1d(self.b,wake)
        if len(wake):
            t0 = time.process_time()
            SAW = self.S[np.ix_(active,wake)]
            if mode == "capsule":
                factor = cho_factor(self.S[np.ix_(wake,wake)],lower=True)
                response = cho_solve(factor,SAW.T)
                load = cho_solve(factor,strain*self.fb[wake])
                Q = self.S[np.ix_(active,active)]-SAW @ response
                rhs = strain*self.fb[active]-SAW @ load
            else:
                Q = self.S[np.ix_(active,active)].copy()
                rhs = strain*self.fb[active]-SAW @ old_u[wake]
            rebuild = time.process_time()-t0
            Q[np.diag_indices_from(Q)] += stiff[active]
            ua = cho_solve(cho_factor(Q,lower=True),rhs)
            u = np.zeros(self.m)
            u[active] = ua
            u[wake] = load-response @ ua if mode == "capsule" else old_u[wake]
        else:
            Q = self.S+np.diag(stiff)
            u = cho_solve(cho_factor(Q,lower=True),strain*self.fb)
            rebuild = 0.0
        return u, None, rebuild

    def advance(self, hist, cut, strain, mode):
        kappa = hist["kappa"].copy()
        old_u = hist["u"].copy()
        rebuild, solves = 0.0, 0
        for it in range(400):
            u,U,rc = self.equilibrium(kappa,cut,strain,mode,old_u)
            rebuild += rc
            solves += 1
            next_k = np.maximum(kappa,np.maximum(2*u,0))
            # Exact completed surface, independent of strain unloading.
            next_k[cut] = kappa[cut]
            diff = float(np.max(np.abs(next_k-kappa)))
            kappa = next_k
            if diff < 1e-12:
                break
        else:
            raise RuntimeError(f"cohesive iteration did not converge: {mode}, {diff}")
        # Recompute once with converged history; explicit cost charged.
        u,U,rc = self.equilibrium(kappa,cut,strain,mode,old_u)
        rebuild += rc
        solves += 1
        stiff = self.stiffness(kappa,cut)
        residual = self.S @ u+stiff*u-strain*self.fb
        residual_norm = float(np.linalg.norm(residual)/max(np.linalg.norm(strain*self.fb),1e-20))
        damage = 1-np.exp(-kappa/self.mat.delta0)
        damage[cut] = 1
        x = kappa/self.mat.delta0
        diss = self.mat.G*(1-(1+x+0.5*x*x)*np.exp(-x))
        diss[cut] = self.mat.G
        cohesive_free = np.sum(0.5*self.mat.K*np.exp(-x)*(~cut)*(2*u)**2*self.mat.h*self.weights)
        # Full specimen bulk energy = twice upper-half energy, including BC terms.
        bulk = float(u @ self.S @ u-2*strain*self.fb @ u+strain**2*self.energy_constant_unit)
        diss_J = float(self.mat.h*self.weights @ diss)
        energy = bulk+float(cohesive_free)+diss_J
        if U is not None:
            global_bulk = float(U @ (self.A @ U)-2*strain*self.f_unit @ U+strain**2*self.top_constant_unit)
            energy_check = abs(global_bulk-bulk)
        else:
            energy_check = None
        strip = (self.strip_R @ u+strain*self.strip_c).reshape(self.nrows,self.m)
        normal_strain = np.abs(np.diff(strip,axis=0)/self.dy)
        # Damage-zone width is a declared strain-threshold proxy, not cell death.
        above = normal_strain > 0.15
        widths = np.sum(above,axis=0)*self.dy
        cut_indices = np.flatnonzero(cut)
        tip_idx = cut_indices.max() if len(cut_indices) else 0
        near_tip = np.abs(self.b-tip_idx) <= 3
        width = float(np.max(widths[near_tip]))
        uncut_damaged_length = float(np.sum(self.weights[(~cut) & (damage>=0.10)]))
        return {"kappa":kappa,"u":u,"damage":damage}, {
            "gap_max_m":float(2*np.max(u[cut])) if np.any(cut) else 0.0,
            "exposed_area_m2":float(self.weights @ (2*u*cut)),
            "mechanical_strain_zone_halfwidth_m":width,
            "zone_censored_at_patch_edge":bool(np.any(above[-1,near_tip])),
            "cohesive_damaged_uncut_length_m":uncut_damaged_length,
            "bulk_elastic_energy_J":bulk,"cohesive_free_energy_J":float(cohesive_free),
            "cohesive_dissipation_J":diss_J,"total_energy_J":energy,
            "full_energy_pullback_check_J":energy_check,
            "residual_relative":residual_norm,"linear_solves":solves,
            "capsule_rebuild_cpu_s":rebuild,
            "history_monotone":bool(np.all(kappa>=hist["kappa"]-1e-12)),
            "damage_min":float(damage.min()),"damage_max":float(damage.max())
        }


def history_plan():
    plan = [(0.0,e) for e in (0.03,0.10,0.03,0.10)]
    plan += [(float(a),0.10) for a in np.linspace(0.002,0.032,16)]
    # Additional deformation of a cut with a completed wake exposes memory loss.
    plan += [(0.032,0.03),(0.032,0.10)]
    return plan


def run_history(model,mode):
    hist = {"kappa":np.zeros(model.m),"u":np.zeros(model.m)}
    rows, gaps, damages = [],[],[]
    t0,w0 = time.process_time(),time.perf_counter()
    for length,strain in history_plan():
        cut = (model.x>=0.008-1e-12) & (model.x<0.008+length-1e-12) & (length>0)
        prevE = rows[-1]["total_energy_J"] if rows else None
        prevL = rows[-1]["length_m"] if rows else 0
        prevStrain = rows[-1]["strain"] if rows else None
        hist,metrics = model.advance(hist,cut,strain,mode)
        force = (metrics["total_energy_J"]-prevE)/(length-prevL) if prevE is not None and length>prevL and strain==prevStrain else None
        metrics.update(length_m=length,strain=strain,blade_advance_energy_force_N=force)
        rows.append(metrics)
        gaps.append(2*hist["u"].copy())
        damages.append(hist["damage"].copy())
    return {"mode":mode,"cpu_s":time.process_time()-t0,"wall_s":time.perf_counter()-w0,
            "global_state_reconstructions":sum(r["linear_solves"] for r in rows) if mode=="full" else 0,
            "global_query_factorizations":sum(r["linear_solves"] for r in rows) if mode=="full" else 0,
            "capsule_rebuild_cpu_s":sum(r["capsule_rebuild_cpu_s"] for r in rows),"rows":rows}, np.array(gaps),np.array(damages)


def compare(ref, candidate, rg,cg,rd,cd):
    energy_scale = max(max(abs(r["total_energy_J"]) for r in ref["rows"]),1e-20)
    metrics = {
        "max_gap_error_m":float(np.max(np.abs(rg-cg))),
        "max_damage_error":float(np.max(np.abs(rd-cd))),
        "energy_relative_error":float(max(abs(r["total_energy_J"]-c["total_energy_J"]) for r,c in zip(ref["rows"],candidate["rows"]))/energy_scale),
        "max_residual_relative":float(max(r["residual_relative"] for r in candidate["rows"])),
        "max_strain_zone_width_error_m":float(max(abs(r["mechanical_strain_zone_halfwidth_m"]-c["mechanical_strain_zone_halfwidth_m"]) for r,c in zip(ref["rows"],candidate["rows"]))),
        "all_history_monotone":all(r["history_monotone"] for r in candidate["rows"])
    }
    metrics["gate"] = "PASS" if metrics["max_gap_error_m"]<=1e-7 and metrics["max_damage_error"]<=1e-6 and metrics["energy_relative_error"]<=1e-5 and metrics["max_residual_relative"]<=1e-8 and metrics["all_history_monotone"] else "FAIL"
    return metrics


def world_test():
    prereg = json.loads((LANE/"PREREG_R1.json").read_text())
    obs = prereg["world_gates"]["goda_2024_PAAm"]["observations"]
    cleavage,tear = obs["cutting_intercept_J_m2"],obs["tear_toughness_J_m2"]
    predicted_slope = 4*obs["uniaxial_failure_energy_density_J_m3"]
    slope = obs["cutting_force_per_thickness_radius_slope_Pa"]
    # These are not individual measured force points: evaluations of published fits.
    radii = np.array([0.125,0.210,0.275,0.390,0.470])*1e-3
    return {
        "single_energy_transfer": {"prediction_J_m2":tear,"observed_intercept_J_m2":cleavage,
            "relative_error":abs(tear-cleavage)/cleavage,"gate":"PASS" if abs(tear-cleavage)/cleavage<=0.20 else "FAIL"},
        "changed_operation":"Two constitutive ports: cleavage surface energy Gamma0 [J/m2] and pristine failure work Wf [J/m3]",
        "separate_assay_slope_prediction": {"prediction_Pa":predicted_slope,"observed_Pa":slope,
            "relative_error":abs(predicted_slope-slope)/slope,"gate":"PASS" if abs(predicted_slope-slope)/slope<=0.20 else "FAIL"},
        "evaluations_of_published_fit_NOT_measurements": {"radius_m":radii.tolist(),"published_fit_force_per_thickness_N_m":(cleavage+slope*radii).tolist(),"two_port_prediction_N_m":(cleavage+predicted_slope*radii).tolist()},
        "unmeasured_FPZ_length_if_energy_closure_Wf_times_length":(tear-cleavage)/obs["uniaxial_failure_energy_density_J_m3"],
        "skin_cross_prediction":"UNKNOWN: no matched skin cleavage/strength/gap/cell-injury dataset",
        "novelty":"Published scaling relation, reused as a discriminating mechanism; no claim of discovering it",
        "review":"PENDING_INDEPENDENT_REVIEW"
    }


def identifiable_ports():
    mat = Material()
    strengths = np.array([50000.,100000.,200000.])
    # The LEFM dimensional process-zone scale is distinct from our scalar model's
    # strain threshold. No equality to measured biological width is asserted.
    return {"G_fixed_J_m2":mat.G,"E_fixed_Pa":mat.E,"strength_Pa":strengths.tolist(),
        "cohesive_delta0_m":(mat.G/(math.e*strengths)).tolist(),
        "dimensional_zone_scale_EG_over_strength_squared_m":(mat.E*mat.G/strengths**2).tolist(),
        "largest_to_smallest_scale_ratio":float((strengths.max()/strengths.min())**2),
        "constant_toughness_does_not_identify_strength_or_zone":True,
        "force_radius_observation_Jacobian_columns_G0_Wf_strength":[[1.0,4*r,0.0] for r in np.array([0.125,0.210,0.275,0.390,0.470])*1e-3],
        "rank":2,"null_direction":[0,0,1],
        "required_new_acquisition":"Cohesive tensile peak/strain map plus perfused vascular geometry; cutting radii alone leave strength unobserved"}


def hx_consumer(area, injury_seed):
    """Reuse actual HX code, with explicit synthetic new coupling ports.

    No platelet/coagulation import: those source scripts mutate private folders
    at import time. Darcy exudate is never relabeled as vascular bleeding.
    """
    q33 = parent_module("q33_readonly",ROOT/"results/BT-HX-Q033/model.py")
    q36 = parent_module("q36_readonly",ROOT/"results/BT-HX-Q036/model.py")
    t0 = time.process_time()
    p33 = q33.DEFAULT_PARAMETERS
    trans = q33.transport_flow(p33)
    flux = (trans.cut_total_flow_m3_s/q33.ellipse_geometry(p33)["exposed_area_m2"])
    p = q36.Parameters()
    from scipy.integrate import solve_ivp
    outputs = []
    # Coupled memory: mechanical injury sets M and q_ext; Q036 then retains
    # oxygen/ATP/Ca/ROS/edema feedback. New port constants are synthetic closures.
    for seed in (0.0,injury_seed):
        sc = q36.Scenario(name="incision_edge_synthetic",pre_ischemia_min=0,ischemia_min=10,reperfusion_min=10,q_ischemia=max(0.01,1-seed))
        state = q36.healthy_state(p)
        state[8] = seed*0.2
        state[9] = seed*0.1
        states=[]
        for lo,hi in ((0,10),(10,20)):
            # Explicit split at the discontinuous external-flow port.
            solve = solve_ivp(lambda t,y:q36.rhs(t,y,p,sc),(lo,hi),state,method="DOP853",rtol=1e-8,atol=1e-10,max_step=0.25)
            if not solve.success:
                raise RuntimeError(solve.message)
            state=solve.y[:,-1]
            states.append({"time_min":hi,"state":q36.unpack(state)})
        outputs.append({"mechanical_damage_seed":seed,"q_ext_during_first_10min":sc.q_ischemia,"snapshots":states})
    return {"Q033_exudate_m3_s":flux*area,"Q033_area_input_m2":area,"Q033_evaporation_multiplier_analog":p33.barrier_damage_multiplier,
        "Darcy_flow_is_not_bleeding":True,"vascular_bleeding_m3_s":None,"hemostasis_time_min":None,
        "Q036_cases":outputs,"coupling":"one-way mechanical initial state -> nonlinear Q036 internal no-reflow/edema feedback; no mechanics feedback claimed",
        "biological_status":"SYNTHETIC_CONSTITUTIVE_CLOSURES", "cpu_s":time.process_time()-t0}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--tag",default="base")
    parser.add_argument("--sizes",default="32,64,96")
    args=parser.parse_args()
    out=LANE/"r1"/args.tag
    out.mkdir(exist_ok=False)
    pre=(LANE/"PREREG_R1.json").read_bytes()
    expected=(LANE/"PREREG_R1.sha256").read_text().split()[0]
    if hashlib.sha256(pre).hexdigest()!=expected:
        raise RuntimeError("prereg hash drift")
    total0=time.process_time()
    all_rows=[]
    for nx in map(int,args.sizes.split(",")):
        model=Membrane(nx,nx//2)
        payload={"nx":nx,"ny":nx//2,"global_free_dofs":model.n,"seam_dofs":model.m,
            "build_cpu_s":model.build_cpu,"build_wall_s":model.build_wall,
            "schur_and_observation_acquisition_cpu_s":model.schur_cpu,"schur_and_observation_acquisition_wall_s":model.schur_wall,
            "operator_bytes":model.acquisition_bytes,"methods":{}}
        ref,rg,rd=run_history(model,"full")
        payload["methods"]["full"]=ref
        np.savez(out/f"raw_n{nx}_full.npz",gap=rg,damage=rd,x=model.x)
        for mode in ("condensed","frozen_wake","capsule"):
            result,cg,cd=run_history(model,mode)
            result["comparison"]=compare(ref,result,rg,cg,rd,cd)
            result["speedup_vs_repeated_full_query_only"]=ref["cpu_s"]/result["cpu_s"]
            setup=model.schur_cpu
            savings=ref["cpu_s"]-result["cpu_s"]
            result["finite_break_even_complete_histories_vs_repeated_full"]=setup/savings if savings>0 and result["comparison"]["gate"]=="PASS" else None
            np.savez(out/f"raw_n{nx}_{mode}.npz",gap=cg,damage=cd,x=model.x)
            payload["methods"][mode]=result
        payload["strongest_conventional_control"]={"algorithm":"Identical fine Schur substructuring / completed-wake elimination","speedup":1.0,"gate":"TIE","reason":"Same acquired operator, nonlinear law, histories and output pullbacks"}
        payload["query_without_setup_counts"]={"condensed_global_solves":0,"capsule_global_solves":0,"history_length":len(history_plan()),"local_front_not_constant_size":"All remaining nonlinear seam DOFs retained; completed wake matrix grows. No fixed-cost arbitrary-cut claim."}
        write_json(out/f"mechanical_n{nx}.json",payload)
        np.savez(out/f"operator_n{nx}.npz",S=model.S,fb=model.fb,strip_R=model.strip_R,strip_c=model.strip_c,x=model.x,energy_constant_unit=model.energy_constant_unit)
        all_rows.append(payload)
        print(json.dumps({"nx":nx,"full_cpu_s":ref["cpu_s"],"methods":{m:{"cpu_s":v["cpu_s"],"comparison":v.get("comparison")} for m,v in payload["methods"].items() if m!="full"}}),flush=True)
    largest=all_rows[-1]["methods"]["condensed"]["rows"][-1]
    world=world_test()
    ident=identifiable_ports()
    # Scalar mechanical damage at a supplied scale is a synthetic initial-state
    # port. It is not an independently measured cellular-injury fraction.
    hx=hx_consumer(largest["exposed_area_m2"],0.5)
    write_json(out/"world_transfer.json",world)
    write_json(out/"identifiability.json",ident)
    write_json(out/"HX_chain.json",hx)
    convergence=[]
    for coarse,fine in zip(all_rows[:-1],all_rows[1:]):
        c=coarse["methods"]["condensed"]["rows"][-1]
        f=fine["methods"]["condensed"]["rows"][-1]
        rel=abs(c["gap_max_m"]-f["gap_max_m"])/max(f["gap_max_m"],1e-20)
        convergence.append({"meshes":[coarse["nx"],fine["nx"]],"gap_relative_difference":rel,"gate":"PASS" if rel<=0.05 else "FAIL","not_continuum_error_certificate":True})
    summary={"prereg_sha256":expected,"mechanical_sizes":all_rows,"mesh_convergence":convergence,
        "world_transfer":world,"identifiability":ident,"HX_chain":hx,"total_cpu_s":time.process_time()-total0,
        "max_rss_MiB":resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,
        "gate":"TIE","empirical_skin_validation":"UNKNOWN","pending_independent_review":True}
    write_json(out/"summary.json",summary)
    print(json.dumps({"output":str(out),"gate":summary["gate"],"total_cpu_s":summary["total_cpu_s"],"max_rss_MiB":summary["max_rss_MiB"]}),flush=True)


if __name__=="__main__":
    main()
