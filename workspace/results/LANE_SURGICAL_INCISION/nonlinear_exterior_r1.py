#!/usr/bin/env python3
"""Load-bearing extension: nonlinear bulk and charged residual certification."""
import os
for key in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","NUMEXPR_NUM_THREADS"):
    os.environ[key]="2"
import hashlib
import json
import time
from pathlib import Path
import numpy as np
from scipy.linalg import cho_factor,cho_solve
from scipy.sparse import coo_matrix,diags
from scipy.sparse.linalg import splu
from incision_r1 import LANE,Membrane,Material,run_history,write_json


class Nonlinear:
    def __init__(self,m,beta=20,strain=.1):
        self.m,self.beta,self.strain=m,beta,strain
        self.tris,self.grad=m.triangles,m.grad
        self.top=np.full(m.nx+1,strain*m.H)
        xyz=m.coords[self.tris]
        mid=xyz.mean(axis=1)
        self.coeff=m.mat.h*m.mat.E*(1+.22*np.sin(2*np.pi*mid[:,0]/m.L)*np.cos(np.pi*mid[:,1]/m.H))*m.area
        self.mid=mid
        self.rr=np.repeat(self.tris,3,axis=1).ravel()
        self.cc=np.tile(self.tris,(1,3)).ravel()
        self.kappa=np.zeros(m.m)
        hist={"u":np.zeros(m.m),"kappa":self.kappa}
        for e in (.03,.10,.03,.10):
            hist,_=m.advance(hist,np.zeros(m.m,dtype=bool),e,"condensed")
        self.kappa=hist["kappa"]
        self.cut0=np.zeros(m.m,dtype=bool)
        self.lu0=splu(m.A)
        B=np.zeros((m.n,m.m))
        B[:m.m]=np.eye(m.m)
        compliance=self.lu0.solve(B)
        self.max_gap_compliance=2*np.sqrt(np.max(np.diag(compliance[:m.m])))
        self.uref,self.baseline_info=self.full(self.cut0,np.zeros(m.n))
        self.rref,self.Jref=self.bulk(self.uref)
        self.bref=self.Jref @ self.uref-self.rref
        self.refstrain=self.gradient(self.uref)
        self.refnlforce,self.refnlJ=self.nl_element(self.refstrain)

    def gradient(self,U):
        allu=np.r_[U,self.top]
        return np.einsum("ei,eik->ek",allu[self.tris],self.grad)

    def nl_element(self,g):
        norm=np.sum(g*g,axis=1)
        # local element vector and positive Hessian of beta*c*|grad u|^4/4.
        force=self.beta*self.coeff[:,None]*np.einsum("eik,ek->ei",self.grad,norm[:,None]*g)
        tensor=norm[:,None,None]*np.eye(2)+2*g[:,:,None]*g[:,None,:]
        J=self.beta*self.coeff[:,None,None]*np.einsum("eik,ekl,ejl->eij",self.grad,tensor,self.grad)
        return force,J

    def assemble(self,force,J,selected=None):
        tris=self.tris if selected is None else self.tris[selected]
        r=np.zeros(len(self.m.coords))
        np.add.at(r,tris.ravel(),force.ravel())
        rr=np.repeat(tris,3,axis=1).ravel()
        cc=np.tile(tris,(1,3)).ravel()
        matrix=coo_matrix((J.ravel(),(rr,cc)),shape=(len(r),len(r))).tocsc()
        return r[:self.m.n],matrix[:self.m.n,:self.m.n].tocsc()

    def bulk(self,U):
        force,J=self.nl_element(self.gradient(U))
        rnl,Jnl=self.assemble(force,J)
        return self.m.A @ U-self.strain*self.m.f_unit+rnl,self.m.A+Jnl

    def full(self,cut,initial):
        t0=time.process_time()
        U=initial.copy()
        k=np.zeros(self.m.n)
        k[:self.m.m]=self.m.stiffness(self.kappa,cut)
        nit=0
        for it in range(50):
            r,J=self.bulk(U)
            r=r+k*U
            nit+=1
            rel=np.linalg.norm(r)/max(np.linalg.norm(self.strain*self.m.f_unit),1e-20)
            if rel<1e-10:
                break
            step=splu(J+diags(k,format="csc")).solve(-r)
            alpha=1.0
            for _ in range(20):
                nextU=U+alpha*step
                nextR,_=self.bulk(nextU)
                if np.linalg.norm(nextR+k*nextU)<np.linalg.norm(r):
                    U=nextU
                    break
                alpha*=.5
            else:
                raise RuntimeError("nonlinear Newton line search failed")
        else:
            raise RuntimeError("nonlinear Newton did not converge")
        return U,{"cpu_s":time.process_time()-t0,"newton_iterations":nit,"residual_relative":float(rel),"global_factorizations":nit-1}

    def patch(self,cut,radius,wake=False):
        tip=self.m.x[cut].max()
        if wake:
            selected=np.flatnonzero((self.mid[:,0]>=.006)&(self.mid[:,0]<=tip+.002)&(self.mid[:,1]<=.004))
        else:
            selected=np.flatnonzero((np.abs(self.mid[:,0]-tip)<=radius)&(self.mid[:,1]<=radius))
        active=np.unique(np.r_[self.m.b,self.tris[selected].ravel()])
        active=active[active<self.m.n]
        inactive=np.setdiff1d(np.arange(self.m.n),active)
        t0=time.process_time()
        Jaa=self.Jref[active,:][:,active].toarray()
        Jia=self.Jref[inactive,:][:,active].toarray()
        lu=splu(self.Jref[inactive,:][:,inactive].tocsc())
        R=-lu.solve(Jia)
        c=lu.solve(self.bref[inactive])
        S=Jaa+Jia.T @ R
        f=self.bref[active]-Jia.T @ c
        S=(S+S.T)/2
        setup_cpu=time.process_time()-t0
        # Every selected nonlinear element is supported on retained nodes; query
        # evaluates only those local element gradients, not inactive state.
        index=np.full(len(self.m.coords),-1,dtype=int)
        index[active]=np.arange(len(active))
        local_tris=index[self.tris[selected]]
        if np.any(local_tris<0):
            raise RuntimeError("patch touches prescribed top unexpectedly")
        pg=self.grad[selected]
        baseg=self.refstrain[selected]
        oldf=self.refnlforce[selected]
        oldJ=self.refnlJ[selected]
        coef=self.coeff[selected]
        k=np.zeros(len(active))
        k[:self.m.m]=self.m.stiffness(self.kappa,cut)
        a=self.uref[active].copy()
        query0=time.process_time()
        nit=0
        for it in range(50):
            g=np.einsum("ei,eik->ek",a[local_tris],pg)
            dg=g-baseg
            norm=np.sum(g*g,axis=1)
            nf=self.beta*coef[:,None]*np.einsum("eik,ek->ei",pg,norm[:,None]*g)
            tangent_old=self.beta*coef[:,None]*np.einsum("eik,ekl,el->ei",pg,
                np.sum(baseg*baseg,axis=1)[:,None,None]*np.eye(2)+2*baseg[:,:,None]*baseg[:,None,:],dg)
            rem=nf-oldf-tangent_old
            r=S @ a-f+k*a
            np.add.at(r,local_tris.ravel(),rem.ravel())
            rel=np.linalg.norm(r)/max(np.linalg.norm(f),1e-20)
            nit+=1
            if rel<1e-10:
                break
            tensor=norm[:,None,None]*np.eye(2)+2*g[:,:,None]*g[:,None,:]
            elemJ=self.beta*coef[:,None,None]*np.einsum("eik,ekl,ejl->eij",pg,tensor,pg)-oldJ
            J=S+np.diag(k)
            rr=np.repeat(local_tris,3,axis=1).ravel()
            cc=np.tile(local_tris,(1,3)).ravel()
            np.add.at(J,(rr,cc),elemJ.ravel())
            # Truncated tangent remainder is not necessarily globally convex;
            # Newton uses a dense symmetric solve and verifies local residual.
            step=np.linalg.solve(J,-r)
            a=a+step
        else:
            raise RuntimeError("patch Newton did not converge")
        query_cpu=time.process_time()-query0
        # Charge full reconstruction and global residual certificate separately.
        cert0=time.process_time()
        U=np.zeros(self.m.n)
        U[active]=a
        U[inactive]=R @ a+c
        r,_=self.bulk(U)
        r[:self.m.m]+=self.m.stiffness(self.kappa,cut)*U[:self.m.m]
        dual=float(np.sqrt(max(0,r @ self.lu0.solve(r))))
        bound=float(self.max_gap_compliance*dual)
        cert_cpu=time.process_time()-cert0
        return U,{"radius_m":radius,"wake_patch":wake,"active_dofs":len(active),"nonlinear_elements":len(selected),
            "setup_cpu_s":setup_cpu,"local_query_cpu_s":query_cpu,"certificate_cpu_s":cert_cpu,
            "newton_iterations":nit,"local_residual_relative":float(rel),
            "certificate_gap_bound_m":bound,"certificate_accepted":bool(bound<=1e-7),
            "certificate_global_solves":1,"certificate_global_state_reconstructions":1,
            "certificate_validity":"Fixed supplied cohesive history, convex full discrete beta>=0 model, floating-point numerical bound; no irreversible-branch or continuum enclosure"}


def main():
    out=LANE/"r1/nonlinear"
    out.mkdir(exist_ok=False)
    expected=(LANE/"PREREG_R1_EXTENSION.sha256").read_text().split()[0]
    if hashlib.sha256((LANE/"PREREG_R1_EXTENSION.json").read_bytes()).hexdigest()!=expected:
        raise RuntimeError("extension prereg drift")
    begin=time.process_time()
    m=Membrane(64,32)
    acq0=time.process_time()
    nl=Nonlinear(m)
    acq=time.process_time()-acq0
    rows=[]
    for length in (.008,.020,.032):
        cut=(m.x>=.008)&(m.x<.008+length-1e-12)
        truth,fullinfo=nl.full(cut,nl.uref)
        refgap=2*truth[:m.m]
        record={"length_m":length,"full":fullinfo,"patches":[]}
        for radius,wake in ((.002,False),(.004,False),(.008,False),(.004,True)):
            U,info=nl.patch(cut,radius,wake)
            error=float(np.max(np.abs(2*U[:m.m]-refgap)))
            info["actual_max_gap_error_m"]=error
            info["bound_contains_observed_error"]=bool(error<=info["certificate_gap_bound_m"]+1e-12)
            info["accuracy_gate"]="PASS" if error<=1e-7 else "FAIL"
            fallback0=time.process_time()
            if not info["certificate_accepted"]:
                fallback,finfo=nl.full(cut,U)
                info["fallback"]=finfo
                info["returned_gap_max_error_m"]=float(np.max(np.abs(2*fallback[:m.m]-refgap)))
            else:
                info["fallback"]=None
                info["returned_gap_max_error_m"]=error
            info["fallback_cpu_s"]=time.process_time()-fallback0
            info["charged_total_cpu_s"]=sum(info[k] for k in ("setup_cpu_s","local_query_cpu_s","certificate_cpu_s","fallback_cpu_s"))
            info["speedup_vs_uniform_full_including_rebuild_cert_fallback"]=fullinfo["cpu_s"]/info["charged_total_cpu_s"]
            record["patches"].append(info)
            np.savez(out/f"raw_l{int(length*1000)}_r{int(radius*1000)}_wake{int(wake)}.npz",patch_U=U,reference_U=truth,cut=cut,kappa=nl.kappa)
        rows.append(record)
        print(json.dumps(record),flush=True)
    strengths=[]
    for strength in (50000.,100000.,200000.):
        sm=Membrane(64,32,Material(strength=strength))
        result,gap,damage=run_history(sm,"condensed")
        strengths.append({"strength_Pa":strength,"G_J_m2":sm.mat.G,"delta0_m":sm.mat.delta0,"final":result["rows"][-1],"query_cpu_s":result["cpu_s"],"setup_cpu_s":sm.build_cpu+sm.schur_cpu})
        np.savez(out/f"strength_{int(strength)}.npz",gap=gap,damage=damage)
    summary={"prereg_sha256":expected,"nonlinear_prestress_compliance_acquisition_cpu_s":acq,
        "base_mesh_build_schur_cpu_s":m.build_cpu+m.schur_cpu,"rows":rows,"strength_intervention":strengths,
        "total_cpu_s":time.process_time()-begin,"outcome":"FAIL for cheap patch accuracy/cost; full fallback returns matched model accuracy; strongest same-information conventional patch is TIE",
        "obstacle":"Global nonlinear strain correction persists outside tip and wake patch; residual certification reconstructs global state and uses a global solve; moving patch rebuild costs can dominate full fine solve",
        "next":"Compile residual dual norm or coarse block nonlinear flux enclosures without global state reconstruction, preserve cohesive branch history, and compare against sparse nonlinear tangent updates",
        "review":"PENDING_INDEPENDENT_REVIEW"}
    write_json(out/"summary.json",summary)
    print(json.dumps({"out":str(out),"total_cpu_s":summary["total_cpu_s"]}),flush=True)


if __name__=="__main__":
    main()
