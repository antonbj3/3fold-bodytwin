#!/usr/bin/env python3
"""Three coupled finite-strain membranes: scenario, not calibrated skin.

Reuse R1's triangular half-domain and energy-release convention. Replace the
scalar bulk by two displacement components, NH matrix and recruited fibers;
couple layers through slip energy. Prescribed incision, no inferred failure law.
"""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[key] = '2'
import argparse
import json
import time
from pathlib import Path
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve

LANE = Path(__file__).resolve().parent

def write_json(path, value):
    with Path(path).open('x') as f:
        json.dump(value, f, indent=2, allow_nan=False)
        f.write('\n')

class Layers:
    def __init__(self, nx=32, ny=16, angle=0., coupling=2e7):
        self.nx, self.ny = nx, ny
        self.L, self.H = .064, .016
        self.x = np.linspace(-self.L/2, self.L/2, nx+1)
        self.xy = np.array([(x,y) for y in np.linspace(0,self.H,ny+1) for x in self.x])
        self.n = len(self.xy)
        self.h = np.array([.0001,.0015,.003])
        E = np.array([1e5,2e5,2e4])
        nu = .3
        self.mu = E/(2*(1+nu))
        self.lame = E*nu/((1+nu)*(1-2*nu))
        self.kfiber = np.array([2e4,1e5,0.])
        theta = np.deg2rad(angle)
        self.fiber = np.array([np.cos(theta),np.sin(theta)])
        self.coupling = coupling
        tri = []
        for j in range(ny):
            for i in range(nx):
                a = j*(nx+1)+i
                tri.extend(((a,a+1,a+nx+2),(a,a+nx+2,a+nx+1)))
        self.tri = np.array(tri)
        xyz = self.xy[self.tri]
        e1, e2 = xyz[:,1]-xyz[:,0],xyz[:,2]-xyz[:,0]
        self.area = (e1[:,0]*e2[:,1]-e1[:,1]*e2[:,0])/2
        self.grad = np.stack((np.stack((xyz[:,1,1]-xyz[:,2,1],xyz[:,2,0]-xyz[:,1,0]),axis=1),
                              np.stack((xyz[:,2,1]-xyz[:,0,1],xyz[:,0,0]-xyz[:,2,0]),axis=1),
                              np.stack((xyz[:,0,1]-xyz[:,1,1],xyz[:,1,0]-xyz[:,0,0]),axis=1)),axis=1)/(2*self.area[:,None,None])
        self.mass = np.zeros(self.n)
        np.add.at(self.mass,self.tri.ravel(),np.repeat(self.area/3,3))
        self.bottom = np.arange(nx+1)
        self.outer = np.flatnonzero((np.abs(self.xy[:,0])>=self.L/2-1e-12)|(self.xy[:,1]>=self.H-1e-12))
        self.seam_w = np.full(nx+1,self.L/nx)
        self.seam_w[[0,-1]] *= .5
        self.dofs = (2*self.tri[:,:,None]+np.arange(2)).reshape(-1,6)
        self.rr = np.repeat(self.dofs,6,axis=1).ravel()
        self.cc = np.tile(self.dofs,(1,6)).ravel()
        # Full specimen energy: two mirror half-domains, including layer slip.
        d = np.arange(2*self.n)
        self.idd = d
        self.slip_diag = np.repeat(2*coupling*self.mass,2)

    def energy(self, U, hessian=False):
        U = U.reshape(3,self.n,2)
        r = np.zeros_like(U)
        parts = []
        rr,cc,dd = [],[],[]
        eye = np.eye(2)
        minJ = 1e30
        for layer in range(3):
            F = eye+np.einsum('eai,eaj->eij',U[layer,self.tri],self.grad)
            J = np.linalg.det(F)
            minJ = min(minJ,float(J.min()))
            if np.any(J<=0):
                return np.inf, r.ravel(), None, {'minJ':minJ}
            invt = np.linalg.inv(F).transpose(0,2,1)
            logJ = np.log(J)
            mu,la,kf = self.mu[layer],self.lame[layer],self.kfiber[layer]
            fa = F@self.fiber
            q = np.maximum(np.sum(fa*fa,axis=1)-1,0)
            W = mu/2*(np.sum(F*F,axis=(1,2))-2-2*logJ)+la/2*logJ**2+kf/2*q**2
            P = mu*(F-invt)+la*logJ[:,None,None]*invt+2*kf*q[:,None,None]*fa[:,:,None]*self.fiber[None,None,:]
            fac = 2*self.h[layer]*self.area
            parts.append(float(fac@W))
            local = fac[:,None,None]*np.einsum('eij,eaj->eai',P,self.grad)
            np.add.at(r[layer],self.tri.ravel(),local.reshape(-1,2))
            if hessian:
                A = mu*np.einsum('ik,jl->ijkl',eye,eye)[None,:,:,:,:]
                A = A+(mu-la*logJ)[:,None,None,None,None]*np.einsum('eil,ekj->eijkl',invt,invt)
                A += la*np.einsum('eij,ekl->eijkl',invt,invt)
                active = (q>0).astype(float)
                A += 4*kf*active[:,None,None,None,None]*np.einsum('ei,ek,j,l->eijkl',fa,fa,self.fiber,self.fiber)
                A += 2*kf*q[:,None,None,None,None]*np.einsum('ik,j,l->ijkl',eye,self.fiber,self.fiber)[None,:,:,:,:]
                locJ = fac[:,None,None,None,None]*np.einsum('eaj,eijkl,ebl->eaibk',self.grad,A,self.grad)
                offset = layer*2*self.n
                rr.append(self.rr+offset); cc.append(self.cc+offset); dd.append(locJ.reshape(-1))
        slip = 0.
        for layer in (0,1):
            diff = U[layer]-U[layer+1]
            slip += float(self.coupling*np.sum(self.mass[:,None]*diff**2))
            force = 2*self.coupling*self.mass[:,None]*diff
            r[layer] += force; r[layer+1] -= force
            if hessian:
                d1 = self.idd+2*self.n*layer
                d2 = d1+2*self.n
                for a,b,sgn in ((d1,d1,1),(d2,d2,1),(d1,d2,-1),(d2,d1,-1)):
                    rr.append(a); cc.append(b); dd.append(sgn*self.slip_diag)
        matrix = coo_matrix((np.concatenate(dd),(np.concatenate(rr),np.concatenate(cc))),shape=(6*self.n,6*self.n)).tocsc() if hessian else None
        return sum(parts)+slip,r.ravel(),matrix,{'layer_bulk_J':parts,'slip_J':slip,'minJ':minJ}

    def solve(self, lengths, stretch=1.21, initial=None):
        t0 = time.process_time()
        fixed = np.zeros((3,self.n,2),dtype=bool)
        fixed[:,self.outer,:] = True
        cut = np.array([np.abs(self.x)<a/2-1e-12 for a in lengths])
        for layer in range(3):
            fixed[layer,self.bottom[~cut[layer]],1] = True
        fixed = fixed.ravel()
        U = np.tile((stretch-1)*self.xy,(3,1,1)).ravel() if initial is None else initial.copy().ravel()
        U[fixed] = np.tile((stretch-1)*self.xy,(3,1,1)).ravel()[fixed]
        U.reshape(3,self.n,2)[:,self.bottom,1][~cut] = 0
        free = np.flatnonzero(~fixed)
        scale = max(1.,float(np.sum(self.h)*np.max(self.mu)*self.H*(stretch-1)))
        failures = 0
        for it in range(60):
            en,r,J,detail = self.energy(U,True)
            rel = np.linalg.norm(r[free])/scale
            if rel<1e-10:
                break
            step = spsolve(J[free,:][:,free],-r[free])
            alpha = 1.
            descent = float(r[free]@step)
            for bt in range(30):
                candidate = U.copy(); candidate[free] += alpha*step
                en1,_,_,_ = self.energy(candidate)
                if en1 <= en+1e-4*alpha*descent+1e-16:
                    U = candidate
                    break
                alpha *= .5
            else:
                failures += 1
                break
        en,r,_,detail = self.energy(U)
        gaps = 2*U.reshape(3,self.n,2)[:,self.bottom,1]
        maxgap = [float(gaps[j,cut[j]].max()) if cut[j].any() else 0. for j in range(3)]
        row = {'lengths_m':list(lengths),'stretch':stretch,'energy_J':float(en),**detail,
               'max_gap_m':maxgap,'area_m2':[float(self.seam_w[cut[j]]@gaps[j,cut[j]]) for j in range(3)],
               'projected_crack_length_m':(cut@self.seam_w).tolist(),
               'iterations':it+1,'residual_relative':float(np.linalg.norm(r[free])/scale),
               'line_search_failures':failures,'cpu_s':time.process_time()-t0,
               'global_solves':it,'empirical_accuracy':'UNKNOWN'}
        return U,gaps,row

def verify(model):
    rng = np.random.default_rng(61)
    U = np.tile(.1*model.xy,(3,1,1)).ravel()
    en,g,J,_ = model.energy(U,True)
    v = rng.normal(size=U.size); v /= np.linalg.norm(v)
    eps = 1e-7
    ep,gp,_,_ = model.energy(U+eps*v)
    em,gm,_,_ = model.energy(U-eps*v)
    graderr = abs((ep-em)/(2*eps)-g@v)/max(abs(g@v),1e-8)
    hesserr = np.linalg.norm((gp-gm)/(2*eps)-J@v)/max(np.linalg.norm(J@v),1e-8)
    rest = np.zeros_like(U)
    angle = .31
    rot = np.array([[np.cos(angle),-np.sin(angle)],[np.sin(angle),np.cos(angle)]])
    rigid = np.tile(model.xy@(rot-np.eye(2)).T,(3,1,1)).ravel()
    rigidE,rigidG,_,_ = model.energy(rigid)
    return {'energy_gradient_directional_relative_error':float(graderr),
            'gradient_hessian_directional_relative_error':float(hesserr),
            'rigid_rotation_energy_J':float(rigidE),'rigid_rotation_gradient_N':float(np.linalg.norm(rigidG)),
            'rest_energy_J':float(model.energy(rest)[0]),
            'scope':'finite-difference/invariance verification of declared constitutive instrument; no skin accuracy certificate'}

def run(out, nx):
    out.mkdir(exist_ok=False)
    start = time.process_time()
    records = []
    for angle in (0,45,90):
        m = Layers(nx,nx//2,angle)
        if angle==0:
            write_json(out/'verification.json',verify(m))
        # A prescribed trajectory retains topology; surface work is accumulated
        # irreversibly once. No crack selection, cohesive evolution or healing.
        U = None
        for a in (0,.008,.012,.016,.020,.024):
            U,g,row = m.solve([a,a,a],initial=U)
            row['angle_deg'] = angle
            row['kind'] = 'full_depth'
            records.append(row)
            with (out/f'full_angle{angle}_a{int(a*1000)}.npz').open('xb') as f:
                np.savez_compressed(f,U=U,gap=g,x=m.x)
        if angle==0:
            for ls in ([.020,0,0],[.020,.020,0]):
                U,g,row = m.solve(ls)
                row.update(angle_deg=angle,kind='partial_depth')
                records.append(row)
            for stretch in (1.,1.10):
                U,g,row = m.solve([.020]*3,stretch=stretch)
                row.update(angle_deg=angle,kind='prestretch_ablation')
                records.append(row)
    for angle in (0,45,90):
        seq = [r for r in records if r['angle_deg']==angle and r['kind']=='full_depth']
        for prev,row in zip(seq,seq[1:]):
            da = row['projected_crack_length_m'][1]-prev['projected_crack_length_m'][1]
            released = (prev['energy_J']-row['energy_J'])/da
            row['all_layer_elastic_release_force_N'] = released
            row['dermis_intrinsic_surface_force_N'] = .0015*150
            row['conditional_blade_force_known_dermis_minus_all_release_N'] = .0015*150-released
            row['other_layer_fracture_force_N'] = None
            row['contact_and_process_force_N'] = None
            row['meaning'] = 'Known-term balance only; missing other layers/contact. Negative term is not a negative measured cutting force.'
    summary = {'scope':'SYNTHETIC_FINITE_STRAIN_COUPLED_LAYER_INSTRUMENT',
               'empirical_force_gap_accuracy':'UNKNOWN','strongest_same_information_equilibrium_control':'TIE',
               'nx':nx,'ny':nx//2,'layers':3,'displacement_components':2,
               'scenario_defaults':'SOURCE_TARGETS_R2.json','records':records,
               'max_residual_relative':max(r['residual_relative'] for r in records),
               'CPU_s':time.process_time()-start}
    write_json(out/'summary.json',summary)
    print(json.dumps({'out':str(out),'CPU_s':summary['CPU_s'],'residual':summary['max_residual_relative'],'cases':len(records)}))

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--nx',type=int,default=32)
    parser.add_argument('--out',type=Path,required=True)
    args = parser.parse_args()
    run(args.out,args.nx)
