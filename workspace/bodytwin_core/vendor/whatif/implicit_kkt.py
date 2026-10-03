"""Differentiate the fixed-support KKT equations after ±1° system rebuilds.

The geometric system derivative is central finite difference. The recruitment
derivative is implicit KKT, including A, b, N and hip-map changes. The N40b
SVD row gauge is aligned by orthogonal Procrustes before differencing.
"""
from pathlib import Path
import json, os
import numpy as np
from scipy.linalg import lstsq, svd

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
BASE=Path(os.environ.get('CXW2_BASE','/media/anton/sdc1-tmp/bodytwin/CX-SOLVER2/data'))
CLOUD=Path(os.environ.get('CXW2_MATRICES','/media/anton/sdc1-tmp/bodytwin/CX-WHATIF2/matrices_out'))
FBASE=Path(os.environ.get('CXW2_FBASE',ROOT/'results/CX-SOLVER2'))
OUTPUT=Path(os.environ.get('CXW2_OUTPUT',HERE/'implicit_kkt_results.json'))

def load(tag,param,delta):
    stem=f'{tag}_{param}_{delta:+g}'.replace('+','p').replace('-','m')
    with np.load(CLOUD/'data'/f'sys_{stem}.npz') as source:
        z={k:source[k] for k in ('A','b','m','N','names')}
    with np.load(CLOUD/'data'/f'leg_sys_{stem}.npz') as source:
        l={k:source[k] for k in ('Whip','dhip')}
    f=np.load(CLOUD/'out'/f'{stem}.npz')['f']
    receipt=json.load(open(CLOUD/'out'/f'{stem}_matrix_receipt.json'))
    return z,l,f,receipt

def align(A0,b0,A1,b1):
    # R A1 ~= A0, with equal rank and row count.
    U,_,Vt=svd(A0@A1.T,full_matrices=False)
    R=U@Vt
    return R@A1,R@b1,float(np.linalg.norm(R@A1-A0)/np.linalg.norm(A0))

def frame(t,z0,l0,f0,zm,lm,fm,zp,lp,fp):
    m0,mm,mp=[int(z['m'][t]) for z in (z0,zm,zp)]
    if not (m0==mm==mp):return dict(reason='row count changed')
    A0=z0['A'][t,:m0];b0=z0['b'][t,:m0]
    Am,bm,errm=align(A0,b0,zm['A'][t,:mm],zm['b'][t,:mm])
    Ap,bp,errp=align(A0,b0,zp['A'][t,:mp],zp['b'][t,:mp])
    dA=(Ap-Am)/2;db=(bp-bm)/2
    N0=z0['N'][:,t];dN=(zp['N'][:,t]-zm['N'][:,t])/2
    names=z0['names']; contact=np.array(['Thorax_ContactReaction' in str(n) for n in names])
    muscle=(N0>0)&~contact
    active=(f0[t]>1e-7)&(muscle|contact)
    am=(fm[t]>1e-7)&(muscle|contact);ap=(fp[t]>1e-7)&(muscle|contact)
    stable=bool(np.array_equal(active,am) and np.array_equal(active,ap))
    if not stable:return dict(reason='active set changed',align_rel=[errm,errp])
    idx=np.flatnonzero(active);A=A0[:,idx];dAs=dA[:,idx];fs=f0[t,idx]
    H=np.zeros(len(idx));dH=np.zeros(len(idx));is_mus=muscle[idx]
    H[is_mus]=2/N0[idx[is_mus]]**2
    dH[is_mus]=-4*dN[idx[is_mus]]/N0[idx[is_mus]]**3
    lam=lstsq(A.T,H*fs,cond=1e-12)[0]
    K=np.block([[np.diag(H),-A.T],[A,np.zeros((m0,m0))]])
    rhs=np.r_[-dH*fs+dAs.T@lam,db-dA@f0[t]]
    sol=lstsq(K,rhs,cond=1e-12)[0]
    df=np.zeros_like(f0[t]);df[idx]=sol[:len(idx)]
    W0=l0['Whip'][t];d0=l0['dhip'][t]
    dW=(lp['Whip'][t]-lm['Whip'][t])/2
    dd=(lp['dhip'][t]-lm['dhip'][t])/2
    v=d0-W0@f0[t]
    pred=float(v@(dd-dW@f0[t]-W0@df)/np.linalg.norm(v))
    vm=lm['dhip'][t]-lm['Whip'][t]@fm[t]
    vp=lp['dhip'][t]-lp['Whip'][t]@fp[t]
    fd=float((np.linalg.norm(vp)-np.linalg.norm(vm))/2)
    residual=float(np.linalg.norm(K@sol-rhs)/max(np.linalg.norm(rhs),1e-12))
    stationarity=float(np.linalg.norm(A.T@lam-H*fs)/max(np.linalg.norm(H*fs),1e-12))
    return dict(implicit_N_per_deg=pred,full_fd_N_per_deg=fd,
                relative_error=abs(pred-fd)/max(abs(fd),1e-9),
                abs_error_N_per_deg=abs(pred-fd),KKT_linear_residual=residual,
                KKT_stationarity_residual=stationarity,align_rel=[errm,errp])

def main():
    out={}
    for tag in ('identity','vsd_z001','vsd_z009'):
        with np.load(BASE/f'sys_{tag}.npz') as source:
            z0={k:source[k] for k in ('A','b','m','N','names')}
        with np.load(BASE/f'leg_sys_{tag}.npz') as source:
            l0={k:source[k] for k in ('Whip','dhip')}
        f0=np.load(FBASE/f'{tag}.npz')['f']
        out[tag]={}
        for param in ('CCD','AV'):
            zm,lm,fm,rm=load(tag,param,-1);zp,lp,fp,rp=load(tag,param,1)
            if rm.get('kkt_pass')!=141 or rp.get('kkt_pass')!=141:
                out[tag][param]=dict(negative_receipt=rm,positive_receipt=rp,
                    status='not analyzed: rebuilt ±1° point failed 141/141 KKT criterion')
                print(tag,param,'SKIP: rebuilt KKT criterion failed',flush=True)
                continue
            rows=[frame(t,z0,l0,f0,zm,lm,fm,zp,lp,fp) for t in range(141)]
            good=[r for r in rows if 'implicit_N_per_deg' in r]
            peak=int(np.argmax(np.linalg.norm(l0['dhip']-np.einsum('tkn,tn->tk',l0['Whip'],f0),axis=1)))
            out[tag][param]=dict(negative_receipt=rm,positive_receipt=rp,
                n_stable=len(good),n_active_switch=sum(r.get('reason')=='active set changed' for r in rows),
                n_row_count_change=sum(r.get('reason')=='row count changed' for r in rows),
                peak_frame=peak,peak=rows[peak],
                median_relative_error=float(np.median([r['relative_error'] for r in good])) if good else None,
                p95_relative_error=float(np.quantile([r['relative_error'] for r in good],.95)) if good else None,
                max_abs_error_N_per_deg=max((r['abs_error_N_per_deg'] for r in good),default=None),
                rows=rows)
            print(tag,param,len(good),'peak',rows[peak],flush=True)
    OUTPUT.write_text(json.dumps(out,indent=2))

if __name__=='__main__':main()
