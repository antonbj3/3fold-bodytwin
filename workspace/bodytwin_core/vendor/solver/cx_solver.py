"""CX-SOLVER: row equilibrated semismooth dual Newton, explicit unilateral contacts.
Source formulation: results/N14b/n14_solve.py; new row scaling and convergence logic here.
"""
import numpy as np
from scipy.linalg import solve, lstsq, LinAlgWarning
from scipy.optimize import linprog
import warnings

def polish_active(A,b,a0,maxit=500):
    """Feasible primal active set. HiGHS supplies feasibility only; objective is solved here."""
    rd=1/np.maximum(np.linalg.norm(A,axis=1),1e-10)
    D=A*rd[:,None]; y=b*rd
    support=a0>1e-9
    for trial in range(2):
        idx=np.where(support)[0]
        lp=linprog(np.zeros(len(idx)),A_eq=D[:,idx],b_eq=y,bounds=(0,None),method='highs')
        if lp.success:break
        support[:]=True
    if not lp.success:raise RuntimeError('no feasible active-set start')
    x=np.zeros(A.shape[1]);x[idx]=lp.x
    for k in range(maxit):
        idx=np.where(support)[0]
        cand=lstsq(D[:,idx],y,cond=1e-13,lapack_driver='gelsy')[0]
        if np.min(cand)<-1e-10:
            delta=cand-x[idx];neg=delta<0
            ratios=-x[idx][neg]/delta[neg]
            alpha=float(np.min(ratios)); x[idx]=x[idx]+alpha*delta
            hit=idx[np.where(neg)[0][int(np.argmin(ratios))]]
            support[hit]=False;x[hit]=0
        else:
            x[:]=0;x[idx]=cand
            lam=lstsq(D[:,idx].T,2*x[idx],cond=1e-13,lapack_driver='gelsy')[0]
            dual=D.T@lam;viol=np.where(~support,dual,0)
            j=int(np.argmax(viol))
            if viol[j]<=1e-11*max(1,float(np.max(2*x))):
                return x,lam*rd,k+1
            support[j]=True
    raise RuntimeError('primal active-set iteration limit')

def inner(A,b,U=np.inf,lam0=None,maxit=200,polish_enabled=True):
    m,n=A.shape
    d=1/np.maximum(np.linalg.norm(A,axis=1),1e-12)
    C=A*d[:,None]; y=b*d
    scale=max(np.linalg.norm(y),1e-12); C=C/scale; y=y/scale
    if np.isfinite(U): U=np.broadcast_to(np.asarray(U,float),(n,))
    else: U=np.full(n,np.inf)
    if lam0 is None:
        lam=2*solve(C@C.T+1e-14*np.eye(m),y,assume_a='pos')
    else: lam=lam0.copy()
    def state(lam):
        s=C.T@lam; a=np.clip(s/2,0,U); g=C@a-y
        obj=.5*np.dot(a,a)-np.dot(lam,g) # dual minimization = lam*y - .5 ||a||² (negative? )
        obj=-np.dot(lam,y)+np.dot(s,a)-np.dot(a,a)
        return s,a,g,obj
    s,a,g,obj=state(lam)
    for k in range(maxit):
        if np.linalg.norm(g)<=2e-12:break
        free=(s>0)&(s<2*U)
        H=(C[:,free]@C[:,free].T)*.5
        # Solve the Newton equation without a fixed regularizer that imposes a residual floor.
        with warnings.catch_warnings():
            warnings.simplefilter('error', LinAlgWarning)
            try: step=solve(H,-g,assume_a='pos')
            except (LinAlgWarning, np.linalg.LinAlgError): step=lstsq(H,-g,cond=1e-15,lapack_driver='gelsy')[0]
        slope=np.dot(g,step); alpha=1.
        for ls in range(45):
            sn,an,gn,on=state(lam+alpha*step)
            if on<=obj+1e-4*alpha*slope or np.linalg.norm(gn)<np.linalg.norm(g)*(1-1e-4*alpha):break
            alpha*=.5
        if ls==44:break
        lam+=alpha*step;s,a,g,obj=sn,an,gn,on
    physical_res=np.linalg.norm(A@a-b)/max(np.linalg.norm(b),1e-300)
    if polish_enabled and physical_res>1e-10 and np.all(np.isinf(U)):
        a,phys_lam,polish_iters=polish_active(A,b,a)
        physical_res=np.linalg.norm(A@a-b)/max(np.linalg.norm(b),1e-300)
        return a,phys_lam,float(physical_res),k+1+polish_iters
    return a,lam*d/scale,float(physical_res),k+1

def solve_step(A,b,N,names,U=np.inf,polish_enabled=True):
    xm=np.array(['Thorax_ContactReaction' in str(v) for v in names]); mus=(N>0)&~xm
    M=A[:,mus]*N[mus]; X=A[:,xm]; nc=X.shape[1]; active=list(range(nc))
    for outer in range(20):
        Xa=X[:,active] if active else np.empty((len(b),0))
        if active:
            q,s,_=np.linalg.svd(Xa,full_matrices=True);rank=sum(s>1e-12*s[0]);Q=q[:,rank:]
        else:Q=np.eye(len(b))
        a,lr,res,it=inner(Q.T@M,Q.T@b,U,polish_enabled=polish_enabled)
        lam=Q@lr
        x=np.linalg.lstsq(Xa,b-M@a,rcond=None)[0] if active else np.empty(0)
        if len(x) and np.min(x)<-1e-9:
            active.pop(int(np.argmin(x)));continue
        cand=[j for j in range(nc) if j not in active and X[:,j]@lam>1e-9]
        if cand:active.append(max(cand,key=lambda j:X[:,j]@lam));continue
        f=np.zeros(A.shape[1]);f[mus]=a*N[mus];f[np.where(xm)[0][active]]=x
        primal=np.linalg.norm(A@f-b)/max(np.linalg.norm(b),1e-300)
        # multiplier stationarity for muscle forces: 2a-M^T lam and contact sign.
        s=M.T@lam
        stat=np.max(np.abs((2*a-s)[(a>0)&(a<U)]),initial=0)
        lo=np.max(np.maximum(s[a<=0],0),initial=0)
        hi=np.max(np.maximum(2*U-s[a>=U-1e-9],0),initial=0) if np.isfinite(U).all() else 0.
        cx=np.max(np.maximum(X.T@lam,0)[[j for j in range(nc) if j not in active]],initial=0)
        dualscale=max(np.max(np.abs(2*a),initial=0),np.max(np.abs(s),initial=0),1)
        kkt=max(primal,stat/dualscale,lo/dualscale,hi/dualscale,cx/dualscale,max(-np.min(a,initial=0),-np.min(x,initial=0)))
        return f,dict(kkt=float(kkt),primal=float(primal),stat=float(stat/dualscale),lo=float(lo/dualscale),hi=float(hi/dualscale),contact=float(cx/dualscale),inner_iters=it,outer_iters=outer+1,active_contact=active)
    raise RuntimeError('contact active set did not converge')
