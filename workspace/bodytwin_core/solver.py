"""A295 exact solver, copied unchanged to vendor/solver/cx_solver.py."""
import numpy as np
from functools import lru_cache
from .vendor.solver import cx_solver
from .vendor.solver.n40_solver_np import solve_poly
from scipy.optimize import linprog

def _inner_n43(A,b,U=np.inf,lam0=None,maxit=200,polish_enabled=True):
    """Exact adapter copied from results/CX-SOLVER2/run_hybrid.py."""
    if np.isfinite(U).any():raise NotImplementedError("upper bounds not used in the collaborator's system")
    r=solve_poly(A[None],b[None],2,tol=1e-12,maxit=maxit)
    sc=float(r["sc"][0]);a=r["a"][0];lam=r["lam"][0]/sc
    resid=float(np.linalg.norm(A@a-b)/np.linalg.norm(b))
    return a,lam,resid,int(r["iters"][0])

cx_solver.inner=_inner_n43
solve_step=cx_solver.solve_step

def solve_bounds(A, b, joint_maps, offsets=None, lower=None, upper=None, *, A_ub=None, b_ub=None, tol=1e-7):
    """Exact LP ranges of affine joint-force *components* per frame.

    A has shape (time, equations, forces), b (time, equations), and each
    joint map (time, components, forces). The reported intervals are for
    ``offset + map @ forces``. A vector norm is not an LP observable.
    Infinite endpoints are preserved; infeasible frames raise ValueError.
    """
    A=np.asarray(A,float); b=np.asarray(b,float)
    if A.ndim==2: A=A[None]
    if b.ndim==1: b=b[None]
    T,m,n=A.shape
    if b.shape!=(T,m): raise ValueError('A/b shape mismatch')
    G=None if A_ub is None else np.asarray(A_ub,float)
    h=None if b_ub is None else np.asarray(b_ub,float)
    if (G is None)!=(h is None): raise ValueError('A_ub and b_ub required together')
    if G is not None:
        if G.ndim==2:G=G[None]
        if h.ndim==1:h=h[None]
        if G.shape[:1]!=(T,) or G.shape[2]!=n or h.shape!=G.shape[:2]:raise ValueError('inequality shape mismatch')
    lo=np.zeros((T,n)) if lower is None else np.broadcast_to(np.asarray(lower,float),(T,n))
    hi=np.full((T,n),np.inf) if upper is None else np.broadcast_to(np.asarray(upper,float),(T,n))
    if np.any(lo>hi): raise ValueError('lower exceeds upper')
    result={}
    for name,raw_map in joint_maps.items():
        W=np.asarray(raw_map,float)
        if W.ndim==2: W=W[None]
        if W.shape[0]!=T or W.shape[2]!=n: raise ValueError(f'{name}: map shape mismatch')
        d=np.zeros((T,W.shape[1])) if offsets is None or name not in offsets else np.broadcast_to(np.asarray(offsets[name],float),(T,W.shape[1]))
        result[name]=np.empty((T,W.shape[1],2))
        for t in range(T):
            scale=1/np.maximum(np.linalg.norm(A[t],axis=1),1e-12)
            eq=A[t]*scale[:,None]; rhs=b[t]*scale
            if G is not None:
                gs=1/np.maximum(np.linalg.norm(G[t],axis=1),1e-12)
                ub=G[t]*gs[:,None]; ub_rhs=h[t]*gs
            else:ub=ub_rhs=None
            bounds=list(zip(lo[t],hi[t]))
            for j,c in enumerate(W[t]):
                for end,sign in enumerate((1,-1)):
                    r=linprog(sign*c,A_eq=eq,b_eq=rhs,A_ub=ub,b_ub=ub_rhs,bounds=bounds,method='highs')
                    if r.status==2:
                        check=linprog(np.zeros(n),A_eq=eq,b_eq=rhs,A_ub=ub,b_ub=ub_rhs,
                                      bounds=bounds,method='highs-ds',options={'presolve':False})
                        if check.success and np.linalg.norm(A[t]@check.x-b[t])/max(1,np.linalg.norm(b[t]))<=tol:
                            raise RuntimeError(f'{name}: inconsistent HiGHS infeasibility at frame {t}')
                        raise ValueError(f'{name}: infeasible frame {t}')
                    if r.status==3: value=-np.inf if end==0 else np.inf
                    elif r.status==4 and np.isinf(hi[t]).any():
                        # HiGHS sometimes reports a solve error on a highly
                        # degenerate unbounded the collaborator LP. Verify a nonnegative
                        # recession direction before labeling it unbounded.
                        recession_bounds=[(0,1) if np.isinf(v) else (0,0) for v in hi[t]]
                        q=linprog(sign*c,A_eq=eq,b_eq=np.zeros(m),A_ub=ub,b_ub=np.zeros(len(ub)) if ub is not None else None,
                                  bounds=recession_bounds,method='highs')
                        q_res=np.linalg.norm(eq@q.x) if q.success else np.inf
                        if q.success and q_res<1e-8 and q.fun < -1e-7*max(1,np.linalg.norm(c)):
                            value=-np.inf if end==0 else np.inf
                        else:
                            r=linprog(sign*c,A_eq=eq,b_eq=rhs,A_ub=ub,b_ub=ub_rhs,bounds=bounds,
                                      method='highs-ds',options={'presolve':False})
                            if not r.success:raise RuntimeError(f'{name}: unresolved HiGHS solve error at frame {t}')
                            residual=np.linalg.norm(A[t]@r.x-b[t])/max(1,np.linalg.norm(b[t]))
                            if residual>tol:raise RuntimeError(f'{name}: fallback LP residual {residual:g} at frame {t}')
                            value=float(c@r.x+d[t,j])
                    elif r.success:
                        residual=np.linalg.norm(A[t]@r.x-b[t])/max(1,np.linalg.norm(b[t]))
                        if residual>tol: raise RuntimeError(f'{name}: LP residual {residual:g} at frame {t}')
                        value=float(c@r.x+d[t,j])
                    else: raise RuntimeError(f'{name}: HiGHS {r.message} at frame {t}')
                    result[name][t,j,end]=value
    return result
from .config import paths

@lru_cache(maxsize=2)
def _xf4_training(path):
    with np.load(path) as data:
        return data['Xtr'], data['NU_train']

def solve_lift(tag='identity', config=None, frames=None, mode='exact'):
    if mode not in ('exact', 'certified'):
        raise ValueError("mode must be 'exact' or 'certified'")
    p=paths(config)
    z=np.load(p['solver_data']/f'sys_{tag}.npz')
    tids=range(len(z['m'])) if frames is None else frames
    forces=[]; receipts=[]
    for t in tids:
        m=int(z['m'][t])
        f,receipt=solve_step(z['A'][t,:m], z['b'][t,:m], z['N'][:,t], z['names'])
        if mode == 'certified':
            # XF4 is trained on K1's four-equation QP. The the collaborator lift has
            # 93-94 equations and unilateral contacts: its safe result is
            # the exact fallback until a matching emulator is calibrated.
            receipt=dict(receipt, accepted=False, fallback=True,
                         fallback_reason='no_emulator_for_john_contact_system')
        forces.append(f); receipts.append(receipt)
    return np.asarray(forces),receipts

def hip_curve(tag, forces, config=None):
    p=paths(config)
    z=np.load(p['solver_data']/f'leg_sys_{tag}.npz')
    return np.linalg.norm(z['dhip']-np.einsum('tkn,tn->tk',z['Whip'],forces),axis=1)


def solve_k1(x, mode='exact', config=None, candidate_hook=None, model=None):
    """Solve the 141-frame XF4/K1 p=2 model; return forces, reactions, receipts.

    ``certified`` uses XF4's nearest training dual to classify active sets,
    six one-sided KKT sweeps, then checks a numerical primal-dual bound.
    A candidate hook supports counter-tests through the previous-frame path.
    """
    if mode not in ('exact', 'certified'):
        raise ValueError("mode must be 'exact' or 'certified'")
    from . import emulator
    from .certified import solve_step as certified_step, certificate, warm_exact
    model = emulator.from_source(config) if model is None else model
    e = model
    delta, sa, ss, r0, r, N = e._geom(np.asarray(x, float))
    _, Cfull, b = e._C_b(r0, r, N, delta)
    forces = np.empty((e.T, e.n)); reactions = np.empty((e.T, 3))
    receipts = []
    mask = None; previous_nu = None
    if mode == 'certified' and candidate_hook is None:
        # XF4's 1-NN dual classifier followed by six one-sided KKT sweeps.
        # The archived training duals are valid only for this K1 model.
        from .config import paths as core_paths
        xtr, nutr = _xf4_training(str(core_paths(config)['xf4_results']/'inputs/BT-N29/n29_data.npz'))
        j = int(np.argmin(np.sum((xtr - np.asarray(x))**2, axis=1)))
        nu_batch = np.array(nutr[j], copy=True)
        for _ in range(6):
            masks = np.einsum('tdf,td->tf', Cfull, nu_batch) > 0
            active = Cfull * masks[:, None, :]
            mat = active @ active.transpose(0, 2, 1) + 1e-14*np.eye(4)
            nu_batch = np.linalg.solve(mat, b[..., None])[..., 0]
        # A final frozen solve is XF4's one-sided reconstruction at K=6.
        masks = np.einsum('tdf,td->tf', Cfull, nu_batch) > 0
        active = Cfull * masks[:, None, :]
        mat = active @ active.transpose(0, 2, 1) + 1e-14*np.eye(4)
        nu_batch = np.linalg.solve(mat, b[..., None])[..., 0]
        a_batch = np.einsum('tdf,td->tf', active, nu_batch)
        scores = np.einsum('tdf,td->tf', Cfull, nu_batch)
        projected = np.maximum(scores, 0)
        primal = 0.5*np.sum(a_batch*a_batch, axis=1)
        lower = np.sum(nu_batch*b, axis=1)-0.5*np.sum(projected*projected, axis=1)
        allowance = 64*np.finfo(float).eps*np.maximum.reduce((np.ones(e.T),np.abs(primal),np.abs(lower)))
        gap = primal-lower
        resid = np.linalg.norm(np.einsum('tdf,tf->td',Cfull,a_batch)-b,axis=1)/np.maximum(1,np.linalg.norm(b,axis=1))
        sign = np.maximum(np.max(np.maximum(-a_batch,0),axis=1),
                          np.max(np.abs(a_batch-projected),axis=1))/np.maximum(1,np.max(np.abs(a_batch),axis=1))
        bound = np.max(N[e.free])*np.sqrt(2*(np.maximum(gap,0)+allowance))/np.maximum(1e-12,np.linalg.norm(a_batch*N[e.free],axis=1))
        batch_ok = np.isfinite(a_batch).all(axis=1)&np.isfinite(nu_batch).all(axis=1)&(resid<=1e-10)&(sign<=1e-10)&(gap>=-allowance)&(bound<=0.005)
    for t in range(e.T):
        C = Cfull[t]; bt = b[t]
        if mode == 'certified':
            if candidate_hook is None:
                a, nu = a_batch[t], nu_batch[t]
                if batch_ok[t]:
                    rec = dict(accepted=True, fallback=False, mask=a > 0,
                               primal_residual=float(resid[t]),dual_residual=float(sign[t]),
                               dual_gap=float(gap[t]),relative_force_bound=float(bound[t]))
                else:
                    a, nu = warm_exact(C, bt, previous_nu)
                    check = certificate(C, bt, a, nu, N[e.free])
                    if check['primal_residual'] > 1e-10 or check['dual_residual'] > 1e-10:
                        raise RuntimeError('exact fallback failed KKT')
                    rec = dict(check,
                               accepted=False, fallback=True,
                               fallback_reason='uncertain', mask=a > 0)
            else:
                a, nu, rec = certified_step(C, bt, lambda C, b: warm_exact(C, b, previous_nu),
                                            mask, N[e.free], candidate_hook=candidate_hook)
        else:
            a, nu = emulator.xf4_model._dual_newton(C, bt)
            rec = dict(certificate(C, bt, a, nu, N[e.free]), accepted=False,
                       fallback=False, mask=a > 0)
        mask = a > 0
        previous_nu = nu
        fi = e.FJ[t].copy(); fi[e.free] = a * N[e.free]
        forces[t] = fi
        reactions[t] = e.RJ[t] + e.g[t].T @ e.FJ[t] - e.g[t].T @ fi
        receipts.append(rec)
    return forces, reactions, receipts
