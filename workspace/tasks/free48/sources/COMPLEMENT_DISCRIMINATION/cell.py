"""Complement source reproduction and explicit convertase-decay fate fork.

No native imports: rhs AST extraction prevents inherited writes/experiments.
Default source_sink preserves old dynamics; recycle is an uncalibrated fork.
Run python3 -B cell.py --output <new.json> [--fate recycle].
"""
import sys
sys.dont_write_bytecode = True
import ast
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp

# PENDING_INDEPENDENT_REVIEW -- path portability only; no model change.
import os
HERE = Path(__file__).resolve().parent
REPO = Path(os.environ.get('BODYTWIN_REPO', 'source_repository'))


def _acquired(rel):
    """Input beside this file if it travelled with the cell, else from the repo."""
    local = HERE / Path(rel).name
    return local if local.exists() else REPO / rel


SOURCE = _acquired('data/complement_cascade/complement_cascade_results.json')
SCRIPT = _acquired('scripts/msk/complement_cascade.py')
DATA = json.loads(SOURCE.read_text())
C = DATA['constants']
BASE = DATA['locked_params']
tree = ast.parse(SCRIPT.read_text())
fn = [n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='rhs']
assert len(fn)==1
env = {}
exec(compile(ast.Module(body=fn,type_ignores=[]),str(SCRIPT),'exec'),env)
original_rhs = env['rhs']

def params(surface='activator', cd55=6., k_mac=None, multipliers=(1.,1.,1.)):
    p = dict(BASE)
    for name,m in zip(['k_f','kcat_amp','Km_amp'],multipliers): p[name]*=m
    if surface=='host':
        p.update(k_decay=C['K_DECAY_ACTIVATOR']*C['FOLD_H']*cd55,
                 k_I=C['K_I_HOST'],k_MAC=C['K_MAC_HOST'])
    if k_mac is not None: p['k_MAC']=k_mac
    p['S_cl']=2.
    return p

def growth(p, fate='source_sink'):
    rho = 1. if fate=='recycle' else p.get('decay_recycle_fraction',0.)
    if not 0<=rho<=1: raise ValueError('decay_recycle_fraction must be in [0,1]')
    b = p['kcat_amp']*C['C3_0_NM']/(p['Km_amp']+C['C3_0_NM'])
    J = np.array([[-p['k_f']-p['k_I'],b+rho*p['k_decay']],
                  [p['k_f'],-p['k_decay']]])
    crit = p['k_f']*b/(p['k_I']+p['k_f']*(1-rho))
    return dict(lambda_per_min=float(np.max(np.linalg.eigvals(J).real)),
                decay_critical_per_min=crit,decay_margin_per_min=crit-p['k_decay'],
                half_life_critical_min=float(np.log(2)/crit),
                half_life_current_min=float(np.log(2)/p['k_decay']),J=J.tolist())

def rhs(t,y,p,fate):
    dy = original_rhs(t,y[:8],p)
    decay = p['k_decay']*max(y[2],0.)
    rho = 1. if fate=='recycle' else p.get('decay_recycle_fraction',0.)
    if not 0<=rho<=1: raise ValueError('decay_recycle_fraction must be in [0,1]')
    dy[1]+=rho*decay
    cb=max(y[1],0.);pool=max(y[2],0.)*cb/(cb+p['K_C3b_for_C5conv'])+p['kappa_cl']*(p['S_cl']>0)
    return dy+[(1-rho)*decay,p['kcat_C5']*pool]

def simulate(p,fate='source_sink',tmax=60.,rtol=1e-10,max_step=.2):
    y0=[C['C3_0_NM'],0.,0.,0.,C['C5_0_NM'],0.,0.,0.,0.,0.]
    sol=solve_ivp(rhs,[0,tmax],y0,args=(p,fate),method='DOP853',rtol=rtol,atol=1e-11,dense_output=True,max_step=max_step)
    assert sol.success,sol.message
    return sol

def terminal_mac(time, c5, k_mac):
    """Exact linear-reservoir evolution for piecewise-constant release flux.

    C5 trajectory comes from an acquired nonlinear history, not query replay.
    The release representation error is measured by refinement, not certified.
    """
    time=np.asarray(time);c5=np.asarray(c5);dt=np.diff(time)
    if k_mac<0 or np.any(dt<=0): raise ValueError('positive time order / nonnegative rate required')
    k=np.asarray(k_mac)
    if np.any(k==0): raise ValueError('use explicit zero-rate limit')
    releases=c5[:-1]-c5[1:]
    # Each interval contributes its integrated constant source weighted to T.
    weights=np.exp(-k*(time[-1]-time[1:]))*(-np.expm1(-k*dt))/(k*dt)
    pending=float(weights@releases)
    mac=float(c5[0]-c5[-1]-pending)
    return dict(MAC_nM=mac,C5b_pending_nM=pending,C5_nM=float(c5[-1]),
                C5_conservation_error_nM=float(abs(c5[-1]+pending+mac-c5[0])),
                full_ODE_calls=0,history_rebuilds=0)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--output',required=True);ap.add_argument('--fate',choices=['source_sink','recycle'],default='source_sink');args=ap.parse_args()
    out={}
    for surf in ['activator','host']:
        p=params(surf);s=simulate(p,args.fate)
        out[surf]=dict(MAC_nM=float(s.y[6,-1]),growth=growth(p,args.fate))
    out.update(activator_over_host_MAC=out['activator']['MAC_nM']/out['host']['MAC_nM'],fate=args.fate,
               CD55_status='FOLD_CD55_ILLUSTRATIVE; no measured clinical threshold',review='PENDING_INDEPENDENT_REVIEW')
    with Path(args.output).open('x') as f:json.dump(out,f,indent=2,allow_nan=False)
