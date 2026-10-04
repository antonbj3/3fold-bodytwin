"""End-to-end conditional femur design evaluation on the the collaborator lift.

Uses D1's immutable geometry and mechanics builders; all new outputs are local.
"""
from __future__ import annotations
import importlib.util
import os
import subprocess
import sys
import time
from pathlib import Path
import numpy as np
from scipy.optimize import linprog
from .solver import solve_step, solve_bounds

ROOT = Path(__file__).resolve().parents[1]
D1 = ROOT / 'results/CX-D1PARITY'
TMP = Path('external_media')
SRC = Path('external_media')
D1DATA = Path('external_media')
sys.path[:0] = [str(D1), str(ROOT / 'results/CX-FIELDSHARE')]
import direct_geometry
import run_d1 as D


def _capacity(A, b, N, names):
    """Minimum uniform strength scale with unconstrained nonnegative reactions."""
    n = len(N)
    contact = np.array(['Thorax_ContactReaction' in str(v) for v in names])
    mus = (N > 0) & ~contact
    cols = mus | contact
    AA = A[:, cols]
    caps = N[cols]
    ismus = mus[cols]
    c = np.zeros(len(caps)+1); c[-1] = 1
    ub = np.zeros((int(ismus.sum()), len(caps)+1))
    ub[np.arange(len(ub)), np.flatnonzero(ismus)] = 1
    ub[:, -1] = -caps[ismus]
    scale = 1 / np.maximum(np.linalg.norm(AA, axis=1), 1e-12)
    r = linprog(c, A_eq=np.column_stack([AA*scale[:,None], np.zeros(len(b))]),
                b_eq=b*scale, A_ub=ub, b_ub=np.zeros(len(ub)),
                bounds=[(0,None)]*(len(caps)+1), method='highs')
    if not r.success:
        return {'status': r.message, 'kappa': None, 'feasible_at_1': False}
    residual = np.linalg.norm(A @ np.pad(r.x[:-1], (0,n-len(caps))) - b) if cols.all() else np.linalg.norm(AA@r.x[:-1]-b)
    return {'status':'optimal', 'kappa':float(r.x[-1]), 'feasible_at_1':bool(r.x[-1] <= 1+1e-7),
            'relative_residual':float(residual/max(1,np.linalg.norm(b)))}


def _component_band(A,b,W,d,N,names):
    """LP component ranges with muscle Fmax; infinity remains explicit."""
    contact=np.array(['Thorax_ContactReaction' in str(v) for v in names])
    upper=np.where(contact,np.inf,N)
    upper=np.where((N>0)|contact,upper,0.)
    try:
        ans=solve_bounds(A,b,{'hip_xyz':-W}, {'hip_xyz':d},upper=upper)['hip_xyz'][0]
        return ans.tolist()
    except (ValueError,RuntimeError) as e:
        return {'status':str(e)}


def design_eval(shape_params, subject='z001', task='john_lift', *, diagnostics=True):
    """Rebuild attachments, mechanics, bone demand and all 141 recruitment steps.

    shape_params: {'lengthening_mm': number in [-10,10]}.
    subject: 'z001' or 'C01RFE' (the latter is a TLEM proxy).
    """
    if task != 'john_lift': raise ValueError('only john_lift is supported')
    if set(shape_params) != {'lengthening_mm'}: raise ValueError('requires lengthening_mm only')
    length=float(shape_params['lengthening_mm'])
    if not np.isfinite(length) or abs(length)>10: raise ValueError('length outside [-10,10] mm')
    if subject not in ('z001','C01RFE'): raise ValueError('unsupported subject')
    tag='vsd_z001' if subject=='z001' else 'identity'
    TMP.mkdir(parents=True,exist_ok=True)
    for name in ('buckle_geo.npy','buckle_reac.npy','cyl_reac.npy','refs_fix3.json'):
        link=TMP/name
        if not link.exists(): link.symlink_to(SRC/name)
    u=np.array([0.,0.,length]); key=f'{subject}_length_{length:+.5f}'.replace('+','p').replace('-','m').replace('.','d')
    stem=f'{tag}_{subject}_{key}'
    timings={}; start=time.perf_counter()
    direct_geometry.OUT=TMP
    direct_geometry.make(subject,tag,u,key)
    timings['attachment_s']=time.perf_counter()-start
    source_morph=TMP/f'morph_{stem}.npz'
    if not source_morph.exists(): raise RuntimeError('attachment output missing')
    env=os.environ.copy()
    env.update(N14W=str(SRC)+'/',N14BW=str(TMP)+'/',N14B_BUCKLE='1',N14B_FIX='1',
               N14B_FIXFILE='refs_fix3.json',N14B_CYL='1',N14_EXCL_CYL='0',N14B_SAVELEG='1',
               N14B_MORPH=str(source_morph),N14SYS=f'sys_{stem}.npz',
               N14TAG=str(TMP/(stem+'_')),N40B_N2B='1',N40B_D1_ALL='1',
               N40B_D1_PROJECTION_NAME=f'projection_{stem}.npz',
               OPENBLAS_NUM_THREADS='2',OMP_NUM_THREADS='2')
    start=time.perf_counter()
    with (TMP/f'build_{stem}.log').open('w') as log:
        subprocess.run([sys.executable,str(D1/'n40b_projection.py')],cwd=TMP,env=env,
                       stdout=log,stderr=subprocess.STDOUT,check=True)
    timings['mechanics_build_s']=time.perf_counter()-start
    start=time.perf_counter()
    shape,head=D.bone(subject); p0=D.fields(shape,np.zeros(3)); p1=D.fields(shape,u)
    timings['mass_properties_s']=time.perf_counter()-start
    with np.load(TMP/f'sys_{stem}.npz') as z, np.load(TMP/f'leg_sys_{stem}.npz') as l, np.load(TMP/f'projection_{stem}.npz') as pr:
        Aall=z['A']; ball=z['b']; mall=z['m']; Nall=z['N']; names=z['names']
        Wall=l['Whip']; dall=l['dhip']; P=pr['P']; eq=list(pr['eq']); PA=pr['A']; pm=pr['m']
    start=time.perf_counter(); forces=[]; profile=[]; kkt=[]; projerr=[]; massdelta=[]
    for t in range(141):
        m=int(mall[t]); A=Aall[t,:m]; b=ball[t,:m].copy()
        if int(pm[t])!=m: raise RuntimeError('projection dimension mismatch')
        projerr.append(float(np.linalg.norm(A-PA[t,:m])/max(1,np.linalg.norm(A))))
        kin=D.kinematics(t)
        w0,F0=D.wrench(p0,head,kin); w1,F1=D.wrench(p1,head,kin)
        delta=np.zeros(P.shape[2]); row=eq.index(kin['segment'])
        delta[6*row:6*row+6]=-(w1-w0)
        b+=P[t,:m]@delta
        f,receipt=solve_step(A,b,Nall[:,t],names)
        forces.append(f); kkt.append(float(receipt['kkt']))
        profile.append(float(np.linalg.norm(dall[t]-(F1-F0)-Wall[t]@f)))
        massdelta.append(float(np.linalg.norm(w1-w0)))
    timings['solve_141_s']=time.perf_counter()-start
    profile=np.asarray(profile); peak_t=int(profile.argmax())
    checks={}
    if diagnostics:
        start=time.perf_counter()
        for t in sorted({7,75,peak_t}):
            m=int(mall[t]); kin=D.kinematics(t)
            w0,F0=D.wrench(p0,head,kin);w1,F1=D.wrench(p1,head,kin)
            delta=np.zeros(P.shape[2]);row=eq.index(kin['segment']);delta[6*row:6*row+6]=-(w1-w0)
            b=ball[t,:m]+P[t,:m]@delta
            d=dall[t]-(F1-F0)
            checks[str(t)]={'capacity':_capacity(Aall[t,:m],b,Nall[:,t],names),
                            'hip_component_band_N':_component_band(Aall[t,:m],b,Wall[t],d,Nall[:,t],names)}
        timings['checks_s']=time.perf_counter()-start
    return {'subject':subject,'tag':tag,'task':task,'lengthening_mm':length,
            'stem':stem,'moved_attachments':int(__import__('json').loads((TMP/f'morph_{stem}.json').read_text())['moved']),
            'bone_mass_kg':float(p1['m']),'bone_com_mm':p1['c'].tolist(),
            'bone_inertia_kgm2':p1['I'].tolist(),'peak_N':float(profile.max()),
            'peak_frame':peak_t,'profile_N':profile.tolist(),
            'max_kkt':max(kkt),'kkt_pass_1e10':sum(v<=1e-10 for v in kkt),
            'projection_error_max':max(projerr),'max_bone_wrench_delta':max(massdelta),
            'checks':checks,'timings_s':timings,
            'parameter_free_law':'not applicable: knee law does not estimate hip contact'}
