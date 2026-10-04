"""Resolved conditional corneal scattering port; no transmission refitting.

Run with python3 -B cell.py --output <new.json>. All biological gates pending.
Reuses the eye lane's Maxwell and ray primitives; angle is never replaced by tau.
"""
import sys
sys.dont_write_bytecode = True
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import j1

# PENDING_INDEPENDENT_REVIEW -- path portability only; no model change.
import os
HERE = Path(__file__).resolve().parent
EYE_DEFAULT = 'results/LANE_EYE_OPTICAL_TWIN'


def _eye_dir():
    """Eye-lane primitives: env override, beside this file, any ancestor workspace,
    then the canonical lane directory. No fixed directory depth assumed."""
    cands = []
    env = os.environ.get('BODYTWIN_EYE_DIR')
    if env:
        cands.append(Path(env))
    cands += [HERE, HERE / 'LANE_EYE_OPTICAL_TWIN']
    cands += [p / 'results/LANE_EYE_OPTICAL_TWIN' for p in HERE.parents]
    cands.append(Path(EYE_DEFAULT))
    for c in cands:
        if (c / 'scatter_port_r1.py').exists():
            return c
    raise FileNotFoundError('eye-lane primitives not found; set BODYTWIN_EYE_DIR')


EYE = _eye_dir()
sys.path.insert(0, str(EYE))
from scatter_port_r1 import ray_angle
from maxwell_cylinder_r1 import density as maxwell_density, coeff

REPO = Path(os.environ.get('BODYTWIN_REPO', 'source_repository/'))
_REL = 'data/corneal_transparency/corneal_transparency_results.json'
_LOCAL = HERE / Path(_REL).name
SOURCE = _LOCAL if _LOCAL.exists() else REPO / _REL
DATA = json.loads(SOURCE.read_text())
PARAMS = DATA['params']

def density(theta, lam=550., sf=.15, mode='inherited', scale=1., poisson=False):
    """Width nm/rad (Maxwell) or calibrated width nm/rad (inherited RGD)."""
    theta = np.asarray(theta)
    a = PARAMS['a_nm'] * scale
    k = 2*np.pi*PARAMS['n_matrix']/lam
    q = 2*k*np.abs(np.sin(theta/2))
    s = np.ones_like(q) if poisson else -np.expm1(-(q*sf*PARAMS['d_hex_nm']*scale)**2)
    if mode == 'maxwell':
        return maxwell_density(theta, lam, a=a)*s
    x = q*a
    f = np.ones_like(x)
    np.divide(2*j1(x), x, out=f, where=x != 0)
    return DATA['C_PREF']*k**3*a**4*(PARAMS['dn']/PARAMS['n_matrix'])**2*f*f*s

def integral(fn, lo, hi):
    return quad(lambda t: float(fn(t)), lo, hi, epsabs=1e-16, epsrel=2e-11)[0]

def attenuation(lam=550., sf=.15, mode='inherited', scale=1., poisson=False):
    return 2*integral(lambda t: density(t,lam,sf,mode,scale,poisson),0,np.pi)*PARAMS['rho2d_per_nm2']/scale**2*1000

def acceptance():
    def accepted(t):
        return bool(ray_angle(np.array([t]), .1, 5.7)[0][0])
    lo, hi = 0., 1.4
    assert accepted(lo) and not accepted(hi)
    for _ in range(45):
        mid = (lo+hi)/2
        if accepted(mid): lo = mid
        else: hi = mid
    cut = (lo+hi)/2
    core = brentq(lambda t: float(ray_angle(np.array([t]),.1,5.7)[2][0])-np.pi/(180*60),0,cut,xtol=1e-14)
    return cut, core

def port(mode='inherited', sf=.15, poisson=False, cut=None, core=None):
    if cut is None: cut,core = acceptance()
    fn = lambda t: density(t,550.,sf,mode,poisson=poisson)
    tot = integral(fn,0,np.pi)
    acc = integral(fn,0,cut)/tot
    cr = integral(fn,0,core)/tot
    tau = 2*tot*PARAMS['rho2d_per_nm2']*1000
    optical_depth = tau*PARAMS['L_nm']/1000
    ballistic = float(np.exp(-optical_depth))
    retinal = ballistic+(1-ballistic)*acc
    tail = (1-ballistic)*(acc-cr)/retinal
    qs = {}
    for p in [.01,.05,.5,.95,.99]:
        t = brentq(lambda t: integral(fn,0,t)/tot-p,0,np.pi,xtol=1e-12)
        qs[str(p)] = t
    return dict(mode=mode,sigma_fraction=sf,poisson=poisson,tau_per_um=tau,
                optical_depth=optical_depth,ballistic=ballistic,accepted_fraction=acc,
                core_fraction=cr,retinal_tail_fraction=tail,theta_quantiles_rad=qs,
                angular_tail_mass_outside_aperture=1-acc,theta_cutoff_rad=cut,
                theta_core_rad=core,closure='conditional on-axis meridional single-scatter; not 3D retinal halo',
                transport_gate='THIN' if optical_depth<=.1 else 'EXTRAPOLATION',
                amplitude_status='PHENOMENOLOGICAL calibrated to transmission' if mode=='inherited' else 'derived Maxwell single cylinder; ensemble closure synthetic')

if __name__ == '__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--mode',choices=['inherited','maxwell'],default='inherited');ap.add_argument('--output',required=True)
    args=ap.parse_args()
    with Path(args.output).open('x') as f: json.dump(port(args.mode),f,indent=2,allow_nan=False)
