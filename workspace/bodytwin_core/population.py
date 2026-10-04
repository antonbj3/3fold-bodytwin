"""A363 batched N2b solver and BT-N50 calibrated hip band."""
import sys
from .vendor.population import n2b_solver_np
sys.modules.setdefault('n2b_solver_np',n2b_solver_np)
from .vendor.population.n2b_core import BE, Problem
from .config import paths
import json

def hip_band(peak_N, config=None):
    p=paths(config)
    result=json.loads((p['n50_results']/'results.json').read_text())
    k=result['hip_r']['full_data']['kappa']
    # N50 calibrates a band from a person's posterior-draw SD, which must be supplied.
    return {'peak_N':float(peak_N),'kappa':k,'z90':result['z90'],
            'rule':'peak ± kappa*z90*posterior_sd_N', 'posterior_sd_N':None}

def band_with_sd(peak_N, posterior_sd_N, config=None):
    band=hip_band(peak_N,config)
    radius=band['kappa']*band['z90']*posterior_sd_N
    band.update(posterior_sd_N=float(posterior_sd_N),lower_N=float(peak_N-radius),upper_N=float(peak_N+radius))
    return band
