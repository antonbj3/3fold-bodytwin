"""Six source-number regression checks. No optional-data skips: missing files fail."""
import json
import numpy as np
import pytest
from bodytwin_core.config import paths
from bodytwin_core import geometry, solver, whatif, emulator, population, nullmodels

P=paths()

def source(key,file):
    return json.loads((P[key]/file).read_text())

def test_geometry_n7c_registry_and_certificate():
    r=geometry.instantiate('z001')
    assert r.accepted
    assert len(r.instance.reg.entities)==216  # N7c/CX-GEOMCERT, exact IDs
    assert r.instance.reg.manifest_hash()  # stable identity manifest

def test_solver_a295_peak_frame():
    f,rec=solver.solve_lift('identity',frames=[75])
    assert rec[0]['kkt']<=1e-10
    peak=float(solver.hip_curve('identity',np.tile(f,(141,1)))[75])
    expected=source('solver_results','summary.json')['population']['individuals']['identity']['hip_r_N']
    assert peak==pytest.approx(expected,abs=0.05) # 0.05 N allows source formulation/rounding

def test_whatif_a369_implicit_peak():
    r=whatif.reported_derivative('identity','CCD',75)
    rows=source('whatif_results','implicit_kkt_results.json')['identity']['CCD']['rows']
    assert r['frame']==75
    assert r['implicit_N_per_deg']==pytest.approx(rows[75]['implicit_N_per_deg'],abs=1e-12)
    assert r['implicit_N_per_deg']==pytest.approx(8.351990217102463,abs=1e-9) # N/degree
    assert r['relative_error']<1e-4

def test_emulator_a303_max_relative_error():
    z=np.load(P['xf4_results']/'xf4_pred.npz')
    e=emulator.rel_error(z['R_xf4'],z['R_oracle'])
    expected=source('xf4_results','results.json')['test_error']['XF4']['max']
    assert float(e.max())==pytest.approx(expected,abs=1e-10)

def test_population_a363_band_and_n50_rule():
    s=source('n2b_results','n2b_results.json')['bands']['ens_s9p5_p2_lm_M1000']['hip_r']
    b=population.band_with_sd(2583.0,100.0)
    k=source('n50_results','results.json')['hip_r']['full_data']['kappa']
    assert b['upper_N']==pytest.approx(2583+100*k*1.6448536269514722,abs=1e-9)
    assert s['peak_p2_5_pct']==pytest.approx(-8.486082649167702,abs=1e-9)
    assert s['peak_p97_5_pct']==pytest.approx(15.877548804245535,abs=1e-9)

def test_nullmodels_b24_loso_and_n1():
    rows=source('b24_results','results.json')['facets']
    row=next(r for r in rows if r['joint']=='hip' and r['sub_group']=='gen1' and r['activity_bucket']=='Walking')
    errors=[r['signed_error_pctBW'] for r in row['error']['patient_loso']]
    rmse=float(np.sqrt(np.mean(np.square(errors))))
    assert rmse==pytest.approx(row['error']['rmse_pctBW'],abs=1e-5) # source rounded 6 decimals
    assert rmse==pytest.approx(70.374640,abs=1e-5)
    assert nullmodels.n1_knee(1000,2.2)==2200
    assert nullmodels.n1g_knee(1000,2.37)==2370
