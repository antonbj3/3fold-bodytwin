import json
from pathlib import Path
import numpy as np
import pytest
from bodytwin_core.determined import joint_load,knee_load,measurement_plan

ROOT=Path(__file__).resolve().parents[2]

def test_joint_load_affine_contacts_and_infeasibility():
    system=dict(A=np.array([[[1.]]]),b=np.array([[5.]]),cj=np.array([[1.]]),
                C0=np.array([2.]),Mkx=np.array([[0.]]),mx=np.array([0.]),d=1.,F0=np.array([10.]))
    got=joint_load(system)
    np.testing.assert_allclose(got['total'],[[7.,7.]])
    np.testing.assert_allclose(got['medial'],[[3.5,3.5]])
    np.testing.assert_allclose(got['lateral'],[[3.5,3.5]])
    system['F0']=np.array([4.])
    assert np.isnan(joint_load(system)['total']).all()
    with pytest.raises(ValueError,match='infeasible'):
        joint_load(system,on_infeasible='raise')

def test_law_envelope_and_missing_set():
    trial=dict(grf=np.array([[0.,0.,700.]]),knee_moment_Nm=np.array([50.]),
               ankle_moment_Nm=np.array([20.]),knee_flex_deg=np.array([45.]),
               measured_N=np.array([1900.]))
    got=knee_load(trial)
    assert got['law_N'][0]==pytest.approx(700+50/.05+.29*20/.045)
    assert got['law_band_N'][0,0]==pytest.approx(700+50/.06+.20*20/.06)
    assert got['law_band_N'][0,1]==pytest.approx(700+50/.03+.40*20/.035)
    assert got['law_contains_measured'][0]
    assert got['set_N'] is None and not got['measurement']['current_data_enough']

def test_frozen_gait_cohort_and_person_unit():
    result=json.loads((ROOT/'results/CX-DETERMINED/results.json').read_text())
    assert result['all']['trials']==108
    assert result['all']['frames']==11521
    assert result['all']['set_frames']==sum(v['set_frames'] for v in result['by_person'].values())
    assert len(result['by_person'])==4
    assert all(v['law_criterion']==(v['law_fraction']>=.5) for v in result['by_person'].values())
    assert measurement_plan()['expected_95pct_width_N']['stance']['cost_46_triple']==[452,403]
    new=json.loads((ROOT/'results/CX-DETERMINED/new_activity_results.json').read_text())
    assert new['trials']==11 and new['frames']==9401
    assert all('lunge' not in row['trial'] for row in new['rows'])
    assert sum(row['law_inside'] for row in new['rows'])==new['law_inside']
