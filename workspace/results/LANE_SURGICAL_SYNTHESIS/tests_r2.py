"""Regressions for diagnosis and acquisition algebra, plus unchanged R1 suite."""
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:
    os.environ[k]='2'
import argparse, copy, json, unittest
from pathlib import Path
import numpy as np
import tests_r1
from diagnose_r2 import ROOT, load_reference, factors, score
from plan_measurements_r2 import contraction, GROUPS
from surgical_chain.chain import run
from surgical_chain.scenario import scenario_config
from nonlinear_information_r2 import kernel_predict

class AttributionRegression(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        with np.load(ROOT/'r1/nominal_complete_v4/fields.npz') as z:cls.z={k:z[k].copy() for k in z.files}
        cls.ref=load_reference()
        cls.report=json.loads((ROOT/'r2/attribution_v1/ATTRIBUTION.json').read_text())
        cls.rows={x['case']:x for x in cls.report['rows']}
    def test_baseline_preserved(self):
        self.assertAlmostEqual(score(self.z['days'],self.z['strength_pct'])['RMSE_pp'],43.703695914540624,places=9)
    def test_isotropic_factorization(self):
        A,m=factors(self.z)
        np.testing.assert_allclose(75*A*m,self.z['strength_pct'],atol=1e-11,rtol=0)
    def test_same_observation_contract_reproduces_R3(self):
        t=self.z['days'];q=75*np.interp(t,self.ref['days'],self.ref['C'])*np.interp(t,self.ref['days'],self.ref['HP_wound_reference_fraction'])
        self.assertAlmostEqual(score(t,q)['RMSE_pp'],8.306367094671637,places=9)
    def test_abundance_dominates_single_swap(self):
        A,m=factors(self.z);t=self.z['days']
        a=score(t,75*np.interp(t,self.ref['days'],self.ref['C'])*m)['RMSE_pp']
        b=score(t,75*A*np.interp(t,self.ref['days'],self.ref['HP_wound_reference_fraction']))['RMSE_pp']
        self.assertAlmostEqual(a,13.000632079497242,places=7)
        self.assertAlmostEqual(b,39.45678148933305,places=7)
        self.assertLess(a,b)
    def test_shapley_accounts_for_interaction(self):
        sh=self.report['two_factor_Shapley_RMSE_reduction_pp']
        self.assertAlmostEqual(sum(sh.values()),43.703695914540624-8.306367094671637,places=9)
        self.assertGreater(sh['abundance'],6*sh['maturity'])
    def test_rate_substitution_exposes_birth_saturation(self):
        self.assertGreater(self.rows['upstream_k_deposit']['RMSE_pp'],self.rows['baseline']['RMSE_pp'])
        self.assertAlmostEqual(self.rows['upstream_k_deposit']['RMSE_pp'],59.82772921327508,places=7)
    def test_no_new_empirical_admission_or_same_chemistry_claim(self):
        self.assertFalse(self.report['scientific_admission'])
        self.assertFalse(self.report['lineage_correction']['same_chemistry_input_claim'])
        self.assertFalse(self.report['lineage_correction']['R3_C_is_measured_collagen'])
    def test_legacy_maturity_has_no_UIM_consumer(self):
        self.assertAlmostEqual(self.rows['upstream_k_mature_legacy']['RMSE_pp'],43.703695914540624,places=9)
    def test_actual_chemistry_intervention_cannot_change_strength(self):
        cfg=scenario_config();a=run(cfg,mode='scenario')
        changed=copy.deepcopy(cfg);changed['parameters']['chemistry_rate']['value']*=2
        b=run(changed,mode='scenario')
        np.testing.assert_array_equal(a['strength_time']['value']['pct'],b['strength_time']['value']['pct'])
        self.assertNotEqual(a['chemistry_D_H_A']['value']['D_H_A'][1][-1],b['chemistry_D_H_A']['value']['D_H_A'][1][-1])

class NoisyMeasurementAlgebra(unittest.TestCase):
    def test_scalar_exact_gaussian_update(self):
        x=np.array([-1.,0.,1.]);v=x.var(ddof=1)
        d=contraction(x[:,None],x[:,None],[2.],[0.])
        self.assertAlmostEqual(d['conditional_linear_SD'][0]**2,v*4/(v+4),places=12)
    def test_noisy_observation_never_increases_conditional_moment_variance(self):
        rng=np.random.default_rng(6122);z=rng.normal(size=(200,3));y=np.column_stack([z[:,0]+z[:,1],z[:,2]**2])
        q=contraction(y,z,[.1,.2,.3],[.05,.05,.05])
        self.assertTrue(np.all(np.array(q['conditional_linear_SD'])<=q['prior_SD']))
    def test_shared_calibration_cannot_be_counted_twice(self):
        x=np.linspace(-1,1,100);y=x[:,None];z=np.column_stack([x,x])
        independent=contraction(y,z,[.1,.1],[0,0])
        common=contraction(y,z,[.1,.1],[.5,.5])
        self.assertGreater(common['conditional_linear_SD'][0],independent['conditional_linear_SD'][0])
    def test_disconnected_measurement_has_zero_implemented_gain(self):
        x=np.arange(50.)[:,None]
        q=contraction(x,x,[.1],[0.],disconnected=True)
        self.assertEqual(q['variance_reduction_fraction'],[0.])
    def test_unresolved_floor_remains_after_precise_measurement(self):
        x=np.arange(50.)[:,None]
        q=contraction(x,x,[1e-6],[0.],floor=[20.])
        self.assertGreaterEqual(q['conditional_linear_SD'][0],20.)
    def test_all_six_R1_rows_and_cost_ranges_present(self):
        self.assertEqual(set(GROUPS),{'M1','M2','M3','M4','M5','M6'})
        d=json.loads((ROOT/'r2/measurement_v1/MEASUREMENT_VALUE.json').read_text())
        self.assertEqual(d['sample_count'],256)
        self.assertTrue(all(v is None for v in d['native_outcomes'].values()))
        for key in ['M4','M6']:
            q=next(x for x in d['variants']['n256_noise1.0_floorFalse'] if x['id']==key)
            self.assertEqual(q['variance_reduction_fraction'],[0.]*6)
    def test_refinement_and_nonlinear_conservation(self):
        d=json.loads((ROOT/'r2/measurement_v1/REFINEMENT.json').read_text())
        self.assertEqual(d['gate'],'PASS')
        v=json.loads((ROOT/'r2/measurement_v1/MEASUREMENT_VALUE.json').read_text())
        self.assertLess(v['max_mass_balance'],1e-10)
        self.assertLess(v['max_oxygen_conservation'],1e-7)

class NonlinearInformationRegression(unittest.TestCase):
    def test_positive_kernel_preserves_observable_bounds(self):
        pred,w=kernel_predict(np.array([[0.,1.,5.],[5.,1.,0.]]),np.array([[0.,120.],[50.,120.],[100.,120.]]),.5)
        np.testing.assert_allclose(w.sum(1),1.,atol=1e-14)
        self.assertTrue(np.all((pred[:,0]>=0)&(pred[:,0]<=100)))
        np.testing.assert_allclose(pred[:,1],120.,atol=1e-12)
    def test_noise_and_fresh_validation_are_separate(self):
        root=ROOT/'r2/nonlinear_information_v1'
        t=json.loads((root/'TRAINING_V2.json').read_text());d=json.loads((root/'DESIGN.json').read_text())
        self.assertTrue(t['trained_without_new64_outputs'])
        self.assertNotEqual(t['noise_seed'],d['measurement_noise_seed'])
        self.assertEqual(d['new_samples'],64)
        v=json.loads((root/'VALIDATION.json').read_text())
        self.assertEqual(v['training_samples'],256);self.assertEqual(v['new_validation_samples'],64)
        self.assertIsNone(v['native_uncertainty'])
    def test_original_linear_hypoxia_failure_is_retained(self):
        d=json.loads((ROOT/'r2/measurement_v1/MEASUREMENT_VALUE.json').read_text())
        q=next(x for x in d['variants']['n256_noise1.0_floorFalse'] if x['id']=='M2')
        self.assertLess(q['heldout_train128_test128_MSE_reduction_fraction'][5],-.5)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    a.out.mkdir(exist_ok=False,parents=True);tests_r1.TEST_OUT=a.out
    loader=unittest.TestLoader();suite=loader.loadTestsFromModule(tests_r1)
    suite.addTests(loader.loadTestsFromTestCase(AttributionRegression));suite.addTests(loader.loadTestsFromTestCase(NoisyMeasurementAlgebra))
    suite.addTests(loader.loadTestsFromTestCase(NonlinearInformationRegression))
    result=unittest.TextTestRunner(verbosity=2).run(suite)
    if not result.wasSuccessful():raise SystemExit(1)
