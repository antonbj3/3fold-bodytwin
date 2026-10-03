"""Retain all seventy old tests; exercise shared-reference physics and output gates."""
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[k]='2'
import argparse,copy,hashlib,json,time,unittest
from pathlib import Path
import numpy as np
from scipy.linalg import expm
import tests_r1,tests_r2
from surgical_chain.shared import mark_rhs,cohort_rhs,chemical_rhs,late_response,VERSION
from surgical_chain.chain import resume,run,write
from surgical_chain.scenario import scenario_config,evidence_config
from surgical_chain.vendor import response_r1 as r
from diagnose_r2 import score
from experiment_r3 import all_r3_config
ROOT=Path(__file__).resolve().parent
DATA=ROOT/'r3/experiment_v2'

class SharedReferenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.report=json.loads((DATA/'EXPERIMENT.json').read_text());cls.rows={v['case']:v for v in cls.report['rows']}
        with np.load(DATA/'all_R3/fields.npz') as z:cls.fields={k:z[k].copy() for k in z.files}
        with np.load(DATA/'all_R3/shared_reference.npz') as z:cls.shared={k:z[k].copy() for k in z.files}
    def test_integrated_strength_RMSE_regression(self):
        actual=score(self.fields['days'],self.fields['strength_pct'])['RMSE_pp']
        self.assertAlmostEqual(actual,21.811181366256044,places=8)
        self.assertLess(actual,43.703695914540624)
    def test_component_target_failure_remains_explicit(self):
        self.assertEqual(self.report['autonomous_component_reproduction_gate'],'FAIL')
        self.assertGreater(self.rows['all_R3']['score']['RMSE_pp'],8.406367094671637)
        self.assertFalse(self.report['external_HP_fed_into_candidate']);self.assertFalse(self.report['late_strength_fit'])
    def test_pointwise_inventory_and_removal_ledgers_all_cases(self):
        for row in self.rows.values():
            q=row['shared_reference']
            self.assertLess(q['max_pointwise_trace_error'],1e-10)
            self.assertLess(q['max_collagen_ledger_error'],1e-10)
            self.assertGreaterEqual(q['min_mark_eigenvalue'],-1e-10)
            self.assertEqual(q['additional_UI_death'],0.)
        np.testing.assert_allclose(self.shared['UIM'].sum(1),self.shared['C'],rtol=0,atol=1e-12)
    def test_birth_reference_single_capacity(self):
        self.assertAlmostEqual(self.rows['all_R3']['C90'],.9983202811882146,places=9)
        self.assertGreater(sum(self.rows['all_R3']['UIM90']),.99)
    def test_chemical_site_ledger_does_not_create_collagen(self):
        z=self.shared['chemical_marks'];s=.102
        np.testing.assert_allclose(z[:,0]+2*z[:,1]+z[:,2],s*self.shared['C'],rtol=0,atol=1e-12)
        self.assertLess(self.rows['all_R3']['shared_reference']['max_site_ledger_error'],1e-10)
    def test_internal_maturation_cancels_without_removal(self):
        Q=np.zeros((3,2,2,3));Q[:,0,0]=np.array([[.1],[.2],[.3]]);Q[:,1,1]=.5*Q[:,0,0]
        dy=mark_rhs(Q,np.array([.1,.2,.3]),np.zeros(3),np.ones(3),.1,.05)
        np.testing.assert_allclose(np.trace(dy.sum(0),axis1=0,axis2=1),[.1,.2,.3],atol=1e-15)
    def test_anisotropic_proportional_loss_preserves_reference(self):
        v=np.array([np.cos(.7),np.sin(.7)]);Q=np.repeat((np.outer(v,v)*np.array([.2,.3,.5])[:,None,None])[:,:,:,None],2,axis=3)
        b=np.array([.04,.05]);loss=np.array([.03,.08]);h=np.ones(2)
        dy=mark_rhs(Q,b,loss,h,.1,.05)
        np.testing.assert_allclose(np.trace(dy.sum(0),axis1=0,axis2=1),b-loss,atol=1e-15)
    def test_no_birth_no_loss_exact_sequential_conversion(self):
        # Independent analytic matrix-exponential limit against midpoint Q.
        A=np.array([[-.1,0,0],[.1,-.05,0],[0,.05,0]])
        q=np.zeros((3,2,2,1));q[0,0,0]=q[0,1,1]=.5
        dt=.001
        for j in range(10000):
            mid=q+.5*dt*mark_rhs(q,np.zeros(1),np.zeros(1),np.ones(1),.1,.05)
            q+=dt*mark_rhs(mid,np.zeros(1),np.zeros(1),np.ones(1),.1,.05)
        expected=expm(10*A)@np.array([1.,0,0])
        np.testing.assert_allclose(np.trace(q,axis1=1,axis2=2).ravel(),expected,atol=2e-9,rtol=0)
    def test_pure_loss_nonzero_initial_mass(self):
        q=np.zeros((3,2,2,1));q[2,0,0]=.6;q[2,1,1]=.4
        dt=.01;loss=np.array([.08])
        for j in range(1000):
            mid=q+.5*dt*mark_rhs(q,np.zeros(1),loss,np.ones(1),0,0)
            q+=dt*mark_rhs(mid,np.zeros(1),loss,np.ones(1),0,0)
        self.assertAlmostEqual(float(np.trace(q.sum(0),axis1=0,axis2=1)[0]),np.exp(-.8),delta=4e-8)
    def test_strongest_cohort_control_matches(self):
        row=self.rows['cohort_control']
        self.assertLess(row['shared_reference']['max_cohort_tensor_error'],1e-12)
        self.assertAlmostEqual(row['score']['RMSE_pp'],self.rows['all_R3']['score']['RMSE_pp'],places=10)
        self.assertEqual(self.report['strongest_control_gate'],'TIE')
    def test_pathway_blockades_consume_same_C(self):
        for key in ['pathway_block_k1','pathway_block_k2']:
            self.assertEqual(self.rows[key]['strength90_pct'],0.)
            self.assertAlmostEqual(self.rows[key]['C90'],self.rows['all_R3']['C90'],places=12)
    def test_noformation_and_noremove(self):
        with np.load(DATA/'noformation.npz') as z:
            self.assertEqual(float(abs(z['C']).max()),0);self.assertEqual(float(abs(z['ledger']).max()),0)
        with np.load(DATA/'noremove.npz') as z:
            self.assertEqual(float(abs(z['ledger'][:,1]).max()),0.)
            self.assertGreaterEqual(float(np.diff(z['C']).min()),-1e-12)
    def test_refinement_is_below_error_budget(self):
        self.assertLess(max(self.report['refinement_max_strength_pp'].values()),.01)
    def test_checkpoint_replays_new_law_and_chemistry(self):
        out=TEST_OUT/'shared_restart';l=resume(DATA/'all_R3/checkpoint21.npz',out)
        mask=self.fields['days']>21+1e-8
        np.testing.assert_allclose(l['strength_pct'],self.fields['strength_pct'][mask],rtol=0,atol=1e-9)
        np.testing.assert_allclose(l['chemical_marks_series'],self.shared['chemical_marks'][mask],rtol=0,atol=1e-12)
        self.assertTrue((out/'shared_reference_resumed.npz').exists())
    def test_old_checkpoint_rejected_for_shared_reference(self):
        cfg=all_r3_config();pa={k:v['value'] for k,v in cfg['parameters'].items()}
        meta=json.loads((DATA/'all_R3/checkpoint21.json').read_text());p=r.Params(**meta['parameters']);g=r.Grid(p)
        with np.load(ROOT/'r1/nominal_complete_v4/checkpoint21.npz') as z:ck={k:z[k].copy() for k in z.files}
        with self.assertRaisesRegex(ValueError,'Legacy collagen reference'):late_response(g,p,{'hazard':ck['hazard'],'state':ck['early_state']},pa,90.,.1,checkpoint=ck)
    def test_corrupted_checkpoint_reference_rejected(self):
        cfg=all_r3_config();pa={k:v['value'] for k,v in cfg['parameters'].items()}
        meta=json.loads((DATA/'all_R3/checkpoint21.json').read_text());p=r.Params(**meta['parameters']);g=r.Grid(p)
        with np.load(DATA/'all_R3/checkpoint21.npz') as z:ck={k:z[k].copy() for k in z.files}
        ck['Q']*=.9
        with self.assertRaisesRegex(ValueError,'ledger mismatch'):late_response(g,p,{'hazard':ck['hazard'],'state':ck['early_state']},pa,90.,.1,checkpoint=ck)
    def test_chemical_budget_prevents_HP_target_even_with_fast_conversion(self):
        q=self.rows['all_R3']['shared_reference'];bound=q['chemical_maturity_upper_bound']
        self.assertAlmostEqual(bound,.28625,places=12)
        self.assertLess(float(self.shared['chemical_strength_pct'].max()),75*bound+1e-10)
        self.assertAlmostEqual(self.rows['all_R3']['chemical_score']['RMSE_pp'],44.34045698668247,places=8)
    def test_stoichiometric_flux_with_new_birth_and_removal(self):
        pa={k:v['value'] for k,v in scenario_config()['parameters'].items()}
        z=np.array([[.04,.06],[.01,.005],[.02,.01]])
        dy=chemical_rhs(z,np.array([.2,.1]),np.array([.03,.04]),np.ones(2),pa,.102)
        np.testing.assert_allclose(dy[0]+2*dy[1]+dy[2],.102*np.array([.2,.1])-np.array([.03,.04])*(z[0]+2*z[1]+z[2]),atol=1e-15)
    def test_native_unknown_and_prereg_frozen(self):
        d=run(evidence_config(),inventory='shared')
        self.assertIsNone(d['native_strength']['value']);self.assertFalse(self.report['scientific_admission'])
        self.assertEqual(hashlib.sha256((ROOT/'PREREG_R3.json').read_bytes()).hexdigest(),(ROOT/'PREREG_R3.sha256').read_text().strip())
    def test_unused_extra_loss_does_not_get_added(self):
        q=self.rows['all_R3']['shared_reference']
        self.assertEqual(q['unused_legacy_ports'],['birth_scale','turnover_U_I'])

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--out',type=Path,required=True);args=a.parse_args();TEST_OUT=args.out
    if not TEST_OUT.resolve().is_relative_to(ROOT):raise ValueError('Only lane writes')
    TEST_OUT.mkdir(exist_ok=False,parents=True);tests_r1.TEST_OUT=TEST_OUT
    loader=unittest.TestLoader();suite=loader.loadTestsFromModule(tests_r1)
    for c in [tests_r2.AttributionRegression,tests_r2.NoisyMeasurementAlgebra,tests_r2.NonlinearInformationRegression,SharedReferenceTests]:suite.addTests(loader.loadTestsFromTestCase(c))
    st=time.perf_counter();result=unittest.TextTestRunner(verbosity=2).run(suite)
    write(TEST_OUT/'verification.json',dict(review_state='PENDING_INDEPENDENT_REVIEW',scientific_admission=False,count=result.testsRun,passed=result.testsRun-len(result.failures)-len(result.errors),gate='PASS' if result.wasSuccessful() else 'FAIL',failures=[{'test':str(t),'trace':s} for t,s in result.failures+result.errors],wall_s=time.perf_counter()-st,scope='Mathematical/source/shared-reference regressions; empirical prediction UNKNOWN; autonomous8.3 target FAIL retained'))
    raise SystemExit(0 if result.wasSuccessful() else 1)
