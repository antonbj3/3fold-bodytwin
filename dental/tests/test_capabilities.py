import copy
import importlib
import json
import math
import unittest
import numpy as np
from dental_release.capabilities import RUNNERS, MUTATIONS, CONTRACTS, check, fixture, clean, package, defs, synthetic_bite, verify_freeze
from dental_release.load import ROOT, kernel

class PortableCapabilities(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        verify_freeze()
        cls.contracts=json.loads(CONTRACTS.read_text())['profiles']
        cls.results={k:clean(f()) for k,f in RUNNERS.items()}

    def test_profiles_match_frozen_criteria(self):
        self.assertEqual(set(RUNNERS),set(self.contracts))
        for k,r in self.results.items():
            with self.subTest(profile=k):check(k,r,self.contracts[k])

    def test_every_profile_rejects_wrong_result(self):
        for k,r in self.results.items():
            bad=copy.deepcopy(r);field,value=MUTATIONS[k];bad[field]=value
            with self.subTest(profile=k),self.assertRaises(AssertionError):check(k,bad,self.contracts[k])

    def test_nonfinite_result_cannot_pass_a_numeric_gate(self):
        for k,field in [('CAD','coordinate_error_mm'),('X87','max_control_error_mm'),('X98','max_bracket_width_mm')]:
            for value in [math.inf,math.nan]:
                with self.subTest(profile=k,value=str(value)),self.assertRaises(AssertionError):
                    check(k,dict(self.results[k],**{field:value}),self.contracts[k])

    def test_exact_summary_identity_does_not_accept_a_tolerance(self):
        for k,r in self.results.items():
            if 'summary_identity_error' in r:
                self.assertEqual(r['summary_identity_error'],0)
                with self.subTest(profile=k),self.assertRaises(AssertionError):
                    check(k,dict(r,summary_identity_error=1e-16),self.contracts[k])
        for k,field in [('CAD','waist_rebound_difference_mm'),('X89','midpoint_difference_rational'),('X96','shape_gap_difference_mm'),('X97','spatial_difference_mm2'),('LAYER_NESTING','changed_area')]:
            self.assertTrue(self.results[k][field])

    def test_crown_failure_is_preserved(self):
        self.assertEqual(self.results['CROWN_PREP']['complete_original_chain'],0)
        self.assertEqual(fixture('crown_diagnosis.json')['complete_form_and_nominal_function'],0)
        self.assertEqual(fixture('crown_diagnosis.json')['requested_designs'],72)

    def test_digital_results_cannot_promote_physical_status(self):
        for k,field in [('X87','physical_margin'),('X88','target_site_risk'),('X89','physical_qualification'),('X98','physical_safety')]:
            self.assertEqual(self.results[k][field],'UNKNOWN')
            with self.subTest(profile=k),self.assertRaises(AssertionError):check(k,dict(self.results[k],**{field:'PASS'}),self.contracts[k])

class OcclusionInputContracts(unittest.TestCase):
    def test_mismatched_pose_frame_case_or_units_reject(self):
        with package('X97','occlusion_module') as m:
            for field,value in [('case_id','other'),('pose_id','other'),('frame_id','other'),('unit','um')]:
                b,c=synthetic_bite();c[field]=value
                with self.subTest(field=field),self.assertRaises(m.ContractError):m.analyze(b,c,{'propose_adjustment':False})

    def test_invalid_roles_facets_and_nonfinite_vertices_reject(self):
        with package('X97','occlusion_module') as m:
            for change in ['roles','indices','vertices']:
                b,c=synthetic_bite()
                if change=='roles':c['face_roles'][0]=3
                elif change=='indices':c['faces'][0,0]=100
                else:c['vertices_mm'][0,0]=np.nan
                with self.subTest(change=change),self.assertRaises(m.ContractError):m.analyze(b,c,{'propose_adjustment':False})

    def test_measured_total_force_requires_calibration(self):
        with package('X97','occlusion_module') as m:
            b,c=synthetic_bite();b.update(total_force_interval_N=[100,110],total_force_evidence={'kind':'MEASURED','locator':'synthetic test'})
            with self.assertRaises(m.ContractError):m.analyze(b,c,{'propose_adjustment':False})

    def test_empty_duplicate_reversed_pair_graph_reject(self):
        with package('X97','occlusion_module') as m:
            for change in ['empty','duplicate','reversed']:
                b,c=synthetic_bite()
                if change=='empty':b['tooth_pairs']=[]
                elif change=='duplicate':b['tooth_pairs'].append(copy.deepcopy(b['tooth_pairs'][0]))
                else:b['tooth_pairs'][0]['gap_interval_mm']=['1','0']
                with self.subTest(change=change),self.assertRaises(m.ContractError):m.analyze(b,c,{'propose_adjustment':False})

    def test_inputs_are_preserved_and_geometry_identity_changes(self):
        with package('X97','occlusion_module') as m:
            b,c=synthetic_bite();before={k:v.copy() for k,v in c.items() if isinstance(v,np.ndarray)};h=m.geometry_sha256(c)
            m.analyze(b,c,{'propose_adjustment':False})
            for k,v in before.items():np.testing.assert_array_equal(c[k],v)
            c['vertices_mm'][0,2]+=.01
            self.assertNotEqual(h,m.geometry_sha256(c))

class ImplantInputContracts(unittest.TestCase):
    def test_axis_and_coordinates_must_be_finite_and_unit(self):
        with package('X98','implant_safety') as m:
            mod=importlib.import_module('implant_safety.module');frame={'unit':'synthetic'}
            for axis in [[0,0,2],[0,0,0],[0,np.nan,1],[1,0]]:
                with self.subTest(axis=str(axis)),self.assertRaises(m.ContractError):mod.validate_pose(dict(frame=frame,entry_zyx_mm=[0,0,0],axis_zyx=axis),frame)
            with self.assertRaises(m.ContractError):mod.validate_pose(dict(frame=frame,entry_zyx_mm=[0,0,np.inf],axis_zyx=[0,0,1]),frame)

    def test_frame_cannot_silently_change_units_or_voxel_box_axes(self):
        with package('X98','implant_safety') as m:
            mod=importlib.import_module('implant_safety.module');frame=dict(order='zyx',units='mm',spacing_mm=[1,1,1],origin_mm=[0,0,0],direction=np.eye(3).tolist())
            mod.validate_frame(frame)
            for field,value in [('units','um'),('order','xyz'),('spacing_mm',[1,0,1]),('direction',[[0,1,0],[1,0,0],[0,0,1]])]:
                with self.subTest(field=field),self.assertRaises(m.ContractError):mod.validate_frame(dict(frame,**{field:value}))

    def test_booleans_and_nonpositive_dimensions_reject(self):
        with package('X98','implant_safety') as m:
            mod=importlib.import_module('implant_safety.module')
            for value in [True,0,-1,math.nan,math.inf,'1']:
                with self.subTest(value=str(value)),self.assertRaises(m.ContractError):mod.finite_number(value,'dimension',True)

class LaboratoryContracts(unittest.TestCase):
    def test_missing_curvature_pose_and_extrapolation_abstain(self):
        with kernel('X89','code/calibration.py') as m:
            p=dict(pose_id='synthetic',force_levels_N=[0,1],tip_displacement_lower_mm=[0,.001],tip_displacement_upper_mm=[0,.002],curvature_bound_mm_N2=None)
            self.assertEqual(m.displacement_bound(p,.5,'synthetic')['status'],'UNKNOWN_NO_NONLINEAR_ENCLOSURE')
            self.assertEqual(m.displacement_bound(dict(p,curvature_bound_mm_N2=0),2,'synthetic')['status'],'UNKNOWN_EXTRAPOLATION')
            self.assertEqual(m.displacement_bound(p,.5,'other')['status'],'UNKNOWN_POSE_MISMATCH')
            with self.assertRaises(ValueError):m.displacement_bound(dict(p,curvature_bound_mm_N2=-1),.5,'synthetic')

    def test_published_counts_cannot_be_fractional_or_negative(self):
        with kernel('X88','code/model.py') as m:
            for k,n in [(-1,10),(11,10),(1.5,10),(1,0)]:
                with self.subTest(k=k,n=n),self.assertRaises(ValueError):m.cp(k,n)

    def test_three_lab_protocols_require_identity_and_prediction_freeze(self):
        rows=fixture('lab_protocols.json')['protocols']
        self.assertEqual([r['id'] for r in rows],['R4P','R4C','R4F'])
        for row in rows:
            self.assertIn('specimen_id',row['required_record_fields'])
            self.assertIn('prediction_hash',row['required_record_fields'])
            self.assertTrue(row['freeze_before_measurement'])
            self.assertIn('UNKNOWN',row['cost'])

    def test_coefficient_correction_preserves_frozen_values_only(self):
        with kernel('X24_CBCT_HU_CALIBRATION','frozen_coefficients_match.py') as m:
            self.assertTrue(m.models_match({'synthetic':{'coefficients':[1.,2.]}},{'synthetic':{'coefficients':[1.,2.+1e-12]}}))
            self.assertFalse(m.models_match({'synthetic':{'coefficients':[1.,2.]}},{'synthetic':{'coefficients':[1.,2.1]}}))

class ReviewBindingContracts(unittest.TestCase):
    def test_all_added_entries_bind_accepted_result_and_review(self):
        from dental_release.review_scope import code_eligible
        records=json.loads((ROOT/'provenance/REVIEW_SCOPE.json').read_text())['records']
        entries=json.loads((ROOT/'demos.json').read_text())['demos']
        added=[r for r in records if r['action']=='ADD']
        scope=json.loads((ROOT/'provenance/RELEASE_SCOPE.json').read_text())
        self.assertEqual(len(added),scope['added_review_receipts'])
        for row in added:
            self.assertTrue(code_eligible(row))
            pair=dict(result_sha256=row['result_sha256'],review_sha256=row['review_sha256'])
            self.assertTrue(any((e.get('reviewed_result_sha256')==row['result_sha256'] and e.get('review_sha256')==row['review_sha256'])
                                or pair in e.get('accepted_receipts',[]) for e in entries))

    def test_referent_only_and_missing_binding_are_not_code_approval(self):
        from dental_release.review_scope import code_eligible
        row=dict(classification='ACCEPTED',action='ADD',result_sha256='a'*64,review_sha256='b'*64)
        self.assertTrue(code_eligible(row))
        for change in [dict(action='SCOPED_REFERENT_ONLY'),dict(classification='UNKNOWN'),dict(result_sha256=None),dict(review_sha256='wrong')]:
            with self.subTest(change=change):self.assertFalse(code_eligible(dict(row,**change)))

    def test_identical_entry_count_can_hide_lost_review_scope(self):
        from dental_release.review_scope import count_sufficiency
        records=json.loads((ROOT/'provenance/REVIEW_SCOPE.json').read_text())['records']
        r=count_sufficiency([r for r in records if r['action']=='ADD'][:9])
        self.assertEqual(r['identity_error'],0)
        self.assertEqual(r['downstream_eligible_counts'],[9,8])
        self.assertEqual(r['downstream_difference'],1)

if __name__=='__main__':unittest.main()
