import copy
import importlib.util
import json
from pathlib import Path
import subprocess
import tempfile
import unittest
import numpy as np
from dental_release.demo import RUNNERS, MUTATIONS, CONTRACTS, check, synthetic_response, verify_freeze
from dental_release.load import ROOT, kernel
from dental_release.paths import expand, MissingInput

spec=importlib.util.spec_from_file_location('scrub',ROOT/'tools/scrub.py')
scrub=importlib.util.module_from_spec(spec);spec.loader.exec_module(scrub)

class FrozenKernelProfiles(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        verify_freeze()
        cls.contracts=json.loads(CONTRACTS.read_text())['profiles']
        cls.results={name:fn() for name,fn in RUNNERS.items()}

    def test_all_declared_profiles_match_frozen_contract(self):
        for name,result in self.results.items():
            with self.subTest(demo=name):check(name,result,self.contracts[name])

    def test_every_profile_rejects_changed_numeric_output(self):
        for name,result in self.results.items():
            changed=copy.deepcopy(result);key,value=MUTATIONS[name];changed[key]=value
            with self.subTest(demo=name),self.assertRaises(AssertionError):check(name,changed,self.contracts[name])

    def test_summaries_preserve_machine_identity_but_change_consumers(self):
        for name in ['X14','X59','X82']:
            with self.subTest(demo=name):self.assertEqual(self.results[name]['summary_identity_error'],0)
        self.assertEqual(self.results['X63']['closed_gap_identity_error'],0)
        self.assertEqual(self.results['X59']['classification_difference'],4)
        self.assertGreater(self.results['X82']['downstream_difference_N'],.49)

    def test_preserves_published_negative_regression_result(self):
        self.assertEqual(self.results['X59']['regression_gate_Dentsply_Sirona'],'FAIL')
        self.assertEqual(self.results['X59']['individual_prediction'],'UNKNOWN')

    def test_colour_result_requires_finite_enclosure(self):
        nominal=self.results['X74']
        for update in [dict(deltaE=500.),dict(deltaE=np.inf),dict(deltaE=np.nan),
                       dict(interval=[-1.,1.]),dict(interval=[6.,5.]),
                       dict(interval=[0.,np.inf]),dict(interval=[0.,9.,10.])]:
            with self.subTest(update=str(update)),self.assertRaises(AssertionError):
                check('X74',dict(nominal,**update),self.contracts['X74'])

class ForceContracts(unittest.TestCase):
    def test_bad_units_and_identity_reject(self):
        with kernel('X82','code/height_force.py') as m:
            s=synthetic_response();edit=m.make_edit(s,0,.02)
            for key,value in [('unit','um'),('case','different'),('geometry_sha256','different')]:
                with self.subTest(key=key),self.assertRaises(ValueError):m.predict(s,dict(edit,**{key:value}))

    def test_uncertainty_validation_precedes_unknown(self):
        with kernel('X82','code/height_force.py') as m:
            s=synthetic_response()
            for key in ['spectral_error_N_per_mm','baseline_error_N']:
                for value in [-1.,np.nan,np.inf,True]:
                    bad=dict(s,complete_active_calibration=False,**{key:value})
                    with self.subTest(key=key,value=str(value)),self.assertRaises(ValueError):m.predict(bad,m.make_edit(bad,0,.02))

    def test_missing_calibration_and_new_load_abstain(self):
        with kernel('X82','code/height_force.py') as m:
            s=synthetic_response();edit=m.make_edit(s,0,.02)
            self.assertEqual(m.predict(dict(s,complete_active_calibration=False),edit)['status'],'UNKNOWN')
            self.assertEqual(m.predict(s,dict(edit,load_N=110))['status'],'UNKNOWN')

    def test_height_domain_rejects(self):
        with kernel('X82','code/height_force.py') as m:
            s=synthetic_response()
            for height in [.051,np.nan,np.inf]:
                with self.subTest(height=str(height)),self.assertRaises(ValueError):m.predict(s,m.make_edit(s,0,height))

    def test_position_rate_history_closure_required(self):
        with kernel('X82','code/height_force.py') as m:
            self.assertEqual(m.project_reaction(1,2,.005,.01,0,0)['status'],'UNKNOWN')
            for err in [-1,np.nan,np.inf]:
                with self.subTest(error=str(err)),self.assertRaises(ValueError):m.project_reaction(1,2,.005,.01,err,0,True,True,True)

class MeasurementContracts(unittest.TestCase):
    def test_missing_local_reference_and_nonreciprocal_calibration(self):
        with kernel('X56','code/motion_port.py') as m:
            z=np.array([-1.,1.]);H=np.column_stack([np.ones(2),z]);W=np.eye(2);C=np.array([[.001,.0001],[.0002,.001]])
            self.assertEqual(m.identify(z,W,H@C,np.zeros((2,2)))['status'],'REJECT_NONRECIPROCAL_CALIBRATION')
            self.assertEqual(m.identify([0,0],W,H@C,np.zeros((2,2)))['status'],'UNKNOWN_RANK_DEFICIENT_MEASUREMENT')

    def test_lab_observation_operator_is_explicit(self):
        with kernel('X85','code/lab_compare.py') as m:
            f=dict(case='synthetic',geometry_sha256='synthetic',frame='synthetic',force_interval_N=[[7,9]])
            y=dict(case='synthetic',geometry_sha256='synthetic',frame='synthetic',unit='percent',matched_load_rate_history=True,reaction_N=[8],reaction_error_N=[.1])
            self.assertEqual(m.compare(f,y)['status'],'UNKNOWN')
            self.assertEqual(m.compare(f,dict(y,unit='N',matched_load_rate_history=False))['status'],'UNKNOWN')
            with self.assertRaises(ValueError):m.compare(f,dict(y,unit='N',reaction_error_N=[-.1]))

    def test_frozen_force_intervals_reject_malformed_endpoints(self):
        with kernel('X85','code/lab_compare.py') as m:
            frozen=dict(case='synthetic',geometry_sha256='synthetic',frame='synthetic',force_interval_N=[[7.,9.]])
            measured=dict(case='synthetic',geometry_sha256='synthetic',frame='synthetic',unit='N',matched_load_rate_history=True,reaction_N=[8.],reaction_error_N=[.1])
            self.assertEqual(m.compare(frozen,measured)['status'],'PASS_FROZEN_OBSERVATION_GATE')
            malformed=[[[7.,9.,999.]],[[-np.inf,np.inf]],[[np.nan,9.]],[[9.,7.]],[],[7.,9.],
                       8.,[[7.]],[[7.,9.],[7.]],[[[7.,9.]]],[[7.,np.inf]],[[np.nan,np.nan]]]
            for intervals in malformed:
                with self.subTest(intervals=str(intervals)),self.assertRaises(ValueError):
                    m.compare(dict(frozen,force_interval_N=intervals),measured)

    def test_fixed_plan_does_not_turn_unknown_shape_into_exact_bound(self):
        with kernel('X51','code/acceptance.py') as m:
            with self.assertRaises(ValueError):m.known_shape_lower([100,110,120],3)
            for n,k in [(0,0),(2,3),(2,-1)]:
                with self.subTest(n=n,k=k),self.assertRaises(ValueError):m.upper(n,k)

    def test_relocation_requires_external_inputs(self):
        self.assertEqual(expand('@DENTAL_IMPLEMENTATIONS@/X82'),str(ROOT/'implementations/X82'))
        with self.assertRaises(MissingInput):expand('@DENTAL_MISSING_TEST_ROOT@/input')

class ScrubFaults(unittest.TestCase):
    def probe(self,name,content,reason):
        with tempfile.TemporaryDirectory(prefix='dental-scrub-') as temp:
            p=Path(temp,name);p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(content)
            r=scrub.scan(temp,verify_manifest=False,inspect_git=False)
            self.assertEqual(r['status'],'FAIL')
            self.assertIn(reason,[x['reason'] for x in r['findings']])

    def test_every_restricted_term_rejects_case_insensitively(self):
        for token in scrub.PATTERNS:
            with self.subTest(pattern_index=scrub.PATTERNS.index(token)):self.probe('README.md',token.swapcase().encode(),'restricted_vocabulary')

    def test_private_path_rejects(self):
        path=bytes.fromhex('2f686f6d652f616e746f6e2f70726976617465')
        self.probe('README.md',path,'private_absolute_path')

    def test_dataset_and_derived_mesh_payloads_reject(self):
        for name in ['patient.stl','case.npz','archive.zip','weights.pt','records.csv']:
            with self.subTest(file=name):self.probe(name,b'synthetic adversarial payload','unapproved_payload_extension')

    def test_large_payload_rejects(self):
        self.probe('payload.txt',b'0'*(scrub.MAX_BYTES+1),'file_over_10_MB')

    def test_binary_payload_rejects(self):
        self.probe('disguised.txt',b'\x00\xff','binary_payload')

    def test_personal_metadata_rejects(self):
        key=bytes.fromhex('70617469656e745f6e616d65').decode()
        self.probe('record.json',json.dumps({key:'synthetic injected record'}).encode(),'personal_metadata')

    def test_email_rejects(self):
        self.probe('README.md',('example'+'@'+'invalid.test').encode(),'embedded_email')

    def test_personal_record_literal_in_source_rejects(self):
        key=bytes.fromhex('70617469656e745f6e616d65').decode()
        self.probe('record.py',('record='+repr({key:'injected record'})).encode(),'personal_metadata_literal')

    def test_symlink_rejects(self):
        with tempfile.TemporaryDirectory() as temp:
            Path(temp,'link.py').symlink_to('/tmp/nonexistent-dental-test')
            self.assertIn('symlink',[x['reason'] for x in scrub.scan(temp,False,False)['findings']])

    def test_unlisted_json_cannot_be_smuggled_in(self):
        with tempfile.TemporaryDirectory() as temp:
            Path(temp,'RELEASE_FILES.json').write_text(json.dumps({'files':{}}))
            Path(temp,'new_record.json').write_text('{}')
            self.assertIn('unlisted_or_changed_release_file',[x['reason'] for x in scrub.scan(temp,True,False)['findings']])

    def test_large_numeric_table_rejects(self):
        self.probe('constants.py',('data='+repr(list(range(130)))).encode(),'large_embedded_numeric_payload')

    def test_git_history_retains_removed_fault(self):
        with tempfile.TemporaryDirectory() as temp:
            def git(*args):return subprocess.run(['git','-C',temp,*args],check=True,capture_output=True)
            git('init','-q');git('config','user.name','Research test');git('config','user.email','test')
            p=Path(temp,'README.md');p.write_text(scrub.PATTERNS[0]);git('add','README.md');git('commit','-qm','fixture')
            p.write_text('safe');r=scrub.scan(temp,False,True)
            self.assertTrue(any(x.get('commit') and x['reason']=='restricted_vocabulary' for x in r['findings']))

    def test_removed_personal_literal_in_history_rejects(self):
        with tempfile.TemporaryDirectory() as temp:
            def git(*args):return subprocess.run(['git','-C',temp,*args],check=True,capture_output=True)
            git('init','-q');git('config','user.name','Research test');git('config','user.email','test')
            p=Path(temp,'record.py')
            p.write_bytes(bytes.fromhex('7265636f72643d7b2770617469656e745f6e616d65273a2773796e74686574696320696e6a6563746564207265636f7264277d0a'))
            git('add','record.py');git('commit','-qm','synthetic fixture')
            p.unlink();git('add','-A');git('commit','-qm','remove fixture')
            result=scrub.scan(temp,False,True)
            self.assertTrue(any(x.get('commit') and x['reason']=='personal_metadata_literal' for x in result['findings']))

    def test_generic_attribution_rejects(self):
        self.probe('README.md',bytes.fromhex('47656e65726174656420627920436861744750540a'),'restricted_vocabulary')

    def test_git_identity_email_is_scanned_and_public_noreply_allowed(self):
        with tempfile.TemporaryDirectory() as temp:
            def git(*args):return subprocess.run(['git','-C',temp,*args],check=True,capture_output=True)
            git('init','-q');git('config','user.name','Research test')
            # An invented private-looking address; no personal contact data.
            git('config','user.email',bytes.fromhex('73796e746865746963406578616d706c652e74657374').decode())
            Path(temp,'README.md').write_text('safe\n')
            git('add','README.md');git('commit','-qm','synthetic fixture')
            result=scrub.scan(temp,False,True)
            self.assertTrue(any(x['reason']=='embedded_email' and x.get('identity_field')=='author_email' for x in result['findings']))
            git('config','user.email',scrub.PUBLIC_IDENTITY_EMAIL)
            git('commit','--amend','--no-edit','--reset-author')
            self.assertEqual(scrub.scan(temp,False,True)['status'],'PASS')

    def test_license_and_notice_are_permitted_plain_text(self):
        with tempfile.TemporaryDirectory() as temp:
            for name in ['LICENSE','NOTICE']:
                Path(temp,name).write_text('License and attribution text\n')
            self.assertEqual(scrub.scan(temp,False,False)['status'],'PASS')

    def test_commit_message_is_checked(self):
        with tempfile.TemporaryDirectory() as temp:
            def git(*args):return subprocess.run(['git','-C',temp,*args],check=True,capture_output=True)
            git('init','-q');git('config','user.name','Research test');git('config','user.email','test')
            Path(temp,'README.md').write_text('safe');git('add','README.md');git('commit','-qm',scrub.PATTERNS[-1])
            r=scrub.scan(temp,False,True)
            self.assertTrue(any(x['file'].startswith('commit ') and x['reason']=='restricted_vocabulary' for x in r['findings']))

if __name__=='__main__':unittest.main()
