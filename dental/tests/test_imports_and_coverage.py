import copy
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from dental_release.coverage import validate
from dental_release.load import ROOT

spec=importlib.util.spec_from_file_location('release_import_check',ROOT/'tools/import_check.py')
imports=importlib.util.module_from_spec(spec)
spec.loader.exec_module(imports)

class AllSourceImports(unittest.TestCase):
    def test_every_implementation_file_with_declared_external_prerequisites(self):
        frozen=json.loads((ROOT/'FROZEN_IMPORTS.json').read_text())
        for path, expected_sha in frozen['files'].items():
            self.assertEqual(hashlib.sha256((ROOT/path).read_bytes()).hexdigest(),expected_sha,path)
        report=imports.check()
        expected=sorted(str(p.relative_to(ROOT)) for p in (ROOT/'implementations').rglob('*.py') if '__pycache__' not in p.parts)
        self.assertEqual(sorted(r['file'] for r in report['rows']),expected)
        self.assertEqual(report['entries'],129)
        self.assertEqual(report['status'],'PASS',[r for r in report['rows'] if r['status']=='FAIL'])
        self.assertEqual(report['failures'],0)

    def probe(self, source, dependency=False):
        with tempfile.TemporaryDirectory(prefix='dental-import-fault-') as temp:
            root=Path(temp)
            (root/'tools').mkdir();(root/'provenance').mkdir()
            (root/'implementations/SYNTHETIC').mkdir(parents=True)
            (root/'dental_release').mkdir()
            (root/'dental_release/__init__.py').write_text('')
            (root/'dental_release/paths.py').write_text('class MissingInput(ValueError): pass\n')
            (root/'tools/import_check.py').write_bytes((ROOT/'tools/import_check.py').read_bytes())
            (root/'provenance/IMPORT_CONTRACTS.json').write_text(json.dumps(dict(paths={},external_code={},optional_packages=[])))
            (root/'implementations/SYNTHETIC/consumer.py').write_text(source)
            if dependency:
                (root/'implementations/SYNTHETIC/designgate').mkdir()
                (root/'implementations/SYNTHETIC/designgate/__init__.py').write_text('')
            result=subprocess.run([sys.executable,'-I','-B',str(root/'tools/import_check.py'),'--worker','SYNTHETIC'],capture_output=True,text=True,check=True)
            return json.loads(result.stdout)

    def test_missing_module_after_input_blocker_is_rejected(self):
        rows=self.probe("raise FileNotFoundError('external input absent')\nimport injected_missing_release_module\n")
        self.assertEqual(rows[0]['status'],'FAIL')
        self.assertEqual(rows[0]['missing_module'],'injected_missing_release_module')

    def test_missing_predecessor_gate_is_rejected(self):
        rows=self.probe('import designgate.gate\n',True)
        failed=[r for r in rows if r['status']=='FAIL']
        self.assertTrue(failed)
        self.assertEqual(failed[0]['missing_module'],'designgate.gate')

class AcceptedDeliveryCoverage(unittest.TestCase):
    def test_independent_frame_has_no_omitted_scope(self):
        self.assertEqual(validate(),dict(deliveries=110,direct=109,superseded=1,missing=0,scope='Reviewed source/document coverage, not physical validity'))

    def test_identical_count_different_delivery_is_rejected(self):
        rows=copy.deepcopy(json.loads((ROOT/'demos.json').read_text())['demos'])
        before=len(rows)
        next(r for r in rows if r['id']=='X34')['id']='injected-different-delivery'
        self.assertEqual(len(rows),before)
        with self.assertRaisesRegex(ValueError,'LANE_X34_STL_DESIGN_GATE'):validate(rows)

    def test_changed_bound_result_is_rejected(self):
        rows=copy.deepcopy(json.loads((ROOT/'demos.json').read_text())['demos'])
        next(r for r in rows if r['id']=='X38')['reviewed_result_sha256']='0'*64
        with self.assertRaises(ValueError):validate(rows)

    def test_five_corrected_exclusions_preserve_original_decisions(self):
        rows=json.loads((ROOT/'provenance/REVIEW_SCOPE.json').read_text())['records']
        corrected=[r for r in rows if r.get('previous_action')=='EXCLUDED_NO_ACCEPTED_SCOPE']
        self.assertEqual(len(corrected),5)
        self.assertEqual(sum(r['action']=='SUPERSEDED_BY_BOUND_ENTRY' for r in corrected),1)
        for r in corrected:
            self.assertTrue(any(str(r.get(k,'')).startswith('ACCEPT') for k in ('decision','graph_decision')) or r.get('sample_class')=='ACCEPTED')

if __name__=='__main__':unittest.main()
