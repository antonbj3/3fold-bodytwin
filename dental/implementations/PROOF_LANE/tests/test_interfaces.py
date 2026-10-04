import unittest, copy, json
import numpy as np
from gencad_bench.io import ROOT, digest, freeze
from gencad_bench.schema import validate_design, validate_task
from gencad_bench.evaluate import evaluate
from gencad_bench.generators.rounded import round_intaglio
from gencad_bench.geometry import dome_profile, planes_from_cap, cap
from gencad_bench.checks import milling

class Interfaces(unittest.TestCase):

    def task(self):
        return json.loads((ROOT / 'data/tasks.json').read_text())[0]

    def test_all_invalid_plugin_keeps_denominator_and_null_shape(self):
        from run_bench import material_shape_summary
        rr = [dict(generator='broken', material='M', shape=None), dict(generator='broken', material='M', shape=None)]
        s = material_shape_summary(rr, 'broken')['M']
        self.assertIsNone(s['RMS_mm'])
        self.assertEqual(s['n_scored'], 0)
        self.assertEqual(s['n_attempted'], 2)

    def test_mesh_plugin_is_scored_unknown_not_falsely_certified(self):
        t = self.task()
        mesh = dict(vertices_mm=[[0, 0, 0], [1, 0, 0], [0, 1, 0]], faces=[[0, 1, 2]])
        d = dict(kind='triangle_mesh_v1', units='mm', frame=t['frame'], task_id=t['task_id'], outer_mesh=mesh, intaglio_mesh=mesh)
        validate_design(d, t)
        r = evaluate(t, d, dict(vertices=np.array(mesh['vertices_mm']), faces=np.array(mesh['faces'])), {})
        self.assertEqual(r['level1']['status'], 'UNKNOWN')
        self.assertAlmostEqual(r['shape']['symmetric_sampled_RMS_mm'], 0)
        d['outer_mesh']['faces'] = []
        self.assertEqual(evaluate(t, d, {}, {})['level1']['status'], 'INVALID')

    def test_units_and_requirements_are_not_generator_options(self):
        t = self.task()
        bad = copy.deepcopy(t)
        bad['requirements']['cement_min_mm'] = -1
        with self.assertRaises(ValueError):
            validate_task(bad)
        bad = copy.deepcopy(t)
        bad['units'] = 'um'
        with self.assertRaises(ValueError):
            validate_task(bad)

    def test_rounding_operation_and_radius_mutation(self):
        p = round_intaglio(dome_profile(1, 3))
        (v, _) = cap(p)
        pl = planes_from_cap(p)
        self.assertEqual(milling.check(pl, v.tolist(), 0.5, 0.05)['status'], 'PASS')
        self.assertEqual(milling.check(pl, v.tolist(), 3.0, 0.05)['status'], 'FAIL')

    def test_json_schema_matches_frozen_tasks_and_rejects_bad_units(self):
        import jsonschema
        schema = json.loads((ROOT / 'data/task.schema.json').read_text())
        for t in json.loads((ROOT / 'data/tasks.json').read_text()):
            jsonschema.validate(t, schema)
        t = self.task()
        t['units'] = 'm'
        with self.assertRaises(jsonschema.ValidationError):
            jsonschema.validate(t, schema)

    def test_frozen_artifact_rejects_changed_payload(self):
        from pathlib import Path
        from tempfile import TemporaryDirectory
        with TemporaryDirectory() as td:
            path = Path(td) / 'freeze.json'
            freeze(path, {'value': 1})
            freeze(path, {'value': 1})
            with self.assertRaises(ValueError):
                freeze(path, {'value': 2})
            stored = json.loads(path.read_text())
            stored['payload']['value'] = 7
            path.write_text(json.dumps(stored))
            with self.assertRaises(ValueError):
                freeze(path, {'value': 1})

    def test_homothetic_control_has_exact_supporting_faces(self):
        from gencad_bench.generators.rolling_ball import profile
        p = profile(dome_profile(1, 3))
        self.assertTrue(planes_from_cap(p))
        bad = copy.deepcopy(p)
        bad['radii_mm'][2][0] *= 2
        with self.assertRaises(ValueError):
            planes_from_cap(bad)
