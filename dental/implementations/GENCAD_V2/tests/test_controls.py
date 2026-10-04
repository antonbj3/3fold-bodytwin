import sys, unittest, copy, tempfile, subprocess, json
from pathlib import Path
import numpy as np
from scipy.sparse import csr_matrix
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / 'gencad_bench_v2/vendor/v1'))
from gencad_bench_v2.checks import all_checks, validate, feasibility
from gencad_bench_v2.generators import submit, optimize
from gencad_bench_v2.geometry import grid, shell, section_weights, area_weights
from gencad_bench_v2.common import dump, freeze, read, digest, DATA
from gencad_bench.checks import cement, milling, cement_max

def task(family='molar_crown'):
    xy = np.array([[-2.0, -2.0], [2.0, -2.0], [-2.0, 2.0], [2.0, 2.0]])
    f = np.array([[0, 1, 2], [1, 3, 2]])
    return dict(task_id='fixture', frame='fixed', family=family, xy=xy, faces=f, A=csr_matrix(np.eye(4)), obstacle_b=np.ones(4) * 2.0, preparation_z=np.zeros(4), requirements=dict(wall_mm=1.0, film_min_mm=0.04, film_max_mm=0.12, clearance_mm=0.02, connector_mm2=4.0, channel_max_deg=25.0, strut_mm=0.35), connector_x=[-0.5, 0.5], implant_axis=[0, 0, 1], minimum_strut_mm=0.5, prior=np.ones(4) * 1.2, weights=np.ones(4))

class Controls(unittest.TestCase):

    def good(self, fam='molar_crown'):
        t = task(fam)
        d = submit(t, np.ones(4) * 1.2)
        self.assertEqual(all_checks(t, d)['validity'], 'PASS')
        return (t, d)

    def test_contract_dimension_units_topology_nan(self):
        (t, d) = self.good()
        mutations = [('units', 'm'), ('frame', 'forged'), ('outer_vertices', [[0, 0]] * 4), ('inner_vertices', [[0, 0, 0, 1]] * 4), ('faces', [[0, 1, -1], [1, 3, 2]]), ('faces', [[0, 1, 2]]), ('reported_wall_mm', 999)]
        for (k, v) in mutations:
            bad = copy.deepcopy(d)
            bad[k] = v
            self.assertEqual(all_checks(t, bad)['validity'], 'INVALID', (k, v))
        for x in [float('nan'), float('inf')]:
            bad = copy.deepcopy(d)
            bad['outer_vertices'][0][2] = x
            self.assertEqual(all_checks(t, bad)['validity'], 'INVALID')
        bad = copy.deepcopy(d)
        bad['outer_vertices'][0][0] += 1
        self.assertEqual(all_checks(t, bad)['validity'], 'INVALID')

    def test_wall_nesting_and_film(self):
        (t, d) = self.good()
        for (z, inner, key) in [(0.9, 0.08, 'wall'), (-1, 0.08, 'nesting'), (1.2, 0.01, 'film_min'), (1.5, 0.2, 'film_max')]:
            b = submit(t, np.full(4, z), np.full(4, inner))
            self.assertEqual(all_checks(t, b)['checks'][key]['status'], 'FAIL')

    def test_antagonist_and_impossibility(self):
        (t, d) = self.good()
        self.assertEqual(all_checks(t, submit(t, np.ones(4) * 3))['checks']['antagonist']['status'], 'FAIL')
        self.assertEqual(feasibility(t)['status'], 'FEASIBLE')
        t['obstacle_b'] = np.ones(4) * 0.9
        self.assertEqual(feasibility(t)['status'], 'INFEASIBLE')
        self.assertEqual(optimize(t)['status'], 'ABSTAIN')

    def test_connector_channel_lattice(self):
        for (fam, key, change) in [('bridge3', 'connector_area', lambda t: t['requirements'].update(connector_mm2=5)), ('implant_crown', 'channel_angle', lambda t: t.update(implant_axis=[0.6, 0, 0.8])), ('lattice_onlay', 'strut_width', lambda t: t.update(minimum_strut_mm=0.2))]:
            (t, d) = self.good(fam)
            change(t)
            self.assertEqual(all_checks(t, d)['checks'][key]['status'], 'FAIL')

    def test_unknown_not_success(self):
        (t, d) = self.good()
        t['A'] = csr_matrix((0, 4))
        t['obstacle_b'] = []
        r = all_checks(t, d)
        self.assertEqual(r['validity'], 'UNKNOWN')
        self.assertIsNone(r['digital_score'])

    def test_section_closed_form(self):
        t = task()
        w = section_weights(t['xy'], t['faces'], 0.25)
        self.assertAlmostEqual(w @ np.full(4, 2.0), 8.0)
        self.assertGreater(abs(w @ np.full(4, 3.0) - 8), 1)

    def test_shell_closed_positive_volume(self):
        import trimesh
        (t, d) = self.good()
        (v, f) = shell(t['xy'], np.ones(4) * 1.2, np.ones(4) * 0.08, t['faces'])
        m = trimesh.Trimesh(v, f, process=False)
        self.assertTrue(m.is_watertight)
        self.assertGreater(m.volume, 0)
        self.assertFalse(trimesh.Trimesh(v, f[:-1], process=False).is_watertight)

    def test_patched_v1_dimension_attacks(self):
        planes = [[1, 0, 0, 2], [-1, 0, 0, 2], [0, 1, 0, 2], [0, -1, 0, 2], [0, 0, 1, 2], [0, 0, -1, 2]]
        self.assertEqual(cement.check([[0, 0, 0]], planes, 0.1)['status'], 'PASS')
        self.assertEqual(cement.check([[0, 0]], planes, 0.1)['status'], 'INVALID')
        openplane = [[0, 0, 1, 0]]
        self.assertTrue(milling.verify_positive(openplane, [[0, 0, 0]], 1, 0, [0, 0, -1], [[0, 0, -1]]))
        self.assertFalse(milling.verify_positive(openplane, [[0, 0, 0]], 1, 0, [0, 0, -1], [[0, 0]]))
        self.assertFalse(milling.verify_positive(openplane, [[0, 0, 0]], 1, 0, [0, 0, 0, 1], [[0, 0, -1]]))
        self.assertTrue(cement_max.verify_positive(planes, [[0, 0, 0]], [[0, 0, 0]], 0.12))
        self.assertFalse(cement_max.verify_positive(planes, [[0, 0, 100]], [], 0.12))

    def test_freeze_rejects_edit(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / 'lock.json'
            freeze(p, {'x': 1})
            freeze(p, {'x': 1})
            with self.assertRaises(ValueError):
                freeze(p, {'x': 2})
if __name__ == '__main__':
    unittest.main()
