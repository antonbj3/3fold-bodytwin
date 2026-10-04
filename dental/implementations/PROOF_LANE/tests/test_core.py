import unittest, copy, math
import numpy as np
from gencad_bench.checks import cement, milling, wall, occlusion, material, inherited
from gencad_bench.checks.exact import triangle_distance, q, sqrt_bounds
from gencad_bench.geometry import dome_profile, cap, planes_from_cap, shell
from gencad_bench.scoring import level1, interval_score, setup_match
from gencad_bench.schema import validate_split, validate_design
BOX = [[1, 0, 0, 1], [-1, 0, 0, 1], [0, 1, 0, 1], [0, -1, 0, 1], [0, 0, 1, 1]]

class ExactCore(unittest.TestCase):

    def test_cement_positive_and_injected_thin_film(self):
        self.assertEqual(cement.check([[0, 0, 0]], BOX, 0.9)['status'], 'PASS')
        self.assertEqual(cement.check([[0, 0, 0.99]], BOX, 0.03)['status'], 'FAIL')
        self.assertEqual(cement.check([], BOX, 0.03)['status'], 'INVALID')

    def test_cement_steep_wall_normal_not_vertical(self):
        self.assertEqual(cement.check([[0, 0, 0]], [[3, 0, 4, 1]], 0.2)['status'], 'PASS')
        self.assertEqual(cement.check([[0, 0, 0]], [[3, 0, 4, 0.9]], 0.2)['status'], 'FAIL')

    def test_mill_positive_and_corrupt_center(self):
        p = [[0, 0, 1], [0.1, 0, 1], [0, 0.1, 1]]
        c = [[0, 0, 0.5], [0.1, 0, 0.5], [0, 0.1, 0.5]]
        self.assertTrue(milling.verify_positive(BOX, p, 0.5, 0, [0, 0, -1], c))
        c[0][2] = 0.8
        self.assertFalse(milling.verify_positive(BOX, p, 0.5, 0, [0, 0, -1], c))

    def test_mill_negative_and_corrupt_dual(self):
        y = [1, 0, 1, 0, 1]
        self.assertTrue(milling.verify_negative(BOX, [1, 1, 1], 0.5, 0.05, y))
        self.assertFalse(milling.verify_negative(BOX, [1, 1, 1], 0.5, 0.05, [0] * 5))
        self.assertFalse(milling.verify_negative(BOX, [1, 1, 1], 0.5, 0.05, [-1, 0, 1, 0, 1]))

    def test_mill_entry_direction_is_not_erosion(self):
        self.assertFalse(milling.verify_positive(BOX, [[0, 0, 1]], 0.5, 0, [0, 0, 1], [[0, 0, 0.5]]))

    def test_mill_smooth_face_impossibility(self):
        self.assertTrue(milling.tangent_obstruction([0, 0, 0], [0, 0, 1], [0.1, 0, 0.1], 0.5))
        self.assertFalse(milling.tangent_obstruction([0, 0, 0], [0, 0, 1], [2, 0, 0], 0.5))

    def test_mill_discovery_matches_exact_cube_control(self):
        good = milling.check(BOX, [[0, 0, 1]], 0.5, 0.05)
        bad = milling.check(BOX, [[1, 1, 1]], 0.5, 0.05)
        self.assertEqual(good['status'], 'PASS')
        self.assertEqual(bad['status'], 'FAIL')
        self.assertAlmostEqual(bad['max_discovery_distance_mm'], math.sqrt(3) / 2, places=7)

    def test_triangle_interior_intersection_not_vertices(self):
        a = [[0, 0, 0], [2, 0, 0], [0, 2, 0]]
        b = [[0.5, 0.5, -1], [0.5, 0.5, 1], [2, 2, 0]]
        self.assertEqual(triangle_distance(a, b)[0], 0)
        self.assertEqual(triangle_distance(a, [[0, 0, 1], [2, 0, 1], [0, 2, 1]])[0], 1)

    def test_triangle_exact_against_known_skew_edge_case(self):
        a = [[0, 0, 0], [1, 0, 0], [0, 1, 0]]
        b = [[0.5, -1, 1], [0.5, 2, 1], [0.5, 2, 2]]
        self.assertEqual(triangle_distance(a, b)[0], 1)

    def test_wall_positive_and_thinned_outer(self):
        inner = dome_profile(2, 3)
        self.assertEqual(wall.check(dome_profile(4, 5), inner, 0.5)['status'], 'PASS')
        self.assertEqual(wall.check(dome_profile(2.1, 3.1), inner, 0.5)['status'], 'FAIL')

    def test_insertion_parent_positive_and_impossible(self):
        positive = inherited.insertion_cone([[1, 0, 1], [-1, 0, 1], [0, 1, 1], [0, -1, 1]])
        negative = inherited.insertion_cone(np.vstack([np.eye(3), -np.eye(3)]))
        self.assertEqual(positive['status'], 'YES_STRICT')
        self.assertEqual(negative['status'], 'NO_NONZERO_DIRECTION')
        self.assertTrue(inherited.verify_cone_negative(np.vstack([np.eye(3), -np.eye(3)]), negative))
        missing = copy.deepcopy(negative)
        missing['certificates'] = []
        self.assertFalse(inherited.verify_cone_negative(np.vstack([np.eye(3), -np.eye(3)]), missing))
        shifted = copy.deepcopy(negative)
        shifted['certificates'][0]['A'][0][0] = '999'
        self.assertFalse(inherited.verify_cone_negative(np.vstack([np.eye(3), -np.eye(3)]), shifted))
        parent = inherited.load('test_parent_cone', 'LANE_NEXT_D_INSERTION_PROOF/code/cone.py')
        self.assertTrue(parent.recheck(negative))
        negative['certificates'][0]['y'] = ['0'] * len(negative['certificates'][0]['y'])
        self.assertFalse(parent.recheck(negative))

    def test_occlusion_known_gap_and_injected_penetration(self):
        lo = [[0, 0, 0], [1, 0, 0], [0, 1, 0]]
        up = [[0, 0, 1], [1, 0, 1], [0, 1, 1]]
        self.assertEqual(occlusion.pair(up, lo)['status'], 'PASS')
        up[0][2] = -0.1
        self.assertEqual(occlusion.pair(up, lo)['status'], 'FAIL')
        self.assertEqual(occlusion.check_scene([up], [lo])['status'], 'UNKNOWN')

    def test_occlusion_parent_full_pair_parity(self):
        U = [[[0, 0, 1], [1, 0, 1], [0, 1, 1]]]
        L = [[[0, 0, 0], [1, 0, 0], [0, 1, 0]]]
        r = inherited.projected_extremum(U, L)
        self.assertLessEqual(r['parity_mm'], 1e-07)
        self.assertEqual(r['candidate']['minimum_gap_mm'], float(occlusion.pair(U[0], L[0])['minimum_gap_mm']))

    def test_material_shape_independent_veto(self):
        self.assertEqual(material.point_conflict([0.9, 0, 0], [1.1, 0, 0], BOX, 0.3)['status'], 'FAIL')
        self.assertEqual(material.point_conflict([0, 0, 0], [0.1, 0, 0], BOX, 0.3)['status'], 'UNKNOWN')
        self.assertEqual(cement.check([[0, 0, 0]], BOX, 0.3)['status'], 'PASS')

    def test_score_veto_missing_and_correct(self):
        required = ['wall', 'cement']
        a = {k: dict(status='PASS') for k in required}
        self.assertEqual(level1(a, required)['score'], 1)
        a['wall']['status'] = 'FAIL'
        self.assertEqual(level1(a, required)['score'], 0)
        del a['wall']
        self.assertIsNone(level1(a, required)['score'])
        self.assertFalse(level1(a, required)['eligible'])

    def test_interval_score_rejects_bad_prediction(self):
        self.assertGreater(interval_score(0, 1, 5), interval_score(4, 6, 5))
        self.assertGreater(interval_score(-1000, 1000, 5), interval_score(4, 6, 5))
        with self.assertRaises(ValueError):
            interval_score(-math.inf, math.inf, 5)

    def test_setup_mismatch_cannot_be_calibrated(self):
        fields = ['material_product', 'specimen', 'thickness_mm', 'support', 'cement', 'load_angle_deg', 'ageing']
        t = {k: 1 for k in fields}
        self.assertTrue(setup_match(t, t))
        other = dict(t)
        other['support'] = 2
        self.assertFalse(setup_match(t, other))
        other['support'] = None
        self.assertFalse(setup_match(t, other))

    def test_sqrt_outward(self):
        for v in [q('0'), q('2'), q('0.00001'), q('9/4')]:
            (lo, hi) = sqrt_bounds(v)
            self.assertLessEqual(lo * lo, v)
            self.assertGreaterEqual(hi * hi, v)

    def test_geometry_topology_and_nan(self):
        a = dome_profile(4, 5)
        b = dome_profile(2, 3)
        m = shell(a, b)
        self.assertTrue(m.is_watertight)
        self.assertTrue(m.is_winding_consistent)
        a['radii_mm'][0][0] = float('nan')
        with self.assertRaises(ValueError):
            cap(a)

    def test_split_patient_leakage(self):
        a = dict(subject_group='a', split='development', public=dict(population_training_groups=['b']))
        b = copy.deepcopy(a)
        b['split'] = 'test'
        with self.assertRaises(ValueError):
            validate_split([a, b])
        b['subject_group'] = 'c'
        self.assertTrue(validate_split([a, b]))
        b['subject_group'] = 'b'
        with self.assertRaises(ValueError):
            validate_split([a, b])
if __name__ == '__main__':
    unittest.main()
