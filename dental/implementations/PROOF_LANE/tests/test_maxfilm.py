import unittest
from gencad_bench.checks.cement_max import verify_negative, verify_positive

class MaximumFilm(unittest.TestCase):

    def test_exact_positive_and_bad_point(self):
        planes = [[1, 0, 0, 1], [-1, 0, 0, 1], [0, 1, 0, 1], [0, -1, 0, 1], [0, 0, 1, 1], [0, 0, -1, 0]]
        self.assertTrue(verify_positive(planes, [[1.08, 0, 0.5]], [[1, 0, 0.5]], 0.12))
        self.assertFalse(verify_positive(planes, [[1.3, 0, 0.5]], [[1, 0, 0.5]], 0.12))
        self.assertFalse(verify_positive(planes, [[1.08, 0, 0.5]], [[1.08, 0, 0.5]], 0.12))

    def test_exact_support_negative_and_tamper(self):
        cube = [[x, y, z] for x in [-1, 1] for y in [-1, 1] for z in [0, 1]]
        self.assertTrue(verify_negative(cube, [1.3, 0, 0.5], [1, 0, 0], 0.12))
        self.assertFalse(verify_negative(cube, [1.3, 0, 0.5], [0, 0, 0], 0.12))
        self.assertFalse(verify_negative(cube, [1.3, 0, 0.5], [-1, 0, 0], 0.12))
