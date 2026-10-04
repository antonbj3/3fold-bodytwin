import sys
import unittest
from pathlib import Path
import numpy as np
from scipy.spatial.transform import Rotation
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'code'))
from metrology import kabsch, transform, validate_transform, summarize

class MetrologyTests(unittest.TestCase):

    def test_known_pose_and_reflection_rejection(self):
        x = np.array([[0, 0, 0], [2, 0, 0], [0, 3, 0], [0, 0, 4]], float)
        T = np.eye(4)
        T[:3, :3] = Rotation.from_euler('xyz', [30, 7, -40], degrees=True).as_matrix()
        T[:3, 3] = [3, 5, -4]
        fitted = kabsch(transform(x, T), x)
        self.assertLess(np.max(abs(transform(transform(x, T), fitted) - x)), 1e-12)
        bad = np.eye(4)
        bad[0, 0] = -1
        with self.assertRaises(ValueError):
            validate_transform(bad)
        with self.assertRaises(ValueError):
            kabsch(x[:3] * [1, 0, 0], x[:3])

    def test_precision_is_not_pairwise_or_pooled_rms(self):
        d = np.array([[0.08, 0.1, 0.12] * 8, [0.09, 0.1, 0.11] * 8, [0.1, 0.1, 0.1] * 8])
        region = np.tile(np.repeat(['marginal', 'intaglio', 'occlusal', 'axial'], 6), 1)
        out = summarize(d, abs(d), region, np.ones(24), 0.3)['all']
        self.assertAlmostEqual(out['rms_pairwise_difference_um'] / np.sqrt(2), out['repeatability_normal_sd_um'], places=12)
        self.assertGreater(out['pooled_normal_rms_um'], 10 * out['repeatability_normal_sd_um'])

    def test_missing_surface_not_hidden_by_rms(self):
        d = np.zeros((2, 24))
        distance = np.zeros_like(d)
        distance[:, 0:6] = 2
        region = np.repeat(['marginal', 'intaglio', 'occlusal', 'axial'], 6)
        r = summarize(d, distance, region, np.ones(24), 0.3)
        self.assertEqual(r['marginal']['status'], 'INSUFFICIENT_COVERAGE')
        self.assertEqual(r['marginal']['retained_sites'], 0)
if __name__ == '__main__':
    unittest.main()
