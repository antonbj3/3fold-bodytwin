import sys, unittest
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'code'))
import numpy as np, trimesh
from collision import *

def cup(aperture=0.7, height=6):
    meshes = []
    for (ext, loc) in [([6, 6, 0.4], [0, 0, -0.2]), ([1, 6, height], [aperture + 0.5, 0, height / 2]), ([1, 6, height], [-aperture - 0.5, 0, height / 2]), ([2 * aperture, 1, height], [0, aperture + 0.5, height / 2]), ([2 * aperture, 1, height], [0, -aperture - 0.5, height / 2])]:
        m = trimesh.creation.box(ext)
        m.apply_translation(loc)
        meshes.append(m)
    m = trimesh.util.concatenate(meshes)
    return Scene(m.vertices, m.faces)

def tool(d=0.6, neck=3, gauge=22):
    return dict(diameter_mm=d, neck_reach_mm=neck, shank_mm=3.0, gauge_mm=gauge, holder_diameter_mm=10.0, holder_length_mm=12.0)

class CollisionTests(unittest.TestCase):

    def test_plane_closed_form(self):
        tri = np.array([[[-10.0, -10, 0], [10, -10, 0], [0, 10, 0]]])
        self.assertAlmostEqual(segment_triangle(np.array([0, 0, 2.0]), np.array([0, 0, 4.0]), tri)[0], 2.0, places=12)
        self.assertEqual(segment_triangle(np.array([0, 0, -2.0]), np.array([0, 0, 4.0]), tri)[0], 0.0)

    def test_continuous_control_enclosure(self):
        rng = np.random.default_rng(53)
        for _ in range(30):
            tri = rng.normal(size=(5, 3, 3))
            (a, b) = rng.normal(size=(2, 3))
            got = segment_triangle(a, b, tri)
            v = tri.reshape(-1, 3)
            f = np.arange(len(v)).reshape(-1, 3)
            mesh = trimesh.Trimesh(v, f, process=False)
            pts = a + (b - a) * np.linspace(0, 1, 401)[:, None]
            (_, dd, _) = trimesh.proximity.closest_point_naive(mesh, pts)
            upper = dd.min()
            lower = max(0.0, upper - np.linalg.norm(b - a) / 800)
            self.assertGreaterEqual(got.min(), lower - 1e-07)
            self.assertLessEqual(got.min(), upper + 1e-07)
            self.assertGreater(got.min() + 0.1, upper + 1e-07)

    def test_larger_tool_can_reach_when_small_shank_cannot(self):
        scene = cup()
        p = np.array([0.0, 0.0, 0.0])
        n = np.array([0.0, 0.0, 1.0])
        small = search(scene, p, n, tool(0.6, 3), 4)
        big = search(scene, p, n, tool(1.0, 16), 4)
        self.assertEqual(small['status'], 'NOT_FOUND_ON_POSE_GRID')
        self.assertEqual(big['status'], 'FOUND')
        self.assertIn('shank', small['tested_rejections'])

    def test_holder_rejects(self):
        scene = cup()
        p = np.zeros(3)
        n = np.array([0.0, 0.0, 1.0])
        good = search(scene, p, n, tool(1.0, 16, 22), 4)
        bad = search(scene, p, n, tool(1.0, 16, 4), 4)
        self.assertEqual(good['status'], 'FOUND')
        self.assertEqual(bad['status'], 'NOT_FOUND_ON_POSE_GRID')
        self.assertIn('holder', bad['tested_rejections'])

    def test_axes_actual_plane_not_z_cone(self):
        d4 = directions(4)
        d5 = directions(5)
        self.assertEqual(float(abs(d4[:, 0]).max()), 0.0)
        self.assertAlmostEqual(abs(d5[:, 0]).max(), np.sin(np.radians(35)), places=14)
        self.assertTrue(np.allclose(d5[:24], d4))
        self.assertAlmostEqual(d4[:, 1].max(), 1.0)

    def test_enclosed_void_impossible_straight_exit(self):
        outer = trimesh.creation.box([6, 6, 6])
        inner = trimesh.creation.icosphere(subdivisions=2, radius=1)
        inner.faces = inner.faces[:, ::-1]
        m = trimesh.util.concatenate([outer, inner])
        s = Scene(m.vertices, m.faces)
        p = inner.triangles_center[0]
        n = inner.face_normals[0]
        self.assertEqual(search(s, p, n, tool(0.1, 3), 5)['status'], 'NOT_FOUND_ON_POSE_GRID')

    def test_union_signed_stock_and_overcut(self):
        p = np.array([[0.0, 0.0, 0.0]])
        n = np.array([[0.0, 0.0, 1.0]])
        self.assertAlmostEqual(first_union_ball_entry(p, n, [[0, 0, 0.55]], [0.5])[0], 0.05, places=14)
        self.assertAlmostEqual(first_union_ball_entry(p, n, [[0, 0, 0.45]], [0.5])[0], -0.05, places=14)

    def test_input_mutations(self):
        with self.assertRaises(ValueError):
            Scene(np.array([[0, 0, 0], [1, 0, 0], [2, 0, 0]]), [[0, 1, 2]])
if __name__ == '__main__':
    unittest.main()
