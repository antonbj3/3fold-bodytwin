"""Independent published controls plus adversarial input tests."""
import copy, hashlib, io, json, math, os, sys, tempfile, unittest, zipfile
from pathlib import Path
import numpy as np
import trimesh
from scipy.stats import binom, binomtest, fisher_exact
from scipy.integrate import quad
from scipy.stats import chi2, norm, t
from statsmodels.stats.power import TTestIndPower
R = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(R / 'code'))
from mesh_transport import *
from exportgate import export, check
from labcount import upper, plan, endpoint_summary, binary_power, fisher_pvalues, continuous_power, compare

class Statistics(unittest.TestCase):

    def test_exact_bound_published_scipy(self):
        differences = []
        for n in range(1, 41):
            for k in range(n + 1):
                control = binomtest(k, n, alternative='less').proportion_ci(0.95, method='exact').high
                differences.append(abs(upper(k, n, 0.05) - control))
        self.assertLess(max(differences), 1e-10)
        self.assertAlmostEqual(upper(0, 3, 0.05), 0.6315968501359613, places=12)
        self.assertEqual(upper(0, 0, 0.05), 1)
        self.assertGreater(abs(0.05 - binomtest(0, 3, alternative='less').proportion_ci(0.95).high), 1e-10)

    def test_minimum_and_acceptance(self):
        for k in [0, 1, 2]:
            for target in [0.01, 0.1, 0.3]:
                q = plan(target, 0.95, k)
                self.assertLess(q['upper_if_acceptance_boundary'], target)
                self.assertGreaterEqual(q['previous_n_upper'], target)
        self.assertEqual(plan(0.1, 0.95)['n_per_condition'], 29)
        self.assertEqual(plan(0.1, 0.95, family_tests=4)['n_per_condition'], 42)
        q = plan(0.1, 0.95, true_risk=0.1)
        self.assertAlmostEqual(q['probability_of_acceptance_at_assumed_risk'], 0.9 ** 29, places=12)
        self.assertLess(q['probability_of_acceptance_at_assumed_risk'], 0.05)
        self.assertGreater(abs(0.8 - 0.9 ** 29), 1e-08)

    def test_fisher_published_control(self):
        errors = []
        for n in [2, 3, 5, 8, 12]:
            pv = fisher_pvalues(n)
            for x in range(n + 1):
                for y in range(n + 1):
                    p = fisher_exact([[x, n - x], [y, n - y]], alternative='two-sided').pvalue
                    errors.append(abs(p - pv[x, y]))
        self.assertLess(max(errors), 1e-10)

    def test_power_independent_count_enumeration(self):
        n = 69
        pa = 0.3
        pb = 0.1
        published = 0.0
        for x in range(n + 1):
            for y in range(n + 1):
                if fisher_exact([[x, n - x], [y, n - y]], alternative='two-sided').pvalue <= 0.05:
                    published += binom.pmf(x, n, pa) * binom.pmf(y, n, pb)
        self.assertAlmostEqual(binary_power(n, pa, pb), published, places=10)
        self.assertGreater(abs(0.95 - published), 1e-08)
        q = compare('binary', p_a=pa, p_b=pb, max_n=100)
        self.assertEqual(q['n_per_group'], 69)
        self.assertLess(q['maximum_previous_power'], 0.8)
        self.assertGreaterEqual(q['achieved_power'], 0.8)
        for p in [0.01, 0.1, 0.3, 0.5, 0.8]:
            self.assertLessEqual(binary_power(29, p, p), 0.05 + 1e-12)

    def test_continuous_power_published_statsmodels(self):
        for n in [2, 10, 26, 80]:
            for d in [0.2, 0.8, 1.5]:
                control = float(TTestIndPower().power(d, n, alpha=0.05))
                if math.isfinite(control):
                    self.assertAlmostEqual(continuous_power(n, d), control, places=10)
                df = 2 * n - 2
                delta = d * math.sqrt(n / 2)
                cut = t.isf(0.025, df)
                (lo, hi) = chi2.ppf([1e-12, 1 - 1e-12], df)
                (value, error) = quad(lambda v: chi2.pdf(v, df) * (norm.cdf(delta - cut * math.sqrt(v / df)) + norm.cdf(-delta - cut * math.sqrt(v / df))), lo, hi, epsabs=1e-10, epsrel=1e-10, limit=200)
                self.assertLess(error, 1e-09)
                self.assertLess(abs(continuous_power(n, d) - value), 1e-08)
                self.assertGreater(abs(0.4 - value), 1e-08)
        q = compare('continuous', effect=0.8)
        self.assertEqual(q['n_per_group'], 26)
        self.assertGreaterEqual(q['achieved_power'], 0.8)
        self.assertLess(q['maximum_previous_power'], 0.8)

    def test_censored_observation_and_duplicates(self):
        rows = json.loads((R / 'inputs/specimens.json').read_text())
        q = endpoint_summary(rows, 5000000.0, 'MUSLA 04020', 100)
        self.assertEqual((q['n'], q['known_failures'], q['complete_survivors'], q['unknown_early_censors']), (3, 0, 3, 0))
        self.assertAlmostEqual(q['population_upper_one_sided'], 0.6315968501359613, places=12)
        with self.assertRaises(ValueError):
            endpoint_summary(rows + [rows[0]], 5000000.0, 'MUSLA 04020', 100)
        base = dict(specimen_id=1, configuration='g', max_load_N=100, cycles=2, cycles_unit='cycles', force_unit='N', failure=0, censored=1, locator='our_own_fixture')
        q = endpoint_summary([base], 10, 'g', 100)
        self.assertEqual(q['unknown_early_censors'], 1)
        self.assertEqual(q['population_upper_one_sided'], 1)
        rr = [dict(base, specimen_id=i, cycles=2 if i < 2 else 10) for i in range(4)]
        q = endpoint_summary(rr, 10, 'g', 100)
        for k in [0, 1, 2]:
            self.assertLessEqual(upper(k, 4, 0.05), q['population_upper_one_sided'])
        self.assertGreater(q['population_upper_one_sided'], upper(0, 4, 0.05))
        with self.assertRaises(ValueError):
            endpoint_summary([dict(base, force_unit=None)], 10, 'g', 100)
        with self.assertRaises(ValueError):
            endpoint_summary([dict(base, failure=1, censored=1)], 10, 'g', 100)

    def test_invalid_statistics_inputs(self):
        for kwargs in [dict(risk=0), dict(confidence=1), dict(risk=float('nan')), dict(max_failures=-1), dict(family_tests=0)]:
            with self.assertRaises(ValueError):
                plan(**kwargs)
        with self.assertRaises(ValueError):
            upper(4, 3, 0.05)
        with self.assertRaises(ValueError):
            compare('binary', p_a=0.1, p_b=0.1)
        with self.assertRaises(ValueError):
            continuous_power(10, float('nan'))
        self.assertEqual(compare('binary', p_a=0.01, p_b=0.011, max_n=3)['status'], 'SEARCH_LIMIT')

def mutate_3mf(source, dst, operation):
    with zipfile.ZipFile(source) as src, zipfile.ZipFile(dst, 'w', zipfile.ZIP_DEFLATED) as out:
        for name in src.namelist():
            b = src.read(name)
            if name.endswith('.model'):
                r = ET.fromstring(b)
                operation(r)
                b = ET.tostring(r)
            out.writestr(name, b)

class Transport(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.m = trimesh.creation.box(extents=[2, 3, 4])
        self.tri = self.m.vertices[self.m.faces]
        self.source = self.root / 'source.stl'
        write_stl(self.source, self.tri)
        self.pred = self.root / 'pred.json'
        dump(self.pred, dict(frozen_before_measurement=True))
        self.reg = self.root / 'regions.json'
        dump(self.reg, dict(labels=['fixture_A' if i % 2 else 'fixture_B' for i in range(12)], semantics='our_own_fixture transport probe'))
        self.out = self.root / 'export'
        export(self.source, self.out, 'mm', self.reg, 'fixture_material', 1.0, 'final_sintered', self.pred)
        self.asset = self.out / 'model.stl'
        self.side = self.out / 'model.stl.json'
        self.h = sha(self.side)
        self.native = self.out / 'model.3mf'
        self.nside = self.out / 'model.3mf.json'
        self.nh = sha(self.nside)

    def tearDown(self):
        self.tmp.cleanup()

    def test_independent_trimesh_read_and_reorder(self):
        for name in ['model.stl', 'model_ascii.stl', 'model.3mf']:
            obj = trimesh.load(self.out / name, process=False)
            m = obj.to_geometry() if isinstance(obj, trimesh.Scene) else obj
            current = np.asarray(m.vertices)[m.faces]
            correspondence(self.tri, current)
        perm = np.arange(12)[::-1]
        path = self.root / 'reordered.stl'
        write_stl(path, np.roll(self.tri[perm], 1, axis=1))
        q = check(path, self.side, self.h, 'mm', True)
        labels = json.loads(self.reg.read_text())['labels']
        self.assertEqual(q['regions_in_imported_face_order'], [labels[i] for i in perm])

    def test_units_and_coordinate_faults(self):
        with self.assertRaises(Refused):
            check(self.asset, self.side, self.h)
        with self.assertRaises(Refused):
            check(self.asset, self.side, self.h, 'inch')
        for tri in [self.tri * 1.01, self.tri + [0.001, 0, 0], self.tri[:, ::-1]]:
            p = self.root / 'wrong.stl'
            write_stl(p, tri)
            with self.assertRaises(Refused):
                check(p, self.side, self.h, 'mm', True)
        matrix = np.array([[1, 0.1, 0], [0, 1, 0], [0, 0, 1.0]])
        p = self.root / 'shear.stl'
        write_stl(p, self.tri @ matrix)
        with self.assertRaises(Refused):
            check(p, self.side, self.h, 'mm', True)

    def test_asset_contract_and_prediction_hash_faults(self):
        p = self.root / 'changed.stl'
        write_stl(p, self.tri[::-1])
        with self.assertRaises(Refused):
            check(p, self.side, self.h, 'mm')
        c = json.loads(self.side.read_text())
        c['payload']['sinter']['factor'] = 2
        dump(self.side, c)
        with self.assertRaises(Refused):
            check(self.asset, self.side, self.h, 'mm')
        with self.assertRaises(Refused):
            check(self.asset, self.side, sha(self.side), 'mm')
        (self.out / 'FROZEN_PREDICTIONS.json').write_text('{}')
        with self.assertRaises(Refused):
            check(self.native, self.nside, self.nh)

    def test_topology_and_truncation(self):
        for tri in [self.tri[:-1], np.r_[self.tri, self.tri[:1]], np.r_[self.tri[:-1], self.tri[:1]], np.zeros((12, 3, 3))]:
            with self.assertRaises(Refused):
                validate_triangles(tri)
        p = self.root / 'truncated.stl'
        p.write_bytes(self.source.read_bytes()[:-1])
        with self.assertRaises((ValueError, UnicodeDecodeError)):
            load_mesh(p, 'mm')

    def test_native_metadata_faults_and_fallback(self):

        def label_swap(r):
            gs = r.findall('.//{' + SETS + '}triangleset')
            (a, b) = (gs[0].get('name'), gs[1].get('name'))
            gs[0].set('name', b)
            gs[1].set('name', a)

        def material(r):
            r.find('.//{' + CORE + '}base').set('name', 'wrong_material')

        def wrong_unit(r):
            r.set('unit', 'inch')

        def transformed(r):
            r.find('.//{' + CORE + '}item').set('transform', '1 0 0 0 1 0 0 0 1 .1 0 0')

        def wrong_contract(r):
            r.find('{' + CORE + '}metadata').text = '0' * 64
        for (i, op) in enumerate([label_swap, material, wrong_unit, transformed, wrong_contract]):
            p = self.root / f'bad{i}.3mf'
            mutate_3mf(self.native, p, op)
            with self.assertRaises(Refused):
                check(p, self.nside, self.nh, converted=True)

        def stripped(r):
            for x in list(r.findall('{' + CORE + '}metadata')):
                r.remove(x)
            mesh = r.find('.//{' + CORE + '}mesh')
            mesh.remove(mesh.find('{' + SETS + '}trianglesets'))
        p = self.root / 'stripped.3mf'
        mutate_3mf(self.native, p, stripped)
        self.assertEqual(check(p, self.nside, self.nh, converted=True)['region_source'], 'geometry-bound sidecar')

        def in_inches(r):
            r.set('unit', 'inch')
            for v in r.findall('.//{' + CORE + '}vertex'):
                for k in ['x', 'y', 'z']:
                    v.set(k, str(float(v.get(k)) / 25.4))
        p = self.root / 'inches.3mf'
        mutate_3mf(self.native, p, in_inches)
        self.assertEqual(check(p, self.nside, self.nh, converted=True)['status'], 'PASS')

    def test_unsupported_objects_and_rewrite(self):

        def multi(r):
            build = r.find('{' + CORE + '}build')
            ET.SubElement(build, '{' + CORE + '}item', {'objectid': '2'})
        p = self.root / 'multi.3mf'
        mutate_3mf(self.native, p, multi)
        with self.assertRaises(Refused):
            check(p, self.nside, self.nh, converted=True)
        with self.assertRaises(Refused):
            export(self.source, self.out, 'mm', None, 'm', 1, 'final_sintered', self.pred)
        for f in [0, float('nan')]:
            with self.assertRaises(Refused):
                export(self.source, self.root / 'invalid', 'mm', None, 'm', f, 'final_sintered', self.pred)

    def test_real_source_refusals_preserve_cutoff(self):
        refused = []
        for p in sorted((R / 'inputs/exports').glob('*/*.stl')):
            (tri, _) = load_mesh(p, 'mm')
            try:
                validate_triangles(tri)
            except Refused:
                refused.append((p.parent.name, p.stem))
        self.assertEqual(refused, [('M1', 'die'), ('M1', 'preparation'), ('M2', 'die'), ('M2', 'preparation')])

    def test_design_gate_source_join(self):
        report = self.root / 'report.json'
        dump(report, dict(inputs={'crown': {'sha256': sha(self.source)}}, verdict='FAIL'))
        out = self.root / 'matched_report'
        export(self.source, out, 'mm', None, 'm', 1, 'final_sintered', self.pred, report)
        q = check(out / 'model.stl', out / 'model.stl.json', sha(out / 'model.stl.json'), 'mm')
        self.assertEqual(q['design_gate']['verdict'], 'FAIL')
        self.assertEqual(q['manufacturing_release'].split(';')[0], 'UNKNOWN')
        with self.assertRaises(Refused):
            export(self.source, self.root / 'wrong_report', 'mm', None, 'm', 1, 'final_sintered', self.pred, R / 'inputs/X34_real_designgate.json')
        self.assertFalse((self.root / 'wrong_report').exists())
if __name__ == '__main__':
    unittest.main()
