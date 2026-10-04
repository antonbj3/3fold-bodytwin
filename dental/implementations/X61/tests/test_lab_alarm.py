import copy, json, math, sys, unittest
from pathlib import Path
import numpy as np
from scipy.integrate import quad
from scipy.special import logsumexp
from scipy.stats import norm, weibull_min
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'code'))
from lab_alarm import *
from synthetic import series

class LabTests(unittest.TestCase):

    @classmethod
    def setUpClass(c):
        c.p = verified(R / 'inputs/DEMO_PROFILE.json')
        c.f = json.loads((R / 'inputs/X1B_FROZEN_PREDICTIONS.json').read_text())
        c.groups = model_groups(c.f)

    def test_external_weibull_cdf_and_density(self):
        for m in (3.0, 6.0, 13.0):
            for f in (200.0, 1800.0, 3000.0, 6000.0):
                expected = quad(lambda z: weibull_min.cdf(f * math.exp(-z), m, scale=3000) * norm.pdf(z, scale=0.01), -0.12, 0.12, epsabs=2e-13)[0]
                self.assertLess(abs(force_cdf(f, m, 3000.0, 0.01) - expected), 1e-10)
                got = math.exp(float(log_force_likelihood(f, np.array([m]), np.array([3000.0]), 0.01)[0])) / f
                ref = quad(lambda z: weibull_min.pdf(f * math.exp(-z), m, scale=3000) * math.exp(-z) * norm.pdf(z, scale=0.01), -0.12, 0.12, epsabs=2e-13)[0]
                self.assertLess(abs(got - ref), 1e-10)

    def test_full_data_posterior_recompute(self):
        rows = series(6100, self.p, self.f)
        post = JointPosterior(self.p, self.groups)
        for r in rows[:18]:
            post.update_force(r)
        z = post.initial.copy()
        for row in rows[:18]:
            key = (row['design_id'], int(row['load_angle_deg']))
            g = self.groups[key]
            factor = g['base_factor'] * np.where(post.form == 1, g['Q'] if g['interaction'] else 1.0, 1.0)
            f = float(row['fracture_force_N'])
            sc = post.lam * factor
            x = math.log(f) - math.sqrt(2) * self.p['null']['force_log_noise_sd'] * _HX if False else None
            from numpy.polynomial.hermite import hermgauss
            (xx, ww) = hermgauss(32)
            ww = ww / math.sqrt(math.pi)
            true = f * np.exp(-math.sqrt(2) * 0.01 * xx)
            dens = weibull_min.pdf(true[None, :], post.m[:, None], scale=sc[:, None]) * true[None, :]
            z += np.log(np.maximum(dens @ ww, 1e-300))
        w = np.exp(z - logsumexp(z))
        self.assertLess(float(np.max(np.abs(post.weights() - w))), 1e-10)

    def test_gaussian_port_closed_form(self):
        p = GaussianPort(0.8, 0.02)
        values = [0.799, 0.801, 0.802]
        sds = [0.001, 0.002, 0.001]
        for (x, s) in zip(values, sds):
            p.update(x, s)
        r = p.report('1')
        precision = 1 / 0.02 ** 2 + sum((1 / s ** 2 for s in sds))
        mu = (0.8 / 0.02 ** 2 + sum((x / s ** 2 for (x, s) in zip(values, sds)))) / precision
        self.assertLess(abs(r['mean'] - mu), 1e-12)
        self.assertLess(abs(r['sd'] - 1 / math.sqrt(precision)), 1e-12)

    def test_categorical_expected_wealth_exact(self):
        K = 8
        vals = (np.arange(K) + 0.5) / K - 0.5
        for stat in [vals, vals ** 2]:
            w = []
            for x in stat:
                a = Alarm(alpha=0.01)
                (wealth, _) = a.update(float(stat.mean()), float(x), [(1 / K, float(v)) for v in stat])
                w.append(wealth)
            self.assertLess(abs(np.mean(w) - 1), 1e-12)

    def test_record_fault_does_not_update_strength(self):
        a = series(6100, self.p, self.f, 'nominal')
        b = series(6100, self.p, self.f, 'one_bad_record')
        (r, post) = analyze(b, self.p, self.f)
        self.assertEqual(r['monitor']['first_alarm_at_specimen']['measurement_record'], 9)
        self.assertEqual(post.n, 71)
        control = JointPosterior(self.p, self.groups)
        for (i, x) in enumerate(a):
            if i != 8:
                control.update_force(x)
        self.assertLess(float(np.max(abs(control.weights() - post.weights()))), 1e-12)

    def test_missing_channels_stay_unknown(self):
        rows = series(6100, self.p, self.f)
        for x in rows:
            x['_sidecar'] = {}
        (out, p) = analyze(rows, self.p, self.f)
        self.assertEqual(out['posterior']['sinter_factor']['observations'], 0)
        self.assertEqual(out['posterior']['cement_film_mean']['observations'], 0)
        self.assertIn('UNKNOWN', out['posterior']['sigma0_material_MPa']['status'])

    def test_mechanism_and_protocol_cannot_fit_away(self):
        rows = series(6100, self.p, self.f)[:2]
        rows[0]['failure_mode'] = 'cement_debonding'
        rows[1]['die_E_MPa'] = '2000'
        (out, p) = analyze(rows, self.p, self.f)
        self.assertEqual(p.n, 0)
        self.assertEqual(out['dropout']['force_updates_withheld'], 2)

    def test_batch_change_does_not_merge(self):
        rows = series(6100, self.p, self.f)[:3]
        rows[1]['material_batch'] = 'NEW'
        (out, p) = analyze(rows, self.p, self.f)
        self.assertEqual(p.n, 1)
        self.assertEqual(p.sinter.n, 1)

    def test_profile_hash_falsifier(self):
        import tempfile
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'p.json'
            path.write_text('{}')
            path.with_name(path.name + '.sha256').write_text('0' * 64)
            with self.assertRaises(ValueError):
                verified(path)
if __name__ == '__main__':
    unittest.main()

class PairedCLI(unittest.TestCase):

    def test_paired_updates_real_sinter_and_scanner_separately(self):
        p = verified(R / 'inputs/DEMO_PAIRED_PROFILE.json')
        f = json.loads((R / 'inputs/X1B_FROZEN_PREDICTIONS.json').read_text())
        rows = series(6100, p, f)
        for row in rows:
            s = row['_sidecar']
            s['post_sinter_ref_length_mm'] = '8.06'
            s['gauge_sinter_ratio'] = '.8'
            s['gauge_sinter_sd'] = '.0005'
        (out, post) = analyze(rows, p, f)
        rep = out['posterior']
        self.assertLess(abs(rep['sinter_factor']['mean'] - 0.8), 1e-05)
        self.assertLess(abs(rep['scanner_bias']['mean'] - 0.006), 1e-05)
        self.assertIn('scanner_bias', out['monitor']['first_alarm_at_specimen'])
        self.assertNotIn('sinter', out['monitor']['first_alarm_at_specimen'])
        self.assertLess(rep['paired_sinter_metrology']['covariance'][0][1], 0)

    def test_paired_without_gauge_cannot_promote_ratio(self):
        p = verified(R / 'inputs/DEMO_PAIRED_PROFILE.json')
        f = json.loads((R / 'inputs/X1B_FROZEN_PREDICTIONS.json').read_text())
        rows = series(6100, p, f)
        (out, post) = analyze(rows, p, f)
        self.assertEqual(out['posterior']['sinter_factor']['observations'], 0)
        self.assertNotIn('sinter', out['monitor']['first_alarm_at_specimen'])

    def test_physical_sigma_requires_exact_bound_operator(self):
        p = verified(R / 'inputs/DEMO_PROFILE.json')
        f = json.loads((R / 'inputs/X1B_FROZEN_PREDICTIONS.json').read_text())
        groups = model_groups(f)
        p['stress_operator'] = {'status': 'UNVALIDATED', 'locator': 'synthetic', 'sha256': '1' * 64, 'm_grid': [], 'sigma0_per_lambda_MPa_per_N': []}
        with self.assertRaises(ValueError):
            JointPosterior(p, groups).report()

    def test_bad_source_sidecar_and_duplicate_rejected(self):
        import tempfile, csv
        from synthetic import save_series
        p = verified(R / 'inputs/DEMO_PROFILE.json')
        f = json.loads((R / 'inputs/X1B_FROZEN_PREDICTIONS.json').read_text())
        rows = series(6100, p, f)
        with tempfile.TemporaryDirectory() as tmp:
            (path, side) = save_series(Path(tmp) / 'a.csv', rows)
            ss = side.read_text().replace(rows[0]['force_trace_sha256'], '0' * 64)
            side.write_text(ss)
            with self.assertRaises(ValueError):
                read_csv(path, side)
            (path, side) = save_series(Path(tmp) / 'b.csv', rows + [rows[0]])
            with self.assertRaises(ValueError):
                read_csv(path, side)
