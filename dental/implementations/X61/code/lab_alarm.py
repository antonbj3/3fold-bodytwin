"""Research K50: fixed-law alarms and separately evolving conditional posteriors.
No physical absolute-force prediction or clinical recommendation is asserted.
"""
from __future__ import annotations
from dental_release.paths import expand as _release_expand
import argparse, csv, hashlib, importlib.util, json, math, re, sys, time
from io import StringIO
from pathlib import Path
from types import SimpleNamespace
import numpy as np
from numpy.polynomial.hermite import hermgauss
from scipy.special import logsumexp
from scipy.stats import norm
R = Path(__file__).resolve().parents[1]
ENGINE = Path(_release_expand('@DENTAL_EXTERNAL_ROOT@/projects/graph_workspace/graph-night-integration-20260923/src/graph_engine'))

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def write(p, o):
    Path(p).write_text(json.dumps(o, indent=2, allow_nan=False, default=lambda x: x.tolist() if isinstance(x, np.ndarray) else x.item()) + '\n')

def csv_rows(p):
    return list(csv.DictReader(StringIO(Path(p).read_text())))

def verified(p):
    p = Path(p)
    expected = p.with_name(p.name + '.sha256')
    if not expected.exists() or sha(p) != expected.read_text().strip():
        raise ValueError(f'Frozen artifact hash drift or missing seal: {p}')
    return json.loads(p.read_text())

def load_engine(name):
    spec = importlib.util.spec_from_file_location('x61_' + name, ENGINE / (name + '.py'))
    m = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = m
    spec.loader.exec_module(m)
    return m
Alarm = load_engine('alarm').Alarm
DC = load_engine('decision_cert')
(_HX, _HW) = hermgauss(16)
_HW = _HW / math.sqrt(math.pi)

def force_cdf(f, m, scale, noise):
    z = m * (np.log(f) - np.log(scale) - math.sqrt(2) * noise * _HX)
    return float(np.dot(_HW, -np.expm1(-np.exp(np.clip(z, -700, 30)))))

def log_force_likelihood(f, m, scale, noise):
    z = m[..., None] * (math.log(f) - np.log(scale)[..., None] - math.sqrt(2) * noise * _HX)
    return logsumexp(np.log(_HW) + np.log(m)[..., None] + z - np.exp(np.clip(z, -700, 700)), axis=-1)

def model_groups(frozen):
    fam = frozen['matched_family']
    ref = next((d for d in fam if d['design'] == 'D1'))
    ar = ref['FE_absolute_diagnostic_N']['30'] / ref['FE_absolute_diagnostic_N']['0']
    out = {}
    for d in fam:
        rr = d['FE_conditional_ratios']
        q = rr['30'] / rr['0']
        for angle in (0, 30):
            out[d['design'], angle] = dict(base_factor=rr['0'] * (ar if angle else 1.0), Q=q, interaction=bool(angle), gap_um=d['gap']['marginal_mean_um'])
    return out

def numeric(row, key, positive=False, nonnegative=False):
    try:
        x = float(row[key])
    except (KeyError, ValueError, TypeError):
        raise ValueError(f"{row.get('specimen_id', '?')}: missing/invalid {key}")
    if not math.isfinite(x) or (positive and x <= 0) or (nonnegative and x < 0):
        raise ValueError(f"{row.get('specimen_id', '?')}: invalid {key}")
    return x

def read_csv(path, sidecar=None):
    rows = csv_rows(path)
    sides = {}
    if not rows:
        raise ValueError('No specimen rows; blank protocol template is not measured data')
    if sidecar:
        for s in csv_rows(sidecar):
            sid = s.get('specimen_id', '')
            if not sid or sid in sides:
                raise ValueError('Empty or duplicate sidecar specimen_id')
            sides[sid] = s
    allowed = set((r['specimen_id'] for r in csv_rows(R / 'inputs/LAB_MEASUREMENTS_TEMPLATE.csv')))
    seen = set()
    for row in rows:
        sid = row.get('specimen_id', '')
        if sid not in allowed or sid in seen:
            raise ValueError(f'Unknown or duplicate specimen_id: {sid}')
        seen.add(sid)
        a = numeric(row, 'load_angle_deg')
        d = row.get('design_id')
        if a not in (0, 30) or d not in ('D1', 'M1', 'M2') or (not sid.startswith(f'{d}_{int(a)}_')):
            raise ValueError('Specimen identity/design/angle mismatch')
        numeric(row, 'fracture_force_N', positive=True)
        numeric(row, 'MG_distance_mean_um', nonnegative=True)
        for key in ('source_scan_sha256', 'force_trace_sha256'):
            if row.get(key) and (not re.fullmatch('[0-9a-f]{64}', row[key])):
                raise ValueError(f'Invalid source hash: {key}')
        s = sides.get(sid)
        if s:
            if s.get('post_sinter_ref_length_mm'):
                numeric(s, 'pre_sinter_ref_length_mm', positive=True)
                numeric(s, 'post_sinter_ref_length_mm', positive=True)
                numeric(s, 'sinter_ratio_sd', positive=True)
            if s.get('assembled_film_mean_um'):
                numeric(s, 'assembled_film_mean_um', nonnegative=True)
                numeric(s, 'assembled_film_sd_um', positive=True)
            if s.get('gauge_sinter_ratio'):
                numeric(s, 'gauge_sinter_ratio', positive=True)
                numeric(s, 'gauge_sinter_sd', positive=True)
                if not s.get('post_sinter_ref_length_mm'):
                    raise ValueError('Independent gauge must be paired with scanner ratio')
            if s.get('trace_peak_force_N'):
                numeric(s, 'trace_peak_force_N', positive=True)
                if s.get('source_trace_sha256') != row.get('force_trace_sha256') or not row.get('force_trace_sha256'):
                    raise ValueError('Force trace provenance mismatch')
        row['_sidecar'] = s or {}
    if set(sides) - seen:
        raise ValueError('Orphan sidecar specimens not present in input CSV')
    return rows

class GaussianPort:

    def __init__(self, mean, sd):
        self.prec = 1 / sd ** 2
        self.h = mean * self.prec
        self.n = 0

    def update(self, x, sd):
        self.prec += 1 / sd ** 2
        self.h += x / sd ** 2
        self.n += 1

    def report(self, unit):
        mu = self.h / self.prec
        sd = 1 / math.sqrt(self.prec)
        return dict(mean=mu, sd=sd, credible_95=[mu - 1.959963984540054 * sd, mu + 1.959963984540054 * sd], observations=self.n, unit=unit, status='UPDATED_CONDITIONAL' if self.n else 'UNKNOWN_NO_MEASUREMENT_PRIOR_ONLY', resolution='POPULATION')

class GapPort:
    """Unsigned distance observation: retain both latent normal signs."""

    def __init__(self, mean, sd):
        self.x = np.linspace(-80, 80, 321)
        self.logw = -0.5 * ((self.x - mean) / sd) ** 2
        self.n = 0

    def update(self, y, g0, sd):
        mu = g0 + self.x
        self.logw += np.logaddexp(norm.logpdf(y, mu, sd), norm.logpdf(-y, mu, sd))
        self.logw -= logsumexp(self.logw)
        self.n += 1

    def report(self, unit):
        w = np.exp(self.logw - logsumexp(self.logw))
        cdf = np.cumsum(w)
        mu = float(w @ self.x)
        return dict(mean=mu, sd=math.sqrt(float(w @ (self.x - mu) ** 2)), credible_95=[float(self.x[min(np.searchsorted(cdf, q), len(self.x) - 1)]) for q in [0.025, 0.975]], observations=self.n, unit=unit, status='CONDITIONAL_FOLDED_NORMAL_GRID' if self.n else 'UNKNOWN_NO_MEASUREMENT_PRIOR_ONLY', resolution='POPULATION', grid_boundary_mass=float(w[0] + w[-1]))

class JointPosterior:

    def __init__(self, p, groups):
        pr = p['priors']
        self.groups = groups
        self.p = p
        mg = np.linspace(*pr['m_grid'])
        lg = p['null']['lambda_reference_N'] * np.linspace(*pr['lambda_ratio_grid'])
        (self.m, self.lam, self.form) = np.meshgrid(mg, lg, np.array([0, 1]), indexing='ij')
        self.m = self.m.ravel()
        self.lam = self.lam.ravel()
        self.form = self.form.ravel()
        fp = pr['model_probabilities']
        self.initial = np.log(np.where(self.form == 1, fp['FE'], fp['separable'])) - math.log(len(mg) * len(lg))
        self.logw = self.initial.copy()
        self.n = 0
        self.gap = GapPort(pr['dry_gap_bias_mean_um'], pr['dry_gap_bias_sd_um'])
        self.sinter = GaussianPort(pr['sinter_mean'], pr['sinter_sd'])
        self.cement = GaussianPort(pr['cement_mean_um'], pr['cement_sd_um'])
        self.paired = None
        if p.get('paired_metrology'):
            from paired_metrology import PairedPort
            self.paired = PairedPort(p['paired_metrology']['prior_mean'], p['paired_metrology']['prior_sd'])
        self.accepted = []

    def factor(self, key):
        g = self.groups[key]
        return g['base_factor'] * np.where(self.form == 1, g['Q'] if g['interaction'] else 1.0, 1.0)

    def update_force(self, row):
        key = (row['design_id'], int(float(row['load_angle_deg'])))
        f = float(row['fracture_force_N'])
        self.logw += log_force_likelihood(f, self.m, self.lam * self.factor(key), self.p['null']['force_log_noise_sd'])
        self.logw -= logsumexp(self.logw)
        self.n += 1
        self.accepted.append((key, f))

    def update_other(self, row):
        g = self.groups[row['design_id'], int(float(row['load_angle_deg']))]
        n = self.p['null']
        self.gap.update(float(row['MG_distance_mean_um']), g['gap_um'], n['dry_gap_sd_um'])
        s = row.get('_sidecar', {})
        if self.paired:
            if s.get('gauge_sinter_ratio'):
                self.paired.update(float(s['post_sinter_ref_length_mm']) / float(s['pre_sinter_ref_length_mm']), float(s['gauge_sinter_ratio']), float(s['sinter_ratio_sd']), float(s['gauge_sinter_sd']), self.p['paired_metrology']['cross_error_correlation'])
        elif s.get('post_sinter_ref_length_mm'):
            self.sinter.update(float(s['post_sinter_ref_length_mm']) / float(s['pre_sinter_ref_length_mm']), float(s['sinter_ratio_sd']))
        if s.get('assembled_film_mean_um'):
            self.cement.update(float(s['assembled_film_mean_um']), float(s['assembled_film_sd_um']))

    def weights(self):
        return np.exp(self.logw - logsumexp(self.logw))

    def quantity(self, x, unit):
        w = self.weights()
        order = np.argsort(x)
        v = x[order]
        cdf = np.cumsum(w[order])
        mu = float(w @ x)
        var = float(w @ (x - mu) ** 2)
        return dict(mean=mu, sd=math.sqrt(var), credible_95=[float(v[min(np.searchsorted(cdf, q), len(v) - 1)]) for q in [0.025, 0.975]], unit=unit, resolution='POPULATION', grid_boundary_mass=float(w[(x == x.min()) | (x == x.max())].sum()))

    def report(self):
        w = self.weights()
        pf = float(w[self.form == 1].sum())
        cov = float(w @ ((self.m - w @ self.m) * (self.lam - w @ self.lam)))
        cert = DC.certify(DC.Decision(name='FE family form supported in this finite model space', kind='chain', nodes=('FE_form',), alpha=0.05), SimpleNamespace(unlock_graph=None, p_holds={'FE_form': pf}))
        out = dict(force_observations=self.n, m=self.quantity(self.m, '1'), sigma0_force_equivalent=self.quantity(self.lam, 'N'), sigma0_material_MPa=dict(status='UNKNOWN_MISSING_M_DEPENDENT_STRESS_HAZARD', resolution='PHENOMENOLOGICAL', replacement_measurement='Independently validated geometry/contact/support stress hazard I_j(m) in reference-area convention, plus matched strength specimens'), model_probabilities={'FE': pf, 'separable': 1 - pf}, model_certificate=cert, joint_m_scale_covariance_N=cov, sinter_factor=self.sinter.report('post/pre length'), cement_film_mean=self.cement.report('um'), dry_marginal_gap_bias=self.gap.report('um'), limitations='Conditional finite prior grid and two selected force forms; physical stress/gap/contact closures unvalidated, no grid-truncation enclosure.')
        if self.paired:
            out['paired_sinter_metrology'] = self.paired.report()
            out['sinter_factor'] = out['paired_sinter_metrology']['sinter_factor']
            out['scanner_bias'] = out['paired_sinter_metrology']['scanner_bias']
        op = self.p.get('stress_operator')
        if op:
            if op.get('status') != 'INDEPENDENTLY_VALIDATED' or not op.get('locator') or (not re.fullmatch('[0-9a-f]{64}', op.get('sha256', ''))):
                raise ValueError('Stress operator is not independently bound')
            ms = np.array(op['m_grid'])
            cs = np.array(op['sigma0_per_lambda_MPa_per_N'])
            uq = np.unique(self.m)
            if not np.array_equal(ms, uq) or cs.shape != ms.shape or (not np.all(np.isfinite(cs) & (cs > 0))):
                raise ValueError('Stress operator must cover exact posterior m grid without unbounded interpolation')
            out['sigma0_material_MPa'] = self.quantity(self.lam * cs[np.searchsorted(ms, self.m)], 'MPa for operator reference area')
        return out

class Monitor:

    def __init__(self, p, groups):
        self.p = p
        self.groups = groups
        self.alarms = {k: Alarm(alpha=a) for (k, a) in p['alpha_by_channel'].items() if k != 'measurement_record'}
        self.t = 0
        self.first_hits = {}
        self.record_faults = []
        K = p['bins']
        self.centers = (np.arange(K) + 0.5) / K

    def score_category(self, name, u, stat='scale', sign=1):
        K = self.p['bins']
        k = min(max(int(u * K), 0), K - 1)
        vals = self.centers - 0.5
        if stat == 'shape':
            vals = vals ** 2
        vals = sign * vals
        mu = float(vals.mean())
        v = float(vals[k])
        outcomes = [(1 / K, float(x)) for x in vals]
        (w, hit) = self.alarms[name].update(mu, v, outcomes)
        if self.alarms[name].hit and name not in self.first_hits:
            self.first_hits[name] = self.t
        return dict(pit=float(u), category=k, wealth=float(w), ever_hit=self.alarms[name].hit, alarm_now=hit)

    def step(self, row):
        self.t += 1
        sid = row['specimen_id']
        g = self.groups[row['design_id'], int(float(row['load_angle_deg']))]
        null = self.p['null']
        s = row.get('_sidecar', {})
        events = []
        scores = {}
        protocol = []
        for (col, expected, tol) in [('die_E_MPa', 18000, 1800), ('indenter_diameter_mm', 5, 0.05), ('crosshead_mm_min', 0.5, 0.025)]:
            try:
                v = float(row.get(col, ''))
            except ValueError:
                v = float('nan')
            if not math.isfinite(v):
                protocol.append(col + ':UNKNOWN')
            elif abs(v - expected) > tol:
                protocol.append(col + ':OUT_OF_PROTOCOL')
        for col in ('material_batch', 'cement_batch', 'die_material_batch'):
            if not row.get(col):
                protocol.append(col + ':UNKNOWN')
        if protocol:
            events.append(dict(kind='PROTOCOL', cause_candidates=['support/contact/acquisition'], details=protocol, false_alarm_claim=False))
        incompatible = row.get('failure_mode') in ('contact_damage', 'cement_debonding', 'die_fracture') or row.get('fracture_origin') in ('contact_zone', 'cement_interface', 'die')
        mechanism_unknown = row.get('failure_mode') != 'crown_tensile_fracture' or row.get('fracture_origin') != 'intaglio_tensile_zone'
        if incompatible:
            events.append(dict(kind='MECHANISM_FALSIFIER', cause_candidates=['contact/cement/support_model'], false_alarm_claim=False, action='Retain specimen and report informative alternate failure; do not fit it as tensile strength'))
        record_bad = False
        if s.get('trace_peak_force_N'):
            f = float(row['fracture_force_N'])
            ft = float(s['trace_peak_force_N'])
            noise = self.p['record_log_noise_floor']
            err = abs(math.log(f / ft))
            pval = float(2 * norm.sf(err / noise))
            budget = self.p['alpha_by_channel']['measurement_record'] / (self.t * (self.t + 1))
            record_bad = err > max(self.p['record_absolute_floor_N'] / ft, noise * norm.isf(budget / 2))
            scores['measurement_record'] = dict(log_error=err, two_sided_p=pval, spending_alpha=budget, rejected=bool(record_bad), scope='Recorded maximum vs source trace; sensor calibration fault remains unidentifiable')
            if record_bad:
                self.record_faults.append(sid)
                self.first_hits.setdefault('measurement_record', self.t)
                events.append(dict(kind='MEASUREMENT_RECORD_ALARM', cause_candidates=['measurement_record'], specimen_id=sid, action='Reconcile CSV with trace; keep raw record; do not update strength with disputed endpoint'))
        if not record_bad:
            factor = g['base_factor'] * (g['Q'] if null['model'] == 'FE' and g['interaction'] else 1.0)
            u = force_cdf(float(row['fracture_force_N']), null['m'], null['lambda_reference_N'] * factor, null['force_log_noise_sd'])
            scores['material_scale'] = self.score_category('material_scale', u)
            scores['material_shape'] = self.score_category('material_shape', u, 'shape')
            if g['interaction'] and g['Q'] != 1:
                scores['model_form'] = self.score_category('model_form', u, 'scale', sign=1 if math.log(g['Q']) > 0 else -1)
        y = float(row['MG_distance_mean_um'])
        mu = g['gap_um'] + null['dry_gap_bias_um']
        sd = null['dry_gap_sd_um']
        u = norm.cdf((y - mu) / sd) - norm.cdf((-y - mu) / sd)
        scores['dry_gap'] = self.score_category('dry_gap', u)
        if self.p.get('paired_metrology'):
            if s.get('gauge_sinter_ratio'):
                val = float(s['post_sinter_ref_length_mm']) / float(s['pre_sinter_ref_length_mm'])
                gauge = float(s['gauge_sinter_ratio'])
                ss = float(s['sinter_ratio_sd'])
                gs = float(s['gauge_sinter_sd'])
                rho = self.p['paired_metrology']['cross_error_correlation']
                ds = math.sqrt(ss ** 2 + gs ** 2 - 2 * rho * ss * gs)
                scores['sinter'] = self.score_category('sinter', norm.cdf((gauge - null['sinter_factor']) / gs))
                scores['scanner_bias'] = self.score_category('scanner_bias', norm.cdf((val - gauge - self.p['paired_metrology']['scanner_bias_null']) / ds))
        elif s.get('post_sinter_ref_length_mm'):
            val = float(s['post_sinter_ref_length_mm']) / float(s['pre_sinter_ref_length_mm'])
            u = norm.cdf((val - null['sinter_factor']) / float(s['sinter_ratio_sd']))
            scores['sinter'] = self.score_category('sinter', u)
        if s.get('assembled_film_mean_um'):
            u = norm.cdf((float(s['assembled_film_mean_um']) - null['cement_film_mean_um']) / float(s['assembled_film_sd_um']))
            scores['cement_film'] = self.score_category('cement_film', u)
        for (k, v) in scores.items():
            if v.get('ever_hit') and self.first_hits[k] == self.t:
                events.append(dict(kind='FROZEN_LAW_ALARM', channel=k, cause_candidates={'material_scale': ['material_strength', 'stress/contact/support_scale', 'force_sensor_bias'], 'material_shape': ['flaw_population', 'mixed_failure_modes', 'measurement_noise'], 'model_form': ['geometry_angle_response', 'contact/support_model'], 'dry_gap': ['seating', 'sinter', 'scanner_bias'], 'sinter': ['sinter_or_length_metrology'] if not self.p.get('paired_metrology') else ['sinter_given_independent_gauge_calibration'], 'scanner_bias': ['scanner_vs_reference_bias'], 'cement_film': ['cementation_or_film_metrology']}[k]))
        return dict(index=self.t, specimen_id=sid, scores=scores, events=events, force_update_allowed=not record_bad and (not protocol) and (not mechanism_unknown))

    def state(self):
        return dict(first_alarm_at_specimen=self.first_hits, any_alarm=bool(self.first_hits), measurement_record_specimens=self.record_faults, channel_states={k: a.state() for (k, a) in self.alarms.items()})

def analyze(rows, p, frozen, posterior=True):
    groups = model_groups(frozen)
    mo = Monitor(p, groups)
    post = JointPosterior(p, groups) if posterior else None
    timeline = []
    batches = set()
    for row in rows:
        batches.add((row.get('material_batch', ''), row.get('cement_batch', ''), row.get('die_material_batch', '')))
        e = mo.step(row)
        if len(batches) > 1:
            e['events'].append(dict(kind='BATCH_CHANGE', cause_candidates=['new_batch'], action='Separate posterior and new frozen profile required'))
            e['force_update_allowed'] = False
        if post:
            if e['force_update_allowed']:
                post.update_force(row)
            if len(batches) == 1:
                post.update_other(row)
            e['model_probability_FE_after'] = post.report()['model_probabilities']['FE'] if mo.first_hits or mo.t == len(rows) else None
        timeline.append(e)
    report = dict(schema='K50-lab-alarm-v1', claim_type='capability', input_specimens=len(rows), resolution='PER_TOOTH', timescale='HANDOVER', monitor=mo.state(), timeline=timeline, posterior=post.report() if post else None, real_lab_alpha_validity='UNKNOWN_UNTIL_INDEPENDENT_NULL_AND_NOISE_CALIBRATION', profile_status=p['status'], localization='COMPETING_CAUSES_UNLESS_ORTHOGONAL_MEASUREMENT_DISCRIMINATES', physical_validation='UNKNOWN_NO_MATCHED_LAB_DATA', dropout=dict(received=len(rows), invalid_rejected=0, force_updates_withheld=sum((not e['force_update_allowed'] for e in timeline)), absent_planned=72 - len(rows)), method='Score frozen prediction before updating posterior; no posterior-driven alarm reset', engine_reuse={name: dict(path=str(ENGINE / (name + '.py')), sha256=sha(ENGINE / (name + '.py'))) for name in ('alarm', 'decision_cert')})
    return (report, post)

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('csv')
    ap.add_argument('--sidecar')
    ap.add_argument('--profile')
    ap.add_argument('--predictions', default=str(R / 'inputs/X1B_FROZEN_PREDICTIONS.json'))
    ap.add_argument('--output', required=True)
    args = ap.parse_args()
    rows = read_csv(args.csv, args.sidecar)
    f = json.loads(Path(args.predictions).read_text())
    if args.profile:
        p = verified(args.profile)
        if sha(args.predictions) != p['predictions_sha256']:
            raise ValueError('Prediction/profile binding drift')
        (out, post) = analyze(rows, p, f)
    else:
        out = dict(schema='K50-lab-alarm-v1', claim_type='capability', input_specimens=len(rows), status='UNKNOWN_NO_FROZEN_LAB_LIKELIHOOD', posterior={'Weibull_m': 'UNKNOWN', 'sigma0_MPa': 'UNKNOWN', 'cement_film': 'UNKNOWN_NO_ASSEMBLED_FILM', 'sinter_factor': 'UNKNOWN_NO_PAIRED_LENGTH'}, next_required='Supply a hash-sealed lab profile fixed before outcomes, independently estimated measurement errors and separate process/trace channels', physical_validation='UNKNOWN', incompatible_specimens=[r['specimen_id'] for r in rows if r.get('failure_mode') in ('contact_damage', 'cement_debonding', 'die_fracture') or r.get('fracture_origin') in ('contact_zone', 'cement_interface', 'die')])
        post = None
    out['input_csv_sha256'] = sha(args.csv)
    out['sidecar_sha256'] = sha(args.sidecar) if args.sidecar else None
    out['predictions_sha256'] = sha(args.predictions)
    dest = Path(args.output)
    if dest.exists():
        raise ValueError('Output already exists; use a new path to preserve previous gates')
    dest.parent.mkdir(parents=True, exist_ok=True)
    write(dest, out)
    if post:
        np.savez_compressed(dest.with_suffix('.posterior.npz'), m=post.m, lambda_reference_N=post.lam, form=post.form, weights=post.weights())
        out['joint_posterior_file'] = str(dest.with_suffix('.posterior.npz'))
        out['joint_posterior_sha256'] = sha(dest.with_suffix('.posterior.npz'))
        write(dest, out)
    print(json.dumps({'output': str(dest), 'specimens': len(rows), 'alarm': out.get('monitor', {}).get('any_alarm'), 'status': out.get('status', out.get('physical_validation'))}))
if __name__ == '__main__':
    try:
        main()
    except (ValueError, KeyError) as e:
        print(str(e), file=sys.stderr)
        sys.exit(2)
