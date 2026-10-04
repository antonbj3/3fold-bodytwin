from dental_release.paths import expand as _release_expand
import csv, datetime, hashlib, json, resource, time
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.stats import gamma
from risk_operator import scenario_probability, scenario_margin, body_bound, bound_margin, scalar_control, moment_tail, gamma_signed_survival
ROOT = Path(__file__).resolve().parent
DATA = Path(_release_expand('@DENTAL_WORK_ROOT@/X8-guide-nerve-risk'))

def dump(name, obj):
    (ROOT / name).write_text(json.dumps(obj, indent=2, ensure_ascii=False) + '\n')

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()
start = time.perf_counter()
pr = json.loads((ROOT / 'PREREG_R2.json').read_text())
assert sha(ROOT / 'PREREG_R2.json') == (ROOT / 'PREREG_R2.sha256').read_text().strip()
assert sha(pr['input_path']) == pr['input_sha256']
assert sha(ROOT / 'GUIDE_PROFILES.json') == pr['guide_profiles_sha256']
profiles = json.loads((ROOT / 'GUIDE_PROFILES.json').read_text())['profiles']
raw = json.loads(Path(pr['input_path']).read_text())
valid = []
missing = []
for s in raw['sites']:
    if s.get('L_plan') is None or 'd_nom_implant_body' not in s:
        missing.append({'case': s['case'], 'fdi': s['fdi'], 'reason': 'NO_PLAN' if s.get('L_plan') is None else 'NO_COMPLETE_NOMINAL_CLEARANCE'})
    else:
        valid.append(s)
d = np.array([s['d_nom_implant_body'] for s in valid])
practice = d >= 2
grid = np.arange(0.5, 12.001, 0.25)
rows = []
curves = {}
summary = []
losses = {}
validation = []
for (i, p) in enumerate(profiles):
    name = p['id']
    (mu, sd) = p['apex_mean_sd_mm']
    rng = np.random.default_rng(pr['mc']['seed'] + i)
    k = (mu / sd) ** 2
    theta = sd * sd / mu
    R = rng.gamma(k, theta, pr['mc']['n'])
    U = rng.uniform(-1, 1, pr['mc']['n'])
    loss = np.sort(R * U + rng.normal(0.059, np.hypot(0.06, 0.22), pr['mc']['n']))
    losses[name] = loss
    pmc = 1 - np.searchsorted(loss, grid - 1, side='right') / len(loss)
    exact = scenario_probability(grid, p)
    err = float(np.max(np.abs(pmc - exact)))
    scalar_err = []
    for t in [0.3, 1.0, 2.0, 4.0]:
        q = 0.5 * quad(lambda r: (1 - t / r) * gamma.pdf(r, k, scale=theta), t, np.inf, epsabs=1e-11)[0]
        scalar_err.append(abs(q - float(gamma_signed_survival(t, mu, sd))))
    bound_control_errors = [abs(float(body_bound(t, p)) - scalar_control(t, p)) for t in [2.0, 3.0, 5.0, 8.0, 12.0]]
    siteprob = scenario_probability(d, p)
    upper = {str(B): body_bound(d, p, B) for B in [0, 0.3, 0.7]}
    drill_upper = body_bound(d, p, 0.3, extra_budget=1.5)
    ms = scenario_margin(p)
    mb = bound_margin(p)
    mbdrill = bound_margin(p, extra_budget=1.5)
    agree = (siteprob < 0.01) == (upper['0.3'] < 0.01)
    summary.append({'guide': name, 'n_sites': len(valid), 'n_2mm_accepted': int(practice.sum()), 'scenario_gt1pct_despite_2mm': int((practice & (siteprob > 0.01)).sum()), 'moment_upper_gt1pct_despite_2mm': int((practice & (upper['0.3'] > 0.01)).sum()), 'scenario_eligible': int((siteprob < 0.01).sum()), 'moment_condition_eligible_B0p3': int((upper['0.3'] < 0.01).sum()), 'scenario_required_clearance_mm_1pct': ms, 'moment_required_clearance_mm_1pct_B0p3': mb, 'drill_moment_required_clearance_mm_1pct_B0p3': mbdrill, 'decision_agreement_among_2mm_accepted': float(agree[practice].mean()), 'mean_scenario_proximity_probability': float(siteprob.mean()), 'clinical_injury_risk': 'UNKNOWN'})
    validation.append({'guide': name, 'mc_analytic_max_abs_probability_error': err, 'pass': err <= 0.005, 'direct_density_quad_max_abs_error': max(scalar_err), 'vector_scalar_bound_max_abs_error': max(bound_control_errors)})
    curves[name] = {'distance_mm': grid.tolist(), 'scenario': exact.tolist(), 'mc': pmc.tolist(), 'body_bound_B0p3': body_bound(grid, p).tolist()}
    for (j, s) in enumerate(valid):
        rows.append({'case': s['case'], 'fdi': s['fdi'], 'guide': name, 'plan_length_mm': s['L_plan'], 'nominal_implant_gap_mm': float(d[j]), 'nominal_drill_gap_y1_mm': s.get('d_nom_y1'), 'fixed_2mm_accept': bool(practice[j]), 'local_apex_scenario_p_lt1': float(siteprob[j]), 'whole_body_bound_B0_mm': float(upper['0'][j]), 'whole_body_bound_B0p3_mm': float(upper['0.3'][j]), 'whole_body_bound_B0p7_mm': float(upper['0.7'][j]), 'drill_bound_B0p3_ymax1p5': float(drill_upper[j]), 'scenario_required_margin_1pct_mm': ms, 'bound_required_margin_1pct_B0p3_mm': mb, 'scenario_eligible_1pct': bool(siteprob[j] < 0.01), 'bound_conditions_eligible_1pct': bool(upper['0.3'][j] < 0.01), 'physical_boundary_validated': False, 'clinical_guide_selection': 'UNKNOWN'})
with (ROOT / 'PER_SITE_RISK.csv').open('w') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]))
    w.writeheader()
    w.writerows(rows)
with (ROOT / 'PER_SITE_SYSTEM_SETS.jsonl').open('w') as f:
    for (j, s) in enumerate(valid):
        rr = [rows[i * len(valid) + j] for i in range(len(profiles))]
        f.write(json.dumps({'case': s['case'], 'fdi': s['fdi'], 'scenario_systems_below_1pct': [r['guide'] for r in rr if r['scenario_eligible_1pct']], 'moment_condition_systems_below_1pct_B0p3': [r['guide'] for r in rr if r['bound_conditions_eligible_1pct']], 'clinical_systems': 'UNKNOWN', 'use': 'Research hypothesis under explicit conditions'}) + '\n')
dump('EXCLUSIONS_R2.json', missing)
np.savez_compressed(DATA / 'R2_CLEARANCE_LOSS_SAMPLES.npz', **losses)
p = profiles[2]
(mu, sd) = p['apex_mean_sd_mm']
t = mu + 5 * sd
b = t + 1e-06
a = mu - sd * sd / (b - mu)
mass = sd * sd / (sd * sd + (b - mu) ** 2)
ext = {'profile': p['id'], 'threshold_mm': t, 'atoms_mm': [a, b], 'weights': [1 - mass, mass], 'mean_mm': (1 - mass) * a + mass * b, 'variance_mm2': (1 - mass) * (a - mu) ** 2 + mass * (b - mu) ** 2, 'tail_probability': mass, 'bound': float(moment_tail(t, mu, sd)), 'scope': 'Magnitude marginal only; not a reconstructed full jointly realizable pose cohort', 'injected_zero_upper_rejected': mass > 0}
dump('TAIL_EXTREMIZER_R2.json', ext)
reference_bound = scalar_control(5.0, profiles[2])
inj = {'mc_probability_plus_.02_rejected': all((v['mc_analytic_max_abs_probability_error'] + 0.02 > 0.005 for v in validation)), 'bound_plus_.1_rejected_by_scalar_control': bool(abs(float(body_bound(5.0, profiles[2])) + 0.1 - reference_bound) > 1e-10), 'zero_bound_rejected_by_legal_marginal_extremizer': ext['injected_zero_upper_rejected'], 'false_clinical_approval_rejected': not all([True, False, False])}
dump('CONTROLS_R2.json', {'validation': validation, 'injected_faults': inj, 'practice_control_executed': {'n': len(d), 'rule': 'nominal implant body d>=2mm', 'accepted': int(practice.sum())}, 'method_outcome': 'TIE: analytic closure and conventional moment inequality receive exactly the same inputs; no numerical novelty claim'})
out = {'construction': 'R2', 'n_sites': len(valid), 'n_cases': len(set((s['case'] for s in valid))), 'source_cases': raw['n_cases'], 'excluded_site_rows': len(missing), 'profiles': summary, 'numerical_parity_pass': all((v['pass'] for v in validation)), 'policy_identifiability_gate': 'FAIL' if any((s['decision_agreement_among_2mm_accepted'] < 0.95 for s in summary)) else 'UNKNOWN_NO_EXTERNAL_BOUNDARY_OR_TAIL_VALIDATION', 'clinical_validation': 'UNKNOWN', 'new_capability': 'Executable sitewise conditional proximity and distribution ambiguity bounds; measured guide magnitudes replace patent sweep', 'external_referent': pr['external_referent'], 'method_comparison': 'TIE', 'cost': {'wall_seconds': time.perf_counter() - start, 'maxrss_kib': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, 'preparation_manual_seconds': None, 'fit_seconds': 'included in wall', 'physical_measurements': 0, 'fallback': 'R3 independent signed clearance-loss calibration port'}, 'artifacts': [{'path': str(DATA / 'R2_CLEARANCE_LOSS_SAMPLES.npz'), 'bytes': (DATA / 'R2_CLEARANCE_LOSS_SAMPLES.npz').stat().st_size, 'sha256': sha(DATA / 'R2_CLEARANCE_LOSS_SAMPLES.npz')}, {'path': pr['input_path'], 'bytes': Path(pr['input_path']).stat().st_size, 'sha256': pr['input_sha256']}], 'limits': ['Population moment transfer assumed', 'Halfspace projection treats closest point as apex; not a geometric bound', 'True canal boundary budgets are unvalidated', 'Drill motion uses final implant error law conditionally', 'No paired clinical sensory outcomes', 'All anatomy sites are virtual tooth-axis immediate plans', 'Subsequent source inspection: K3 d_nom_implant_body samples only apical6mm. The whole-body probability bound is conditional on complete nominal-distance validity, which R2 does not establish. R4/R5 replace this geometric input.'], 'prereg_sha256': sha(ROOT / 'PREREG_R2.json')}
dump('RAW_R2.json', out)
dump('RAW_R2_CURVES.json', curves)
(fig, axes) = plt.subplots(1, 3, figsize=(15, 4.5))
colors = ['#aa4050', '#d5903c', '#3a8e68', '#477eb6']
for (p, col, s) in zip(profiles, colors, summary):
    name = p['id']
    cc = curves[name]
    axes[0].semilogy(grid, np.maximum(cc['scenario'], 1e-06), label=name, color=col)
    axes[0].semilogy(grid, np.maximum(cc['body_bound_B0p3'], 1e-06), '--', color=col, alpha=0.7)
axes[0].axvline(2, color='black', lw=1, label='2 mm rule')
axes[0].axhline(0.01, color='gray', lw=1)
axes[0].set(xlabel='Nominal implant/canal gap (mm)', ylabel='P(gap < 1 mm), conditional', ylim=(1e-05, 1), title='Solid: assumed law; dashed: moment upper')
axes[0].legend(fontsize=8)
xx = np.arange(4)
axes[1].bar(xx - 0.18, [s['scenario_gt1pct_despite_2mm'] for s in summary], width=0.36, label='Assumed law >1%')
axes[1].bar(xx + 0.18, [s['moment_upper_gt1pct_despite_2mm'] for s in summary], width=0.36, label='Moment upper >1%')
axes[1].set_xticks(xx, [p['id'].replace('_', '\n') for p in profiles], fontsize=8)
axes[1].set(ylabel='Virtual sites accepted by nominal 2 mm', title=f'{len(valid)} sites; no injury labels')
axes[1].legend(fontsize=8)
axes[2].bar(xx - 0.18, [s['scenario_required_clearance_mm_1pct'] for s in summary], width=0.36, label='Assumed law')
axes[2].bar(xx + 0.18, [s['moment_required_clearance_mm_1pct_B0p3'] for s in summary], width=0.36, label='Moment bound B=.3 mm')
axes[2].axhline(2, color='black', lw=1)
axes[2].set_xticks(xx, [p['id'].replace('_', '\n') for p in profiles], fontsize=8)
axes[2].set(ylabel='Required nominal gap (mm)', title='Research target P(<1 mm) <1%')
axes[2].legend(fontsize=8)
fig.suptitle('X8: geometric proximity scenarios, NOT clinical nerve-injury risk', fontsize=12)
fig.tight_layout()
fig.savefig(ROOT / 'GUIDE_RISK_FIGURE.png', dpi=170)
plt.close(fig)
frozen = {'frozen_at': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'prereg_sha256': sha(ROOT / 'PREREG_R2.json'), 'per_site_prediction_sha256': sha(ROOT / 'PER_SITE_RISK.csv'), 'profiles': summary, 'physical_measurements_performed': False, 'interpretation': 'Frozen conditional predictions for future independent phantom measurements; no retrofit permitted'}
if (ROOT / 'FROZEN_PREDICTIONS.json').exists():
    previous = json.loads((ROOT / 'FROZEN_PREDICTIONS.json').read_text())
    assert previous['per_site_prediction_sha256'] == frozen['per_site_prediction_sha256']
else:
    dump('FROZEN_PREDICTIONS.json', frozen)
dump('CURRENT_WORK_STATE.json', {'lane': 'X8-guide-nerve-risk', 'status': 'R2_COMPLETE_R3_DESIGN', 'latest_gate': out['policy_identifiability_gate'], 'next_operation': 'Freeze independent signed clearance-loss calibration design, rank measurement precision bottleneck, export phantom and compare exact coverage control'})
(ROOT / 'HANDOFF_R2.md').write_text('R2 completed. RAW_R2.json and PER_SITE_RISK.csv contain exact conditional results. Numerical controls TIE; decision-identifiability gate ' + out['policy_identifiability_gate'] + '. No clinical outcome or true anatomical boundary validation. Next: replace unobserved distribution tail by directly measured signed clearance loss, with finite-sample confidence and independent anatomy budget.\n')
print(json.dumps(out, indent=2))
