"""Offline biological/material observation ports. No clinical decision engine."""
from pathlib import Path
import json, math, hashlib, time, datetime, resource, itertools, copy
import numpy as np
from scipy.stats import t as student_t
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from source_io import Sources, ROOT, IDS, nums, text
START = time.perf_counter()

def stamp():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def write(p, obj):
    p = ROOT / p
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(obj, ensure_ascii=False, indent=2, allow_nan=False) + '\n')

def state(round, last, next):
    write('CURRENT_WORK_STATE.json', dict(lane='X80-biocompatibility', status='RUNNING', milestone=round, last_gate=last, current_operation=next, next_operation=next, updated_at=stamp(), threads=4, intermediate_limit_bytes=3000000000))

def referent(id, locator, quantity, refutes):
    m = next((m for m in S.manifest if m['id'] == id))
    return dict(kind='independent_measurement', locator='https://doi.org/' + m['doi'] + '; ' + locator, local_path=m['path'], source_sha256=m['sha256'], compared_quantity=quantity, refutes_us=bool(refutes))

def finish_round(r, obj, narrative, next):
    obj['claim_type'] = 'information_link'
    obj['review_state'] = 'PENDING_INDEPENDENT_REVIEW'
    write('rounds/' + r + '/results.json', obj)
    (ROOT / 'rounds' / r / 'HANDOFF.md').write_text(narrative + '\n\nNext construction: ' + next + '\n')
    state(r, obj['outcome'], next)
S = Sources()
M = S.extract()
PREREGS = ['PREREG_R1_TITANIUM.json', 'PREREG_R2_RESIN_KINETICS_v2.json', 'PREREG_R3_SPATIAL_PROCESS_v2.json', 'PREREG_R4_MATERIAL_DECISION.json']
contracts = {p: json.loads((ROOT / p).read_text()) for p in PREREGS}
assert all((c['claim_type'] == 'information_link' for c in contracts.values()))
assert (ROOT / 'DECOMPOSITION.json').is_file()
state('R1', 'PREREG_FROZEN', 'Run titanium source, tissue and recovery comparison')
mass = M['retained Ti mass']['values']
recovery = M['measurement-recovery']['values']
curr = M['corrosion current']['values']
org = M['organ Ti at30d']['values']
blood = M['blood Ti at30d']['values']
reductions = {k: 1 - v['mean'] / mass['Control +']['mean'] for (k, v) in mass.items() if k not in ['Control +', 'Control −']}
organ_delta = {k: v['experimental_median'] - v['control_median'] for (k, v) in org.items()}
blood_delta = blood['after']['experimental_median'] - blood['after']['control_median']
k = math.log(2) / 10.0
total_a = sum([1.0])
total_b = sum([1.0])
end_a = 1.0 * math.exp(-k * 30)
end_b = 1.0
assert abs(end_a + (1 - end_a) - total_a) < 1e-15
r1 = dict(outcome='PROXY_REFUTED_EXPOSURE_UNKNOWN', resolution='POPULATION', current_ratio_H_over_R=curr['H']['mean'] / curr['R']['mean'], current_proxy_order='H>R', independent_release_order='R>H', current_proxy_pass=False, release_upper_bound_ug_in100ml_30d=6 * 0.1, release_average_upper_bound_ug_per_day=6 * 0.1 / 30, release_rate_caveat='Endpoint upper bound divided by duration; not instantaneous rate or monthly in vivo source.', release_current_unit_conflict='Table5 µA/cm² versus abstract µA/mm². No absolute Faraday dose produced.', implantoplasty_retained_mass_ug=mass, retained_mass_relative_reductions=reductions, aas_mass_recovery_fraction=[x / 19000 for x in recovery], aas_caveat='No recovery correction extrapolated to particles; source digestion destroys size/speciation.', blood_difference_ng_g=blood_delta, organ_differences_ng_g=organ_delta, blood_proxy_pass=all((np.sign(blood_delta) == np.sign(x) for x in organ_delta.values())), tissue_concentration='UNKNOWN', inflammatory_assay={'locator': 'https://doi.org/10.3390/ijms241411644; Figure7 and cytokine Results', 'species': 'TiO2; <100nm nanoparticles or <5µm microparticles', 'tested_doses_ug_ml': [20, 100], 'duration_h': 24, 'cells': 'THP1-derived M0 macrophages', 'observation': '100µg/mL increased IL1β,IL8,IL18; TNFα did not significantly change. This is an observed-effect dose, not a universal threshold or cytotoxicity dose.', 'resolution': 'POPULATION'}, sufficiency={'summary': 'total injected elemental mass', 'identity_error': abs(total_a - total_b), 'bit_identical': total_a.hex() == total_b.hex(), 'downstream_difference': end_b - end_a, 'downstream_unit': 'fixture normalized inventory', 'early_final': end_a, 'late_final': end_b, 'assumed_half_life_days': 10, 'resolution': 'PHENOMENOLOGICAL', 'minimum_extension': 'Time-weighted retained inventory for a single terminal query; species inventories and calibrated transport for future queries', 'external_referent': {'kind': 'our_own_fixture', 'locator': 'code/run_demo.py:R1 pulse mass balance', 'compared_quantity': 'inventory under first-order removal, not dental tissue measurement', 'refutes_us': True}}, upstream_links={'X56': 'Rejected as absolute-release input: bone/system micromotion not local implant–abutment slip', 'X65': 'Accepted cycle clock metadata only; no tribochemical coefficient or wet history', 'microgap': 'UNKNOWN source; gap geometry alone supplies no release rate'}, external_referent=referent('PMC10223376', 'Table5, Results Figure5 description', 'current versus ion-release ranking', True), external_referents=[referent('PMC10415425', 'Tables2–3', 'retained Ti µg and recovery', False), referent('PMC10087269', 'Tables2–3', 'blood vs organ Ti difference at30d', True)], control='Direct table reading reproduces all mass and current inputs exactly; exact scalar propagator matches vector inventory. No algorithm advantage claimed.', next_measurement='Same wet implant–abutment assembly: measure relative slip, collected fractionated Ti/Al/V and particle size over known cycles; multiple collection times, blanks and recovery standards; measure local volume and retention before tissue prediction.')
finish_round('R1', r1, 'A measured release ordering and rat tissue endpoints refute two source proxies. Source means are population observations. The patient dose chain fails at speciation, source rate, local volume and clearance; no tissue threshold is evaluated. Retained particle mass and AAS recovery are replayed from original tables. Source unit discrepancy remains quarantined.', 'R2: test whether a conditional DC time model and leave-one-study-out can support an unvisited recipe')
dc = M['post-cure DC']['values']
y = np.array(dc['dc_percent'])
tt = np.array(dc['time_min'])

def fit_dc(values):
    (y0, y15, y45) = [float(v) for v in values]
    ratio = (y45 - y0) / (y15 - y0)
    q = (-1 + math.sqrt(4 * ratio - 3)) / 2
    if not 0 < q < 1:
        raise ValueError('No positive saturating first-order solution')
    rate = -math.log(q) / 15
    infinity = y0 + (y15 - y0) / (1 - q)
    if not y0 < infinity <= 100:
        raise ValueError('Invalid conversion plateau')
    return (y0, infinity, rate)

def evaluate(p, t):
    (a, b, k) = p
    return a + (b - a) * -math.expm1(-k * t)
p = fit_dc(y)
replay = [evaluate(p, float(t)) for t in tt]
rounding = np.array([evaluate(fit_dc(y + np.array(e)), 30) for e in itertools.product([-0.05, 0.05], repeat=3)])
pred = evaluate(p, 30)
frozen = {'schema': 'x80-frozen-prediction-v1', 'prediction_type': 'model_conditional_not_validated', 'recipe': {'source_resin': 'PMC10892052 model resin', 'oven': 'Triad2000', 'time_min': 30, 'temperature_C': None}, 'quantity': 'DC_percent', 'resolution': 'POPULATION', 'predicted_value': pred, 'rounding_scenarios_min_max': [float(rounding.min()), float(rounding.max())], 'statistical_prediction_interval': None, 'model_error': 'UNKNOWN', 'cytotoxicity_percent': 'UNKNOWN', 'parameters': {'DC0': p[0], 'DCinfinity': p[1], 'k_per_min': p[2]}, 'input_sha256': sha(ROOT / 'raw/measurements.json'), 'code_sha256': sha(Path(__file__)), 'prereg_sha256': sha(ROOT / PREREGS[1]), 'before_physical_measurement': True, 'measurement_performed': False}
payload = json.dumps(frozen, sort_keys=True, ensure_ascii=False).encode()
frozen['payload_sha256'] = hashlib.sha256(payload).hexdigest()
frozen['frozen_at'] = stamp()
predpath = ROOT / 'FROZEN_PREDICTIONS.json'
if predpath.exists():
    old = json.loads(predpath.read_text())
    assert old['input_sha256'] == frozen['input_sha256'] and old['prereg_sha256'] == frozen['prereg_sha256'] and (old['predicted_value'] == frozen['predicted_value']), 'Frozen input or prediction changed: version before rerun'
    assert old['code_sha256'] == sha(ROOT / 'rounds/R4/first_execution/run_demo.py'), 'Initial prediction code snapshot hash mismatch'
else:
    predpath.write_text(json.dumps(frozen, indent=2, ensure_ascii=False) + '\n')
studies = [{'source': 'PMC10892052', 'compatible_time_series': True, 'reason': 'Only3 prose points; temperature/SD missing'}, {'source': 'PMC10780983', 'compatible_time_series': False, 'reason': 'Different polyurethane aligner; Table1 has no quantitative time series or same oven'}, {'source': 'PMC10182118', 'compatible_time_series': False, 'reason': 'Experimental composites; Fig4 curve not numerically available locally, post-cure flashes not time points'}]
folds = []
for study in studies:
    train = [s for s in studies if s != study and s['compatible_time_series']]
    folds.append({'heldout': study['source'], 'training_compatible_studies': len(train), 'prediction': None, 'absolute_error': None, 'outcome': 'UNKNOWN', 'reason': 'No compatible independent training time curve' if not train else 'Held-out formulation/oven/time observation incompatible; rate transfer unsupported'})
summary_a = 75.0
summary_b = 75.0
free_a = 0.01
free_b = 0.25
mon = M['normalized monomer elution']['values']
composition = {'TEGDMA': 0.12, 'Bis-GMA': 0.35, 'Bis-EMA': 0.53}
weighted = [sum((composition[s] * mon[s][j] for s in composition)) for j in range(6)]
weighted_error = max((abs(a - b) for (a, b) in zip(weighted, mon['Total'])))
assert weighted_error <= 1
r2 = dict(outcome='CONDITIONAL_TIME_FIT_TRANSFER_NOT_IDENTIFIED', resolution='POPULATION', parameters=frozen['parameters'], replay_absolute_error_pp=float(np.max(np.abs(y - replay))), replay_pass=bool(np.max(np.abs(y - replay)) <= 1e-08), unvisited_30min_DC=pred, rounding_scenarios_min_max=frozen['rounding_scenarios_min_max'], statistical_interval='UNKNOWN: source SD unavailable; scenarios do not bound model error', temperature_rate_fit='UNKNOWN: no calibrated temperature variation; Ea not identifiable', kinetics_n_points=3, kinetics_residual_degrees_freedom=0, loso_folds=folds, compatible_studies_count=1, loso_gate_pass=False, viability_prediction='UNKNOWN: no paired conversion/extract/cell endpoint', monomer_total_ppm=mon['Total'], monomer_composition_weighted_replay_max_error_ppm=weighted_error, monomer_elution_self_over_print_composite_ethanol=mon['Total'][0] / mon['Total'][2], monomer_elution_print_resin_over_composite_ethanol=mon['Total'][4] / mon['Total'][2], sufficiency={'summary': 'FTIR DC (double-bond conversion), composition held fixed', 'identity_error': summary_a - summary_b, 'bit_identical': summary_a.hex() == summary_b.hex(), 'dc_percent': summary_a, 'free_unbound_monomer_fraction_A': free_a, 'free_unbound_monomer_fraction_B': free_b, 'downstream_elutable_fraction_difference': free_b - free_a, 'resolution': 'PHENOMENOLOGICAL', 'minimum_extension': 'Free extractable monomer inventory by species plus local transport/extraction and assay conditions; DC is not monomer mass', 'external_referent': {'kind': 'our_own_fixture', 'locator': 'code/run_demo.py:R2 network-bound versus free monomer construction', 'compared_quantity': 'possible free monomer inventory at identical spectral DC', 'refutes_us': True}}, external_referent=referent('PMC10892052', 'Results/Discussion; three explicit prose DC means', 'conditional DC(t)', False), external_referents=[referent('PMC10182118', 'Table2', 'measured species release with composition-specific ppm denominators', False)], control='Exact scalar saturation fit reproduces same points; direct monomer-weighted arithmetic matches source within integer rounding.', next_measurement='Same batch and oven: measured T and irradiance spectrum, FTIR at0/15/30/45min and depth, LC-MS/HPLC extract species and ISO10993-5 matched cells/geometry; at least3 independent matched studies for transfer. Frozen30min prediction has not been measured.')
finish_round('R2', r2, 'A30min model-conditional DC prediction is frozen. Three points fit three parameters exactly: no predictive validation. All3 leave-one-study-out folds return UNKNOWN under the prespecified matching rule. Temperature activation energy and cytotoxicity are unidentifiable. Monomer ppm arithmetic is composition weighted, not a sum of different denominators.', 'R3: change the representation to measured depth profiles and illumination direction; test a mechanical process query')
depth = M['internal hardness']['values']
depth_rows = []
for row in depth:
    a = np.array(row['mean'])
    threshold = 0.8 * a[0]
    ratios = a[1:] / threshold
    depth_rows.append(dict(time_min=row['time_min'], temperature_C=row['temperature_C'], temperature_label=row['temperature_label'], surface_H=float(a[0]), threshold_H=float(threshold), bottom_H=float(a[-1]), ratios_to80pct=[float(x) for x in ratios], all_sampled_depths_pass=bool(np.all(ratios >= 1)), three_mm_pass=bool(np.all(ratios[:3] >= 1)), uncertainty='Mean response only; SD retained in measurements; no specimen or unsampled-depth guarantee'))
positive_time = [r for r in depth_rows if r['time_min'] > 0]
falseaccept = sum((not r['all_sampled_depths_pass'] for r in positive_time))
bilat = M['directional light response']['values']
bottom2 = bilat['Bottom_2mm']
top2 = bilat['Top_2mm']
equal_cycles = [{'total_cycles': n, 'one_side': one, 'both_sides': both, 'bottom_one': bottom2[one][0], 'bottom_both': bottom2[both][0], 'bottom_ratio': bottom2[both][0] / bottom2[one][0], 'top_one': top2[one][0], 'top_both': top2[both][0]} for (n, one, both) in [(2, 'T2', 'T1B1'), (4, 'T4', 'T2B2'), (8, 'T8', 'T4B4')]]
a = np.array([32.0, 32.0, 32.0, 32.0, 16.0, 16.0])
b = np.array([32.0, 24.0, 24.0, 24.0, 24.0, 32.0])
assert a.sum() == b.sum() and a[0] == b[0]
r3 = dict(outcome='SAMPLED_DEPTH_AND_DIRECTION_QUERY_DELIVERED', resolution='POPULATION; each published hardness observation retains PER_POINT depth', recipes=depth_rows, positive_time_recipes=len(positive_time), surface_uniform_proxy_false_accept_count=falseaccept, surface_uniform_proxy_false_accept_fraction=falseaccept / len(positive_time), untreated_relative_uniform_profiles='0min data are not cytotoxicity evidence or validated acceptable cure; excluded from process acceptance tally', all5mm_sampled_means_pass_count=sum((r['all_sampled_depths_pass'] for r in positive_time)), three_mm_minimum_measured_time_by_T={t: min((r['time_min'] for r in positive_time if r['temperature_label'] == t and r['three_mm_pass']), default=None) for t in ['RT', '40', '60', '80']}, equal_total_cycles=equal_cycles, sufficiency={'summary': ['surface hardness', 'mean depth hardness'], 'surface_identity_error': float(a[0] - b[0]), 'mean_identity_error': float(a.mean() - b.mean()), 'surface_bit_identical': a[0].hex() == b[0].hex(), 'mean_bit_identical': a.mean().hex() == b.mean().hex(), 'bottom_difference_H': float(b[-1] - a[-1]), 'bottom_80pct_pass': [bool(a[-1] >= 0.8 * a[0]), bool(b[-1] >= 0.8 * b[0])], 'profiles': [a.tolist(), b.tolist()], 'minimum_extension': 'Bottom value suffices for a bottom-only query; full regional minimum or depth field needed for all-depth query', 'resolution': 'PHENOMENOLOGICAL', 'external_referent': {'kind': 'our_own_fixture', 'locator': 'code/run_demo.py:R3 depth permutation', 'compared_quantity': 'bottom hardness at identical surface and mean; no measured sample claimed', 'refutes_us': True}}, external_referent=referent('PMC10829619', 'Table1', 'depth-resolved hardness under24 recipes', True), external_referents=[referent('PMC11010625', 'Table1', 'lower-surface hardness at equal total LCU cycles with direction changed', True)], control='Direct table min evaluates same sampled-depth answers to0 error. This information link does not claim an optimization algorithm advantage.', limitations=['80%relative hardness is a mechanical proxy, not a viability or clinical limit', 'Spatial continuum between sampled points and individual parts is UNKNOWN', 'No interpolation or linear sensitivity is claimed', 'Minimum measured time is among observed discrete recipes; not an optimum over unvisited recipes'], next_measurement='Split-direction and one-direction coupons of same thickness, total fluence and batch; depth FTIR, species eluate and matched viability under frozen recipe; validate the biological relevance of the internal cure deficit.')
finish_round('R3', r3, 'Measured internal hardness answers a question the surface-only recipe cannot. Equal total light cycles change lower-face hardness substantially in an independent study. The all-depth80% query is mechanical and sample-mean based; the biological edge stays UNKNOWN.', 'R4: preserve the endpoint vector and paired covariance uncertainty for titanium versus zirconia abutments')
clinical = M['six-month paired abutment endpoints']['values']
decisions = {}
q = float(student_t.ppf(0.975, 14))
for (key, v) in clinical.items():
    diff = v['Zr_mean'] - v['Ti_mean']
    semax = (v['Zr_sd'] + v['Ti_sd']) / math.sqrt(v['n_pairs'])
    semin = abs(v['Zr_sd'] - v['Ti_sd']) / math.sqrt(v['n_pairs'])
    lo = diff - q * semax
    hi = diff + q * semax
    decisions[key] = dict(difference_Zr_minus_Ti=diff, unit=v['unit'], SE_all_correlations_min_max=[semin, semax], worst_correlation_95pct_t_interval=[lo, hi], positive_for_every_correlation=lo > 0, negative_for_every_correlation=hi < 0, spans_zero=lo <= 0 <= hi, published_p=v['published_p'], resolution='POPULATION', interpretation='conditional paired-normal aggregate interval; not equivalence or patient risk')
scores_a = np.array([1.0, 0.0])
scores_b = np.array([0.0, 1.0])
score_a = sum(scores_a)
score_b = sum(scores_b)
r4 = dict(outcome='MARGINAL_SUMMARIES_CANNOT_RESOLVE_MATERIAL_ENDPOINTS', resolution='POPULATION', population='15 partly edentulous patients,30 implants,30–50y; early-loaded customized zirconia/Ti abutments;6months; paired correlation unpublished', endpoints=decisions, paired_t_critical=q, rigorous_enclosure='Cauchy-Schwarz: every paired correlation in[-1,1] gives SE≤(sZ+sTi)/sqrt(n). The t confidence construction remains conditional on paired-normal sampling; floating rounding is not formally interval-certified.', biological_superiority='UNKNOWN: SBI,PD,CBL intervals include0; nonsignificance is not equivalence', patient_material_selection='UNKNOWN: no patient-specific mechanical, allergy-test-accuracy or long-term comparative risk model', TiZr={'source': 'PMC10004271 Tables1–4', 'laboratory_electrochemistry': 'Available for alloy discs inPBS', 'patient_outcome': 'UNKNOWN; laboratory corrosion, implant diameter orpapilla studies do not provide a matched Ti-Zr vsTi/Zr biocompatibility endpoint'}, diagnostic_accuracy={'prevalence': 'UNKNOWN from included primary local sources', 'PPV': 'UNKNOWN: Se/Sp and gold-standard disease labels absent', 'formula': 'Se*p/(Se*p+(1-Sp)*(1-p))', 'review_provenance': 'PMC10670842 Table1 references selected case studies and1500-patient Sicilia cohort; denominator alone does not establish diagnostic ascertainment', 'review_disagreement': 'Some reports endorse MELISA; reviews report inconsistent validity/reliability. Selected symptomatic groups, reagent solubility, cell stimulation and lack of independent reference diagnosis prevent combining positivity as prevalence.'}, sufficiency={'summary': 'equal-weight material score', 'identity_error': float(score_a - score_b), 'bit_identical': score_a.hex() == score_b.hex(), 'downstream_esthetics_difference': float(scores_a[0] - scores_b[0]), 'resolution': 'PHENOMENOLOGICAL', 'minimum_extension': 'Endpoint vector; keep component and indication identity', 'external_referent': {'kind': 'our_own_fixture', 'locator': 'code/run_demo.py:R4 endpoint-vector sum', 'compared_quantity': 'two different endpoint vectors with identical sum', 'refutes_us': True}}, external_referent=referent('PMC10896305', 'Tables1–4', 'six-month paired abutment endpoint vector', False), control='Recomputed primary means; every1001-correlation grid SE lies inside analytic bounds. No method novelty claimed.', next_measurement='Read paired patient endpoint covariance; matched fixture/abutment component mechanical test and fractionated release; independent diagnostic reference labels before an allergy-conditioned comparison.')
for (key, v) in clinical.items():
    gridrho = np.linspace(-1, 1, 1001)
    ses = np.sqrt((v['Zr_sd'] ** 2 + v['Ti_sd'] ** 2 - 2 * gridrho * v['Zr_sd'] * v['Ti_sd']) / v['n_pairs'])
    assert np.max(ses) <= decisions[key]['SE_all_correlations_min_max'][1] + 1e-15
finish_round('R4', r4, 'The all-correlation aggregate bound does not resolve even the esthetic endpoint; the PES interval spans0. The original first-execution label was inconsistent with this computed gate and is preserved under first_execution. Biology endpoints do not resolve superiority; titanium-fixture zirconia-abutment data cannot choose an all-ceramic implant or a patient with suspected allergy. No scalar winner is emitted.', 'Next construction: matched wet-cycle source + fractionated eluate + cell assay, with frozen depth-resolved resin recipes')
checks = []

def checked(name, good, bad, validator):
    assert validator(good), name + ' valid input rejected'
    assert not validator(bad), name + ' corruption accepted'
    checks.append({'name': name, 'valid_pass': True, 'injected_bad_rejected': True})
checked('Table2 implantoplasty mass replay', mass['Control +']['mean'], mass['Control +']['mean'] * 100, lambda x: abs(x - 2265.9) <= 1e-09)
checked('Current Table5 exact mean', curr['H']['mean'], curr['H']['mean'] * 100, lambda x: abs(x - 0.069) <= 1e-12)
checked('Blood Table2 median', blood['after']['experimental_median'], blood['after']['experimental_median'] + 100, lambda x: abs(x - 145.4) <= 1e-09)
checked('Brain Table3 median', org['Brain']['experimental_median'], org['Brain']['experimental_median'] / 1000, lambda x: abs(x - 565.6) <= 1e-09)
checked('Monomer total normalization', weighted, [sum((mon[s][j] for s in composition)) for j in range(6)], lambda x: max((abs(a - b) for (a, b) in zip(x, mon['Total']))) <= 1)
checked('DC exact replay', replay, [v + 10 for v in replay], lambda x: max((abs(a - b) for (a, b) in zip(x, y))) <= 1e-08)
checked('Depth original-table first row', depth[0]['mean'], [v * 1000 for v in depth[0]['mean']], lambda x: max((abs(a - b) for (a, b) in zip(x, [13.3, 13.33, 13.19, 13.33, 13.31, 13.17]))) <= 1e-09)
checked('Bilateral direction at2mm2cycles', bottom2['T1B1'][0], bottom2['T2'][0], lambda x: abs(x - 20.01) <= 1e-09)
checked('PES difference sign', decisions['PES']['difference_Zr_minus_Ti'], -decisions['PES']['difference_Zr_minus_Ti'], lambda x: abs(x - (11.4 - 10.3333)) <= 1e-09)
assay = dict(species='TiO2_NP_lt100nm', unit='µg/mL', cells='THP1_M0', duration_h=24, region='culture_supernatant', locator='PMC10381089 Figure7')
for (field, wrong) in [('species', 'total_digested_Ti'), ('unit', 'ng/g'), ('cells', 'human_bone'), ('duration_h', 720), ('region', 'patient_blood'), ('locator', '')]:
    bad = copy.deepcopy(assay)
    bad[field] = wrong
    checked('Exposure operator ' + field, assay, bad, lambda x: all((x.get(k) == v for (k, v) in assay.items())))
write('raw/control_challenges.json', checks)
candidates = json.loads((ROOT / 'raw/candidates.json').read_text())
selected = set(IDS)
screen = []
for c in candidates:
    id = Path(c['path']).stem
    keep = id in selected
    reason = 'included_primary_or_explicit_secondary_context' if keep else 'XML_parse_failure' if 'error' in c else 'not_selected_by_direct_primary_endpoint_and_table_scope; title/metadata screening only, not a fulltext exclusion'
    screen.append(dict(path=c['path'], title=c.get('title'), article_type=c.get('type'), selected=keep, reason=reason))
write('raw/screening.json', screen)
rejected = sum((not s['selected'] for s in screen))
(fig, axs) = plt.subplots(2, 2, figsize=(12, 8))
fig.subplots_adjust(hspace=0.5, wspace=0.35)
ax = axs[0, 0]
names = ['Control +', 'Rubber dam', 'DAP', 'Bone wax']
ax.bar(range(4), [mass[n]['mean'] for n in names], color=['#777', '#377eb8', '#ff7f00', '#4daf4a'])
ax.set_xticks(range(4))
ax.set_xticklabels(names)
ax.set_ylabel('Retained elemental Ti (µg)')
ax.set_title('Implantoplasty collection: source means\nAAS recovery 53–58%; not tissue dose', fontsize=10)
ax = axs[0, 1]
times = np.linspace(0, 60, 121)
ax.plot(times, [evaluate(p, float(t)) for t in times], label='conditional first-order fit')
ax.scatter(tt, y, color='black', label='3 prose measurements')
ax.scatter([30], [pred], marker='x', s=80, color='#d62728', label='frozen, unmeasured30min')
ax.set_xlabel('Post-cure time (min)')
ax.set_ylabel('FTIR conversion (%)')
ax.set_title('One resin/oven; zero fit residual degrees of freedom\nTemperature and cytotoxicity UNKNOWN', fontsize=10)
ax.legend(fontsize=8)
ax = axs[1, 0]
for (time_, T_) in [(0, 60), (15, 60), (120, 60)]:
    row = next((r for r in depth if r['time_min'] == time_ and r['temperature_C'] == T_))
    ratio = np.array(row['mean']) / row['mean'][0]
    ax.plot(range(6), ratio, marker='o', label=str(time_) + 'min,60°C')
ax.axhline(0.8, color='black', ls='--', label='source80%mechanical proxy')
ax.set_xlabel('Measured depth (mm)')
ax.set_ylabel('Hardness / surface hardness')
ax.set_ylim(0.45, 1.07)
ax.set_title('Depth matters: sampled population means\nHardness is not a cytotoxicity assay', fontsize=10)
ax.legend(fontsize=8)
ax = axs[1, 1]
keys = ['PD', 'CBL']
ax.errorbar([decisions[k]['difference_Zr_minus_Ti'] for k in keys], [0, 1], xerr=[q * decisions[k]['SE_all_correlations_min_max'][1] for k in keys], fmt='o', capsize=5)
ax.axvline(0, color='black', ls='--')
ax.set_yticks([0, 1])
ax.set_yticklabels(['Probing depth', 'Crestal loss'])
ax.set_xlabel('Zr–Ti difference (mm)')
ax.set_title('Paired 15-patient abutment endpoints at6months\n95%t interval enveloped over all correlations', fontsize=10)
fig.savefig(ROOT / 'figure.png', dpi=150)
fig.savefig(ROOT / 'figure.svg')
plt.close(fig)
result = dict(schema='x80-biocompatibility-results-v1', claim_type='information_link', lane='X80-biocompatibility', review_state='PENDING_INDEPENDENT_REVIEW', outcome='SOURCE_AND_SPATIAL_DECISIONS_DELIVERED_TISSUE_AND_TRANSFER_UNKNOWN', created_at=stamp(), external_referent=referent('PMC10829619', 'Table1', 'sampled internal cure versus uniform-depth process proxy', True), rounds={'R1': r1, 'R2': r2, 'R3': r3, 'R4': r4}, preregistrations={p: sha(ROOT / p) for p in PREREGS}, frozen_predictions_sha256=sha(predpath), sources_sha256=sha(ROOT / 'raw/source_manifest.json'), measurements_sha256=sha(ROOT / 'raw/measurements.json'), screening={'candidate_paths': len(screen), 'retained_paths_in_initial_search': len(screen) - rejected, 'rejected_paths': rejected, 'rejected_fraction': rejected / len(screen), 'additional_targeted_sources': [id for id in IDS if id not in {Path(c['path']).stem for c in candidates}], 'scope': 'Initial regex348 paths, title/type/table-count metadata; no claim of complete systematic review'}, controls={'challenges': checks, 'all_injected_bad_rejected': all((c['injected_bad_rejected'] for c in checks)), 'same_information_outcome': 'Direct tables and exact scalar operations give identical answers. No algorithm claim.'}, cost={'compute_wall_seconds': time.perf_counter() - START, 'process_cpu_seconds': time.process_time(), 'peak_rss_kib': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, 'source_bytes_read': sum((s['bytes'] for s in S.manifest)), 'fit_DC_n': 3, 'fit_parameters': 3, 'physical_measurements_performed': 0, 'new_datasets_downloaded': 0, 'questions': 0, 'gpu': False, 'threads_max': 4, 'discovery_time_tokens': 'UNKNOWN; session timestamps retained', 'fallback': 'UNKNOWN plus paired physical measurement specification'}, large_arrays=[], phenomenological_debts=['tissue source R_s(x,t) and retention/clearance/volume', 'free species inventory versus FTIR DC', 'T/irradiance-dependent rate and spatial fit', 'cytotoxicity mapping and interaction with biofilm', 'validated diagnostic sensitivity/specificity and individual baseline risk'])
write('results.json', result)
state('R4', 'COMPLETED_FOUR_CONSTRUCTIONS', 'Paired measurements required before tissue-dose and temperature/viability predictions; graph feedback pending')
print(json.dumps({'outcome': result['outcome'], 'R1_blood_proxy': r1['blood_proxy_pass'], 'R2_predicted_DC30min': pred, 'R2_loso': 'UNKNOWN all3folds', 'R3_false_accepts': str(falseaccept) + '/' + str(len(positive_time)), 'R4_PES_interval': decisions['PES']['worst_correlation_95pct_t_interval'], 'control_challenges_passed': len(checks), 'cost': result['cost']}, indent=2))
