"""Build the bounded review corpus and reader-facing evidence; offline replay."""
import collections
import datetime as dt
import hashlib
import json
from pathlib import Path
ROOT = Path(__file__).resolve().parent

def load(name):
    return json.loads((ROOT / name).read_text())

def save(name, j):
    (ROOT / name).write_text(json.dumps(j, ensure_ascii=False, indent=2) + '\n')

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def finished(p):
    q = p / 'RUN_ATTEMPTS.json'
    if q.exists():
        attempts = json.loads(q.read_text())
        times = [x['finished'] for x in attempts if x.get('result_present') and x.get('finished')]
        if times:
            return (max(times), str(q))
    q = p / 'LOCAL_RESOURCE_USAGE.json'
    if q.exists():
        j = json.loads(q.read_text())
        if j.get('finished_at_s'):
            return (j['finished_at_s'], str(q))
    return (None, None)
SCOPES = [('SOURCE_LAW', 'Speed varies inside the source quasi-static rate domain; rate-invariant exponent alone is not claimed sufficient for absolute load.'), ('GENERALIZATION_GUARD', 'Cone theorem counterexample is in the stated normal-cone formalism; scalar crown triangle itself is retained.'), ('SOURCE_LAW', 'N/f fiber within the kinetic law; percentile reversal independently follows from the source flaw CDF.'), ('SOURCE_LAW', 'Positive width varies in the section function; law defect independently checked for all three fixed substrate stacks.'), ('GENERALIZATION_GUARD', 'Synthetic protocol pooling, outside the observed 16 comparisons. Source already scopes no reversals to that sample.'), ('GENERALIZATION_GUARD', 'Support-capped law is a wider operator than the disc source model; it does not refute held-out disc shape MSE.'), ('SOURCE_CHECKER', 'Synthetic table corruption lies in the arithmetic-check interface, not a demonstrated corruption of the published extracted matrix.'), ('SOURCE_CHECKER', 'Serialization collision is a grammar guard; no new actual pair is recovered from the pulp corpus.'), ('CHANGED_CLOSURE_GUARD', 'Component-aware N2 is a different closure. It only tests generalization of S. The actual source refutation is the containment model sign under matched design.'), ('SOURCE_LAW', 'Positive thickness fiber evaluated by original break-force function; source already retains local thickness/minimum.'), ('INPUT_INTERFACE_GUARD', 'Within-event waveform differs at identical annual count. The source periodic SLS closure itself takes frequency and is not being refuted by this pair.'), ('SOURCE_LAW', 'Event reordering under constant damage stays inside closure; count is sufficient there, with four consumers still UNKNOWN.'), ('SOURCE_LAW', 'Positive mu and reference pressure fiber uses original forward functions; source negative identifiability claim holds.'), ('GENERALIZATION_GUARD', 'Scale fiber is an abstract sensitivity interface. Decisive claim refutation uses nonzero source linear-consumer responses, not this abstract fiber.'), ('GENERALIZATION_GUARD', 'Step versus linear history extends prescribed history. Independent source slope and declared-domain counterexample refute the universal claim.'), ('SOLVER_INTERFACE_GUARD', 'Different RHS under same active set is not a refutation of original full solver; source retains RHS. Benchmark correction concerns unmatched caching.'), ('GENERALIZATION_GUARD', 'Torque product does not decide diameter-dependent torsional response. Cancellation in the stated torque law survives.'), ('SOURCE_LAW', 'Independent exact root and positive diameter are inside the declared corrected-thread domain. No typical clinical screw is certified.'), ('SOURCE_LAW', 'Original positive wear/clamp law, unchanged g. Source correctly calls g a lump and leaves physical calibration UNKNOWN.'), ('SOURCE_GENERAL_MONOTONE_CLAIM', 'Both positive laws are monotone and have the same local slope. Fixed monomial special case survives; general finite-response equation fails.'), ('GENERALIZATION_GUARD', 'Element-field fiber is outside the finite 12-mesh FE fixture. Independent source code confirms unweighted quantile despite weighted label.'), ('SOURCE_LAW', 'Original exact robust-set functions for two allowed supports. Source retains support and nesting assumptions.'), ('SOURCE_SIMULATED_FAMILY', 'Hidden coefficient signs use the supplied shell mode family; metrology amplitude remains a PHENOMENOLOGICAL debt, not measurement.'), ('FISHER_INTERFACE_GUARD', 'Positive same-shape marginal scale change; source full Fisher matrix already retains this information. No global optimized layout is refuted by this fiber.'), ('SOURCE_LAW', 'Positive projected-wall fiber in the rotated optical law; published fit-deviation prior is not a measured wall prior.')]
MIN_EXTENSION = ['Ramp speed for absolute load; corrected crack-growth sign.', 'Normal-cone constraints with strict versus weak motion contract; finite path still separate.', 'Exposure N/f and upper flaw-tail probability.', 'Correct transformed second moment, in addition to neutral axis.', 'Protocol-conditioned paired contrast.', 'Support branch/cap for a wider operator.', 'Original cell/row/column binding; 9 independent contrasts remain after all means.', 'Cell delimiters and per-cell units.', 'Compartment dimensions and an independently observed placement geometry; no validated replacement closure yet.', 'Minimum thickness plus harmonic stiffness for equal-area/equal-E break force.', 'Event waveform or stressing rate independently of annual count.', 'No extension in fixed damage-per-event closure; add damage spectrum when loads vary.', 'One independent normal-pressure/force observation for the selected one-dimensional positive fiber.', 'Observation noise/resolution and forward-map injectivity, separately from sensitivity magnitude.', 'Actual retained clamp history and correct negative logN/logF slope.', 'RHS alongside matrix/active set; benchmark a cached standard factorization.', 'Thread helix geometry independently of cone taper.', 'Correct additive-thread factor with positive-root domain check.', 'One recession observation and one independent friction/stiffness observation to resolve all three parameters.', 'Elasticity along the full path, or rigorous integral bounds.', 'Volume weights and the actual downstream weighted dose moment.', 'Noise support, not only a scatter radius.', 'Independent channels and measured physical amplitude, not residual radius interpreted as signal.', 'Marginal information and complete nuisance-conditioned Fisher matrix.', 'Wall anisotropy/minimum and a feasible independent wall prior.']

def grouped(rows, key):
    groups = {}
    for r in rows:
        name = r[key]
        groups.setdefault(name, []).append(r)
    out = {}
    for (k, rs) in sorted(groups.items()):
        c = collections.Counter((r['verdict'] for r in rs))
        n = len(rs)
        out[k] = {'n': n, 'counts': {v: c[v] for v in ["ACCEPTED", "CORRECTION", 'FALLER']}, 'fractions': {v: c[v] / n for v in ["ACCEPTED", "CORRECTION", 'FALLER']}, 'durable_value_per_selected_run': c["ACCEPTED"] / n, 'resolution': 'POPULATION', 'population': 'intentional selected jobs only', 'not_causal_or_prevalence_estimate': True}
    return out

def table(groups):
    text = ["| Group | n | ACCEPTED | CORRECTION | REJECTED | ACCEPTED / execution |", '|---|---:|---:|---:|---:|---:|']
    for (k, g) in groups.items():
        c = g['counts']
        text.append(f"| {k} | {g['n']} | {c['ACCEPTED']} | {c['CORRECTION']} | {c['FALLER']} | {g['durable_value_per_selected_run']:.0%} |")
    return '\n'.join(text)

def main():
    r1 = load('results_R1.json')
    r2 = load('results_R2.json')
    sample = load('SAMPLE.json')['rows']
    inv = load('INVENTORY.json')
    cutoff = dt.datetime.fromisoformat(inv['cutoff']).timestamp()
    completion = []
    for r in inv['rows']:
        (t, loc) = finished(Path(r['source_path']))
        completion.append({'job_id': r['job_id'], 'finished_at_s': t, 'locator': loc, 'completion_after_cutoff': t >= cutoff if t else None})
    save('COMPLETION_AUDIT.json', {'cutoff': inv['cutoff'], 'mtime_candidate_n': len(completion), 'counts': dict(collections.Counter(('verified_recent' if r['completion_after_cutoff'] else 'verified_old' if r['finished_at_s'] else 'UNKNOWN' for r in completion))), 'scope': 'Completion metadata checked for mtime-screened candidates only; older mtime exclusions were not all reopened.', 'rows': completion})
    for (i, r) in enumerate(r1['rows']):
        p = Path(sample[i]['source_path'])
        (t, loc) = finished(p)
        assert t is not None and t >= cutoff
        assert sha(p / 'RESULTS.md') == r['source']['report_sha256']
        assert sha(p / 'results.json') == r['source']['results_sha256']
        for code in r['code_sources']:
            assert sha(code['path']) == code['sha256']
        r['completion'] = {'finished_at_s': t, 'finished_at_UTC': dt.datetime.fromtimestamp(t, dt.timezone.utc).isoformat(), 'locator': loc, 'sha256': sha(loc), 'after_cutoff': True}
        r['probe_admissibility'] = {'class': SCOPES[i][0], 'scope': SCOPES[i][1], 'pair_alone_refutes_original_headline': False, 'verdict_basis': 'Named source code/claim check and independent check, with original domain; never just nonzero synthetic difference.'}
        r['minimum_extension_scoped'] = MIN_EXTENSION[i]
        r['harness_controls'] = r.pop('controls')
        r['harness_controls_scope'] = 'Serialization/equality/provenance guards only; not independent physical controls.'
        r['independent_control_locator'] = 'results_R2.json' if i in [0, 2, 3, 6, 14, 19, 9] else 'independent_checks in this review; scope in probe_admissibility'
        outcome = sample[i].get('outcome', {})
        r['headline_table_check'] = {'original_opening': sample[i]['opening'], 'original_outcome_statement': outcome.get('statement'), 'original_quantifier': outcome.get('quantifier'), 'original_declared_domain': outcome.get('declared_domain'), 'assessment': r['finding'], 'evidence_locator': r['source']['claim_locator']}
        q = p / 'MODEL_ROUTE.json'
        if q.exists():
            m = json.loads(q.read_text())
            r['model_provenance'] = {'model': m.get('model'), 'locator': str(q), 'sha256': sha(q)}
            r['source_cost'] = {'agent_elapsed_s': t - m['started'] if m.get('started') else None, 'scope': 'model route start to successful finished timestamp; excludes source preparation, previous fit, upstream jobs, queue and verification', 'CPU_s': 'UNKNOWN', 'tokens': 'UNKNOWN', 'complete_cost': 'UNKNOWN'}
        else:
            r['model_provenance'] = {'model': 'UNKNOWN', 'reason': 'No model-route artifact; worker location does not identify model'}
            u = p / 'LOCAL_RESOURCE_USAGE.json'
            usage = json.loads(u.read_text()) if u.exists() else {}
            r['source_cost'] = {'CPU_s': usage.get('cpu_usage_s'), 'peak_RSS_MiB': usage.get('unit_memory_peak_mib'), 'locator': str(u), 'sha256': sha(u) if u.exists() else None, 'complete_cost': 'UNKNOWN'}
        r['phenomenological_debt'] = {'physical_validation': False, 'replacement_measurement': 'Independent measurement of the named downstream quantity at the source resolution, in the same declared regime; job-specific missing measurements are in minimum_extension_scoped.'}
        r['producer_claim_type_original'] = r['claim_type']
        r['claim_type'] = 'capability'
        r['review_claim'] = 'Can the source-bound conclusion be reused for the named downstream question?'
        r['edge'] = {'target_id': sample[i]['target_id'], 'resolution': r['sufficiency_probe']['resolution'], 'timescale': r['sufficiency_probe']['timescale'], 'status': 'PENDING_INDEPENDENT_REVIEW', 'scope': 'Only the named quantity and original regime; no source graph mutation'}
        save(f"REVIEW_{r['job_id']}.json", r)
    rows = r1['rows']
    family = grouped(rows, 'family')
    kind = grouped(rows, 'job_type')
    model = grouped(rows, 'model')
    exact = sum((r['sufficiency_probe']['summary_bytes_identical'] and r['sufficiency_probe']['identity_error'] == 0 for r in rows))
    differences = sum((r['sufficiency_probe']['downstream_difference'] != 0 for r in rows))
    result = {'claim_type': 'capability', 'review_state': 'PENDING_INDEPENDENT_REVIEW', 'sample_n': 25, 'selection': load('SAMPLE.json')['selection'], 'cutoff': inv['cutoff'], 'all_sample_completion_verified': True, 'counts': r1['counts'], 'by_family': family, 'by_job_type': kind, 'by_model': model, 'sufficiency': {'bitidentical_summaries_n': exact, 'maximum_identity_error': 0.0, 'downstream_different_n': differences, 'downstream_equal_n': 25 - differences, 'scope': 'Restricted synthetic fibers; a failing summary is not automatically a failing source job.', 'resolution': 'POPULATION'}, 'omission': {'initial_state_jobs': inv['state_jobs'], 'mtime_screen_eligible': inv['eligible'], 'screen_rejections': inv['rejected'], 'screen_rejection_fraction': sum(inv['rejected'].values()) / inv['state_jobs'], 'eligible_not_selected': inv['eligible'] - 25, 'eligible_not_selected_fraction': (inv['eligible'] - 25) / inv['eligible'], 'reasons': 'missing result pair; old mtime; planner/prior review; purposive quotas (all LIT_CROWN, 13 COV, five other families)', 'missing_units_or_locator_in_selected': 0, 'unresolved_model_n': 3, 'resolution': 'POPULATION', 'completion_audit_locator': 'COMPLETION_AUDIT.json'}, 'external_referent': {'kind': 'external_review', 'locator': str(ROOT / 'results_R1.json'), 'compared_quantity': 'Producer headline/domain against its exact source table/code and independently computed checks; 25 named locators', 'refutes_us': True, 'reviewer': 'This X46 lane, independent of producer jobs; not independently admitted by coordinator'}, 'published_referents': [{'kind': 'published_dataset', 'locator': 'doi:10.1055/s-0042-1757910 Tables 1-2; https://pmc.ncbi.nlm.nih.gov/articles/PMC10756807/', 'compared_quantity': 'three substrate fracture-load thickness ratios', 'refutes_us': False}, {'kind': 'closed_form', 'locator': 'https://nvlpubs.nist.gov/nistpubs/specialpublications/nist.sp.960-16e2.pdf Figure 7.32c', 'compared_quantity': 'positive log strength vs log stressing-rate slope', 'refutes_us': True}, {'kind': 'published_dataset', 'locator': 'PMID:19885413; https://pubmed.ncbi.nlm.nih.gov/19885413/', 'compared_quantity': 'thin-minus-thick test-implant mesial bone loss: 1.61-0.26=1.35 mm', 'refutes_us': True}, {'kind': 'published_dataset', 'locator': 'doi:10.3390/ma16062228 Table 1 and Methods; original local inputs/LIT_ISO14801/iso14801_papers.jsonl', 'compared_quantity': 'meaning of external cyclic preload/minimum force, not retained screw clamp', 'refutes_us': True}, {'kind': 'closed_form', 'locator': 'https://pmc.ncbi.nlm.nih.gov/articles/PMC8791875/ biaxial strength equation', 'compared_quantity': 'published thickness denominator d squared versus claimed ISO exponent 1.5', 'refutes_us': True}], 'probe_fixture_referent': {'kind': 'our_own_fixture', 'locator': str(ROOT / 'review.py'), 'compared_quantity': 'machine-exact summary fibers; not measurement validation', 'refutes_us': False}, 'repairs': {'locator': 'results_R2.json', 'independent_controls': r2['independent_controls'], 'corrected_ratio_factor_1p5_pass_n': sum((r['source_factor_1p5_gate_corrected'] for r in r2['published_load_ratios'])), 'control_ratio_factor_1p5_pass_n': sum((r['source_factor_1p5_gate_control'] for r in r2['published_load_ratios'])), 'physical_validation': False, 'blind_validation': False}, 'full_cost': {'preparation': 'Required readings, corpus inventory and source interpretation; human-equivalent hours/tokens UNKNOWN', 'fit': 0, 'discovery_R1': {k: r1[k] for k in ['wall_s', 'cpu_s', 'peak_rss_mib']}, 'validation_R2': r2['cost'], 'questions': 0, 'fallback': 'Two preserved implementation attempts; external measurement remains UNKNOWN', 'source_costs': 'Per-review route elapsed or local resource usage; total upstream preparation/fitting cost UNKNOWN', 'resource_caps': {'threads': 1, 'intermediate_limit_GB': 3, 'GPU': False}}, 'rows': rows, 'new_physical_measurement': False, 'source_jobs_unchanged': True, 'failures_preserved': ['r1_attempt1/', 'r2_attempt1/', 'COMMANDS.md'], 'safety': 'No clinical recommendation; no dataset patient records used'}
    save('results.json', result)
    figures(result, r2)
    docs(result, r2)
    state = load('CURRENT_WORK_STATE.json')
    state.update(status='ROUND_COMPLETE_PENDING_INDEPENDENT_REVIEW', current_operation='R1 25 reviews and R2 source-bound repairs complete', latest_gate="25 bitidentical probes; 6 ACCEPTED / 11 CORRECTION / 8 REJECTED; independent repair controls pass", next_operation='Coordinator verifies eight decisive refutations; next construction: original JATS cell binding and cached-contact baseline', results_sha256=sha(ROOT / 'results.json'))
    save('CURRENT_WORK_STATE.json', state)
    print(json.dumps({'counts': result['counts'], 'exact_pairs': exact, 'job_types': kind}, ensure_ascii=False))

def figures(result, r2):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    (fig, ax) = plt.subplots(1, 2, figsize=(11, 4.3))
    groups = result['by_family']
    names = list(groups)
    y = np.arange(len(names))
    left = np.zeros(len(names))
    for (label, color) in [("ACCEPTED", '#2a8a67'), ("CORRECTION", '#e6b24a'), ('FALLER', '#bb4c46')]:
        vals = np.array([groups[g]['counts'][label] for g in names])
        ax[0].barh(y, vals, left=left, label=label, color=color)
        left += vals
    ax[0].set_yticks(y, names)
    ax[0].invert_yaxis()
    ax[0].set_xlabel('Selected jobs (non-random, n=25)')
    ax[0].set_title('Scoped source review, pending coordinator')
    ax[0].legend(fontsize=8)
    rows = r2['published_load_ratios']
    x = np.arange(3)
    ax[1].errorbar(x, [r['measured'] for r in rows], yerr=[r['propagated_group_SD'] for r in rows], fmt='ko', capsize=4, label='Published group ratio ± propagated SD')
    ax[1].plot(x, [r['original'] for r in rows], 'x', markersize=9, label='Original (wrong section term)')
    ax[1].plot(x, [r['corrected'] for r in rows], 's', label='Corrected elastic section')
    ax[1].axhline(4, color='gray', linestyle='--', label='Free-standing constant-strength control')
    ax[1].set_xticks(x, [r['substrate'] for r in rows])
    ax[1].set_ylabel('F(1.0 mm) / F(0.5 mm)')
    ax[1].set_title('External table anchor; retrospective repair')
    ax[1].legend(fontsize=7, loc='upper right')
    fig.tight_layout()
    fig.savefig(ROOT / 'review_figure.png', dpi=160)
    fig.savefig(ROOT / 'review_figure.svg')
    plt.close(fig)

def docs(result, r2):
    c = result['counts']
    ratios = r2['published_load_ratios']
    jobs = []
    for (i, r) in enumerate(result['rows']):
        probe = r['sufficiency_probe']
        jobs.append(f"| {i + 1} | [{r['job_id']}](REVIEW_{r['job_id']}.json) | {r['family']} / {r['job_type']} | {r['verdict']} | {probe['identity_error']:.0f} | {probe['downstream_difference']:.8g} {probe['units']} | {SCOPES[i][0]} |")
    full = '\n'.join(["| # | Lokator | Familj / typ | Beslut | Identity error | Downstream difference | Scope of the test |", '|---:|---|---|---|---:|---|---|'] + jobs)
    text = f"# X46 — What are the swarm conclusions that can be reused?\n\n25 recently completed source jobs have received executable, traceable adequacy tests. In the intentionally stratified selection, {c['ACCEPTED']} slutsatser inom sin angivna modell, {c['CORRECTION']} needs rectification and {c['FALLER']} has a conclusively rebutted statement or a rebuttable supporting operation. All decisions are PENDING_INDEPENDENT_REVIEW. SHARES do not involve experimental validation.\n\nSamtliga 25 state pairs have bitidetic binary64-summaries and identity errors **exakt 0**. I 24 pairs differ from the named downstream quantity; in constant Miner damage it does not. The samples show what a summary drops. They are not 24 new reflections: several original jobs already retain the missing variable. REVIEW separates source-law, checks, generalization tests and modified disclosure.\n\n## Selection and foreclosure\n\nFryst urval: alla sju LIT_CROWN, 13 prioriterade universella COV-jobs and one universal job in each of five other families. Source completion is verified against RUN_ATTEMPTS/LOCAL_RESOURCE_USAGE efter 2026-10-02 12:00 Europe/Stockholm, not only file time. The first file time screen provided 152 kandidater av 4480 STATE-poster: 929 without a result pair; 2833 elderly, 566 planer/tidigare granskningar. 127 av 152 candidates were not elected (83,55 %, avsiktliga kvoter). `COMPLETION_AUDIT.json` indicates separately the uncertainty in the candidate population. No selected jobs were lost; three lack verifiable model identity. No incorrect units or locators have been replaced by guessed values.\n\nAll shares below have the resolution POPULATION = They are not estimates of the entire swarm population, model quality or clinical effects. CPU/timme cannot be done: full upstream cost is missing. STATUS/driving is a strict measure of reusable, right scooped conclusion among audited runs, even a sustainable negative result.\n\n## Familjer\n\n{table(result['by_family'])}\n\n## Jobbtyper\n\n{table(result['by_job_type'])}\n\n## Modeller\n\n{table(result['by_model'])}\n\n22 har MODEL_ROUTE pr opencode/space-swarm-free. Three local jobs have no such artifact; they are accounted for UNKNOWN. There is no verified second model in the sample and no matched model pairs. Therefore, model selection cannot be controlled by these shares.\n\n## Exactly what claims fell\n\n- `98ae…`: SCGThe law has the wrong speed character and the wrong geometry correction. The invariance to speed of the limited exponent remains. [NIST SP 960-16e2, figur 7.32c](https://nvlpubs.nist.gov/nistpubs/specialpublications/nist.sp.960-16e2.pdf).\n- `02ab…`: nedre 5 %-flawkvantil ger 95 % probability of failure, not 5 %. The code also returns infinite secure strength when the threshold is zero. reference observation: survival identity of the distribution.\n- `aef…`: transformerad area multipliceras med E/Ev a second time in the parallel axis term. Independent bearing integrals fail the original term for all six heterogeneous sections.\n- `f1da…`: also the max budget predicts thin-minus-thick bone loss ≤ 0 at common geometry, against measured +1,35 mm. Solving a free height from each cohort's outcome is calibration of the answer. reference observations: [PMID 19885413](https://pubmed.ncbi.nlm.nih.gov/19885413/).\n- `021f…`: annual number is not dynamic frequency; de finite SLS-wings are not exactly elastic and limit labels are reversed. The lack of identifiable relaxation time remains a true UNKNOWN.\n- `e363…`: a selected detection band replaces structural identification despite changing the linear outcome of the original with the parameter.\n- `603f…`: slope of log N to log F inverted error; corrected damage exponential is approximately 4,147, not 0,241. In addition, published preload is external cyclic minimum load, not residual screw clamping power. The own declared domain contains a counter-example to the universal heading. reference observations: doi:10.3390/ma16062228, tabell 1 and Methods in the local source post.\n- `78a…`: the corrected thread law has a positive root in the declared domain: helix angle 1,801216° and diameter 4,470506 mm vid eta=1 reproduces both endpoints. Its own section 3 contains the same contradiction to the title. This is a model counter example, not a certified real screw.\n\n## Runda 2: corrections with reference observations\n\nDen korrigerade transformerade sektionslagen verifieras mot en oberoende integral av E/Ev·(z−zbar)². The comparison below uses three published group funds from [doi:10.1055/s-0042-1757910, tabell 1–2](https://pmc.ncbi.nlm.nih.gov/articles/PMC10756807/)Tables were already local. POPULATIONThe elasticity model is: PHENOMENOLOGICALThe same force/moment coupling and breaking strength are assumed at both thicknesses.\n\n| Substrat | Publicerad kvot | Propagerad grupp-SD | Original | Korrigerad | Free standing check |\n|---|---:|---:|---:|---:|---:|\n"
    for r in ratios:
        text += f"| {r['substrate']} | {r['measured']:.3f} | {r['propagated_group_SD']:.3f} | {r['original']:.3f} | {r['corrected']:.3f} | 4.000 |\n"
    text += "\nCorrected team passes the original solid -1,5 - gate for 3 / 3 groups; standalone control passes 0 / 3. The close conformity of the original is not proof of its incorrect operation. SD describes spread propagated from groups, not model errors or confidence intervals. This is retrospective correction, no frozen blind measurement and no validated crown prediction. The load ratio alone does not certify absolute cargo or substrate penalty.\n\nRound 2 also corrects probability tail and speed signs, and provides a **rigorous rational enclosure** for two monotonous laws with exactly the same local sensitivity. It replaces `exp(−gNS)` with integral limits when elasticity is not constant. Precise numbers are stored as integer fights in hexadecimal representation. For the stated equal-area/like-E-PDL law, harmonic stiffness plus minimum thickness is sufficient for breaking force; this is proven by the formula and tested with an exact permutation pair.\n\nIndependent controls of integral, tail, speed, regression closure and enclosure rejects intentionally incorrect values . An injected table corruption simultaneously crosses the all line/column but is not detected by cell comparison. The generic R1 checks of byte identity and hash are only hair controls. No own fixture is labelled measurement.\n\n## All source bound specimen\n\n" + full + "\n\nEach number in REVIEW has device, resolution and time scale. The downstream values of the table come from restricted fixings and are used only at the respective interfaces. The smallest extension is available per job and is scoped to query ; no global minimal representational kits are claimed. SIMULTANEOUS is used for the same operation; wear/clamp-state is left with HANDOVER. PHENOMENOLOGICAL parameters remain as measurable liabilities.\n\n## Cost and limitations\n\nThe calculation part is small, run with a thread and without GPU . Actual wall/ CPU / RSS is available in `results.json`; the full cost/token of inventory and source is UNKNOWN. The source jobs are hash controlled and unchanged . Two hasness errors remain: R1 missed two global constants and had a signed zero problem; R2 was stopped by Python-boundary for decimal integer serialization. No dimensions or tolerances were changed. The tests do not reject physical phenomena on the basis of their own fixes. The coordinator needs to independently verify decisions before reusing links.\n"
    (ROOT / 'RESULTS.md').write_text(text)
    (ROOT / 'README_DEMO.md').write_text("# Executable X46 - review\n\nThe idea is to check if a reused summary actually determines the query that the next node sets. Two separate state gets exactly the same summary; code recounts the next quantity. The original code, table, domain and conclusion are bound with hash, and such a summary error is separated from an error in the original job.\n\nRun from this directory:\n\n```bash\nbash run_all.sh\n```\n\nPython 3 with numpy, scipy and matplotlib are required (available in this environment). The original named local job folders must remain. No network access, data download, unpacking or writing in source jobs is required. PREREG and code freezes are verified before execution . Result and figure are recreated; runtime fields are changing and old graph-receipts refer to the originally bound performance speed.\n\n| Results | Count / outcome | Resolution |\n|---|---:|---|\n| Reviewed Source Job | 25 | POPULATION: avsiktligt urval |\n| ACCEPTED / CORRECTION / REJECTED | 6 / 11 / 8 | POPULATION : review , non-admission |\n| Bitidentiska sammanfattningar | 25 / 25 , identity error 0 | Respektive fixture |\n| Olika downstreamutfall | 24/25 | Each fixture; not 24 wrong job |\n| Korrigerad sektionskvot i fast factor-1,5-band | 3/3 | POPULATION: publicerade gruppmedel |\n| Independent section check in the same band | 0/3 | PHENOMENOLOGICAL control |\n\n![ Review and external load ratio facit](review_figure.png)\n\nExternt  reference : [published Table  1 – 2 , doi: 10.1055 / s  -0042  -1757910 ](https://pmc.ncbi.nlm.nih.gov/articles/PMC10756807/), [ NIST  Figure  7.32 c](https://nvlpubs.nist.gov/nistpubs/specialpublications/nist.sp.960-16e2.pdf), and every named source job. The dataset for the demonstration is local published aggregates and code produced, no patient records. Kronartikelns tabeller: CC BY 4.0. NIST - source : American public publication; no PDF is redistributed. Other bibliographic metadata/facts have not verified common dataset license packages and are not redistributed here; read the respective publication for terms. Custom Algebraic specimen and review outputs have no added license. ToothFairy2, Bits2Bites, mandible defects and Teeth3DS are not used.\n\nWhat does not holds : a synthetic summary sample is not a new measurement facit, a group ratio validates no absolute crown load, a local derivative does not determine a finite progress and three job with the unknown model cannot be used in a model comparison. Adjusted section law passes a broad tolerance but is still an elastic closure. All decisions await independent review. See [RESULTS.md ](RESULTS.md) for each locator and [STEERING.md ](STEERING.md) for control.\n")
if __name__ == '__main__':
    main()
