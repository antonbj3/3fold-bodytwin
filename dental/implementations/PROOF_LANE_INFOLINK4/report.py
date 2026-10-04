from pathlib import Path
import json, csv, hashlib, datetime, time, resource
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
H = Path(__file__).resolve().parent
R = H.parent.parent
j = json.loads((H / 'results.json').read_text())
fac = json.loads((H / 'FACIT.json').read_text())
cs = json.loads((H / 'LINK_CONTRACTS_V2.json').read_text())
it = json.loads((H / 'INTERACTION_RESULTS.json').read_text())
su = json.loads((H / 'SUFFICIENCY_TESTS.json').read_text())
pm = json.loads((H / 'raw/PORE_METRICS.json').read_text())

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def save(n, x):
    (H / n).write_text(json.dumps(x, indent=2, ensure_ascii=False) + '\n')
plt.rcParams.update({'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False})
(fig, ax) = plt.subplots(2, 2, figsize=(11, 8))
colors = ['#00796B', '#E07A22', '#5B62A5']
for (ang, col) in zip([0, 45, 90], colors):
    rr = [r for r in fac['L01'] if r['condition_key']['material'] == 'NextDent C&B' and r['condition_key']['orientation_deg'] == ang]
    ax[0, 0].errorbar([r['condition_key']['postcure_min'] for r in rr], [r['original_value'] for r in rr], yerr=[r['reported_sd'] for r in rr], color=col, marker='o', capsize=3, label=f'{ang}°')
ax[0, 0].set(title='NextDent : ranking depends on post curing', xlabel='Post curing ( min )', ylabel='Failure load (N), mean +/- SD')
ax[0, 0].legend(frameon=False, loc='lower right')
ax[0, 0].text(0.02, 0.96, 'POPULATION · n=10/condition\nPMC10097162, table 3', transform=ax[0, 0].transAxes, va='top', fontsize=8)
rr = fac['L11']
ax[0, 1].bar(range(len(rr)), [r['original_value'] for r in rr], yerr=[r['reported_sd'] for r in rr], capsize=4, color='#00796B')
ax[0, 1].set_xticks(range(4), ['0', '16', '32', '60'])
ax[0, 1].set(title='Resilab Temp: state after hardening', xlabel='Post curing ( min )', ylabel='Bench strength (MPa ), mean ± SD')
ax[0, 1].text(0.03, 0.96, 'POPULATION · n=9/villkor\nPMC11012777, table 1', transform=ax[0, 1].transAxes, va='top', fontsize=8)
ax[0, 1].set_ylim(0, 108)
rr = fac['L07']
ax[1, 0].bar(range(3), [r['original_value'] for r in rr], yerr=[r['reported_sd'] for r in rr], capsize=4, color=['#777777', '#E07A22', '#00796B'])
ax[1, 0].set_xticks(range(3), ['Cylinder', 'Implant', 'Tandreplika\nartificiell PDL'])
ax[1, 0].set(title='Support construction and LDS wear', ylabel='Volume loss (mm³ ), mean ± SD')
ax[1, 0].text(0.03, 0.96, 'POPULATION · n=15/villkor\nPMC10455685, table 3', transform=ax[1, 0].transAxes, va='top', fontsize=8)
ax[1, 0].set_ylim(0, 0.53)
pids = [37, 38, 40, 36, 35, 11]
ax[1, 1].bar([str(x) for x in pids], [pm[str(x)]['count_d_ge_0_5mm'] for x in pids], color=['#E07A22'] * 2 + ['#777777'] * 4)
ax[1, 1].set(title='Ti64: svansantal ur ursprungliga porkataloger', xlabel='Process set / catalogue ID', ylabel='Count of pores with diameter >=0.5 mm')
ax[1, 1].text(0.03, 0.96, 'POPULATION per catalogue\nInputs PER_POINT · Luo 2022', transform=ax[1, 1].transAxes, va='top', fontsize=8)
ax[1, 1].set_ylim(0, 48)
fig.suptitle('New information links: the measurement value retains its condition', fontsize=14)
fig.text(0.02, 0.012, 'Published observations and recalculated catalogue number. SD is spreading, not predictive interval. No clinical transfer.', fontsize=8)
fig.tight_layout(rect=[0, 0.045, 1, 0.955])
for ext in ['png', 'pdf', 'svg']:
    fig.savefig(H / ('information_links_figure.' + ext), dpi=160)
plt.close(fig)
with (H / 'FACIT.csv').open('w') as f:
    w = csv.writer(f)
    w.writerow(['link', 'condition_or_pore', 'original_value', 'reported_sd', 'unit', 'resolution', 'time_scale', 'source_locator', 'source_sha256'])
    for (lid, rs) in fac.items():
        for r in rs:
            w.writerow([lid, json.dumps(r.get('condition_key', {'process_set': r.get('process_set'), 'pore_id': r.get('pore_id')}), ensure_ascii=False), r['original_value'], r.get('reported_sd'), r['unit'], r['resolution_level'], r['time_scale'], r['source_locator'], r['source_sha256']])
rows = []
for c in cs:
    rs = fac[c['id']]
    s = next((z for z in su if z['id'] == c['id']))
    a = s['query_A']
    b = s['query_B']
    unit = rs[0]['unit']
    rows.append(f"| {c['id']} | {c['title']} | {a:g} / {b:g} {unit} | POPULATION | {c['time_scale']} | [{c['source']}](https://pmc.ncbi.nlm.nih.gov/articles/{c['source']}/) |")
rows.append('| L12 | Pork catalog → measured upper tail | 36 / 35 pores >=0.5 mm in catalogues 37/38 | PER_POINT → consumer catalogue count | SIMULTANEOUS | Luo 2022; original CSV-rader i FACIT.csv |')
table = '\n'.join(rows)
head = '| Link | Question | Source values | Resolution | Tidsskala | Reference |\n|---|---|---|---|---|---|'
inst = json.loads((H / 'GRAPH_INSTALL_CHECK.json').read_text()) if (H / 'GRAPH_INSTALL_CHECK.json').exists() else None
if inst:
    assert sha(Path(inst['created_file'])) == inst['created_file_sha256'] == sha(H / 'information_links.preview.jsonl')
    j['graph_install'] = inst
j['additional_claims'] = [{'claim_type': 'capability', 'id': 'R3', 'result': 'INTERACTION_RESULTS.json', 'summary': 'Same marginal summaries exactly, different queried angle choice; one interaction contrast reconstructs the chosen2x2 exactly.', 'physical_prediction': 'NOT_CLAIMED'}]
j['phenomenological_debts'] = [{'quantity': 'pore diameter query threshold', 'value': 0.5, 'unit': 'mm', 'resolution_level': 'PHENOMENOLOGICAL', 'meaning': 'Inherited descriptive cutoff, not a validated dental acceptance threshold', 'replacement_measurement': 'Matched XCT defects and fatigue outcomes at controlled geometry, load and process; calibrate decision-specific threshold prospectively.'}, {'quantity': 'legacy CV family-branch threshold', 'value': 0.5, 'unit': '1', 'resolution_level': 'PHENOMENOLOGICAL', 'meaning': 'Classifier heuristic from old code; no physical mechanism inference', 'replacement_measurement': 'Independent labeled microstructural populations and held-out classification validation.'}]
j['source_semantic_review'] = {'R1': '2 time-scale assignments failed self-review; original files preserved', 'V2': 'L06/L07 corrected to SIMULTANEOUS; all other numeric thresholds unchanged', 'independent_review': 'PENDING'}
j['fault_injections_rejected_total'] = j['fault_injections_rejected'] + 2
j['cost']['final_report_peak_rss_kib'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
j['artifact_bytes'] = sum((p.stat().st_size for p in H.rglob('*') if p.is_file()))
j['artifacts'] = {name: sha(H / name) for name in ['FACIT.csv', 'INTERACTION_RESULTS.json', 'information_links.preview.jsonl', 'information_links_figure.png', 'SEMANTIC_REVIEW_R1.json']}
save('results.json', j)
readme = f"""# Conditional measurement values for dental materials and structures\n\nEleven new information links make published test results searchable with the right materials, manufacturing conditions and measurement definition. Reference is ten local original articles and six original Ti64por catalogs. The links are drafts waiting for independent review .\n\nThe idea is to give the consumer the full measurement condition. Example: for NextDent -treled bridges, 45° has higher published mean fracture load than 0° after 30 min curing ( 1307,32 against 980,72 N ), but the order is turned after 120 min ( 1507,19 against 1683,56 N ). These are **POPULATION**-means with n=10 per condition; SD is available in reference. They decide which table line fits a lab question, not which clinical choice is best. Source : [original table 3 ](https://pmc.ncbi.nlm.nih.gov/articles/PMC10097162/#polymers-15-01737-t003).\n\n| Link | What Can Be Asked | Two source values; the conditions are in reference | Resolution | Tidsskala | External reference |\n|---|---|---|---|---|---|\n{table}\n\nRun from this directory:\n\n```sh\n./run_all.sh\n```\n\nThe command rereads the originals, recounts the control, tests the information loss and creates table/figure. It uses existing `/usr/bin/python3`, NumPy, lxml and Matplotlib with user packages turned off. The original files need to be on the hashbound local paths in SOURCE_MANIFEST.json . No data is downloaded. The graph is only changed by the separate installation command.\n\nA concrete query can then be set:\n\n```sh\nPYTHONNOUSERSITE=1 /usr/bin/python3 query.py L01 '{{"material":"NextDent C&B","orientation_deg":0,"postcure_min":120}}'\n```\n\nThe answer is 1683,56 N , reported SD 207,57 N , n= 10 , POPULATION , and exactly original cell. A unknown condition gives `UNKNOWN_SOURCE_CONDITION` ; the code does not interpolate new material properties.\n\n![Published Averages and Catalog Number](information_links_figure.png)\n\nThe adequacy test goes beyond comparing different averages. Two designed 2 × 2 tables have exactly the same row and column sum, even when floating point control: identity error **0**. Nevertheless, the direction that has the highest mean at both 30 and 120 min is changed, and the answer for a given box distinguishes **400 N**. A single signed interaction contrast is the smallest extension that restores all four squares. Table A comes from the article; Table B is a mathematical counter-example, not new test measurements. Separate positive pore constructions have identical number, mean and other elements but different tail numbers: 0 versus 1 . Full values, identity error and reconstruction are available in INTERACTION_RESULTS.json and SUFFICIENCY_TESTS.json.\n\nRestrictions governing use:\n\n- The group host still has the POPULATION level. They don't become local material fields or patient predictions. SD is not a confidence or prediction interval. Unmeasured covariance is not filled in.\n- The module change in L05 is N cm, not preload loss in N. L03/L04 is the remaining strength after load history, not a fatigue limit. L07 uses artificial PDL.\n- Five source conflicts follow: among other things HL group SD51 N in table against 43 N in abstract, the drilling study's conflicting inference and Resilab's printer/age status. SOURCE_CONFLICTS.json is part of the substrate. Optic values are published TP, not an independent reconstruction of the problem printed formula.\n- Ti64 cell's narrow hypothesis of CV < 0,5 failed already before and fails also here. Katalog37 / 38 provides CV {pm['37']['legacy_family']['cv']:.3f}/{pm['38']['legacy_family']['cv']:.3f}. The classification label is `lognormal`; it proves no physical pore mechanism. No BIC mixture analysis or association to fatigue life is retested here.\n- Control without the new information lacks the conditional answers. A separate direct spread with the same original data matches all 85 group values and is used as extraction control. No algorithmic speed or prediction gain is claimed. All 132 injected faults rejected.\n\nLoss: **10 / 21 candidates ( 47,62% )** was rejected or excluded: six already established priority links, an already reused crown study, an incorrect flexure purpose and two older cells without located local primary facit. Bland15 other reviewed proposals are the loss of 4/15 ( 26,67%). This is not a review grade for all 196 index entries. De85 selected publication tables were approved; 71 tail entries come from:{j['attrition']['pore_record_attrition']['read_rows']} read pore posts. Items below the tail threshold are retained in the raw register and are not lost.\n\nLicenses: PMC10097162 , PMC10007339 , PMC10569846 , PMC10743638 , PMC10051685 , PMC10455685 , PMC10453045 and PMC11012777 specify CC BY4.0 in the originals permissions. PMC10350602 states CC BY-NC-ND4.0; PMC10088447 CC BY-NC-SA4.0. These measurement values and locators are treated here local . Ti64 -Luo catalogues license is **UNKNOWN** in reviewed local file; the original files are not redistributed. No patient volumes or journal data are used. No crude originals are copied to the demonstration.\n\nAll numbers in the figure are source measurements or deterministic catalogue calculations. The constructed counter-examples are derivations. No affinity susceptibility approximation, physical simulation or clinical recommendation is reported.\n"""
readme = readme.replace('and 120', 'and 120').replace('mot1', 'against 1').replace('all 85', 'all 85').replace('Totalt132', 'Total 132').replace('Bland15', 'Bland 15').replace('attritionet4', 'attritionet 4').replace('all 196', 'all 196').replace('De85', 'De 85').replace('CC BY4', 'CC BY 4').replace('CC BY-NC-ND4', 'CC BY-NC-ND 4').replace('CC BY-NC-SA4', 'CC BY-NC-SA 4').replace('Katalog37/38', 'Catalogue 37/38').replace('SD51', 'SD 51').replace('mot43', 'against 43')
(H / 'README_DEMO.md').write_text(readme)
(H / 'RESULTS.md').write_text(f"# PROOF_LANE -infolink 4 — eleven new conditional information links\n\nEleven new connections to dental nodes/chains are source bound : 85 published averages as well as 71 individual tail spore items. Reference are ten primary articles and Luo2022 directories. [ README_DEMO.md ](README_DEMO.md) has result , figure, licenses and a run command.\n\n{head}\n{table}\n\nBortfall10 / 21 = 47,62% ; unique proposal 4 / 15 = 26,67% . See ATTRITION.json and REJECTIONS.json . No clause about all 196 cells or reduction of the old number 112/158 without a new chain review.\n\nAll 85 agents and reported SD are reproduced by two separate XML readings and a manually controlled facit list. 130 error injections within the link contracts plus 2 in the interaction consumer are rejected. The check without the information link gives UNKNOWN at the specific source query. Direct spreads with the same source are only an extraction check, no method contest.\n\nSufficiency : 11 tests with identity error 0 in the summary and various downstream responses. In addition, R3 displays exactly identical row/column numbers in two 2 × 2 tables, 400 N grid difference and turn angle order. An interaction contrast restores exactly. The permits are information structures; A in R3 is observed source agent, B synthetic. No physical experimental intervention performed.\n\nPreserved negative : first semantic review failed two time scales (L06 , L07 ); they were corrected to SIMULTANEOUS in a separate frozen V2 . No numerical tolerances were changed. Previously result remains with suffix . R1 . G1335 : s CV -hypotheses are still REFUTED . Five source conflicts are preserved. L02 incorrectly called flexur got a new explicit drag purpose L15 in frozen R2 . Status is PENDING_INDEPENDENT_REVIEW .\n\nGraph: {('installation verified; all 11 ID available in work view; see GRAPH_INSTALL_CHECK.json' if inst else 'installation remaining')}. No personal scientific admission. Graph feedback is held local according to the specific one-file boundary.\n\nFull cost: demo before figure {j['cost']['demo_wall_seconds_before_plot']:.3f}s , peak RSS before figure {j['cost']['peak_rss_kib']}KiB; source exchange{j['cost']['source_bytes_read_and_hashed']}. Figurprocessens RSS{j['cost']['final_report_peak_rss_kib']}KiB. Cumulative command times are available in runes/*/ commands.json . Preparation, model reasoning and previous data collection are not individually timed; no full cost advantage is claimed. Fit0 ; new measurements 0 ; network issues 0 . No GPU , no big midway.\n")
feedback = {'target_id': 'DENT-MAT-RESTORATIVE-PSP', 'additional_targets': sorted({x['consumer'] for x in cs if x['target_kind'] == 'draft_node'} | {'DENT-MAT-TI64-LPBF'}), 'chain_targets': ['K34', 'IM05'], 'result_file': 'results/PROOF_LANE_INFOLINK4/results.json', 'sha256': sha(H / 'results.json'), 'measured_quantity': '85 source-conditioned means and71 pore records in11 new information links; exact summary counterexamples', 'units': 'Explicit per observation; N, MPa, N cm, mm^3, degC, CIE TP, mm', 'resolution_levels': ['POPULATION', 'PER_POINT'], 'time_scales': ['SIMULTANEOUS', 'HANDOVER'], 'uncertainty': 'Source SD and identified conflicts retained; no clinical/part prediction or confidence claim', 'population_regime': 'Laboratory dental materials, assembly/crown tests, external Ti64-XCT catalogs', 'preregistered_gate': 'PREREG_INFOLINK4_R1/R2 and PREREG_INTERACTION_R3; exact source/semantic fields and exact sufficiency test', 'baseline': 'Missing source-conditioned input; legacy summary for pore query. Not algorithm competition.', 'outcome': '11 ESTABLISHED_DRAFT;10/21 excluded;132 injected faults rejected', 'negative_result': True, 'review_state': 'PENDING_INDEPENDENT_REVIEW', 'review_artifact': '', 'review_scope': '', 'review_sha256': '', 'submission_status': 'LOCAL_ONLY_SPECIFIC_WRITE_BOUNDARY', 'reason_not_dispatched': 'User specifies one new shared expansion file, no other file changes, plus working build. tasks/graph_runs writes would exceed that scope.'}
save('GRAPH_FEEDBACK.json', feedback)
save('CURRENT_WORK_STATE.json', {'lane': 'PROOF_LANE-infolink4', 'phase': 'COMPLETE_PENDING_INDEPENDENT_REVIEW' if inst else 'READY_TO_INSTALL', 'last_gate': '11 source links,132 injected faults rejected; exact joint-interaction sufficiency passed', 'next_operation': 'Independent semantic/source review; then matched manufacture of NextDent0/45deg at30/120min with frozen out-of-sample predictions', 'updated_at_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'graph_installed': bool(inst), 'handoff': 'HANDOFF.md'})
print('Report/figure/ reference /query complete; feedback hash:', feedback['sha256'])
