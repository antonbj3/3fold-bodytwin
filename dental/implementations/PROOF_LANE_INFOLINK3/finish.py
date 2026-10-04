from pathlib import Path
import json, hashlib, datetime
H = Path(__file__).resolve().parent
r = json.loads((H / 'results.json').read_text())
g = json.loads((H / 'GRAPH_INSTALL_CHECK.json').read_text())
m = r['metrics']
valid = json.loads((H / 'VALIDATION.json').read_text())
ex = json.loads((H / 'REJECTIONS.json').read_text())
if not any((x['id'] == 'L11' for x in ex)):
    ex.append({'id': 'L11', 'candidate': 'PLA Resumen modulus to lattice design', 'reason': 'Primary summary does not declare original modulus unit. Relative benchmark retained locally; no absolute modulus edge admitted.', 'phase': 'AFTER_ORIGINAL_REVIEW', 'artifact': 'REJECTED_L11.json'})
(H / 'REJECTIONS.json').write_text(json.dumps(ex, indent=2, ensure_ascii=False))
hand = f"""# Handoff — PROOF_LANE-infolink3

Eleven new information links installed in `notes/expansion/information_links_astra3_20261003.jsonl`. `GRAPH_INSTALL_CHECK.json` verifies all eleven IDs in the working view. Build: 768 -> {g['build_receipt']['draft_nodes']} drafts, 6 conflicts and 49 gaps unchanged. All other expansion files have identical hashes to the start of the round. Native statuses unchanged. Exact expansion hash: `{g['created_file_sha256']}`.

Read `README_DEMO.md` for idea, table, figure and limitations. One command: `./run_all.sh`. It reproduces analysis without graph writes. `python3 install_and_check.py` is the separately executed, explicitly authorized installation/build; it refuses to overwrite a differing file. Result control hashes in `VALIDATION.json`. `GRAPH_FEEDBACK.json` is local and PENDING_INDEPENDENT_REVIEW: the specific one-file boundary excludes dispatch/feedback to `tasks/graph_runs`. Coordinator may bind the result after review. Suitable targets exist; no invented coverage target required.

Outcome: L01-L10 and L12 established as scoped drafts. L11 failed for missing explicit modulus unit in original table. Six prioritized brief links already existed in PROOF_LANE2 and were not repeated. Four stronger interpretations rejected: certified grain diameter transform, sphericity -> Kt/fatigue initiation, c-axis components -> true grains, elastic projection deviation -> nonidentifiability theorem. Attrition 11/22 = 50.00% including duplicates; 5/16 = 31.25% for unique proposals. This describes reviewed proposals, not all 196 index records.

Principal raw results:
- Ley fatigue: 138 characterized fracture sites, 136 test regions, 99 physical specimens. 30 fracture site rows have retests after 10 million initial cycles. Primary key is BuildID x SpecimenPosition x TestRegion x InitDefectID, not independent row numbers. Prior loading must follow the specimen.
- L02 partial rank correlation {m['L02']['by_mode']['all']['partial_rank_correlation']:.6f}; first site per physical specimen gives {m['L02']['first_site_per_specimen_partial']:.6f}. L03 location association descriptive; no causal surface factor established.
- Plastic curves: 84 curves. Hollomon RMSE {m['L04']['holl_rmse_median']:.6f} MPa fails older <15 MPa gate despite original cell exit code 0. Frozen gate retained.
- EBSD: pixel-weighted lambda1 {m['L05']['maps'][0]['eigenvalues'][-1]:.9f}/{m['L05']['maps'][1]['eigenvalues'][-1]:.9f}. No conclusion about E anisotropy magnitude. C-axis component diameter PHENOMENOLOGICAL.
- 33/35 pore files accepted. Files 5 and 16 lack diameter units and have `sorted` as fifth column; earlier cells assumed sphericity. New 33-file results separate from preserved old 35-file replays. L07 process-set-adjusted correlation {m['L07']['partial_max_given_porosity']:.6f}; bootstrap interval {m['L07']['bootstrap_process_set_CI95']}. Same-coupon pairing unproven.
- Grain observables differ by factor {m['L10']['mean_table_to_feature_ratio']:.6f}. Observed definition gap; original word stereology does not prove these data types are 2D and 3D.
- PLA first peak: 35 curves, median absolute error {m['L12']['median_abs_projection_residual']:.6f} for retrospective elastic projection. First peak is not independently MEASURED buckling.

Before work, prioritized older links were expected missing; they already existed. Useful addition instead preserves specimen history, size+location, raw curves, orientation, pore shape and measurement definitions. Surprising correction: apparently reproducible cell values can still carry wrong specimen identity, unidentified columns or excessive physical interpretation. Successful replay alone must not become a physical edge.

Next construction: hierarchical fatigue consumer for L01-L03, treating initial runout and retest as two loading stages sharing history. Separate physical specimen key, test region and fracture site; preserve specimen key across validation split. Freeze predictions before comparison. First quantify how incorrect independent-row assumptions and missing retest history change uncertainty. Do not combine Ley coupon stresses with external ISO14801-RAW36 assembly loads; geometry/load transfer model absent. Next physical reference needs finished dental specimens with same postprocessing, 3D defect measurements and observed fracture initiation.

Independent review should first check CSV identity (99 physical specimens is our derivation from source keys), pore file exceptions, grain measurement definitions and UNKNOWN licensing. External EBSD deposition locator not recovered; original CTF and unchanged cells locally hash-bound. No license assumed or original data redistributed.

Resources: one thread, about {r['cost']['wall_seconds']:.1f} s per demo, {r['cost']['max_rss_kib'] / 1024:.1f} MiB peak RSS, {r['cost']['data_bytes'] / 1000000.0:.2f} MB extracted intermediates, no GPU jobs/downloads. {valid['source_values_passed']} original value checks and {valid['fault_injections_rejected']} injected faults checked; certify extraction/contracts, not new prediction accuracy.
"""
(H / 'HANDOFF.md').write_text(hand)
state = {'lane': 'PROOF_LANE-infolink3', 'phase': 'ROUND_COMPLETE_PENDING_INDEPENDENT_REVIEW', 'updated_at_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'latest_gate': '11 links installed; all imported; other expansion files unchanged; 20 exact source witnesses and 68 mutations checked', 'next_operation': 'Independent source/semantics review, then specimen-grouped fatigue-history consumer described in HANDOFF.md', 'established_ids': g['imported_ids'], 'attrition': r['candidate_attrition'], 'graph_feedback': 'LOCAL_ONLY; specific write boundary', 'expansion_sha256': g['created_file_sha256']}
(H / 'CURRENT_WORK_STATE.json').write_text(json.dumps(state, indent=2, ensure_ascii=False))
final = 'Eleven new information links established as drafts, with original values, locators, resolution and timescale.\nThey connect specimen history, defect location, plastic curves and pore shape to dental consumers.\nGraph rebuilt: all eleven links imported; other expansion files unchanged.\nAttrition: 11/22 proposals (50%), six already established; excluding duplicates 5/16 (31.25%).\nPrincipal correction: 138 fatigue fracture site rows correspond to 99 physical specimens.\nLargest obstacle remains transfer from coupon measurements to finished dental constructions.\nRun `./run_all.sh`; references, figure and limitations in README_DEMO.md. Independent review pending.\nNext construction is a fatigue consumer preserving specimen identity and retest history.'
(H / 'SOL_NIGHT_FINAL_PROOF_LANE-infolink3.md').write_text(final + '\n')
print('Handoff and state saved')
