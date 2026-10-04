"""Build reader-facing artifacts from preserved raw, frozen results."""
from dental_release.paths import expand as _release_expand
import os
for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[k] = '2'
import json, hashlib, datetime, resource, re
from pathlib import Path
import numpy as np
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent / 'vendor_plotting'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from benchmark import P, write, sha, state
from parser import FIELDS

def run():
    r1 = json.load(open(P / 'raw/RESULTS_R1.json'))
    r2 = json.load(open(P / 'raw/RESULTS_R2.json'))
    r3 = json.load(open(P / 'raw/RESULTS_R3.json'))
    r4 = json.load(open(P / 'raw/RESULTS_R4.json'))
    joint = json.load(open(P / 'raw/JOINT_DIAGNOSTIC.json'))
    m = json.load(open(P / 'raw/DATA_MANIFEST.json'))
    d = json.load(open(P / 'raw/REPORT_DISAGREEMENT.json'))
    pr = json.load(open(P / 'raw/PARSER_REVIEW_CORRECTED.json'))['summary']
    hold = json.load(open(P / 'raw/PARSER_HOLDOUT_REVIEW.json'))['summary']
    g = [json.loads(l) for l in (P / 'raw/geometry.jsonl').read_text().splitlines()]
    fit = json.load(open(P / 'raw/FIT_COST.json'))
    names = ['Overbite', 'Overjet', 'Korsbett', 'Mittlinjer', 'Spee', 'Molar right', 'Molar left', 'Cornerand right', 'Corner-left']
    (fig, axs) = plt.subplots(1, 3, figsize=(16, 6.8), layout='constrained')
    y = np.arange(9)
    for (offset, meth, label, col) in [(-0.22, 'majority', 'Training majority practice proxy', '#bdbdbd'), (0, 'threshold', 'Simple metric thresholds', '#e3a558'), (0.22, 'point', '31-feature model', '#467ba9')]:
        axs[0].barh(y + offset, [r1['metrics'][f][meth]['balanced_accuracy'] for f in FIELDS], 0.21, label=label, color=col)
    axs[0].set_yticks(y, names)
    axs[0].invert_yaxis()
    axs[0].set(xlim=(0, 1), xlabel='Balanced accuracy vs first clinician report', title='Held-out patients: categorical agreement')
    axs[0].legend(fontsize=8, loc='upper center', bbox_to_anchor=(0.5, -0.14))
    for (meth, label, col) in [('point', 'Forced single category', '#e3a558'), ('first_reference_conformal', 'First-reference conformal', '#9ea6b2'), ('candidate', 'All-report category set', '#467ba9')]:
        axs[1].plot([r2['metrics'][f][meth]['all_report_patient_coverage'] for f in FIELDS], y, 'o-', label=label, color=col)
    axs[1].axvline(0.85, color='black', linestyle='--', linewidth=1, label='Frozen coverage gate')
    axs[1].set_yticks(y, names)
    axs[1].invert_yaxis()
    axs[1].set(xlim=(0, 1), xlabel='Fraction of patients: all readable reports covered', title='Uncertainty retains report disagreement')
    axs[1].legend(fontsize=8, loc='upper center', bbox_to_anchor=(0.5, -0.14))
    axs[2].bar(['Published axes', 'Vector frame'], [r3['max_published_error_mm'], r3['max_candidate_error_mm']], color=['#b86a54', '#467ba9'])
    axs[2].axhline(0.05, color='black', linestyle='--', linewidth=1)
    axs[2].set(ylabel='Maximum rigid-rotation change (mm)', title='3 real pairs × 3 common rotations')
    axs[2].text(0.95, 0.95, 'Anatomical landmark accuracy: UNKNOWN\nNumerical frame diagnostic only', transform=axs[2].transAxes, ha='right', va='top', fontsize=8)
    fig.savefig(P / 'figures/bite_geometry_vs_reports.png', dpi=170)
    fig.savefig(P / 'figures/bite_geometry_vs_reports.pdf')
    plt.close(fig)
    (fig, ax) = plt.subplots(1, 2, figsize=(9, 4.5), layout='constrained')
    ax[0].bar(['Marginal sets R2', 'Whole-bite sets R4'], [joint['all_available_fields_all_reports_covered'], r4['joint_patient_coverage']], color=['#467ba9', '#b86a54'])
    ax[0].axhline(0.85, color='black', linestyle='--')
    ax[0].set(ylim=(0, 1), ylabel='All readable fields/reports covered per patient', title='Joint coverage vs frozen 0.85 gate')
    ax[1].bar(['Marginal sets R2', 'Whole-bite sets R4'], [r2['macro']['mean_size'], r4['mean_field_set_size']], color=['#467ba9', '#b86a54'])
    ax[1].axhline(2, color='black', linestyle='--')
    ax[1].set(ylabel='Mean categories per field', title='Information vs frozen size <=2 gate')
    fig.suptitle('Whole-bite calibration buys coverage with uninformative sets; R4 FAIL', fontsize=11)
    fig.savefig(P / 'figures/joint_bite_tradeoff.png', dpi=170)
    fig.savefig(P / 'figures/joint_bite_tradeoff.pdf')
    plt.close(fig)
    table = '| Quantity | n test | Model accuracy | Balanced accuracy | Kappa | Enkla regler BA | Report–rapport* | All reports in amount | Average Size |\n|---|---:|---:|---:|---:|---:|---:|---:|---:|\n'
    for (f, name) in zip(FIELDS, names):
        a = r1['metrics'][f]['point']
        v = r2['metrics'][f]['candidate']
        agree = d[f]['first_second_agreement']
        table += f"| {name} | {a['n']} | {a['accuracy']:.3f} | {a['balanced_accuracy']:.3f} | {a['kappa']:.3f} | {r1['metrics'][f]['threshold']['balanced_accuracy']:.3f} | {agree:.3f} | {v['all_report_patient_coverage']:.3f} | {v['mean_size']:.2f} |\n"
    err = [x for x in g if x['error']]
    cost = dict(preparation_seconds=sum((x['wall_seconds'] for x in g)), preparation_n=len(g), geometry_errors=len(err), fit=fit, frozen_prediction_seconds=json.load(open(P / 'FROZEN_PREDICTIONS.json'))['query_wall_seconds'], validation_r1_seconds=r1['validation_wall_seconds'], calibration_validation_r2_seconds=r2['wall_seconds'], orientation_diagnostic_seconds=r3['wall_seconds'], joint_calibration_validation_seconds=r4['wall_seconds'], discovery_code_source_search_seconds=None, manual_review_seconds=None, physical_measurement_cost=None, fallback_cost=None, comment='Known computation phases instrumented. No measured full cost or 10x claim; inherited upstream source development UNKNOWN.')
    cost['strong_hgb_control'] = json.load(open(P / 'raw/STRONG_CONTROL_R1.json'))
    log = (P / 'raw/extract.log').read_text()
    rss = re.findall('peakRSS (\\d+)', log)
    cost['geometry_peak_rss_kib'] = int(rss[-1]) if rss else None
    control_runs = []
    for logname in ['raw/strong_control.log', 'raw/run_all_02.log', 'raw/run_all_03_FINAL.log', 'raw/run_all_04_R4_FINAL.log']:
        lp = P / logname
        if lp.exists():
            for seconds in re.findall('Actual conventional HGB refit: all9fields TIE ([0-9.]+)', lp.read_text()):
                control_runs.append(dict(log=logname, fit_query_seconds=float(seconds)))
    cost['explicit_control_replays'] = control_runs
    cost['known_explicit_control_seconds_sum'] = sum((v['fit_query_seconds'] for v in control_runs))
    cost['full_project_cost_status'] = 'UNKNOWN: coding/source inspection/manual reading, plotting/dependency repairs and some replay overhead not completely instrumented. Repeated measured comparator fits are charged explicitly, not discarded.'
    write('COSTS.json', cost)
    refs = dict(r1=r1['external_referent'], r2=r2['external_referent'], r3=r3['external_referent'], r4=r4['external_referent'])
    result = dict(schema='X7-results-v1', lane='X7-bite2text', review_state='PENDING_INDEPENDENT_REVIEW', external_referent=r1['external_referent'], external_referents=refs, rounds={'R1': r1, 'R2': r2, 'R3': {k: v for (k, v) in r3.items() if k != 'rows'}, 'R4': {k: v for (k, v) in r4.items() if k != 'rows'}}, population=dict(paired_arches=m['paired_cases'], eligible=m['eligible_cases'], train=len(m['train']), calibration=len(m['calibration']), test=len(m['test']), geometry_errors=len(err)), parser_manual_review=pr, parser_holdout_review=hold, executed_matched_control='raw/STRONG_CONTROL_R1.json', full_clinical_bite_diagnosis='UNKNOWN_NOT_ESTABLISHED', numeric_landmark_accuracy='UNKNOWN_NO_INDEPENDENT_MM_GROUND_TRUTH', contact='Sampled point-distance upper-bound candidates; no exact surface/loaded contact or force prediction', large_arrays=[], intermediate_budget_bytes=3000000000, costs=cost, scientific_gain='Geometry adds patient-specific categorical information; calibrated sets retain independent report disagreement; numerical orientation repaired. Strong matched algorithms TIE; no methodological novelty or10x claim.', source_zip_manifest='raw/DATA_MANIFEST.json', frozen_predictions=['FROZEN_PREDICTIONS.json', 'FROZEN_PREDICTIONS_R2.json', 'FROZEN_PREDICTIONS_R4.json'])
    write('results.json', result)
    intro = _release_expand(f"Two scans in registered bite contains information about the bite relationship. Demot measures arc profiles and provides categories with a variety of possible clinics assessments. The angle fields for molar/corner are statistical text proxies: the actual hive, the furrow and the dental identity are still unknown .\n\nRun from this directory: `./run_all.sh` . For a couple: `python3 code/diagnose.py --upper upper.stl --lower lower.stl --out assessment.json` . For a local facit case: `python3 code/diagnose.py --case @DENTAL_CASE_ID@ --out assessment.json` . JSON specifies the category set, supporting measures, uncertainty and when review is needed. One-point measurements from a lab can be read with `--landmarks` ; the format is available in LANDMARK_SCHEMA.json . Python3 with NumPy , SciPy , scikit-learn, trimesh and matplotlib needed; the current environment is verified. No network is needed at rerun.\n\nReference is the English IOS reports of clinics in the local Bite2Text -zippen:`{m['zip']}::*/reports_ios_en/*.txt`. It's independent human category assessments, not a numerical reference for millimeters. Dataset source: [Bite2Text ](https://ditto.ing.unimore.it/bite2text/). 994 scan pairs exist; 989 has English IOS reports. Urvalet fryser {len(m['train'])} training, {len(m['calibration'])} calibration and {len(m['test'])} test patients without overlap. The number of reports must not be counted as more independent patients.\n\n{table}\n*Report report applies to first and second readable report per patient in the entire available population, with different n per field ( raw/REPORT_DISAGREEMENT.json ). It is a description of the same case, not a universal model accuracy limit. BA weighs categories equally; kappa corrects for marginal frequencies.\n\nR1 : s frozen improvement gate for majority training: **{r1['primary_gate']}**, macro improvement {r1['macro_balanced_lift']:.3f} BA against the requirement 0,10. The majority rule is a comparison to contentless standard report, not a measured clinician's practice. The strong, equally informed gradient model is TIE .\n\nR2 : s insecure gate: **{r2['primary_gate']}**. All reports are covered in {r2['macro']['all_report_patient_coverage']:.3f} of patients in field agents; the average size of the category is: {r2['macro']['mean_size']:.2f}, andelen enstaka kategorier {r2['macro']['singleton_fraction']:.3f}. Frozen requirements: coverage ≥ 0,85 , size ≤ 2 and single categories ≥ 0,20 . Common split conformal with the same patient maximum gives the exact same response (TIE). Under different centres/protocols is the coverage guarantee UNKNOWN .\n\nR4 tests the entire bite simultaneously. R2 : s separate amounts cover all available quantities /reports together in just **{joint['all_available_fields_all_reports_covered']:.3f}** of the test cases. Joint patient calibration increases this to **{r4['joint_patient_coverage']:.3f}**, but gives**{r4['mean_field_set_size']:.2f}**categories per quantity and only**{r4['singleton_fraction']:.3f}** enstaka kategorier. R4 : s frozen information gate is therefore **FAIL**. The same test population has already been seen in previous rounds: this is a sequential descriptive replay, not new independent confirmation. Figure: [whole bet balance](figures/joint_bite_tradeoff.png).\n\nA concrete remaining error is @DENTAL_CASE_ID@ : the amount for left molar contains only Class I, while the clinic text indicates Class III . The geometric AP proxy is 0,329 mm and does not identify the actual M1 concentration/sheep. A single category is therefore not a certificate of correctness.\n\nR3 : the published axis based frame changes the dimensions by up to {r3['max_published_error_mm']:.3f} mm in joint rotation of three real pairs. A continuous vector frame provides the highest {r3['max_candidate_error_mm']:.3g} mm. This is a numerical rotation check, not anatomical mm - reference. The R1 model is trained on the published frame; the vector frame is not yet used for its category predictions.\n\nWhat does not holds : no independent toothland marks, no physical scanner error budget and no loaded contact forces are available. Proximity from sampled vertex points is exported as candidates ; an upper distance limit can never exclude contact. Laterality follows published framework convention without independent FDI verification in this archive. Spee and midlines can follow the wrong landmarks; missing/mixed teeth breaks the arc simplification. Parserns utvecklingsurval gav {pr['ours']['n_correct']}/{pr['ours']['n_emitted']} correct statements and {pr['ours']['coverage']:.3f} coverage. A second post-selection after code freezing gave {hold['precision']:.3f} precision, {hold['coverage']:.3f} coverage; reverse Spee was misinterpreted. None of these controls is an independent orthodontist review. Twelve test texts were seen by the producer of parseraudit before statistical prediction freezing; no rules or fit settings were changed based on them. Denna exponering redovisas i raw/MANUAL_GOLD_HOLDOUT.json.\n\nWhere the comparisons differ and what can actually be explained is found in raw/ERROR_ANALYSIS.json: match against another report, missing geometric measure, or continued UNKNOWN between landmark, model and parser errors. The port to X2 / NEXT_P is located in OCCLUSION_PORT.json.\n\nThree injected errors rejects controls: error category, negative calibration threshold and non-ortogonal frame . Details are available in VERIFICATION.json and raw results. The next measurement is two assessors' independent cusp/sheep/embrassure points in withheld pairs, see MEASUREMENT_SPEC.md . capability – per tooth bite diagnosis with physical measurement uncertainty – is still open.\n\nLicense: no README / LICENSE files are available in the local Bite2Text -zippen. Local briefen probably indicates CC BY - NC - SA ; the official public data page does not explicitly license. **Exact Bite2Text Terms UNKNOWN**, i.e. no license is relabelled as verified. Only private local research here. Pinned published code has no LICENSE file in the source tree; use/re-distribution conditions UNKNOWN . No second datasets are used and no photos/personal data are exported. Method source: [George et al., arXiv: 2609.13237 ](https://arxiv.org/abs/2609.13237), with exact codehashes in sources/upstream_manifest.json ; no method news is claimed.\n\nFigure: [PNG](figures/bite_geometry_vs_reports.png), [PDF](figures/bite_geometry_vs_reports.pdf). Cost and limitations: COSTS.json . Separate refit of the equally informed standard model is available in raw/STRONG_CONTROL_R1.json. Status: PENDING_INDEPENDENT_REVIEW.\n")
    (P / 'README_DEMO.md').write_text(intro)
    (P / 'RESULTS.md').write_text(intro)
    hand = _release_expand(f"""X7-bite2text — PENDING_INDEPENDENT_REVIEW\n\nExecutable capability : local ZIP → geometric measure → categories and report-disagreement sets → clinic report comparison. Run ./run_all.sh . All 989 eligible pairs processed, {len(err)} geometrifel; 593/198/198 patientdelning.\n\nR1 {r1['primary_gate']}: makro BA-lift {r1['macro_balanced_lift']:.6f} against 0,10 ; the strong model TIE . R4 {r4['primary_gate']}: joint coverage {r4['joint_patient_coverage']:.6f}, size {r4['mean_field_set_size']:.6f}, singleton {r4['singleton_fraction']:.6f}. R2 {r2['primary_gate']}: makro all-report coverage {r2['macro']['all_report_patient_coverage']:.6f}, size {r2['macro']['mean_size']:.6f}, singleton {r2['macro']['singleton_fraction']:.6f}; lika informerad patient-max conformal TIE. R3 rigid-rotation diagnostic {r3['candidate_pass']}: max {r3['max_candidate_error_mm']:.8g} mm vs published {r3['max_published_error_mm']:.6f} mm ; no anatomical validation.\n\nWhat failed : informative uncertainty for the entire bite ( R4 ), and the published frame is rotational sensitive; numerical Angle-landmarks missing; reverse-Spee-parser error remains. Manual source transcription @DENTAL_CASE_ID@ midlines corrected from incorrectly deviated to null, original preserved. Precise millimeter fidelity and clinical diagnosis UNKNOWN. Category amounts must not be presented as scanner error intervals or safe patient diagnosis. Contact files are point proximity candidates, not R6 -triangelgap or force.\n\nNext essential construction: read external tooth-identified point measurement via code/diagnose.py --landmarks and freezer comparisons before the new markup run according to MEASUREMENT_SPEC.md . Auto-segmentation can transfer Teeth3DS morphology but is a new test, no loose pre-requisite. For rotational repair in the prediction model: full vector frame extraction + NY patient partition/gate, not readjustment of these test results. Common vector frame is strong control and must get the same info.\n\nGraf: DENT - VAL - OCCLUSAL - REGIONAL - FORCES package read but its quantity is power share, which does not fit category assessments. Missing coverage proposed DENT - VAL - BITE2TEXT - REPORT - AGREEMENT (proposal only). No dispatch/feedback-receipt or source graphed mission created. Native/generated graph untouched.\n\nExactly reference : local zip members and reporthashar in DATA_MANIFEST , ALL_STRUCTURED_REPORTS , TEST_LABELS . PREREG / R1 / R2 / R3 / R4 , FROZEN_PREDICTIONS , VERIFICATION , COSTS and raw results preserve each outcome. Ingen10x and no clinical recommendation.\n""")
    (P / 'HANDOFF.md').write_text(hand)
    write('GRAPH_COVERAGE_PROPOSAL.json', dict(proposed_id='DENT-VAL-BITE2TEXT-REPORT-AGREEMENT', status='PROPOSAL_ONLY', parent_consumer='DENT-VAL-OCCLUSAL-REGIONAL-FORCES', claim='Categorical bite geometry vs clinician IOS narrative, with patient-disjoint predictions and report-disagreement sets', missing_dependencies=['FDI cusp/groove/embrasure landmarks', 'scanner error', 'operator frame uncertainty'], result_file='results.json', native_mutation=False))
    write('GRAPH_FEEDBACK.json', dict(review_state='PENDING_INDEPENDENT_REVIEW', status='PENDING_INDEPENDENT_REVIEW', target_id=None, coverage_proposal='GRAPH_COVERAGE_PROPOSAL.json', parent_context='DENT-VAL-OCCLUSAL-REGIONAL-FORCES', result_path=str(P / 'results.json'), result_sha256=sha(P / 'results.json'), measured_quantity='Per-field categorical agreement and simultaneous readable report coverage', units='dimensionless', uncertainty='Report disagreement measured; numeric/physical uncertainty UNKNOWN', population_regime='989 local IOS paired cases, 198 held-out test patients', frozen_gate='PREREG_R1/R2/R3/R4.json', baseline='Training-majority, simple thresholds, equal-info HGB and patient-max split conformal; methodological TIE', outcome=dict(R1=r1['primary_gate'], R2=r2['primary_gate'], R3_diagnostic=r3['candidate_pass'], R4=r4['primary_gate']), negative_result=True, dispatch_receipt=None, reason_no_binding='Existing regional-force target is a different measured quantity; proposing missing coverage as instructed'))
    files = [p for p in P.rglob('*') if p.is_file() and p.name != 'ARTIFACT_MANIFEST.json' and ('__pycache__' not in p.parts)]
    write('ARTIFACT_MANIFEST.json', dict(files=[dict(path=str(p.relative_to(P)), bytes=p.stat().st_size, sha256=sha(p)) for p in sorted(files)], large_arrays=[], total_bytes=sum((p.stat().st_size for p in files))))
    state('ROUND_HANDOFF_READY', 'R1 ' + r1['primary_gate'] + '; R2 ' + r2['primary_gate'] + '; R4 joint informative diagnosis ' + r4['primary_gate'] + '; vector-frame numerical PASS; matched controls TIE', 'Next construction: independent per-tooth landmark measurements, or full vector-frame retraining with fresh holdout. Parent capability remains OPEN.')
    print('Packaged reader-facing demo, result hashes and graph coverage proposal', flush=True)
if __name__ == '__main__':
    run()
