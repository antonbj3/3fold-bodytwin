"""One-command, offline, uncertainty-aware research planning demonstration."""
from dental_release.paths import expand as _release_expand
import argparse
from datetime import datetime, timezone
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import resource
import sys
import time
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import scipy
import bs4
import controls
from model import r1_result, r2_result
ROOT = Path(__file__).resolve().parents[1]
LABELS = {'normal': 'Normal bild', 'widening': "Light dilation", 'lesion': 'Lesion'}

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def dump(path, x):

    def native(v):
        if isinstance(v, np.generic):
            return v.item()
        raise TypeError(type(v).__name__)
    path.write_text(json.dumps(x, indent=2, ensure_ascii=False, allow_nan=False, default=native) + '\n')

def verify_frozen():
    for name in ['PREREG_R1', 'PREREG_R2', 'PREREG_R2_EXTERNAL_CODE_CONTROL', 'FROZEN_PREDICTIONS']:
        assert sha(ROOT / (name + '.json')) == (ROOT / (name + '.sha256')).read_text().split()[0], name
    manifest = json.loads((ROOT / 'INPUT_MANIFEST.json').read_text())
    for r in manifest:
        assert sha(ROOT / r['file']) == r['sha256'], r['file']

def make_figure(r1, r2):
    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False})
    (fig, axs) = plt.subplots(1, 2, figsize=(11, 4.9))
    x = np.arange(3)
    p = np.array([r['p'] for r in r1['rows']]) * 100
    lo = np.array([r['CI95'][0] for r in r1['rows']]) * 100
    hi = np.array([r['CI95'][1] for r in r1['rows']]) * 100
    axs[0].errorbar(x, p, yerr=[p - lo, hi - p], fmt='o', color='#156f77', capsize=5, markersize=8, label='Patel: CBCT-category frequency + 95% CP')
    axs[0].axhline(100 * 25 / 86, color='#777777', linestyle='--', label='CBCT-blind cohort: 25/86')
    for (j, r) in enumerate(r1['rows']):
        axs[0].text(j, r['p'] * 100 + 3, f"{r['k']}/{r['n']}", ha='center')
    axs[0].set(xticks=x, xticklabels=['Normal', 'Widening', 'Lesion'], ylim=(0, 100), ylabel='Intraoperative conversion (%)', title='Published study frequencies · POPULATION')
    axs[0].legend(loc='upper left', fontsize=8)
    a = r1['pooled']
    b = r1['external_comparisons'][0]
    for (j, r) in enumerate([a, b]):
        axs[1].errorbar(j, r['p'] * 100, yerr=[[100 * (r['p'] - r['CI95'][0])], [100 * (r['CI95'][1] - r['p'])]], fmt='o', color=['#156f77', '#aa5939'][j], capsize=5, markersize=8)
        axs[1].text(j, r['p'] * 100 + 3, f"{r['k']}/{r['n']}", ha='center')
    axs[1].set(xticks=[0, 1], xticklabels=['Patel\nconfirmed RCT conversion', 'Baranwal\nbleeding-related nonreceipt'], ylim=(0, 60), ylabel='Bleeding-related pathway event (%)', title='Independent study challenge · POPULATION')
    axs[1].text(0.5, 49, 'Different selection/case mix;\nRCT receipt not reported in Baranwal.\nNo external CBCT calibration.', ha='center', fontsize=9)
    fig.suptitle('Preoperative CBCT opens a cohort-bound planning question', fontsize=13)
    fig.text(0.5, 0.015, 'Individual 95% Clopper–Pearson intervals; IID assumption. 1/86 missing CBCT retained separately. No treatment recommendation.', ha='center', fontsize=8)
    fig.tight_layout(rect=(0, 0.055, 1, 0.95))
    fig.savefig(ROOT / 'figures/CBCT_PATHWAY.png', dpi=180)
    fig.savefig(ROOT / 'figures/CBCT_PATHWAY.pdf')
    plt.close(fig)

def make_html(r1):
    data = json.dumps(r1['rows'], ensure_ascii=False)
    template = (ROOT / 'code/planning_template.html').read_text()
    (ROOT / 'planning_demo.html').write_text(template.replace('__DATA__', data))

def main():
    started = time.perf_counter()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--category', choices=['all', 'normal', 'widening', 'lesion', 'unknown'], default='all')
    parser.add_argument('--base-min', nargs=2, type=float, default=[60, 60], metavar=('LOW', 'HIGH'))
    parser.add_argument('--extra-min', nargs=2, type=float, default=[30, 30], metavar=('LOW', 'HIGH'))
    parser.add_argument('--base-cost', nargs=2, type=float, default=[1000, 1000], metavar=('LOW', 'HIGH'))
    parser.add_argument('--extra-cost', nargs=2, type=float, default=[2000, 2000], metavar=('LOW', 'HIGH'))
    args = parser.parse_args()
    for values in [args.base_min, args.extra_min, args.base_cost, args.extra_cost]:
        if not all((np.isfinite(v) and 0 <= v <= 1000000000.0 for v in values)) or values[0] > values[1]:
            parser.error('Resource boxes require finite ordered nonnegative endpoints <=1e9.')
    verify_frozen()
    source = json.loads((ROOT / 'raw/PRIMARY_INPUTS.json').read_text())
    external = json.loads((ROOT / 'raw/INDEPENDENT_INPUTS.json').read_text())
    atlas = json.loads((ROOT / 'raw/ATLAS_CONTEXT.json').read_text())
    r1 = r1_result(source, external)
    r2 = r2_result(r1, args.base_min, args.extra_min, args.base_cost, args.extra_cost)
    checks = controls.run(source, external, (ROOT / 'sources/PRIMARY_TABLE4.html').read_text(), (ROOT / 'sources/PMC9274703.xml').read_text(), atlas, r1, r2, args.base_min, args.extra_min, args.base_cost, args.extra_cost)
    make_figure(r1, r2)
    make_html(r1)
    result = dict(lane=_release_expand('X88'), round_tag='X88-pulpotomy-decision', claim_type='information_link', review_state='PENDING_INDEPENDENT_REVIEW', scientific_admission=False, outcome='SOURCE_COHORT_RISK_AND_CONDITIONAL_PLANNING_DEMO; EXTERNAL_CBCT_CALIBRATION_UNKNOWN', external_referent=dict(kind='independent_measurement', locator='https://pmc.ncbi.nlm.nih.gov/articles/PMC11629050/#iej14144-tbl-0004', compared_quantity='Assigned-pulpotomy CBCT-category x confirmed intraoperative RCT conversion counts', refutes_us=False), external_referents=[dict(kind='independent_measurement', locator='https://pmc.ncbi.nlm.nih.gov/articles/PMC9274703/#F1', compared_quantity='Full-pulpotomy bleeding-related nonreceipt, 1/33, with1 partial-necrosis exclusion; clinical RCT conversion not reported', refutes_us=True, refuted_claim='Unqualified transport of fixed pooled29.1% to all full-pulpotomy cohorts; descriptive challenge, unmatched selection prevents matched gate.'), next((x['external_referent'] for x in checks['checks'] if x['name'] == 'published_scipy_PMF'))], rounds=[r1, r2], controls=checks, atlas_context={k: v for (k, v) in atlas.items() if k != 'rows'}, information_comparator=dict(name='CBCT-blind pooled source probability25/86', status='Explicit information ablation; actual current-practice behavior not measured.', source_frequency_change_pp={r['category']: 100 * (r['p'] - 25 / 86) for r in r1['rows']}), uncertainty=dict(probability='Source-cohort individual CP95; simultaneous Bonferroni over3 categories and missing-case envelope; IID conditional', target_site='UNKNOWN', empirical_CBCT_times='UNKNOWN', clinical_cost='UNKNOWN', affine_enclosure='Exact corner enclosure over declared resource and computed probability boxes; no rigorous beta-quantile rounding enclosure.', reserve='Independent future conversions conditional on source-compatible scenario; sampling confidence and future95% tail are distinct.'), resolution=dict(study_counts_and_risks='POPULATION', preop_input='PER_TOOTH', atlas_geometry='PER_TOOTH', resource_and_schedule_scenarios='PHENOMENOLOGICAL'), edges=[dict(producer='Patel Table4 CBCT+assigned pathway', consumer='preoperative pathway-planning research', resolution='POPULATION', input_observation_resolution='PER_TOOTH', time_scale='HANDOVER', status='SOURCE_BOUND'), dict(producer='X12 corrected atlas R8', consumer='paired geometry+preop/intraop acquisition contract', resolution='PER_TOOTH', time_scale='SIMULTANEOUS', status='ANATOMY_CONTEXT_ONLY; no disease/clinical join')], dropout=dict(primary_CBCT_missing_fraction=1 / 86, independent_study_screen=json.loads((ROOT / 'sources/SEARCH_SCREEN.json').read_text())['dropout'], atlas_domain_truncated_fraction=90 / 3539, atlas_clinical_risk_join_rejected_fraction=1.0), scenarios=dict(time_unit='min', cost_unit='user-selected currency units', base_time=args.base_min, extra_time=args.extra_min, base_cost=args.base_cost, extra_cost=args.extra_cost, evidence_type='PHENOMENOLOGICAL user inputs; defaults illustrative, not measurements'), patient_information_example="In the published study, some planned pulmonaryotomies were changed to root treatment when the bleeding could not be stopped according to the protocol. The frequency of the image group is a study task; your individual probability is not yet calibrated. Later treatment, cost and times require local basis.", full_cost=dict(replay_seconds_before_serialization=time.perf_counter() - started, peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, preparation='Directed readings, reused reviewed primary counts; model/human effort not instrumented', fit='No learned model; finite binomial estimates', discovery='Own calculation frozen07:15UTC before5 primary-study screens; see SEARCH_SCREEN', validation='Source parser, CDF inversion, exact256-outcome enumeration, published code, nine fault controls', query='Single-cohort category lookup and rational corner arithmetic', fallback='Unknown empirical category time/cost/site risk; no atlas-inferred disease labels', physical_acquisition='NOT_RUN', new_imaging_download_bytes=0, new_primary_publication_snapshots='XML, selected table, one CONSORT figure only', token_cost='UNKNOWN', threads_max=4, GPU=False, heavy_job_needed=False), environment=dict(python=sys.version, numpy=np.__version__, scipy=scipy.__version__, matplotlib=matplotlib.__version__, beautifulsoup=bs4.__version__), frozen_predictions=dict(file='FROZEN_PREDICTIONS.json', sha256=sha(ROOT / 'FROZEN_PREDICTIONS.json'), nature='Source-derived transport hypothesis frozen before independent publication retrieval; no new measurement made'), input_manifest=json.loads((ROOT / 'INPUT_MANIFEST.json').read_text()), code_manifest=[dict(file=str(p.relative_to(ROOT)), sha256=sha(p)) for p in sorted((ROOT / 'code').glob('*')) if p.is_file()] + [dict(file='run_all.sh', sha256=sha(ROOT / 'run_all.sh'))])
    stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S_%fZ')
    run = ROOT / 'replay' / stamp
    run.mkdir(parents=True)
    dump(run / 'results.json', result)
    dump(ROOT / 'results.json', result)
    dump(ROOT / 'rounds/R2_RESULTS.json', r2)
    dump(ROOT / 'raw/CONTROLS.json', checks)
    selection = dict(requested_category=args.category, clinical_patient_probability=None, clinical_treatment_choice=None, empirical_time_min=None, source_rows=r1['rows'] if args.category == 'all' else [r for r in r1['rows'] if r['category'] == args.category], resource_scenarios=r2['scenarios'] if args.category == 'all' else [r for r in r2['scenarios'] if r['category'] == args.category], unknown_category_source_context=r1['pooled'] if args.category == 'unknown' else None)
    dump(run / 'query.json', selection)
    dump(ROOT / 'raw/QUERY_LATEST.json', selection)
    dump(ROOT / 'CURRENT_WORK_STATE.json', dict(round_tag='X88-pulpotomy-decision', stage='R2_COMPLETE_DEMO_VERIFIED', timestamp_utc=datetime.now(timezone.utc).isoformat(), latest_gate=r2['gates'], next_operation='Prospective same-tooth local panel and held-out-site CBCT calibration; no local clinical measurement available.', latest_result_sha256=sha(ROOT / 'results.json'), latest_replay=str(run), no_running_jobs=True))
    print(json.dumps(dict(outcome=result['outcome'], source_rows=[{k: r[k] for k in ['category', 'k', 'n', 'p', 'CI95']} for r in r1['rows']], checks_passed=checks['n_valid_pass'], mutants_rejected=checks['n_mutants_rejected'], panel=r2['panel'], results='results.json', figure='figures/CBCT_PATHWAY.png', interactive='planning_demo.html', replay=str(run)), indent=2, ensure_ascii=False))
if __name__ == '__main__':
    main()
