"""One-command local reconstruction, figure, facit and complete result index."""
import csv, json, resource, time
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from extract_sources import ROOT, save, sha, source_manifest
import analyze_r2, analyze_r3, analyze_r4, verify_demo

def figure(r2, r3):
    rows = json.loads((ROOT / 'REGIONAL_METHOD_CONTRASTS.json').read_text())
    (fig, axes) = plt.subplots(2, 2, figsize=(13, 9), layout='constrained')
    colors = {'Amann': '#2463a8', 'Dentsply_Sirona': '#d47712'}
    for (sys, col) in colors.items():
        rr = [r for r in rows if r['system'] == sys]
        x = np.array([(float(r['replica_mean_um']) + float(r['ct_mean_um'])) / 2 for r in rr])
        y = np.array([float(r['bias_replica_minus_ct_um']) for r in rr])
        axes[0, 0].plot(x, y, 'o-', label=sys.replace('_', ' '), color=col)
        for (xx, yy, r) in zip(x, y, rr):
            axes[0, 0].annotate(r['region'], (xx, yy), xytext=(4, 7), textcoords='offset points', fontsize=8)
    axes[0, 0].axhline(0, color='gray', lw=0.8)
    axes[0, 0].set(xlabel='Mean of two reported cohort means (µm)', ylabel='Replica minus dry CT cohort mean (µm)', title='Published regional contrasts; specimen LoA UNKNOWN')
    axes[0, 0].legend(fontsize=9)
    categories = ['Marginal ≤120', 'Marginal ≤100', 'Axial 50–100', 'Occlusal 50–100']
    for (sys, col, pos) in [('Amann', colors['Amann'], -0.18), ('Dentsply_Sirona', colors['Dentsply_Sirona'], 0.18)]:
        sc = r2['x13']['conditional_scenarios']
        order = [('marginal', 120), ('marginal', 100), ('axial', 100), ('occlusal', 100)]
        vals = [next((x['conditional_flips'] for x in sc if x['donor_system'] == sys and x['region'] == reg and (x['upper_um'] == thr))) for (reg, thr) in order]
        bars = axes[0, 1].bar(np.arange(4) + pos, vals, width=0.34, color=col, label=sys.replace('_', ' '))
        axes[0, 1].bar_label(bars, padding=3)
    axes[0, 1].set_xticks(np.arange(4), categories, rotation=14, fontsize=9)
    axes[0, 1].set(ylim=(0, 10), ylabel='Conditional cell flips (count)', title='Unvalidated offset scenarios: 9 cells per region')
    axes[0, 1].text(0.02, 0.97, '27 / 148 cells in scenario scope; 121 rejected\nAll 148 lack validated target correction\n50–100 µm is an X13 engineering scenario', transform=axes[0, 1].transAxes, va='top', fontsize=9, bbox=dict(facecolor='white', alpha=0.9, edgecolor='none'))
    s = r2['sufficiency_test']['states']
    y = [0, 1]
    for (j, st) in enumerate(s):
        (lo, hi) = st['nominal_normal_loa_um']
        axes[1, 0].plot([lo, hi], [j, j], lw=5, color=['#368554', '#b94646'][j])
        axes[1, 0].plot(st['bias_um'], j, 'ko')
    axes[1, 0].set_yticks(y, ['Matched pairing', 'Reversed pairing'])
    axes[1, 0].set(xlabel='Nominal normal sample agreement limits (µm)', title='Proof witness: identical marginal mean + SD')
    axes[1, 0].text(0.03, 0.97, 'Summary identity error = 0.0\n120 µm classification error: 0/4 → 4/4\nConstructed witness; no measured crowns', transform=axes[1, 0].transAxes, va='top', fontsize=9)
    axes[1, 0].set_ylim(-0.5, 1.7)
    device = json.loads((ROOT / 'SAME_FILM_DEVICE_CONTRASTS.json').read_text())
    labels = [r['tooth_type'].split()[0] + ' ' + r['region'] for r in device]
    values = [float(r['phone_minus_microscope_mean_um']) for r in device]
    axes[1, 1].barh(np.arange(len(device)), values, color='#5c7499')
    axes[1, 1].set_yticks(np.arange(len(device)), labels, fontsize=8)
    axes[1, 1].invert_yaxis()
    axes[1, 1].set(xlabel='Smartphone minus microscope cohort mean (µm)', title='Same film: detector contrast isolated in source')
    axes[1, 1].set_xlim(0, 9)
    fig.suptitle('X59 — observation state is required for a replica ↔ micro-CT bridge', fontsize=15)
    for ax in axes.flat:
        ax.spines[['top', 'right']].set_visible(False)
    fig.savefig(ROOT / 'METHOD_BRIDGE_DEMO.png', dpi=180)
    fig.savefig(ROOT / 'METHOD_BRIDGE_DEMO.pdf')
    plt.close(fig)

def facit():
    rr = []
    for fn in ['raw/CUNALI_REGION_SUMMARIES.json', 'raw/PASHA_REGION_SUMMARIES.json']:
        for r in json.loads((ROOT / fn).read_text()):
            rr.append(dict(source=r['source'], system=r['system'], region=r['region'], method_1='silicone replica', method_2='CT', mean_1_um=r['replica_mean_um'], sd_1_um=r['replica_sd_um'], mean_2_um=r['ct_mean_um'], sd_2_um=r['ct_sd_um'], n_copings=r['n_specimens'], pair_level=r['pair_level'], resolution=r['resolution'], locator=r['source_locator'], source_path=r['source_path'], source_sha256=r['source_sha256']))
    for r in json.loads((ROOT / 'raw/SMARTPHONE_REGION_SUMMARIES.json').read_text()):
        rr.append(dict(source=r['source'], system=r['tooth_type'], region=r['region'], method_1='smartphone', method_2='microscope', mean_1_um=r['phone_mean_um'], sd_1_um=r['phone_sd_um'], mean_2_um=r['microscope_mean_um'], sd_2_um=r['microscope_sd_um'], n_copings='10 per tooth type; reported SD/CI pooling unit unclear', pair_level=r['pair_level'], resolution=r['resolution'], locator=r['source_locator'], source_path=r['source_path'], source_sha256=r['source_sha256']))
    with (ROOT / 'FACIT.csv').open('w') as f:
        w = csv.DictWriter(f, fieldnames=list(rr[0]))
        w.writeheader()
        w.writerows(rr)
    return len(rr)

def run():
    start = time.perf_counter()
    before = resource.getrusage(resource.RUSAGE_SELF)
    r2 = analyze_r2.run()
    r3 = analyze_r3.run()
    r4 = analyze_r4.run()
    checks = verify_demo.run()
    nfacit = facit()
    figure(r2, r3)
    source_manifest()
    after = resource.getrusage(resource.RUSAGE_SELF)
    child = resource.getrusage(resource.RUSAGE_CHILDREN)
    accounting = dict(wall_seconds=time.perf_counter() - start, cpu_seconds=after.ru_utime + after.ru_stime - before.ru_utime - before.ru_stime, child_cpu_seconds=child.ru_utime + child.ru_stime, peak_self_rss_MiB=after.ru_maxrss / 1024, peak_child_rss_MiB=child.ru_maxrss / 1024, threads_max=4, gpu_used=False, intermediate_limit_bytes=3000000000, large_arrays_over_50MB=[], full_cost=dict(preparation='Replay includes local table extraction/source hash validation; original metadata/file download/search commands are in COMMANDS.md and raw manifests.', fit='No specimen calibration fit exists; ecological OLS is a source diagnostic only.', discovery='Human reasoning/model token cost and historical harvest UNKNOWN. Search command wall times were observed but not aggregated into CPU total.', validation='External cohort tables/figure plus exact algebra; physical measurements and target specimen validation unavailable.', queries='All 148 X13 cells enumerated in replay.', fallback='Physical lab/phantom costs and instrument acquisition time UNKNOWN, not zero.'))
    save('RUN_ACCOUNTING.json', accounting)
    screen = json.loads((ROOT / 'SOURCE_SCREEN.json').read_text())
    sourcebias = {s: {r['region']: r['bias_replica_minus_ct_um'] for r in json.loads((ROOT / 'REGIONAL_METHOD_CONTRASTS.json').read_text()) if r['system'] == s} for s in ['Amann', 'Dentsply_Sirona']}
    index = dict(lane='X59-replica-microct', claim_type='information_link', review_state='PENDING_INDEPENDENT_REVIEW', outcome='REGIONAL_METHOD_STATE_CONTRASTS_MEASURED_IN_SOURCE; SPECIMEN_AGREEMENT_AND_X13_TRANSFER_UNKNOWN', external_referent=dict(kind='independent_measurement', locator='https://doi.org/10.1590/0103-6440201601531 Table1 p470; https://doi.org/10.7759/cureus.40020 TAB1/TAB2; https://doi.org/10.1038/s41598-024-55711-4 Tab2', compared_quantity='independently reported matched-cohort region means, NOT specimen pair observations or absolute true gaps', refutes_us=True), figshare=dict(locator='https://doi.org/10.6084/m9.figshare.5669308.v1', license='CC BY 4.0', declared_total_bytes=39306, files=5, downloaded=0, status='ALL_DECLARED_FILE_URLS_HTTP_404_ORIGINAL_PUBLICATION_PDF_AVAILABLE'), source_bias_replica_minus_ct_um=sourcebias, source_bias_resolution='PER_SURFACE_REGION', individual_bland_altman=dict(bias='cohort mean difference identified under equal cohort weighting', limits_of_agreement=None, reason='No within-region specimen covariance/registered pairs; ecological region correlation is not specimen covariance', nominal_sample_identified_sets='REGIONAL_METHOD_CONTRASTS.json'), x13=r2['x13'], facit_summary_pairs=nfacit, facit_raw_specimen_pairs=0, source_exclusions=dict(keyword_candidate_files=screen['initial_replica_keyword_files'], admitted=2, rejected=screen['rejected_sources'], rejected_fraction=screen['rejected_fraction'], pasha_missing_rows=2, pasha_candidate_rows=12), sufficiency_tests=dict(marginal_summary=r2['sufficiency_test'], original_two_observations=r4['sufficiency_test']), frozen_gates=dict(R1='FAIL_PAIRED_COVARIANCE_AND_REGISTRATION', R2=dict(Amann_OLS='PASS', Dentsply_OLS='FAIL_intercept_error_0.005404_um_exceeds_frozen_0.005_um'), R3='ALL_6_OUT_OF_DOMAIN_STRESS_ERRORS_FAIL_20UM; DEVICE_SOURCE_CONTRAST_IDENTIFIED', R4='EXACT_CONTRAST_RANK_PASS_CONDITIONAL_ON_CT_STATE_CALIBRATION_PHYSICAL_VALIDATION_UNKNOWN'), rigorous_enclosures=dict(source_statistics='Ordinary rounding boxes and Cauchy-Schwarz SAMPLE SD bounds', affine_source_OLS='Exact rational coefficient enclosure over continuous rounding box; not a specimen predictor', target_error=None, empirical_population_coverage=None), graph_edges_proposed=[dict(source='Cunali Table1 region/system cohort', consumer='X13 matching region cells -> K43 / DENT-IF-RESTORATION-CEMENT-TOOTH context', resolution='PER_SURFACE_REGION', timescale='SIMULTANEOUS', state_semantics='Static comparison of alternative seating states; no physical same-state measurement equivalence', status='PENDING_INDEPENDENT_REVIEW_SCOPE_CONDITIONAL'), dict(source='PMC10912410 same-film optical device contrast', consumer='Replica imaging sub-operator in lab observation port', resolution='PER_SURFACE_REGION', timescale='SIMULTANEOUS', state_semantics='Same sectioned film; apparatus contrast', status='PENDING_INDEPENDENT_REVIEW')], full_cost=accounting, checks=dict(status='PASS', count=len(checks), scope='Software contracts and exact algebra; no scientific admission'), artifacts=dict(facit='FACIT.csv', figure='METHOD_BRIDGE_DEMO.png', decisions='X13_DECISIONS_R2.json', lab_api='code/lab_bridge.py', raw_source_hashes=json.loads((ROOT / 'SOURCE_MANIFEST.json').read_text())), next_construction='Actual registered same-coping four-observation protocol with traceable CT contrast phantom, >=10 independent copings per source protocol/region plus held-out batch; freeze/score via code/lab_bridge.py')
    save('results.json', index)
    freeze = ROOT / 'FROZEN_PREDICTIONS_R3.json'
    save('FROZEN_PREDICTIONS.json', dict(type='INDEX_OF_ORIGINAL_RETROSPECTIVE_FREEZE', file=freeze.name, sha256=sha(freeze), created_utc=json.loads(freeze.read_text())['created_utc'], prospective=False, future_actual_lab_predictions='Not generated: no actual calibration measurements supplied; lab_bridge.py creates non-overwritable prospective freezes from actual CSV.'))
    (ROOT / 'FROZEN_PREDICTIONS.sha256').write_text(sha(ROOT / 'FROZEN_PREDICTIONS.json') + '\n')
    print(json.dumps(dict(status='DEMO_COMPLETE_WITH_SCIENTIFIC_UNKNOWNS', facit_summary_pairs=nfacit, wall_seconds=accounting['wall_seconds'], peak_rss_MiB=accounting['peak_self_rss_MiB'])))
if __name__ == '__main__':
    run()
