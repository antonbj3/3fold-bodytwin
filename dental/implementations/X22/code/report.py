from dental_release.paths import expand as _release_expand
import csv
import json
import resource
import re
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from common import *

def figure(r1, r2, r3):
    (fig, ax) = plt.subplots(2, 3, figsize=(16, 10), layout='constrained')
    case = r3['cases'][0]
    les = next((r for r in case['lesions'] if r['site'] == 'proximal' and r['class'] == 'D2'))
    view = next((r for r in les['views'] if r['angle_deg'] == 0))
    q = np.load(view['path'])
    h = 0.15
    origin = np.array(case['frame']['scene_origin_local_mm'])
    shape = q['image'].shape
    im = ax[0, 0].imshow(q['image'].T, origin='lower', cmap='gray_r', vmin=0, vmax=1.1, extent=[origin[0], origin[0] + shape[0] * h, origin[2], origin[2] + shape[1] * h])
    ax[0, 0].set(title='A. Same-scan multi-tooth DRR / simulated', xlabel='mesiodistal (mm)', ylabel='crown-axis coordinate (mm)')
    fig.colorbar(im, ax=ax[0, 0], label='flat-field normalized counts')
    labels = CLASSES[1:]
    pos = np.arange(5)
    for (name, color, label) in [('reference_logistic', '#317eae', 'Logistic reference'), ('equal_information_trees', '#d78f30', 'Same-input trees')]:
        v = [r1['detectors'][name]['by_depth'][c]['sensitivity'] for c in labels]
        ax[0, 1].plot(pos, v, 'o-', color=color, label=label)
    for i in range(5):
        (lo, hi) = (0.66, 0.75) if i < 2 else (0.8, 0.87)
        ax[0, 1].fill_between([i - 0.25, i + 0.25], [lo, lo], [hi, hi], color='#578c66', alpha=0.25)
    ax[0, 1].set(xticks=pos, xticklabels=labels, ylim=(-0.05, 1.05), ylabel='Sensitivity at frozen validation threshold', title='B. Held-out geometry detection / failed transfer')
    ax[0, 1].legend(fontsize=8)
    ax[0, 1].text(0.02, 0.98, 'Green: published proximal aggregate CI\nnot matched populations or a fitting target', transform=ax[0, 1].transAxes, va='top', fontsize=8)
    lac = r1['gates']['G4_external_attenuation']
    means = np.array(lac['external_mean_sd_cm_inv'])[:, 0]
    sd = np.array(lac['external_mean_sd_cm_inv'])[:, 1]
    ax[0, 2].errorbar([0, 1], means, yerr=2 * sd, fmt='o', capsize=6, label='Measured 70 kVp, ±2SD')
    ax[0, 2].plot([0, 1], lac['model_cm_inv'], 's', label='40 keV assigned model')
    ax[0, 2].set(xticks=[0, 1], xticklabels=['Enamel', 'Dentin'], ylabel='Linear attenuation (cm⁻¹)', title='C. External coefficient check / broad compatibility')
    ax[0, 2].legend(fontsize=8)
    ax[0, 2].text(0.05, 0.03, 'Spectrum unmatched; mineral-loss contrast UNKNOWN', transform=ax[0, 2].transAxes, fontsize=8)
    pair = r2['pairs'][0]
    q = np.load(pair['pair_file'])
    z = int(q['center'][1])
    h = float(q['spacing_mm'])
    for (col, key, title) in [(0, 'shallow_loss', 'D. E2-only constructed mineral deficit'), (1, 'deep_loss', 'E. D1-reaching deficit / same 0° image')]:
        arr = q[key][:, :, z].T
        bg = q['labels'][:, :, z].T
        ax[1, col].imshow(bg, origin='lower', cmap='Greys', alpha=0.35, extent=[0, bg.shape[1] * h, 0, bg.shape[0] * h])
        im = ax[1, col].imshow(np.ma.masked_where(arr <= 0, arr), origin='lower', cmap='magma', vmin=0, vmax=0.5, extent=[0, arr.shape[1] * h, 0, arr.shape[0] * h])
        ax[1, col].set(title=title, xlabel='mesiodistal grid (mm)', ylabel='beam direction grid (mm)')
        fig.colorbar(im, ax=ax[1, col], label='assigned fractional mineral loss')
    dp = [p['scores']['20000']['additional_10deg_view']['dprime_gaussian_approx'] for p in r2['pairs']]
    ax[1, 2].bar(np.arange(7), dp, color=['#538e6d' if v >= 2 else '#ba704f' for v in dp])
    ax[1, 2].axhline(2, color='k', linestyle='--', label='Frozen criterion 2')
    ax[1, 2].set(xticks=np.arange(7), xticklabels=[p['tid'].split('_')[0] for p in r2['pairs']], ylabel='Ideal-observer d-prime, 20 000 photons', title='F. Fixed 10° extra view separates 6/7 pairs')
    ax[1, 2].legend(fontsize=8)
    fig.suptitle('Known synthetic voxel truth; conditional projection ambiguity; clinical realism unvalidated', fontsize=15)
    fn = ROOT / 'figures/benchmark_and_ambiguity.png'
    fig.savefig(fn, dpi=160)
    fig.savefig(ROOT / 'figures/benchmark_and_ambiguity.pdf')
    plt.close(fig)
    return fn

def report():
    r1 = json.loads((ROOT / 'raw/R1_results.json').read_text())
    r2 = json.loads((ROOT / 'raw/R2_results.json').read_text())
    r3 = json.loads((ROOT / 'raw/R3_results.json').read_text())
    r4 = json.loads((ROOT / 'raw/R4_results.json').read_text())
    r5 = json.loads((ROOT / 'raw/R5_results.json').read_text())
    fn = figure(r1, r2, r3)
    logistic = r1['detectors']['reference_logistic']
    trees = r1['detectors']['equal_information_trees']
    max_tau = max((p['projector_nullspace_error'] for p in r2['pairs']))
    max_mean = max((p['mean_transmission_difference'] for p in r2['pairs']))
    acc = [p['scores']['20000']['additional_10deg_view']['sampled_accuracy'] for p in r2['pairs']]
    rows = []
    for c in CLASSES[1:]:
        rows.append([c, logistic['by_depth'][c]['n'], logistic['by_depth'][c]['sensitivity'], trees['by_depth'][c]['sensitivity']])
    with open(ROOT / 'raw/DEPTH_SCORES.csv', 'w') as f:
        w = csv.writer(f)
        w.writerow(['depth_class', 'test_n', 'reference_sensitivity', 'equal_input_trees_sensitivity'])
        w.writerows(rows)
    with open(ROOT / 'raw/AMBIGUITY_SCORES.csv', 'w') as f:
        w = csv.writer(f)
        w.writerow(['geometry', 'zero_angle_tau_error', 'one_view_Bayes_bound', 'ten_deg_dprime_20000', 'sampled_accuracy_20000'])
        for p in r2['pairs']:
            w.writerow([p['tid'], p['projector_nullspace_error'], 0.5, p['scores']['20000']['additional_10deg_view']['dprime_gaussian_approx'], p['scores']['20000']['additional_10deg_view']['sampled_accuracy']])
    manifests = r1['array_manifest'] + r2['array_manifest'] + r3['array_manifest'] + r4['array_manifest'] + r5['array_manifest']
    disk = sum((p.stat().st_size for p in DATA.rglob('*') if p.is_file()))
    overall = {'lane': 'X22-caries-synthetic', 'claim_type': 'capability', 'review_state': 'PENDING_INDEPENDENT_REVIEW', 'scientific_admission': False, 'capability_delivered': 'Known-voxel synthetic caries truth, depth-stratified ROI scoring, same-scan anatomy context and deliberate nonidentifiable tissue-depth pairs', 'n_benchmark_images': r1['n_images'], 'n_sts_geometries': r1['n_real_geometries'], 'n_native_context_views': 60, 'n_ambiguity_pairs': len(r2['pairs']), 'results': {'reference': logistic, 'equal_information_control': trees, 'R1': r1['gates'], 'R2': r2['gates'], 'R3': r3['gates'], 'R4': r4['gates'], 'R5': r5['gates'], 'old_attenuation_pairs_fail_new_material': r5['old_R2_transfer_failure_count'], 'max_R5_areal_mineral_error_g_cm2': max((r['mass_projection_error_g_cm2'] for r in r5['rows'])), 'max_identical_pair_optical_depth_difference': max_tau, 'max_identical_pair_mean_transmission_difference': max_mean, 'one_view_balanced_Bayes_accuracy_bound': 0.5, 'extra_10deg_sampled_accuracy_range_20000': [min(acc), max(acc)]}, 'external_referent': r5['external_referent'], 'external_referents_additional': [r1['external_referent'], *r1['external_referents_additional'], r2['external_referent'], r3['external_referent'], r4['external_referent']], 'empirical_realism': 'UNKNOWN', 'clinical_performance_resemblance': 'FAILED', 'physical_measurement_performed': False, 'source_lineage': {'dataset': 'STS-Tooth CC BY4.0', 'predecessors': [_release_expand('X9'), 'CROWN_pop', 'LANE_X12_PULPY3D semanticwarning only'], 'Pulpy3D_and_ToothFairy2_used': False, 'reason': 'ExpertDEJ absent; no needtoduplicatepairedatlas. Synthetic tissueclosures explicitly replaceunverified anatomy.'}, 'frozen_prediction_file': 'FROZEN_PREDICTIONS.json', 'frozen_prediction_sha256': sha(ROOT / 'FROZEN_PREDICTIONS.json'), 'artifact_hashes': {str(p.relative_to(ROOT)): sha(p) for p in [*ROOT.glob('PREREG_*.json'), *ROOT.glob('DECOMPOSITION_*.json'), ROOT / 'sources/EXTERNAL_REFERENTS.json', fn]}, 'code_manifest': {str(p.relative_to(ROOT)): sha(p) for p in [*(ROOT / 'code').glob('*.py'), ROOT / 'run_all.sh']}, 'array_manifest': manifests, 'disk_bytes': disk, 'intermediate_limit_bytes': 3000000000, 'cost': {'executed_rounds': [r1['cost'], r2['cost'], r3['cost'], r4['cost'], r5['cost']], 'sum_current_replay_wall_s': sum((r['cost']['wall_s'] for r in [r1, r2, r3, r4, r5])), 'peak_round_rss_MiB': max((r['cost']['peak_rss_MiB'] for r in [r1, r2, r3, r4, r5])), 'failed_and_superseded_runs': 'Preserved R1attempt1/2 andR3attempt1 logs/numericfiles; costsareadditional, seeCOMMANDS.md andATTEMPTS.json', 'preparation_search_token_cost': 'UNKNOWN', 'lab_acquisition_cost': 'NOT_MEASURED', 'GPU_s': 0}, 'negative_results': ['Unpaired logistic detector heldout sensitivity0% at frozen validationthreshold', 'Clinical aggregate sensitivity resemblance fails', 'Matched mineral-loss versusclinicalcontrast calibrationabsent', 'Additional10degview fails dprime2 in1/7cases', 'Expert enamel/dentin/pulp boundariesabsentfromusedsource masks'], 'next_construction': 'Replace tissue and beam closures with registered same-specimen microCT/TMR plus calibratedbitewing; prospectivepredictionssealed. Benchmarkdeploy/training efficacy or10x benefit remains unproved.'}
    write(ROOT / 'results.json', overall)
    result_hash = sha(ROOT / 'results.json')
    table = '\n'.join((f'| {c} | {n} | {se:.1%} | {ctrl:.1%} |' for (c, n, se, ctrl) in rows))
    text = f"""# Known voxelvis lesion truth and testable projection loss\n\nBenchmark points to any detector's caries integrity in a given surface, separately for E1 / E2 / D1 / D2 / D3 , mineral loss, angle and noise. It contains **{r1['n_images']:,} images from seven real STS -molar geometry**, as well as 60 local projections with neighboring and opposite jaw teeth in original CT - frame . Mineral loss is known per voxel **within the synthetic tissue model**. Simple detectors generalize poorly to the held out geometry.\n\nSeven designed lesion pairs have different tissue depths and identical expected 0° image. Towards reference in [Kak and Slaney, Chapter 3 , equation 2/3 ](https://engineering.purdue.edu/~malcolm/pct/CTI_Ch03.pdf) is their greatest difference in optical depth {max_tau:.3g}. During the model's common noise, the same image distribution gives a Baye limit of 50 % for balanced two-class classification. A fixed extra 10° view can handle the frozen ideal observer measure in 6/7 geometries. R2 : s actual Poisson and read noise simulation provides {min(acc):.1%}–{max(acc):.1%} pair accuracy at 20 000 fotoner/pixel. The observer knows both forward models; these numbers are simulation, not clinical AI performance.\n\n| Delivery / gate | Results | Giltighet |\n|---|---:|---|\n| Voxel classes and saved mineral fields | All 252 tissue/ site -/class locations passed | Synthetic DEJ and pulp |\n| Projection and noise | Analytical slab, column sum and Poisson torque passed | Deklarerad observation |\n| Ursprunglig anatomi | Target mask overlap with source mask 98,01 %/98,26 % | Two original scan masks |\n| R2:s djuppar, en vy | 7/7 identical bildmedel | Effektiv 40 keV-modell |\n| Extra 10°-vy | 6 / 7 reaches d-prime ≥ 2 | Template-sponsored ideal observer |\n| Clinic and mineral loss contrast | UNKNOWN; prestandalikhet FAIL | Matched measurement missing |\n\n| Djupklass | Held out Image Case | Logistical reference: sensitivity | Equally informed trees: sensitivity |\n|---|---:|---:|---:|\n{table}\n\nThe specificity of the reference is: {logistic['specificity']:.1%}, AUC {logistic['auc']:.3f}The specificity of the tree control is: {trees['specificity']:.1%}, AUC {trees['auc']:.3f}. Both were run with the same image field, split and validation rule. The threshold was chosen on a separate original volume for at least 90 % specificity; the anatomical distribution of the test volumes alters the outcome. No algorithm advantage or 10 × win is shown. Wilson range is available in results.json, but the simulations share geometry and are not independent clinical specimen.\n\n![Figure: anatomy, negative detector result, external reference and depth loss](figures/benchmark_and_ambiguity.png)\n\nR4 changes the mineral operation: removed mineral pulp leaves a fixed control volume element and its volume is filled with water. Organic matrix is not scaled away. The fields contain mineral pulp (g/cm³), pore volume and water inventory. Volume and mass balance, actual projection and a separate mass fraction control against [NIST : s mixture stroke](https://physics.nist.gov/PhysRefData/Xcom/Text/chap3.html) pass on the 14 case.\n\nThe old attaching match pairs lose exact similarity in **7/7 case** when R4 is used. R5 therefore changes representation to **removed mineral pulp per beam**. The new pairs provide the same expected picture at both mineral density 2,99 and 3,15 g/cm³; the largest mineral area error is {overall['results']['max_R5_areal_mineral_error_g_cm2']:.3g} g/cm². The derived common mineral/water model also gives similarity between the ends. The extra 10° view passes again in the 6/7 case at both ends. Natural lesion form and clinical imaging system are not verified. R1 : s benchmark images retain their original material model; R4 / R5 are separate structures.\n\nThe external attaching factor is [Chiu  1996 ,  PMID   8819355 ](https://pubmed.ncbi.nlm.nih.gov/8819355/):  2   mm  snitt,  enamel   2,97   ±   0,71  and  dentin   2,12   ±   0,92   cm⁻¹  at  70  kVp. The model is within the wide ± 2 SD band, but 40 keV and 70 kVp are different observation contracts. It does not validate mineral loss contrast. The water coefficient 0,2683 cm² /g at 40 keV comes from [NIST ](https://physics.nist.gov/PhysRefData/XrayMassCoef/ComTab/water.html). Mineral and matrix coefficients are provisional or assigned Closures, not measurements of these teeth ; source status is found in sources/EXTERNAL_REFERENTS.json .\n\n76 / 91 % is verified in [Abbott 2025 , DOI 10.1016/j.jebdp.2024.102077 ](https://doi.org/10.1016/j.jebdp.2024.102077), including five X-ray and two clinical imaging studies. X-ray sub-group is 73 / 91 %; numbers are nothing E1 – D3 - reference . For bitewing, [Ammer 2024 ](https://doi.org/10.1016/j.jdsr.2024.02.001) provides aggregated enamel sensitivity 0,71 [0,66 ; 0,75 ], dentin 0,84 [ 0,80 ; 0,87 ] and total specificity 0,89 [0,75 ; 0,96 ]. Reference frozen comparison with these bands fails . The data/populations differ and the straps must not be used to re-adjust the simulator. Matched published contrast values per mineral loss and five separate depth classes clinical interval is still UNKNOWN .\n\nWhat does not holds : STS lacks verified enamel -/ dentin -/pulp limits. R1 uses 0,8/1,2/1,6 mm enamel shell and a central pulp model. R1 classes follow a surface beam; R2 / R5 : s tissue depth follows distance to the nearest outer bound . The definitions are kept separate. The beams are parallel, efficient energy 40 keV, noise/PSF /spreading are scenarios and R3 : s other teeths are homogeneous. Natural lesion growth, polychromatics, bone and soft tissue anatomy are not calibrated. Precise synthetic truth does not replace an independent physical measurement.\n\nRun all with `./run_all.sh`. No network access or data set download is required. NumPy , SciPy , scikit-learn, nibabel and matplotlib are required; this machine's installed matplotlib is specified in the script. Data skrivs under `{DATA}` ({disk / 1024 ** 2:.1f} MiB), not more than 3 GB and four wires, without GPU . Precise frozen files, raw datahashar, error injections and costs are available in results.json and COMMANDS.md. The five steps take about {overall['cost']['sum_current_replay_wall_s']:.0f} s on this shared machine. Search and agent cost is UNKNOWN.\n\nOptional detector leaves a JSON list with `{{"sample_id":"…","probability":0.4}}` for each line in`{DATA}/public/index.json`. Full images can be found in `public/<geometri>.npz` , the key `images` ; the index indicates the array index and the candidate area pixel center. Private reference is located separately in `private/truth_index.json` and the voxel bases. The detector gets the picture and candidate surface, never healthy template or lesion mask.\n\n```bash\npython3 code/score_predictions.py --predictions MY_PREDICTIONS.json --threshold 0.5 --output raw/my_scores.json\n```\n\nThe All ID must be unique and complete; the rejects missing/ unknown ID and invalid probabilities. Raw data contains real CLI attempts with missing ID, double ID and probability 2; all were rejected.\n\nFor a segmented research stand: NPZ with `labels` (0 outside, 1 enamel , 2 dentin , 3 pulp ), `spacing_mm` and `origin_mm` , in isotrop physical scale. A mask without verified tissues requires explicit `--shell-mm 1.2`.\n\n```bash\npython3 code/make_case_demo.py --input MY_SEGMENTED_TOOTH.npz --output MY_OUTPUT\n```\n\nDataset/licence: only [STS -Tooth, Wang etc.](https://doi.org/10.1038/s41597-024-04306-9), [Zenodo 10.5281/zenodo.10597292 ](https://zenodo.org/records/10597292), **CC BY 4.0**, is used. Previously CROWN_pop / X9 derivatives are reused with hashia. ToothFairy2 (CC BY - SA ), Pulpy3D , Bits2Bites , mandible defects and Teeth3DS are not used. No personal data are read. Reference publications are not redistributed here.\n\nThe next design is the same extracted tooth with registered micro-CT / TMR and calibrated bitewing. `FROZEN_PREDICTIONS.json` with time stamp and hash freezes future physical contrast comparison after this digital study; no physically specimen has been measured. `LAB_VALIDATION.md` indicates how a 20 %-residual can the reject model. Status **PENDING_INDEPENDENT_REVIEW**. Missing diagnosis coverage in the graph is a limited proposal in HANDOFF.md .\n"""
    (ROOT / 'README_DEMO.md').write_text(text)
    (ROOT / 'RESULTS.md').write_text(text)
    feedback = {'schema': 'dental-research-result-binding-v1', 'lane': 'X22-caries-synthetic', 'review_state': 'PENDING_INDEPENDENT_REVIEW', 'target_id': None, 'coverage_state': 'MISSING_SCOPED_WORKING_TARGET', 'proposed_target_id': 'DENT-DIAG-SYNTHETIC-CARIES-BENCHMARK', 'result_file': 'results/LANE_X22_CARIES_SYNTHETIC/results.json', 'sha256': result_hash, 'claim_type': 'capability', 'measured_quantity': '18144 conditionalsyntheticimages;7nullspacepairs;60same-scan multi-toothviews;clinicalcontrast UNKNOWN', 'units': 'imagecounts,opticaldepth,mm,conditional sensitivity', 'uncertainty': 'Unknown physicalDEJ, spectrum, PSF/scatter; exactdigitaltruth underclosures only', 'population_regime': 'Seven STS derivedmolaroutergeometries; synthetic tissueandlesions; source-volume splits', 'preregistered_gate': 'PREREG_R1 through R5 JSON; frozen hashes in results.json', 'baseline': 'ActualLogistic+ExtraTrees, directaxis-sum, analyticalslab, Poissonmoments andpublishedclinicalcomparison', 'outcome': 'CONDITIONAL_BENCHMARK_AND_AMBIGUITY_DELIVERED; EMPIRICAL_REALISM_UNKNOWN; REFERENCE_GENERALIZATION_FAILED', 'negative_result': True, 'dispatch': 'NOT_CREATED_NO_MATCHING_SCOPED_TARGET', 'reason': 'Directedworking-view/expansionsearchfoundonlysecondarycariesinterface, noappropriate diagnostic goal; missingpacketcall preserved'}
    write(ROOT / 'GRAPH_FEEDBACK.json', feedback)
    assert disk < 3000000000
    print('Demo report:', ROOT / 'README_DEMO.md', '| dataMiB', round(disk / 1024 ** 2, 1))
if __name__ == '__main__':
    report()
