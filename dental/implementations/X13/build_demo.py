from dental_release.paths import expand as _release_expand
from pathlib import Path
import json, csv, re, datetime, time, xml.etree.ElementTree as ET
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from predict_r1 import sha, dump
P = Path(__file__).resolve().parent

def licenses():
    out = []
    for source in json.loads((P / 'SOURCE_MANIFEST.json').read_text())['sources']:
        r = ET.parse(source['path']).getroot()
        p = ' '.join(' '.join(r.find('.//permissions').itertext()).split())
        if 'by-nc-nd/4.0' in p:
            lic = 'CC BY-NC-ND 4.0'
        elif 'Attribution-NonCommercial-ShareAlike 4.0' in p or 'by-nc-sa/4.0' in p:
            lic = 'CC BY-NC-SA 4.0'
        elif 'by-nc/4.0' in p:
            lic = 'CC BY-NC 4.0'
        elif 'by/3.0' in p:
            lic = 'CC BY 3.0'
        elif 'by/4.0' in p:
            lic = 'CC BY 4.0'
        else:
            lic = 'CC BY (version unspecified in local permissions)'
        out.append(dict(study=source['study'], license=lic, permissions_xpath='.//permissions', source_sha256=source['sha256'], doi=source['doi'], pmc_url='https://pmc.ncbi.nlm.nih.gov/articles/' + source['study'] + '/', corrects_frozen_metadata=source['license'] != lic))
    dump(P / 'LICENSES.json', out)
    return out

def figure(rows, r1, r2, r3, r4):
    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False})
    (fig, axes) = plt.subplots(2, 2, figsize=(13, 9))
    ax = axes[0, 0]
    regions = ['marginal', 'axial', 'occlusal']
    x = np.arange(3)
    width = 0.22
    ax.bar(x - width, [r2['metrics'][k]['R1_same_support_rmse_um'] for k in regions], width, label='R1: same held other arms')
    ax.bar(x, [r2['metrics'][k]['candidate_rmse_um'] for k in regions], width, label='R2: one local group')
    ax.bar(x + width, [r2['metrics'][k]['control_rmse_um']['constant_anchor_control_um'] for k in regions], width, label='R2 constant-anchor control')
    ax.axhline(20, color='firebrick', ls='--', label='Frozen 20 µm gate')
    ax.set_xticks(x, regions)
    ax.set_ylabel('Study-balanced group-mean RMSE (µm)')
    ax.set_title('A. Transfer and one-anchor calibration fail')
    ax.legend(fontsize=8)
    ax = axes[0, 1]
    colors = ['#0065A8', '#DD6500']
    for (study, col) in zip(['PMC10246932', 'PMC10721348'], colors):
        d = sorted([d for d in rows if d['study'] == study and d['region'] == 'marginal'], key=lambda d: d['internal_spacer_um'])
        m = next((m for m in json.loads((P / 'FROZEN_PREDICTIONS_R3.json').read_text())['models'] if m['study'] == study))
        s = np.linspace(d[0]['internal_spacer_um'], d[-1]['internal_spacer_um'], 150)
        mu = m['b_um'] + m['a_um_squared'] / s
        ax.plot(s, mu, color=col, label=study + ' local reciprocal')
        ax.scatter([q['internal_spacer_um'] for q in d[:2]], [q['measured_mean_um'] for q in d[:2]], color=col, marker='s', s=50)
        ax.scatter([q['internal_spacer_um'] for q in d[2:]], [q['measured_mean_um'] for q in d[2:]], edgecolor=col, facecolor='white', s=65)
        pred = [q for q in r3['cells'] if q['row_id'].startswith(study)]
        ax.scatter([q['spacer_um'] for q in pred], [q['prediction_um'] for q in pred], color=col, marker='x', s=65)
    ax.set_xlabel('Internal CAD setting (µm)')
    ax.set_ylabel('Dry vertical marginal group mean (µm)')
    ax.set_title('B. Two-dose local pilot: 2.40 µm RMSE')
    ax.legend(fontsize=8)
    ax.text(0.03, 0.03, 'Squares: anchors; circles: held measurement; x: prediction\nOnly 2 studies / 3 held cells; reciprocal control ties.', transform=ax.transAxes, fontsize=8)
    ax = axes[1, 0]
    n = 64
    lo = r4['low_um']
    hi = r4['high_um']
    series = np.where(np.indices((n, n))[1] < n // 2, lo, hi)
    parallel = np.where(np.indices((n, n))[0] < n // 2, lo, hi)
    im = ax.imshow(np.hstack([series, np.full((n, 4), np.nan), parallel]), cmap='viridis', origin='lower', vmin=lo, vmax=hi)
    ax.set_xticks([n / 2, n + n / 2 + 4], ['Series along flow', 'Parallel across flow'])
    ax.set_yticks([])
    ax.set_title('C. Constructed fields, identical histogram + mean')
    fig.colorbar(im, ax=ax, label='Gap (µm)', shrink=0.75)
    ax.set_xlabel('Flow from left to right in EACH field')
    ax = axes[1, 1]
    G = [r4['runs'][-2]['conductance_m3_Pa_s'], r4['uniform_mean_conductance_m3_Pa_s'], r4['runs'][-1]['conductance_m3_Pa_s']]
    ax.bar(['Series', 'Uniform mean', 'Parallel'], np.array(G) * 1e+16, color=['#0065A8', '#888888', '#DD6500'])
    ax.set_ylabel('Conductance (10⁻¹⁶ m³ / Pa s)')
    ax.set_title('D. Same mean, 4.11× different conductance')
    ax.text(0.03, 0.86, 'Existing BTE1 operator vs exact strip integration\nConstructed thin-film probe; no measured seating time.', transform=ax.transAxes, fontsize=8)
    fig.suptitle('CAD spacer ≠ observed gap: empirical calibration and the physical information limit', fontsize=14)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(P / 'CEMENT_GAP_DEMO.png', dpi=180)
    fig.savefig(P / 'CEMENT_GAP_DEMO.pdf')
    plt.close(fig)

def reports():
    rows = json.loads((P / 'measurements.json').read_text())
    r1 = json.loads((P / 'RESULTS_R1.json').read_text())
    r2 = json.loads((P / 'RESULTS_R2.json').read_text())
    r3 = json.loads((P / 'RESULTS_R3.json').read_text())
    r4 = json.loads((P / 'RESULTS_R4.json').read_text())
    dec = json.loads((P / 'INVERSE_DECISIONS.json').read_text())
    lic = licenses()
    figure(rows, r1, r2, r3, r4)
    with (P / 'FACIT.csv').open('w') as f:
        keys = ['row_id', 'study', 'source_family', 'arm', 'region', 'internal_spacer_um', 'marginal_spacer_um', 'measured_mean_um', 'reported_sd_um', 'n_specimens', 'state', 'method', 'doi', 'source_table', 'source_row', 'source_column', 'source_path', 'source_sha256']
        w = csv.DictWriter(f, fieldnames=keys, extrasaction='ignore')
        w.writeheader()
        w.writerows([d for d in rows if d['region'] in ['marginal', 'axial', 'occlusal']])
    result = dict(lane='X13-cement-gap', overall_outcome='REGIONAL_TRANSFER_FAIL_LOCAL_MARGINAL_PILOT_PASS_FULL_INVERSE_UNKNOWN', review_state='PENDING_INDEPENDENT_REVIEW', clinical_recommendation=False, measured_data={'rows': len(rows), 'primary_rows': 94, 'papers': 11, 'conservative_source_families': 10, 'selection': 'auditable purposive local subset; not a systematic review', 'file': 'measurements.csv', 'facit': 'FACIT.csv', 'missing_main_R04_source': 'PMC12237415 absent from prescribed local sources'}, rounds={'R1': {'outcome': r1['prediction_outcome'], 'rmse_um': {k: v['candidate_rmse_um'] for (k, v) in r1['metrics'].items()}, 'control': 'matched GLS ties; simpler controls often better', 'inverse': 'UNKNOWN'}, 'R2': {'outcome': r2['outcome'], 'rmse_um': {k: v['candidate_rmse_um'] for (k, v) in r2['metrics'].items()}, 'anchor_values': r2['anchor_values'], 'constant_anchor_control_rmse_um': {k: v['control_rmse_um']['constant_anchor_control_um'] for (k, v) in r2['metrics'].items()}, 'inverse': 'UNKNOWN'}, 'R3': {'pilot_outcome': r3['pilot_outcome'], 'rmse_um': r3['candidate_rmse_um'], 'max_error_um': r3['maximum_absolute_error_um'], 'affine_control_rmse_um': r3['control_rmse_um']['affine_control_um'], 'constant_control_rmse_um': r3['control_rmse_um']['constant_control_um'], 'strongest_control': 'reciprocal interpolation TIE', 'held_cells': 3, 'studies': 2, 'full_inverse': 'UNKNOWN'}, 'R4': {'outcome': r4['outcome'], 'conductance_ratio_same_mean': r4['aliasing_conductance_ratio'], 'numerical_max_relative_error': max((d['relative_error'] for d in r4['runs'])), 'constructed_fields': True, 'physical_validation': False}}, external_referent={'kind': 'independent_measurement', 'locator': str(P / 'SOURCE_MANIFEST.json') + ' (per-cell DOI or PMCID, local XML, table/row/column and SHA256 in FACIT.csv)', 'compared_quantity': 'regional reported group-mean gap, held whole source families (R1), other arms (R2), higher local doses (R3)', 'refutes_us': True}, external_referents=[r3['external_referent'], r4['external_referent']], inverse_answer={'regional_CAD_setting_um': None, 'status': 'UNKNOWN', 'reason': 'No admitted axial/occlusal within-system spacer dose studies; insufficient marginal independent dose studies and cross-study prediction failure', 'observed_joint_groups': len(dec['observed_joint']), 'observed_joint_groups_with_all_means_in_engineering_target': sum((d['all_group_means_in_target'] for d in dec['observed_joint']))}, physical_consumer={'existing_operator': r4['operator_path'], 'sha256': r4['operator_sha256'], 'port_status': 'UNKNOWN_MISSING_SPATIAL_FILM', 'seating_time_prediction': None, 'reason': 'Regional group means do not identify hydraulic connectivity and many observations already include seating'}, complete_cost={'preparation_parser_seconds': json.loads((P / 'SOURCE_MANIFEST.json').read_text())['seconds'], 'source_collection_search_review_cost': 'UNKNOWN, shared prefix not omitted', 'R1_first_fit_prediction_seconds': json.loads((P / 'FROZEN_PREDICTIONS_R1.json').read_text())['total_seconds'], 'R1_peak_rss_mib': json.loads((P / 'FROZEN_PREDICTIONS_R1.json').read_text())['peak_rss_mib'], 'matched_GLS_cost': 'same variance fitting and source preparation plus independently timed whitening; not whitening alone', 'R2_acquisition': '23 externally reported regional group means from 11 studies; physical lab cost UNKNOWN', 'R3_acquisition': '4 dose group means, n=10 per dose; same teeth reused in endocrown study; not independent patient counts', 'queries': 'INVERSE_DECISIONS.json timing', 'validation': 'all R1-R4 scorers, controls and corruption probes; RUN_ACCOUNTING.json for current full replay', 'reasoning_discovery_cost': 'UNKNOWN; includes failed R1/R2 constructs, not only successful R3', 'fallback_fabrication_metrology_cost': 'UNKNOWN, contract in MEASUREMENT_CONTRACT.md'}, license_manifest='LICENSES.json', license_metadata_correction='Frozen extraction preserved; PMC10557992 correctly CC BY-NC-SA 4.0 and PMC11010717 CC BY-NC-ND 4.0 in LICENSES.json, superseding ambiguous original license strings only.', data_arrays_over_50MB=[], large_data_location=_release_expand('@DENTAL_WORK_ROOT@/X13_cement_gap (unused; no large arrays written)'), graph_binding='BLOCKED_OLD_RESULT_HASH_MISMATCH; local pending feedback retained', limits=['Retrospective group-mean evaluation; no new physical experiment', 'R3 is only a two-study marginal extrapolation pilot and classical control ties', 'Gaussian transfer closures fail accuracy and coverage', 'Purposive local subset, no claim all local papers extracted', 'No spatial field, clinically validated setting or validated cement seating time', 'Reported specimen n and source lineage do not establish independence across sites/arms'])
    dump(P / 'results.json', result)
    freezeindex = dict(kind='index_of_prediction_files_frozen_before_their_scorers', prospective_new_measurement=False, records=[dict(file='FROZEN_PREDICTIONS_' + tag + '.json', sha256=sha(P / ('FROZEN_PREDICTIONS_' + tag + '.json')), created_utc=json.loads((P / ('FROZEN_PREDICTIONS_' + tag + '.json')).read_text())['created_utc'], prereg_sha256=sha(P / ('PREREG_' + tag + '.json'))) for tag in ['R1', 'R2', 'R3']])
    dump(P / 'FROZEN_PREDICTIONS.json', freezeindex)
    (P / 'FROZEN_PREDICTIONS.sha256').write_text(sha(P / 'FROZEN_PREDICTIONS.json') + '\n')
    table = '| Construction | Marginal | Axial | Occlusal | Outcome |\n|---|---:|---:|---:|---|\n'
    for (tag, r) in [('R1', r1), ('R2', r2)]:
        table += f"| {tag}: {('held out study' if tag == 'R1' else 'a local group')} | " + ' | '.join((f"{r['metrics'][k]['candidate_rmse_um']:.1f} µm" for k in ['marginal', 'axial', 'occlusal'])) + ' | FAIL |\n'
    table += f"| R3 : two local doses | {r3['candidate_rmse_um']:.2f} µm | UNKNOWN | UNKNOWN | Local pilot PASS |\n"
    liclines = '| Local article data | local fulltext license |\n|---|---|\n' + ''.join((f"| [{s['study']}]({s['pmc_url']}) | {s['license']} |\n" for s in lic))
    text = _release_expand(f"""# Cement gap as calibrated field\n\nThe CAD setting is a design value. The real gap after manufacturing and seating is a regional observation. This demo replaces the identity statement with published measurements and tests what a local calibration can predict.\n\nRun from this directory: `./run_all.sh` . All data is read local ; no dataset is retrieved. The script uses Python/NumPy/SciPy/Matplotlib in existing `@DENTAL_EXTERNAL_ROOT@/projects/cad-to-simulation/.venv-newton`, maximum four threads and no GPU. Figures are available in [CEMENT_GAP_DEMO.png ](CEMENT_GAP_DEMO.png) and [PDF ](CEMENT_GAP_DEMO.pdf).\n\n{table}\n\nwas frozen tolerance to 20 µm per region. R1 holds complete studies, and the two articles with possible common measurement line are held together. R2 receives a regional calibration group from the held out lab; the other groups are reference. The anchor is a reported agent over n specimen, not a measured individual crown. The simple constant anchor control is better than R2 : 41,9 / 16,5 / 180,3 µm . R1 : s interval coverage is also too low. No thresholds was lowered afterwards.\n\nR3 observes the two lowest local doses and predicates three higher doses in two studies. The maximum error is: {r3['maximum_absolute_error_um']:.1f} µm. Affin control provides {r3['control_rmse_um']['affine_control_um']:.2f} µm , constant second anchor {r3['control_rmse_um']['constant_control_um']:.2f} µm and common reciproc interpolation same {r3['candidate_rmse_um']:.2f} µm as the candidate. It is a small retrospective calibration pilot; not algorithm news or validated full regional inversion.\n\nReference is measured group agent, with n, method and DOI/PMCID /table/cell/hash in [FACIT.csv ](FACIT.csv).  Held out   R3 - measurements  is  15,94  and  13,10   µm  at  CAD   120 / 160  i [ PMC10246932 ](https://doi.org/10.7759/cureus.38688), and  22,90   µm  at  CAD   70  i [ PMC10721348 ](https://doi.org/10.1155/2023/6698453). Prediktionerna var 13,58/9,50/21,41 µm. The originals of the predictions have UTC time and hash and were saved before the respective scooters; the source tables were already read, so this is not a blind or prospective physical trial.\n\n## What does not holds\n\nIt is not yet possible to answer with a CAD transmissible setting that gives 50 – 100 µm axial and occlusal and maximum 120 µm margin. This is a declared engineering scenario, not a clinical advice. The base has two within-system dose studies marginally and zero axial/ occlusal ; the original gate requires at least three per region. Therefore, the API [inverse_query.py](inverse_query.py) returns `UNKNOWN`. The local margin curve can only leave an explicitly unertified query within its anchor range. {len(dec['observed_joint'])} observed groups have all three regions, and {sum((d['all_group_means_in_target'] for d in dec['observed_joint']))} all three has group funds within this particular goal. It describes the read selection and is not an impossibility for second designers.\n\nThe actual [BTE1 solver](../LANE_NEXT_E_CEMENT_SQUEEZE/squeeze.py) is called on two constructed fields with the same two thicknesses and mean. Flow Capacity Distinguishes {r4['aliasing_conductance_ratio']:.2f}× when narrow parts are in series and parallel respectively. The numerical solution matches exactly thin film integration. The fields are constructed, not measured crowns ; the result shows that regional means does not determine a hydraulic flow. [Hydrogen](DECOMPOSITION_R4.json) uses [Reynolds thin film type](https://doi.org/10.1098/rstl.1886.0005). No cement viscosity, curing curve or setting time is empirically calibrated here. [CEMENT_FIELD_PORT.json ](CEMENT_FIELD_PORT.json) retains this as `UNKNOWN_MISSING_SPATIAL_FILM` and avoids double counting of seating .\n\nTable-text conflict 23,22 against 22,22 µm in PMC10721348 is preserved (the table was used). R04 : s main reference PMC12237415 is not available in the prescribed local XML/LIT sources and was not used numerically. Pulplt, axio-occlusal, total internal and absolute margin deviation is kept separate from the primary regions. The selection is targeted, no systematic literature review. N and SD are not used as if all points or groups were independent.\n\n## Try a lab.\n\n[ MEASUREMENT_CONTRACT.md ](MEASUREMENT_CONTRACT.md) describes two calibration settings and a held out setting/batch, three regional observations and paraspatial column field. [ lab_selfcal.py ](lab_selfcal.py) freezes predictions from real calibration CSVs without reading future validation measurements. Such lab data is not available in this shipment.\n\n## Data and licences\n\n148 measurement cells from 11 local primary articles, 94 in primary regions. The LIT package is used for source matching; its duplicates are not counted as more evidence. No personal data or geometry data sets are used. The two previously incompletely classified license strings are fixed separately in [LICENSES.json ](LICENSES.json ; frozen numeric inputs are preserved.\n\n{liclines}\n\nThe whole experiment and all controls can be dropped by injected errors; [DELIVERY_CHECK.json ](DELIVERY_CHECK.json) shows the performed probes. It's reproductive and contract control, not independent scientific review. Status: **PENDING_INDEPENDENT_REVIEW**. Calculation and full cost accounting are available in [results.json ](results.json) and [RUN_ACCOUNTING.json ](RUN_ACCOUNTING.json); previous data collection, research work and physical lab cost is UNKNOWN . Disk size is small, without arrays over 50 MB. Graph dispatch is blocked by an old resulthash conflict; local result and feedback is left for the coordinator.\n""")
    (P / 'README_DEMO.md').write_text(text)
    (P / 'RESULTS.md').write_text(text.replace('# Cement gap as calibrated field', '# X13 — four designs against real cement gap', 1))
    return result
if __name__ == '__main__':
    reports()
