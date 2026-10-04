"""Reader-facing table, external comparisons, figure, frozen lab port, provenance."""
from dental_release.paths import expand as _release_expand
import csv
import datetime
import hashlib
import json
from pathlib import Path
import re
import time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from measure import HERE, ROOT, DATA, sha, dump
from capacity_map import load_local, body_points

def readlines(n):
    return [json.loads(l) for l in (HERE / n).read_text().splitlines()]

def resources(n):
    path = HERE / f'run_{n}_resources.log'
    if not path.exists():
        return {'status': 'unmeasured'}
    text = path.read_text()
    rss = re.search('Maximum resident set size \\(kbytes\\): (\\d+)', text)
    wall = re.search('Elapsed \\(wall clock\\) time \\(h:mm:ss or m:ss\\): (\\S+)', text)
    return {'max_rss_kbytes': int(rss.group(1)) if rss else None, 'wall_hms': wall.group(1) if wall else None}

def main():
    rows = readlines('COMBINED_SITES.jsonl')
    maps = readlines('RAW_R3_CAPACITY_MAPS.jsonl')
    edges = readlines('RAW_R4_LOCAL_EDGES.jsonl')
    sums = {f'R{n}': json.loads((HERE / f'SUMMARY_R{n}.json').read_text()) for n in (1, 2, 3, 4)}
    refs = json.loads((HERE / 'EXTERNAL_REFERENTS.json').read_text())
    valid = [r for r in rows if r['valid']]
    field_counts = {}
    field_evaluations = 0
    field_calls = 0
    for m in maps:
        for p in m['positions']:
            for c in p['candidates']:
                if 'margins_mm' not in c:
                    continue
                pair = (c['D_mm'], c['L_mm'])
                if pair not in field_counts:
                    (_, d) = body_points(np.zeros(3), np.array([1, 0, 0]), np.array([0, 1, 0]), np.array([0, 0, 1]), *pair)
                    field_counts[pair] = 2 * (int((d >= 2 - 1e-09).sum()) + 2 * len(d))
                field_evaluations += field_counts[pair]
                field_calls += 6
    r1 = [r for r in valid if r['construction'] == 'R1']
    passing = [r for r in r1 if r['cross_read_pass']]
    disagreement = lambda r: r['scalar_interval'] == 'height_augmentation_geometry' and r['rule5_flag'] == 'no_lift_trigger' or (r['scalar_interval'] == 'short_geometry' and r['rule5_flag'] == 'lift_trigger')
    external = {}
    for (key, stratum) in (('premolar', 'PM'), ('molar', 'M')):
        group = [r for r in valid if r['stratum'] == stratum]
        external[key] = {'n': len(group), 'our_height_lt8_fraction': sum((r['bone_height_mm'] < 8 for r in group)) / len(group) if group else None, 'our_width1_lt6_fraction': sum((r['widths_mm']['1'] < 6 for r in group)) / len(group) if group else None, 'reference': refs['anatomy'][key], 'gate': 'UNKNOWN_MINIMUM_N_AND_DIFFERENT_COHORT', 'definition_limit': 'Published buccopalatal crest width vs our1mm-subcrestal width; external cohort cannot validate local numerical/physical accuracy.'}
    n_pm = external['premolar']['n']
    n_m = external['molar']['n']
    counts = {k: sum((r.get('scalar_interval') == k for r in rows)) for k in ('short_geometry', 'height_augmentation_geometry', 'uncertain')}
    result = {'schema': 'X20-research-results-v1', 'lane': 'X20-short-implant-sinus', 'claim_type': ['information_link', 'capability'], 'outcome': 'RUNNABLE_ANATOMICAL_DIAGNOSTIC;SCALAR_SUFFICIENCY_REFUTED;COVERAGE_LOCATION_NUMERICAL_GATES_FAILED;PHYSICAL_AND_CLINICAL_UNKNOWN', 'external_referent': refs['anatomy'], 'external_referents': refs, 'rounds': sums, 'frequency_answer': {'isolated_gap_preliminary_disagreements': sum((disagreement(r) for r in r1)), 'isolated_gap_n': len(r1), 'both_disagreements_fail_frozen_numeric_reread': all((not r['cross_read_pass'] for r in r1 if disagreement(r))), 'numerically_accepted_isolated_gap_n': len(passing), 'numerically_accepted_disagreements': sum((disagreement(r) for r in passing)), 'all_queries_preliminary_disagreements': sum((disagreement(r) for r in valid)), 'all_queries_evaluable_n': len(valid), 'extrapolation': 'NONE; mirror location failed and convenience sample is not a population. 5mm baseline deliberately simplified.'}, 'combined_query_summary': {'census_n': 480, 'sinus_labeled_n': 48, 'n_selected_unique': len(rows), 'n_evaluable_labels': len(valid), 'n_excluded': len(rows) - len(valid), 'scalar_interval_counts': counts, 'clinical_unknown_n': len(rows), 'full_body_origin_no_catalog_fit_n': sums['R3']['origin_counts']['no_catalog_fit_geometry'], 'full_body_origin_unknown_n': sums['R3']['origin_counts']['unknown'], 'no_clinical_augmentation_proof': True}, 'independent_anatomy_comparison': external, 'validation': json.loads((HERE / 'VERIFICATION.json').read_text()), 'full_cost': {'preparation': 'Source orientation/literature/preregistration, not inferred from numerical runtimes; interactive agent time not instrumented', 'fit': 'No trained predictive model; mirror per-target noniterative geometry fit', 'discovery': {'R1_sites': 17, 'R2_queries': 24, 'R2_heldout_targets': 296, 'R3_candidate_queries': sums['R3']['candidate_queries'], 'R3_sample_points': sums['R3']['point_queries'], 'R4_profiles': 10}, 'validation': 'Actual voxel reads, heldout label pose, field interpolation, image gradient control,9corruptions, source hashes; external physical validation unavailable', 'R3_actual_field_value_evaluations_including_control': field_evaluations, 'R3_field_lookup_calls_including_control': field_calls, 'questions': 0, 'fallback': 'All clinical choicesUNKNOWN; numerical/selection exclusions retained', 'threads': 1, 'gpu': 0, 'resources': {k: resources(k) for k in ('r1', 'r2', 'r3', 'r4', 'measure', 'mirror', 'capacity_map', 'local_edges')}, 'data_bytes': sum((p.stat().st_size for p in DATA.glob('*') if p.is_file())), 'intermediate_limit_bytes': 3000000000, 'agent_token_cost': 'UNKNOWN'}, 'datasets': [{'name': 'ToothFairy2', 'license': 'CC BY-SA4.0', 'locator': 'https://ditto.ing.unimore.it/toothfairy2/', 'local_archive': str(__import__('tf2_io').ZIP)}], 'honest_unknowns': ['Physical local edge maximum', 'Bony sinus floor vs sinus air/mucosa', 'Edentulous prosthetic pose', 'Bone quality and primary stability', 'Loading and survival', 'Independent exchangeable cohort prevalence'], 'next_construction': 'Fiducial-aligned local edentulous crest/sinus surface measurement and exact cubical voxel distance query replacing EDT surrogate; preserve current failures.'}
    for optional in ('INTERFACE_VERIFICATION.json', 'PHYSICAL_PORT_EMPTY_TEST.json'):
        if (HERE / optional).exists():
            result[optional] = json.loads((HERE / optional).read_text())
    manifest = {'files': []}
    for n in ('DATA_MANIFEST.json', 'DATA_MANIFEST_R2.json', 'DATA_MANIFEST_R4.json'):
        result[n] = json.loads((HERE / n).read_text())
    for n in ('cells/geometry/tf2_io.py', 'results/BONE_HU/tf2_census.jsonl', 'results/DESIGN_implant/catalog.json', 'results/GEOM_UNC/field_model.json'):
        result.setdefault('source_files', []).append({'path': str(ROOT / n), 'sha256': sha(ROOT / n)})
    dump(HERE / 'results.json', result)
    table = []
    with open(HERE / 'SITE_DECISIONS.csv', 'w') as f:
        cols = ['case', 'fdi', 'construction', 'height_mm', 'height_lower_scenario_mm', 'height_upper_scenario_mm', 'width1_mm', 'width_lower_scenario_mm', 'width_upper_scenario_mm', 'rule5', 'scalar_screen', 'numeric_reread_pass', 'full_body_origin', 'shortest_finite_body_catalog_length_mm', 'clinical_choice']
        w = csv.DictWriter(f, cols)
        w.writeheader()
        for r in rows:
            m = next((m for m in maps if (m['case'], m['fdi']) == (r['case'], r['fdi'])))
            h = r.get('height_interval_mm') or [None, None]
            v = r.get('width1_interval_mm') or [None, None]
            z = dict(case=r['case'], fdi=r['fdi'], construction=r['construction'], height_mm=r.get('bone_height_mm'), height_lower_scenario_mm=h[0], height_upper_scenario_mm=h[1], width1_mm=r.get('widths_mm', {}).get('1'), width_lower_scenario_mm=v[0], width_upper_scenario_mm=v[1], rule5=r.get('rule5_flag', 'unknown'), scalar_screen=r.get('scalar_interval', 'uncertain'), numeric_reread_pass=r.get('cross_read_pass'), full_body_origin=m.get('origin_status', 'unknown'), shortest_finite_body_catalog_length_mm=None, clinical_choice='UNKNOWN')
            w.writerow(z)
            if r['valid']:
                table.append(f"| {r['case'].replace('ToothFairy2', '')} / {r['fdi']} | {r['construction']} | {r['bone_height_mm']:.2f} ±0.60 | {r['widths_mm']['1']:.2f} ±0.60 | {r['scalar_interval']} | {('PASS' if r['cross_read_pass'] else 'FAIL')} |")
    plt.rcParams.update({'font.size': 9})
    (fig, axes) = plt.subplots(2, 2, figsize=(13, 10), constrained_layout=True)
    ax = axes[0, 0]
    for r in valid:
        color = '#136a70' if r['construction'] == 'R1' else '#b27318'
        ax.errorbar(r['bone_height_mm'], r['widths_mm']['1'], xerr=0.6, yerr=0.6, fmt='o', color=color, alpha=0.8)
        ax.annotate(f"{r['case'].split('_')[-1]}/{r['fdi']}", (r['bone_height_mm'], r['widths_mm']['1']), fontsize=7, xytext=(3, 3), textcoords='offset points')
    ax.axvline(5, color='#a44', ls='--', label='Simplified5mm trigger')
    ax.axvline(7, color='k', ls=':', label='6mm body +1mm margin')
    ax.axhline(6, color='k', ls=':')
    ax.set(xlabel='Contiguous labeled bone height (mm)', ylabel='Width at1mm below crest (mm)', title='10 exploratory sites; ±0.6mm sensitivity intervals')
    ax.legend(loc='upper left', fontsize=8)
    ax = axes[0, 1]
    ref = refs['anatomy']
    groups = ['PM height<8', 'M height<8', 'PM width<6', 'M width<6']
    observed = [external['premolar']['our_height_lt8_fraction'], external['molar']['our_height_lt8_fraction'], external['premolar']['our_width1_lt6_fraction'], external['molar']['our_width1_lt6_fraction']]
    expected = [ref['premolar']['height_lt8'], ref['molar']['height_lt8'], ref['premolar']['width_lt6'], ref['molar']['width_lt6']]
    x = np.arange(4)
    ax.bar(x - 0.16, observed, 0.32, label=f'TF2 exploratory n={n_pm}PM/{n_m}M')
    ax.bar(x + 0.16, expected, 0.32, label='Independent study349sites')
    ax.set_xticks(x, groups, rotation=20)
    ax.set(ylim=(0, 1), ylabel='Fraction', title='External cohort comparison: gate UNKNOWN (small/different sample)')
    ax.legend(fontsize=8)
    r = next((r for r in valid if r['case'] == 'ToothFairy2F_004' and r['fdi'] == 17))
    (local, anchor, axis, normal, tangent) = load_local(r)
    p = np.asarray(r['crest_mm'])
    u = np.linspace(-7, 7, 141)
    v = np.linspace(-3, 17, 201)
    (U, V) = np.meshgrid(u, v)
    pts = p + U[..., None] * normal + V[..., None] * axis
    bone = local.sample('bone', pts)
    sin = local.sample('sinus', pts)
    image = np.load(DATA / f"{r['case']}_{r['fdi']}_image_crop.npz")
    q = (pts / image['spacing'] - image['lo']).reshape(-1, 3).T
    from scipy import ndimage as ndi
    gray = ndi.map_coordinates(image['image'], q, order=1, mode='constant', cval=0).reshape(U.shape)
    ax = axes[1, 0]
    ax.imshow(gray, cmap='gray', origin='lower', extent=(-7, 7, -3, 17), vmin=np.percentile(gray, 10), vmax=np.percentile(gray, 98), aspect='equal')
    ax.contour(U, V, bone, levels=[0.5], colors=['#35e675'], linewidths=0.8)
    ax.contour(U, V, sin, levels=[0.5], colors=['#31c9fa'], linewidths=0.8)
    ax.plot([-2, 2, 2, -2, -2], [0, 0, 6, 6, 0], color='#edb331', lw=2)
    ax.set(xlabel='Buccopalatal displacement (mm)', ylabel='Depth along query axis (mm)', title=_release_expand('Actual CBCT @DENTAL_CASE_ID@/17;4×6mm test body, no certified fit'))
    ax = axes[1, 1]
    for e in edges:
        if e['valid']:
            z = (np.asarray(e['intensity']) - e['outside_gray']) / (e['inside_gray'] - e['outside_gray'])
            ax.plot(e['profile_t_mm'], z, label=f"{e['case'].split('_')[-1]}/{e['fdi']} δ={e['image_minus_label_offset_mm']:.2f}mm")
    ax.axvspan(-0.3, 0.3, color='#bbb', alpha=0.3)
    ax.axhline(0.5, color='k', ls=':')
    ax.axvline(0, color='k', ls='--')
    ax.set(xlim=(-1.5, 1.5), ylim=(-0.5, 1.8), xlabel='Position relative to label boundary (mm; positive into bone)', ylabel='Normalized raw image profile', title='6 QC-valid local edges;1 exceeds assumed0.3mm box')
    ax.legend(fontsize=7)
    fig.suptitle('Residual maxillary anatomy: label measurements and failed validation gates', fontsize=14)
    fig.savefig(HERE / 'ANATOMY_DECISIONS.png', dpi=170)
    fig.savefig(HERE / 'ANATOMY_DECISIONS.pdf')
    plt.close(fig)
    lines = '\n'.join(table)
    text = f'# Actual posterior maxillary anatomy — research demo\n\nRead actual bone and sinus geometry at missing-tooth queries, then test whether a full short implant body fits. The code now exposes failures hidden by a single residual-height threshold. It delivers measurements and diagnostic maps; physical and clinical decisions remain UNKNOWN.\n\nRun from this directory:\n\n```sh\n./run_all.sh\n```\n\nLocal dependencies: `/usr/bin/python3` with system NumPy, SciPy and Matplotlib; existing TF2 zip and NV1 landmarks. No network access, training, GPU, or whole-archive extraction during the demo. About3minutes per full replay, one thread, observed peak0.44GB,46MB data. The command retains four frozen protocols and their negative gates.\n\n| Quantity | Result |\n|---|---|\n| Census / sinus-labeled scans | 480 /48; independently checked against NV1 |\n| Isolated real gap candidates / label-measurable | 17 /6; coverage gate FAIL (<8) |\n| Preliminary conflict with simplified5mm trigger | 2/6; **both fail numerical reread tolerance** |\n| Numerically accepted isolated measurements / accepted conflicts | 3 /0 |\n| Mirror heldout targets / P95 transverse error | 296 /3.75mm; location gate FAIL (>2mm) |\n| Additional mirror queries / measurable | 11 /4; coverage gate FAIL (<8new) |\n| Combined exploratory queries / label-measurable | 28 /10; 18 exclusions |\n| Scalar interval screens | 3short /4height-deficit /21uncertain; not clinical choices |\n| Full-body origin / feasible region | 0short certificates; 7model no-fit /21unknown; map gate FAIL |\n| Scalar nominal positives rejected by whole body | 5; scalar sufficiency REFUTED |\n| Nearest-grid vs trilinear candidate control | 457/6963 exceed frozen0.3696153mm tolerance; FAIL |\n| Local raw-image edges | 6/10QC-valid;5/6within0.3mm; minimum8 gateUNKNOWN |\n| Corruption controls | 9/9reject injected errors |\n| Clinical short-versus-lift choices | UNKNOWN at all28queries |\n\n![Figure](ANATOMY_DECISIONS.png)\n\n| Query | Locator construction | Bone height mm (scenario) | Width1mm (scenario) | Scalar interval screen | Numerical reread |\n|---|---|---:|---:|---|---|\n{lines}\n\n`SITE_DECISIONS.csv` includes all28queries and exclusions. ±0.60mm intervals are sensitivity boxes from±0.30mm per edge, **not calibrated confidence intervals**. A shortest6mm implant is an actual K37 catalog pair only at compatible diameters; scalar passes do not certify its complete volume. R3 uses the actual paired catalog, not a Cartesian product. No full-body candidate passes the frozen conservative budget, so the final shortest certified length is UNKNOWN.\n\nThe simplified baseline uses <5mm as a lift trigger and >=5mm as no trigger. [Pal2012](https://doi.org/10.4103/0975-5950.102148) describes Misch height classes for choosing augmentation procedures; it does **not** establish a full binary clinical short-versus-lift rule. Our rule discrepancy cannot be presented as contemporary clinicians being wrong. R2 adds a third preliminary conflict (3/10exploratory), but mirror location failed and that number has no population interpretation.\n\nThe independent [CBCT cohort](https://doi.org/10.1007/s12663-019-01236-7) has349edentulous sites. Height<8mm occurred in44.86%PM and67.83%M; width<6mm in54.42%PM and55.45%M. Our{n_pm}PM/{n_m}M sample is below the frozen20-per-stratum minimum, and our1mm-subcrestal width is not identical to its crest measurement: external comparison UNKNOWN. Cohort proportions cannot validate a local physical edge.\n\nThe [5-year randomized trial](https://doi.org/10.1111/jcpe.13025) enrolled101patients at5–7mm residual height: 6mm short versus11–15mm grafted implants, reported patient-level survival98.5%/100% and mean marginal bone levels0.54/0.46mm. The [10-year followup](https://doi.org/10.1111/jcpe.13954) of the same cohort reported96%/100% patient-level survival. These results permit a research comparison and refute automatic extrapolation from our restrictive1mm-clearance geometry to clinical necessity. They provide no TF2 patient-outcome truth, and non-significance is not proven equivalence.\n\nWhat does not hold: reliable reflected edentulous pose; calibrated physical local error; claim that a central height and one width ensure full-body clearance; a validated population frequency; a clinical augmentation certificate. Neighboring roots and local ridge shape can reject an apparently adequate scalar height. Absence of a catalog fit under our frozen proxy does not establish that sinus lifting is the right procedure. The CBCT half-level edge is an observation, not histology or physical surface truth.\n\n`MINIMUM_MEASUREMENT_PORT.md`, `FROZEN_LAB_PREDICTIONS.json` and `compare_measurements.py` let a lab preregister paired local measurements before supplying any values. [Veyre-Goulet2008](https://doi.org/10.1111/j.1708-8208.2008.00083.x) provides an independent precedent with fiducial-aligned sectioning and caliper; it cannot calibrate our local image-label offsets. No physical measurement has been performed here.\n\nData license: ToothFairy2 CC BY-SA4.0, from its local `dataset.json`, source https://ditto.ing.unimore.it/toothfairy2/. Derived crops/figures carry the same dataset attribution and share-alike requirement. FDA GUDID catalog is inherited from K37; local catalog and source code hashes are in `results.json`. No other anatomical dataset, patient records or personal identifiers are used.\n'
    (HERE / 'README_DEMO.md').write_text(text)
    (HERE / 'RESULTS.md').write_text(text + '\nDetailed frozen gates and exact raw counts: PREREG_R1–R4, SUMMARY_R1–R4, RAW_R1/R2/R3/R4, VERIFICATION. All producer findings PENDING_INDEPENDENT_REVIEW.\n')
    labpath = HERE / 'FROZEN_LAB_PREDICTIONS.json'
    lab = {'frozen_at_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'claim_type': 'capability', 'source_measurement_sha256': sha(HERE / 'COMBINED_SITES.jsonl'), 'physical_data_observed': False, 'predictions': [{'case': r['case'], 'fdi': r['fdi'], 'crest_coordinate_mm_zyx': r['crest_mm'], 'axis_zyx': r['axis'], 'height_mm': r['bone_height_mm'], 'height_interval_mm': r['height_interval_mm'], 'width1_mm': r['widths_mm']['1'], 'width1_interval_mm': r['width1_interval_mm'], 'coordinate_validity': 'QUERY_POSITION; independent fiducial matching required'} for r in valid], 'external_referent': {'kind': 'independent_measurement', 'locator': 'https://doi.org/10.1111/j.1708-8208.2008.00083.x', 'compared_quantity': 'Registered posterior maxillary height/width CBCT versus physical section/caliper', 'refutes_us': True}, 'falsifier': 'One independently registered physical height or width outside frozen±0.60mm interval rejects universal box assertion; no after-fit repair', 'physical_port_status': 'UNKNOWN_NO_PHYSICAL_VALUES'}
    if not labpath.exists():
        dump(labpath, lab)
        (HERE / 'FROZEN_LAB_PREDICTIONS.sha256').write_text(sha(labpath) + '\n')
    else:
        assert sha(labpath) == (HERE / 'FROZEN_LAB_PREDICTIONS.sha256').read_text().strip()
    attempts = []
    for i in (1, 2, 3, 4):
        pr = json.loads((HERE / f'PREREG_R{i}.json').read_text())
        attempts.append({'round': f'R{i}', 'claim_type': pr['claim_type'], 'changed_operation': pr['changed_operation'], 'source_prereg_sha256': sha(HERE / f'PREREG_R{i}.json'), 'result': sums[f'R{i}'], 'negative_result': True})
    dump(HERE / 'ATTEMPTS.json', attempts)
    handoff = _release_expand('# X20 handoff\n\nRun `./run_all.sh`; reader-facing idea, table and figure are README_DEMO.md and ANATOMY_DECISIONS.png. Four frozen constructions executed, all negative/UNKNOWN gates preserved. R1 actual isolated gaps:17selected/6measurable; preliminary2/6height-rule disagreements both fail0.15mm voxel reread. No numerically accepted conflict (0/3). R2 mirrored pose:296held-out targets, P95transverse3.751284mm vs2mm FAIL, medianaxis6.692922deg PASS;4new measurable sites, coverageFAIL. R3 full solid:10maps/81positions;6963candidates/189534126samplepoints,0feasible regions;5scalar nominal origins rejected;457lookup discrepancies exceed0.3696153mm. R4 rawlocalimage:6/10QC,5/6within0.3mm, oneoffset−0.442214mm, minimum8UNKNOWN; physicaltruthabsent.\n\nCombined28uniquequeries/10label-measurable;3scalar short/4height deficit/21uncertain; final clinical choices28UNKNOWN. Only48/480census scans have sinus labels≥500voxels; independently matched NV1 census. PopulationfrequencyUNKNOWN. External CBCT source DOI10.1007/s12663-019-01236-7; RCTDOIs10.1111/jcpe.13025 and10.1111/jcpe.13954 (samecohortfollowup), physicalmeasurementprecedent DOI10.1111/j.1708-8208.2008.00083.x. External data are not paired to our local surfaces. No10×orclinicalsuccessclaim.\n\nNext construction: exact cubical complement distance instead of EDT-halfvoxel surrogate, plus an independently marked prosthetic query and bony sinus floor in3registered specimens. Fit/calibrate onlyafter frozen comparison; first run `compare_measurements.py --csv ...` against FROZEN_LAB_PREDICTIONS without altering its hash. MINIMUM_MEASUREMENT_PORT names reference methods, coordinates, fields and fail criteria. A printed labelphantom would be our_own_fixture and cannot supply realanatomicaltruth.\n\nGraph experimentdispatch denied numerical_dispatch_not_allowed/unresolved_prerequisites. A scoped definition/prototype dispatch to DENT-GEOM-NERVE-VESSEL-CANALS succeeded; PENDING_INDEPENDENT_REVIEW. Missingcoverage proposal: edentulous_query_pose + physically observed_bony_sinus_floor + localdirectionalerror port feeding IM14/K37/K42. No native/generated graph modifications. Sourceheight/rule/uncertainty conventions explicit.\n\nKnown limitations retained:1voxel ROIguard is insufficient for a full physicalmargincertificate near truncation; future exactcube implementation must require full needed exteriorcoverage. No positive fullbodycertificate exists in thisround. Frozen controltolfailed, never relaxed. Labels/imagehalfedges and unvalidatedmirror do not justify clinicalnecessity. Actual baseline is simplified5mmtrigger, not complete clinicalpractice.\n\nData46MB in @DENTAL_WORK_ROOT@/X20_short_implant_sinus, hashes in results.json; one thread/GPU0, peak0.44GB. Reproducibility, corruptions and figure checks are executable; final coordinator review remains external. No persondata, mail, publishing, push, subagents or othercheckouts touched.\n')
    (HERE / 'HANDOFF.md').write_text(handoff)
    (HERE / 'HANDOFF_R4.md').write_text(handoff)
    for p in sorted(HERE.glob('*')):
        if p.is_file() and p.suffix in ('.py', '.json', '.jsonl', '.csv', '.png', '.pdf', '.md', '.sh', '.sha256') and (p.name not in ('ARTIFACT_MANIFEST.json', 'CURRENT_WORK_STATE.json', 'GRAPH_FEEDBACK.json', 'GRAPH_FEEDBACK_RECEIPT.json')):
            manifest['files'].append({'path': str(p), 'sha256': sha(p), 'bytes': p.stat().st_size})
    dump(HERE / 'ARTIFACT_MANIFEST.json', manifest)
    print('Reader artifacts and frozen physical measurement port written.')
if __name__ == '__main__':
    main()
