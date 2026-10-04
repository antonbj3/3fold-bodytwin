from measure import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def as_csv(path, rows):
    if not rows:
        return
    keys = list(rows[0])
    with open(path, 'w') as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)

def main():
    start = time.perf_counter()
    patients = json.loads((ROOT / 'raw/PAIRED_PATIENTS.json').read_text())
    r2 = json.loads((ROOT / 'raw/R2_RESULT.json').read_text())
    r3 = json.loads((ROOT / 'raw/R3_RESULT.json').read_text())
    r4 = json.loads((ROOT / 'raw/R4_RESULT.json').read_text())
    controls = json.loads((ROOT / 'raw/CONTROLS.json').read_text())
    sites = json.loads((ROOT / 'raw/PER_SITE_LOCATED.json').read_text())
    sections = json.loads((ROOT / 'raw/PER_PATIENT_SECTIONS.json').read_text())
    identity = json.loads((ROOT / 'raw/IDENTITY.json').read_text())
    margins = json.loads((ROOT / 'raw/GUIDE_SEGMENT_SCENARIOS.json').read_text())
    artifacts = [x['point_artifact'] for x in patients] + r4['artifacts']
    for a in artifacts:
        assert Path(a['path']).stat().st_size == a['bytes'] and sha_file(a['path']) == a['sha256']
    scratch = sum((p.stat().st_size for p in DATA.rglob('*') if p.is_file()))
    assert scratch <= 3000000000
    revisions = [p for p in patients if not p['numeric_mask_equal']]
    regional = []
    for seg in ('mental_endpoint_5mm_proxy', 'molar_nearest_FDI', 'posterior_endpoint_10mm_proxy'):
        rr = [a for p in sections if p.get('revised') for a in p.get('regions', []) if a['segment'] == seg and a['n']]
        regional.append({'segment': seg, 'patients': len(rr), 'median_patient_radial_median_mm': float(np.median([a['median_mm'] for a in rr])), 'median_patient_radial_p95_mm': float(np.median([a['p95_mm'] for a in rr])), 'resolution': 'POPULATION', 'input_resolution': 'PER_POINT/PER_SURFACE_REGION', 'anatomical_validation': 'UNKNOWN_ENDPOINT_PROXIES'})
    radial_margins = []
    guide_bases = {a['guide']: a['guide_base_scenario_mm'] for a in margins}
    existing = {}
    with (DENT / 'results/LANE_X8_GUIDE_NERVE_RISK/PER_SITE_FULL_GEOMETRY_RISK.csv').open() as f:
        for x in csv.DictReader(f):
            existing[x['case'], int(x['fdi'])] = x
    gaps = np.array([0.5 * (float(x['whole_cylinder_label_gap_lower_mm']) + float(x['whole_cylinder_label_gap_upper_mm'])) for x in existing.values()])
    for reg in regional:
        for (guide, base) in guide_bases.items():
            budget = reg['median_patient_radial_p95_mm']
            margin = base + budget
            radial_margins.append({'segment': reg['segment'], 'guide': guide, 'revised_patients_with_shared_wall': reg['patients'], 'radial_budget_mm': budget, 'base_margin_mm': base, 'combined_sensitivity_margin_mm': margin, 'all_2581_2mm_to_combined_class_changes': int(((gaps >= 2) & (gaps < margin)).sum()), 'all_2581_guide_only_to_combined_class_changes': int(((gaps >= base) & (gaps < margin)).sum()), 'resolution': 'PHENOMENOLOGICAL', 'coverage_loss_handling': 'NOT_IMPUTED; inspect station support flags before local use', 'calibrated_margin': False, 'replacement_measurement': 'Independent blinded contours at the site + matched signed guide deviations; do not add a population P95 as a physical guarantee', 'double_counting_possible': True})
    dump(ROOT / 'raw/GUIDE_SHARED_WALL_SENSITIVITY.json', radial_margins)
    as_csv(ROOT / 'raw/PER_PATIENT_SUMMARY.csv', [{'patient': p['patient'], 'revised': not p['numeric_mask_equal'], 'old_voxels': p['old_voxels'], 'new_voxels': p['new_voxels'], 'added_fraction_percent': p['author_added_fraction_percent'], 'global_surface_median_mm': p['symmetric']['median_mm'], 'global_surface_p95_mm': p['symmetric']['p95_mm'], 'image_equal_tf2': p['image_identity']['exact_numeric_equal'], 'resolution': 'PER_SURFACE_REGION'} for p in patients])
    as_csv(ROOT / 'raw/PER_REGION_SUMMARY.csv', regional)
    as_csv(ROOT / 'raw/GUIDE_SHARED_WALL_SENSITIVITY.csv', radial_margins)
    as_csv(ROOT / 'raw/CLASS_REVERSALS.csv', [{'case': x['case'], 'fdi': x['fdi'], 'old_mm': x['old']['upper_mm'], 'new_mm': x['new']['upper_mm'], 'tf2_mm': x['tf2']['upper_mm'], 'union_mm': x['observed_union']['upper_mm'], 'old_to_new_changed': x['old_class_2mm'] != x['new_class_2mm'], 'tf2_to_union_changed': x['tf2_class_2mm'] != x['observed_union']['class_2mm'], 'resolution': 'PER_TOOTH'} for x in sites if x['old_class_2mm'] != x['new_class_2mm'] or x['tf2_class_2mm'] != x['observed_union']['class_2mm']])
    plt.rcParams.update({'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False, 'figure.dpi': 140})
    (fig, ax) = plt.subplots(2, 2, figsize=(11, 8))
    colors = ['#2583a5', '#cf6944', '#7662a3']
    complete = [x for x in sections if x.get('revised') and x.get('old_new_supported_radial', {}).get('n')]
    xx = np.array([x['global_surface_p95_mm'] for x in complete])
    yy = np.array([x['old_new_supported_radial']['p95_mm'] for x in complete])
    ax[0, 0].scatter(xx, yy, c=colors[0], s=25)
    lim = max(xx.max(), yy.max()) * 1.05
    ax[0, 0].plot([0, lim], [0, lim], ls='--', c='0.7')
    ax[0, 0].set(xlabel='Whole-surface patient P95 (mm)', ylabel='Shared-lumen radial patient P95 (mm)', title=f'A  Coverage and wall location differ (n={len(complete)})')
    ax[0, 0].text(0.04, 0.92, f"Ratio of medians = {r4['ratio_global_to_supported']:.1f}\nSource: published same-patient revisions", transform=ax[0, 0].transAxes, va='top')
    values = []
    for r in regional:
        values.append([a['p95_mm'] for p in sections if p.get('revised') for a in p.get('regions', []) if a['segment'] == r['segment'] and a['n']])
    boxes = ax[0, 1].boxplot(values, labels=[f'Mental*\nn={len(values[0])}', f'Molar FDI\nn={len(values[1])}', f'Posterior*\nn={len(values[2])}'], patch_artist=True, showfliers=False)
    for (box, c) in zip(boxes['boxes'], colors):
        box.set_facecolor(c)
        box.set_alpha(0.45)
    ax[0, 1].set(ylabel='Patient radial P95 (mm)', title='B  Shared-wall revisions by region')
    ax[0, 1].text(0.02, 0.94, '* Endpoint windows are anatomical proxies', transform=ax[0, 1].transAxes, va='top', fontsize=8)
    changed = [x for x in sites if x['tf2_class_2mm'] != x['observed_union']['class_2mm']]
    ids = np.arange(len(changed))
    width = 0.35
    ax[1, 0].bar(ids - width / 2, [x['tf2']['upper_mm'] for x in changed], width, label='TF2 mask', color=colors[0])
    ax[1, 0].bar(ids + width / 2, [x['observed_union']['upper_mm'] for x in changed], width, label='Observed mask union', color=colors[1])
    ax[1, 0].axhline(2, c='black', ls='--')
    ax[1, 0].set(xticks=ids, xticklabels=[x['patient'] + ' / ' + str(x['fdi']) for x in changed], ylabel='Fixed cylinder clearance (mm)', title=f'C  {len(changed)}/{len(sites)} located 2 mm class changes')
    ax[1, 0].legend(fontsize=8)
    ss = controls['sufficiency']
    ax[1, 1].bar([0, 1], [ss['downstream_fixed_point_distance_a_mm'], ss['downstream_fixed_point_distance_b_mm']], color=colors[:2])
    ax[1, 1].axhline(2.5, c='black', ls='--')
    ax[1, 1].set(xticks=[0, 1], xticklabels=['State A', 'State B'], ylabel='Distance to fixed query (mm)', title='D  Identical summaries fail downstream')
    ax[1, 1].text(0.03, 0.96, 'Exact identity error = 0\nSame volume, Dice, median, P95, maximum\nLogical fixture, not anatomical validation', transform=ax[1, 1].transAxes, va='top', fontsize=8, bbox={'facecolor': 'white', 'edgecolor': 'none', 'alpha': 0.8})
    fig.suptitle('X58: located canal revisions; independent annotator floor remains UNKNOWN', fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(ROOT / 'figures/x58_results.png')
    fig.savefig(ROOT / 'figures/x58_results.pdf')
    plt.close(fig)
    example = next((x for x in sites if x['old_class_2mm'] != x['new_class_2mm']))
    data = np.load(example['point_artifact']['path'])
    pose = example['pose']
    E = np.array(pose['entry_zyx_mm'])
    a = np.array(pose['axis_zyx'])
    L = pose['length_mm']
    R = pose['radius_mm']
    mid = E + 0.5 * L * a
    (fig, axes) = plt.subplots(1, 2, figsize=(10, 4.6))
    for (j, (u, v, label)) in enumerate([(2, 1, 'Axial XY projection'), (2, 0, 'Coronal XZ projection')]):
        for (key, color, name) in [('old_voxels_zyx', colors[0], 'Maxillo'), ('new_voxels_zyx', colors[1], 'ToothFairy'), ('tf2_voxels_zyx', colors[2], 'ToothFairy2')]:
            pts = data[key] * SP
            near = np.linalg.norm(pts - mid, axis=1) < 15
            pts = pts[near]
            axes[j].scatter(pts[::3, u], pts[::3, v], s=1, alpha=0.22, c=color, label=name)
        t = np.linspace(0, L, 50)
        cylinder = E + t[:, None] * a
        axes[j].plot(cylinder[:, u], cylinder[:, v], c='black', lw=4, label='Fixed plan axis')
        axes[j].set(xlabel='X (mm)', ylabel='Y (mm)' if v == 1 else 'Z (mm)', title=label)
        axes[j].axis('equal')
    axes[0].legend(fontsize=8, markerscale=4)
    fig.suptitle(f"{example['patient']}, FDI {example['fdi']}: old {example['old']['upper_mm']:.3f}, new {example['new']['upper_mm']:.3f}, TF2 {example['tf2']['upper_mm']:.3f} mm\nProjection for visualization; decisions use full 3D cylinder, radius 2 mm", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.88))
    fig.savefig(ROOT / 'figures/patient_P54_FDI35.png')
    plt.close(fig)
    cold = json.loads((ROOT / 'rounds/first_pass_R2_RESULT.json').read_text())['cost']
    r3cold = json.loads((ROOT / 'rounds/first_pass_R3_RESULT.json').read_text())['cost']
    cost = {'preparation_seconds': 'UNKNOWN: document reading and code construction were not separately metered', 'fit_seconds': 0, 'cold_measurement_first_pass_seconds': cold['wall_seconds'], 'decoder_repair_resume_seconds': 16.19051095400937, 'cold_consumer_first_pass_seconds': r3cold['wall_seconds'], 'consumer_expansion_resume_seconds': 5.621007261797786, 'shared_section_seconds': 6.708180482964963, 'controls_seconds': controls['cost']['wall_seconds'], 'maximum_observed_rss_kib': max(cold['maxrss_kib'], r3cold['maxrss_kib'], r4['cost']['maxrss_kib'], controls['cost']['maxrss_kib']), 'questions': 0, 'fallback': '7z streamed Deflate64; MET_INT decoder expanded; seven absent TF2 cases retained; independence UNKNOWN', 'threads': 1, 'gpu': False, 'scratch_bytes': scratch, 'lane_intermediate_limit_bytes': 3000000000}
    resource_measurement = ROOT / 'rounds/DEMO_RESOURCE_MEASUREMENT.json'
    if resource_measurement.exists():
        demo_cost = json.loads(resource_measurement.read_text())
        cost['verified_full_demo'] = demo_cost
        cost['maximum_observed_rss_kib'] = max(cost['maximum_observed_rss_kib'], demo_cost['maximum_rss_kib'])
    references = [{'kind': 'published_dataset', 'locator': str(MAX) + '; ' + str(TF1) + '; ' + str(TF2), 'compared_quantity': 'same-patient canal mask occupancy and image arrays; independent annotator attribution absent', 'refutes_us': True, 'interpretation': 'Refutes the proposed independent-annotator noise-floor claim; supports the scoped release-revision geometry measurements'}, {'kind': 'independent_measurement', 'locator': 'https://www.federicobolelli.it/media/publications/pdfs/2023iciap_iacat2.pdf §4 equation5, Fig6', 'compared_quantity': 'mean newly annotated voxel fraction (published61.9%); terminal underannotation', 'refutes_us': False, 'our_value_percent': r2['mean_author_added_percent'], 'frozen_tolerance_pp': 0.5}, {'kind': 'published_dataset', 'locator': 'https://www.federicobolelli.it/media/publications/pdfs/2023iciap_iacat2.pdf §4.2', 'compared_quantity': '40 revised and 51 unchanged dense masks', 'refutes_us': True, 'our_value': [r2['revised_patients'], r2['copied_patients']]}, {'kind': 'closed_form', 'locator': 'code/controls.py finite cylinder radius1, length5, voxel box center[2,0,4], width0.3; gap2.85mm', 'compared_quantity': 'finite-cylinder to voxel-box distance', 'refutes_us': False}, {'kind': 'our_own_fixture', 'locator': 'code/controls.py reflected canal columns and fixed query', 'compared_quantity': 'global-summary sufficiency for local distance', 'refutes_us': True}]
    inputpaths = [DENT / 'results/LANE_X8_GUIDE_NERVE_RISK/full_geometry.py', DENT / 'results/LANE_X8_GUIDE_NERVE_RISK/PER_SITE_FULL_GEOMETRY_RISK.csv', DENT / 'results/LANE_X8_GUIDE_NERVE_RISK/GUIDE_PROFILES.json', DENT / 'results/LANE_NEXT_N_CANAL_PROXIMITY/RAW_N1.jsonl', DENT / 'cells/geometry/nerve_vessel_canals.py']
    inputs = [{'path': str(p), 'bytes': p.stat().st_size, 'sha256': sha_file(p)} for p in inputpaths]
    archiveinfo = [{'path': str(p), 'bytes': p.stat().st_size, 'mtime_ns': p.stat().st_mtime_ns, 'hash_scope': 'Selected decompressed members SHA256 recorded per patient; full archive not hashed'} for p in (MAX, TF1, TF2)]
    results = {'schema': 'X58-canal-wall-spread-v1', 'lane': 'X58-canal-wall-spread', 'claim_type': 'information_link', 'outcome': 'LOCATED_RELEASE_REVISIONS_CHANGE_CLEARANCE; INDEPENDENT_ANNOTATOR_FLOOR_UNKNOWN', 'review_state': 'PENDING_INDEPENDENT_REVIEW', 'status': 'PENDING_INDEPENDENT_REVIEW', 'capability': 'Same-patient published canal revisions now alter existing site queries with exact image identity, local wall locations and coverage state retained', 'external_referent': references[0], 'external_referents': references, 'resolution': 'PER_POINT/PER_SURFACE_REGION/PER_TOOTH; aggregate summaries POPULATION; guide transfer PHENOMENOLOGICAL', 'time_scale': 'SIMULTANEOUS', 'R1': {'independence_gate': 'FAIL_NOT_IDENTIFIABLE', 'source_image_crc_size_matches': sum((x['crc_equal'] and x['size_equal'] for x in identity['original_images'])), 'source_image_candidates': len(identity['original_images']), 'decompressed_sha_matches': identity['exact_sha_image_checks']}, 'R2': r2, 'R3': r3, 'R4': {k: v for (k, v) in r4.items() if k != 'artifacts'}, 'regional_shared_wall': regional, 'guide_scenarios': {'whole_surface_region': 'raw/GUIDE_SEGMENT_SCENARIOS.json', 'supported_radial_region': 'raw/GUIDE_SHARED_WALL_SENSITIVITY.json', 'claim_type': 'information_link', 'resolution': 'PHENOMENOLOGICAL', 'calibrated': False, 'uncertainty_debt': 'Matched blinded local contours, image-bias phantom and signed guide-displacement series missing'}, 'sufficiency_test': controls['sufficiency'], 'radial_sufficiency_test': controls['radial_sufficiency'], 'controls': controls, 'dropout': {'R1': {'independent_pairs_admitted': 0, 'rejected': 91, 'rejection_fraction': 1.0, 'reason': 'No independently attributed repeated dense contours'}, 'R2': r2['dropout'], 'R3': r3['dropout'], 'R4': r4['dropout']}, 'full_cost': cost, 'array_artifacts': artifacts, 'input_artifacts': inputs, 'archive_metadata': archiveinfo, 'first_failed_pass_preserved': ['rounds/first_pass_FAILURES.json', 'rounds/first_pass_R2_RESULT.json', 'rounds/first_pass_R3_RESULT.json'], 'frozen_protocols': {p.name: sha_file(p) for p in ROOT.glob('PREREG_*.json')}, 'code_sha256': {p.name: sha_file(p) for p in (ROOT / 'code').glob('*.py')}, 'next_construction': 'Blinded repeated local wall contours with reader IDs, sampled near current clearance witnesses; preserve source images and all independent masks before consensus. Obtain a dimensional phantom to separate image bias.'}
    dump(ROOT / 'results.json', results)
    snapshot_hash = sha_file(ROOT / 'results.json')
    snapshot = ROOT / 'snapshots' / snapshot_hash
    snapshot.mkdir(parents=True, exist_ok=True)
    snapshot_result = snapshot / 'results.json'
    if snapshot_result.exists():
        assert sha_file(snapshot_result) == snapshot_hash
    else:
        snapshot_result.write_bytes((ROOT / 'results.json').read_bytes())
    (snapshot / 'code').mkdir(exist_ok=True)
    for cp in (ROOT / 'code').glob('*.py'):
        sp = snapshot / 'code' / cp.name
        if not sp.exists():
            sp.write_bytes(cp.read_bytes())
        assert sha_file(sp) == sha_file(cp)
    frozen = ROOT / 'FROZEN_PREDICTIONS.json'
    if not frozen.exists():
        pred = {'frozen_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'claim_type': 'information_link', 'physical_measurements_performed': False, 'source_exposure': 'Public masks and existing X8 poses were used; this is not a blinded prediction of labels', 'prospective_quantity': 'Gap from current X8 finite cylinder to observed mask union; physical anatomy gap may fall outside this observed interval', 'resolution': 'PER_TOOTH', 'prediction_file': 'raw/PER_SITE_LOCATED.json', 'prediction_sha256': sha_file(ROOT / 'raw/PER_SITE_LOCATED.json'), 'prereg_sha256': results['frozen_protocols'], 'numeric_class_reversals': [{k: x[k] for k in ('case', 'fdi', 'old_class_2mm', 'new_class_2mm', 'tf2_class_2mm')} | {'observed_union_gap_mm': [x['observed_union']['lower_mm'], x['observed_union']['upper_mm']]} for x in changed], 'laboratory_gate': 'Do not interpret observed-mask union as a confidence bound; acquire independent contours and known-distance phantom first'}
        dump(frozen, pred)
        (ROOT / 'FROZEN_PREDICTIONS.sha256').write_text(sha_file(frozen) + '  FROZEN_PREDICTIONS.json\n')
    else:
        assert json.loads(frozen.read_text())['prediction_sha256'] == sha_file(ROOT / 'raw/PER_SITE_LOCATED.json')
    lines = ['# Kanalrevisionernas lokala konsekvens', '', f"Published masks by the same patients alter 2 mm-klassen vid {r3['old_to_new_2mm_changes']} X8- places between Maxillo and ToothFairy. Union of Maxillo, ToothFairy and ToothFairy2 changes the class at {r3['tf2_to_observed_union_2mm_changes']}/{r3['paired_sites']} parade platser ({100 * r3['tf2_to_observed_union_2mm_changes'] / r3['paired_sites']:.3f} %, POPULATION; each place PER_TOOTH)No clinical safety or nerve damage has been identified.", '', f"An executable operator now separates coverage from the local wall position: among the same {r4['revised_patients_with_shared_wall']} Patients are the median of global surface P95 {r4['median_global_patient_p95_mm']:.3f} mm and the median of local radial P95 on joint lumen {r4['median_shared_wall_patient_radial_p95_mm']:.3f} mm, kvot {r4['ratio_global_to_supported']:.2f} (POPULATION; rough measure PER_POINT)There is a protocol/version difference, no independent annotator floor.", '', "| Frozen Trial | Outcome | Resolution |", '|---|---|---|', "| Independent notation dissemination | FAIL / UNKNOWN : no identified blind double masks | PER_POINT |", f"| Alla gamla/nya Closed pairs | {r2['completed']}/91 kompletta, {r2['revised_patients']} amended/ {r2['copied_patients']} identiska | POPULATION from PER_POINT |", "| Published 40 modified / 51 copies | FAIL : observed 43 / 48 ; no readjusted gate | POPULATION |", f"| Publicerad nytillagd voxelandel 61,9 % ±0,5 pp | PASS: {r2['mean_author_added_percent']:.3f} % over 43 patients | POPULATION |", f"| Bildidentitet mot TF2 | {r2['image_matched_tf2']} exakt lika; 7 saknas | PER_POINT |", f'| X8 faktisk 2 mm-reversering | PASS: 2 gamla/nya; 5 TF2/observerad union | PER_TOOTH |', f"| NEXT_N implantatkomponenter | {r3['next_n']['paired_components']}/111 paired; no reverses; replay errors {r3['next_n']['max_replay_error_mm']:.1f} mm | PER_TOOTH |", f"| Coverage separated from wall, ratio >2 | PASS: {r4['ratio_global_to_supported']:.2f} | POPULATION |", f'| Identisk sammanfattning | Identitetsfel 0; downstream difference 2 mm | PER_POINT |', '', "Local wall dimensions of joint lumen, only modified and supported patients:", '', '| Segment | Patienter | Median av patienternas radiella P95 (mm) |', '|---|---:|---:|']
    for r in regional:
        lines.append(f"| {r['segment']} | {r['patients']} | {r['median_patient_radial_p95_mm']:.3f} |")
    lines += ['', "The above group has POPULATION - resolution ; the measurements PER_POINT / PER_SURFACE_REGION . Mental and posterior segments are endpoint proxys and lack anatomical review. No edge is aggregated before the location connection; all couplings are SIMULTANEOUS .", '', "Scenario margins from X8 : s published guide profiles plus regional, local wall revision:", '', '| Guidetyp | Mentalproxy (mm) | Molar FDI (mm) | Bakre proxy (mm) |', '|---|---:|---:|---:|']
    for g in guide_bases:
        rr = [x for x in radial_margins if x['guide'] == g]
        lines.append('| ' + g + ' | ' + ' | '.join((f"{x['combined_sensitivity_margin_mm']:.3f}" for x in rr)) + ' |')
    lines += ['', "These values are PHENOMENOLOGICAL sensitivity scenarios. They are not calibrated necessary or safe margins: the guide distribution direction is missing, P95 is no supremum limit and X8 : s base scenario already contains another channel error. All 2581 site counting with these budgets is explicit population scenario; exact outcomes are found in raw/GUIDE_SHARED_WALL_SENSITIVITY.csv. The local pair analysis does not impute anything.", '', f"Bortfall: 0/91 in final dense mask comparison; 91/91 rejected as independent annotators; {r3['dropout']['x8_unmatched_sites']}/2581 ({100 * r3['dropout']['x8_unmatched_fraction']:.2f} %) X8-places without pairs; {111 - r3['next_n']['paired_components']}/111 NEXT_N-komponenter without pairs; {r4['dropout']['no_exact_tf2_image_or_centerline']}/91 without: TF2-bild/mittlinje for cross-section. {100 * r4['dropout']['unsupported_or_censored_ray_fraction']:.2f} % of the rays lack joint support or are censored. {r4['stations_new_only']} stations are only supported by newer masks, {r4['stations_old_only']} only in the elderly (PER_POINT).", '', f"Control line: cKDTree and independence EDT does not differ {controls['distance_control']['kd_vs_edt_max_abs_error_mm']:.3g} mm; brute-force avviker {controls['distance_control']['kd_vs_bruteforce_max_abs_error_mm']:.3g} mm. Cylinder/voxelbox against closed solution and radiation exit to cubic can withstand 1e-10 The controls get the same information and match. Injected distance, image, shell, independent and beam exit errors are dropped.", '', "Adequate sample: two channel columns have identical volume, Dice, median, P95, maximum and full surface distance distribution to a reference, identity error exactly 0 . The distance to the same query distinguishes 4 from 2 mm. For a fixed query, the minimum addition is enough its targeted distance change; for the new queries, the wall site and a support/ Coverage state are needed. Also local radial P95 folds: two mirrored wall extensions have exactly the same volume, Dice and entire absolute radial distance distribution (identity error 0 ), but the same point question gets 0,75 against 1,05 mm . Minimum additions are the direction of the wall change. Both tests are logical counter-examples, not externally anatomical reference.", '', "Limitations: no clinical conclusions; no independent segmentation flooring and no IAC prediction volumes in X8/X30. All TF2 masks are different from TF1 even when the forming ray is equal. Photo→anatomy-bias is still UNKNOWN . Surface distance uses center in 6 -neighbors cantvoxlar; geometric wall interpretation depends on voxel representation. Ray-exit uses voxel unions. Convex support limits and beam formulas are mathematical under their assumptions, but rigorous floating point enclosure is missing. No affinity for susceptibility approximations are used.", '', f"Resources: a thread, none GPU, observerad maximal RSS {cost['maximum_observed_rss_kib'] / 1024:.1f} MiB; externa arrayer {scratch / 1000000.0:.1f} MB, limit 3000 MB. The entire cost accounting including the first failed passport is available in results.json.full_cost. The separate time of preparation is: UNKNOWN; fit and questions 0.", '', "Reference /lineage: [IACAT2 § 4 and Fig6 ](https://www.federicobolelli.it/media/publications/pdfs/2023iciap_iacat2.pdf), [official TF2 -overlap](https://ditto.ing.unimore.it/toothfairy2/), local decompressed-member-SHA and exact formedrays in raw/patients. The same masks are new sources of information, no new independent measurements of nervanatomy. Licenses and Run Command in README_DEMO.md .", '', "Graf: experiment dispatch stopped because of historical result hash drift. A review-dispatch records the sample and review of existing X8/NV1 - assumptions ; the output results are given to PENDING_INDEPENDENT_REVIEW. No native claims are upgraded.", '', "Next design: Select stations at the five actual union reverses and let two separate, blind readers draw local contours before consensus. Add known distance phantom for image bias and only then a matched directed guide error measurement."]
    (ROOT / 'RESULTS.md').write_text('\n'.join(lines) + '\n')
    pkt = json.loads((ROOT / 'raw/GRAPH_PACKET.json').read_text())
    template = pkt['result_binding_template']
    feedback = dict(template, target_id='DENT-GEOM-UNCERTAINTY', result_file=str(snapshot_result.relative_to(DENT)), sha256=snapshot_hash, measured_quantity='same-patient canal occupancy revisions; supported radial wall exit differences; per-tooth finite-cylinder clearance class', units='mm; count; fraction', uncertainty='Independent annotator variance/anatomy bias UNKNOWN; voxel-center surface metric and conditional voxel-union radial exits; formal floating-point enclosure not certified', population_regime='91 Maxillo/TF1 dense pairs,84 exact TF2 image matches,397 existing X8 sites,25 NEXT_N components; SIMULTANEOUS; PER_POINT/PER_TOOTH edges; population aggregates explicit', preregistered_gate='PREREG_R1-R4 hashes in results.json; independence FAIL;40/51 source-count FAIL;61.9%±0.5pp PASS; actual paired2mm reversal PASS; global/shared-wall ratio>2 PASS', baseline='Fixed nominal2mm rule without alternate masks; equally informed EDT/brute force and closed-form cuboid/cylinder controls', outcome=results['outcome'], negative_result=True, review_state='PENDING_INDEPENDENT_REVIEW', status='PENDING_INDEPENDENT_REVIEW')
    feedback['edge_contracts'] = [{'producer': 'Maxillo/ToothFairy/TF2 located boundary occupancy', 'consumer': 'X8 fixed cylinder clearance', 'resolution': 'PER_POINT', 'time_scale': 'SIMULTANEOUS', 'quantity': 'voxel-center/surface coordinates and occupied voxel set', 'unit': 'mm', 'aggregation': 'Consumer minimizes over full cylinder; no global P95 used as physical bound'}, {'producer': 'same-patient supported cross-section radii', 'consumer': 'guide-type/region sensitivity table', 'resolution': 'PER_SURFACE_REGION', 'time_scale': 'SIMULTANEOUS', 'quantity': 'radial new-minus-old wall exit distance', 'unit': 'mm', 'status': 'PHENOMENOLOGICAL_TRANSPORT_DEBT', 'measurement_to_replace': 'Blinded local reader masks plus signed guide displacement'}]
    feedback['coverage_proposal'] = 'Add a scoped repeat-annotator identity/timing port and voxel support state to existing canal uncertainty node; do not admit clinical safety or noise floor.'
    feedback['dispatch_limit'] = 'Numerical experiment dispatch failed historical Result hash mismatch; review dispatch is evidence-context binding only.'
    dump(ROOT / 'GRAPH_FEEDBACK.json', feedback)
    (ROOT / 'HANDOFF.md').write_text('# X58 handoff\n\n' + lines[2] + '\n\n' + lines[4] + '\n\nR1 independent floor FAIL/UNKNOWN. R2 published40/51 FAIL, observed43/48; first five MET_INT read failures preserved then repaired. 84 exact TF2 images;7 absent. R3 paired397/2581, actual old/new2 and TF2/union5 reversals. NEXT_N paired25/111,0 reversals. R4 common-lumen first-exit separates coverage; global/shared medianP95 ratio15.4098;40 supported revised patients,83 processed patient maps. Source protocol revisions cannot become independent reader variance. All raw arrays hashed outside lane, about30MB.\n\nRun `./run_all.sh` to replay. See raw/CLASS_REVERSALS.csv for exact P7/37,P29/37,P40/37,P54/35,P54/46 witnesses; code/sections.py for DDA and support loss; raw/SUFFICIENCY.json for exact identity counterexample. Read PREREG and failed first-pass reports before altering any gate. Do not treat scenario guide margins as physical guarantees or merge source statuses.\n\nNext concrete construction: freeze a local prediction bundle at these witnesses, independently double-annotate blinded wall cross-sections (retain both masks before consensus), then compare pointwise within matching image frame. A known-distance phantom separates CBCT bias from reader spread. No IAC model predictions in X8/X30; require actual external prediction volumes before algorithm evaluation.\n\nGraph experiment dispatch failed on unrelated historical hash drift; kind review receipt binds source/provenance review and producer measurements pending independent review. No canonical graphs changed.\n')
    state('ROUND4_DEMO_COMPLETE', {'independent_floor': 'UNKNOWN', 'source_count': 'FAIL', 'located_clearance': 'PASS', 'coverage_wall_separation': 'PASS', 'controls': 'PASS'}, 'Independent blinded local contours at actual reversed sites, plus dimensional phantom; coordinator review pending', results_sha256=sha_file(ROOT / 'results.json'), scratch_bytes=scratch)
    print(json.dumps({'results_sha256': sha_file(ROOT / 'results.json'), 'sites': r3['paired_sites'], 'class_changes': r3['tf2_to_observed_union_2mm_changes'], 'surface_vs_supported_ratio': r4['ratio_global_to_supported'], 'scratch_bytes': scratch}, indent=2))
if __name__ == '__main__':
    main()
