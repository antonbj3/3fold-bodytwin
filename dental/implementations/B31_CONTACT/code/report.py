from dental_release.paths import expand as _release_expand
from common import *
from geometry import load_case
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
import time, platform, scipy, shapely

def figures(r1, r2):
    out = ROOT / 'figures'
    out.mkdir(exist_ok=True)
    rows = r1['rows']
    bounds = r2['rows']
    labels = [f"{i + 1}. {x['family'].split('_')[0]}" for (i, x) in enumerate(rows)]
    x = np.arange(len(rows))
    (fig, ax) = plt.subplots(1, 2, figsize=(13, 4.7), gridspec_kw={'width_ratios': [1.8, 1]})
    for (off, name, label, color) in [(-0.25, 'uniform', 'Best uniform (numerical upper endpoint)', '#8796a5'), (0, 'regional', 'Regional relief', '#d48937'), (0.25, 'ordinary_mesh', 'Ordinary mesh LP', '#5e91ad')]:
        vals = [r['outputs'][name]['contact']['symdiff_mm2'] for r in rows]
        ax[0].bar(x + off, vals, 0.25, label=label, color=color)
    floors = [r['nominal']['unavoidable_symdiff_floor_mm2'] for r in bounds]
    ax[0].plot(x, floors, 'ko', label='Lower bound: EVERY frozen-family relief')
    ax[0].set_xticks(x, labels, rotation=60, ha='right')
    ax[0].set_ylabel('Spatial mask mismatch [mm²; PER_SURFACE_REGION]')
    ax[0].legend(fontsize=8)
    ax[0].set_title('R1: earlier wall-cap proxy; fixed supplied pose')
    addition = [r['nominal']['minimum_added_height_max_mm'] * 1000 for r in bounds]
    ax[1].barh(np.arange(len(rows)), addition, color='#a36191')
    ax[1].set_yticks(x, labels)
    ax[1].invert_yaxis()
    ax[1].set_xlabel('Required positive height at source patch [µm; PER_POINT]')
    ax[1].set_title('Lower bound, new manufacture\nNo added height was executed')
    fig.suptitle('X95 R1/R2: earlier proxy and source-height witness — geometry only', fontsize=12)
    fig.tight_layout()
    fig.savefig(out / 'summary.png', dpi=160)
    fig.savefig(out / 'summary.svg')
    plt.close(fig)
    (q, _) = parents()
    chosen = [1, 2, 4]
    (fig, axes) = plt.subplots(2, len(chosen), figsize=(13, 8))
    pr = read(ROOT / 'PREREG_R1.json')
    for (col, i) in enumerate(chosen):
        row = rows[i]
        r = pr['cohort'][i]
        t = load_case(r)
        ff = t['grid_faces']
        xy = t['xy']
        source = t['reference_gap']
        for (ia, name) in enumerate(['uniform', 'regional']):
            ax = axes[ia, col]
            with np.load(row['outputs'][name]['path']) as a:
                pred = a['predicted_gap']
            for (field, color, label, alpha) in [(source, '#268b9a', 'Source-native target', 0.65), (pred, '#e39639', name, 0.6)]:
                (ps, _, _, _) = q.polygons(xy, ff, field)
                verts = [p @ xy[ff[j]] for (j, p) in ps.items()]
                ax.add_collection(PolyCollection(verts, facecolors=color, edgecolors='none', alpha=alpha, label=label))
            (rr, _, _, _) = q.polygons(xy, ff, source)
            missing = []
            for (j, p) in rr.items():
                ids = ff[j]
                pp = q.Q.clip(p, t['original_gap'][ids] - 0.1)
                if len(pp) >= 3:
                    missing.append(pp @ xy[ids])
            ax.add_collection(PolyCollection(missing, facecolors='none', edgecolors='#87315c', linewidths=0.6, label='Source patch needs height addition'))
            ax.set_xlim(xy[:, 0].min(), xy[:, 0].max())
            ax.set_ylim(xy[:, 1].min(), xy[:, 1].max())
            ax.set_aspect('equal')
            ax.set_xlabel('Local x [mm]')
            ax.set_ylabel('Local y [mm]')
            ax.set_title(f"{i + 1}. {r['family']} / {name}\nΔmask={row['outputs'][name]['contact']['symdiff_mm2']:.6g}mm²", fontsize=10)
            if col == 0:
                ax.legend(fontsize=7, loc='upper right')
    fig.suptitle('Spatial source patches and proposals; axes in unchanged site frame\nPER_SURFACE_REGION, SIMULTANEOUS; 0≤gap≤100µm is an operational geometry band', fontsize=11)
    fig.tight_layout()
    fig.savefig(out / 'spatial_masks.png', dpi=180)
    plt.close(fig)

def run():
    r1 = read(ROOT / 'raw/ROUND1.json')
    r2 = read(ROOT / 'raw/ROUND2.json')
    r3 = read(ROOT / 'raw/ROUND3.json') if (ROOT / 'raw/ROUND3.json').exists() else None
    pr = read(ROOT / 'PREREG_R1.json')
    ctrl1 = read(ROOT / 'raw/CONTROLS_R1.json')
    ctrl2 = read(ROOT / 'raw/CONTROLS_R2.json')
    sw = read(ROOT / 'raw/SUFFICIENCY_AND_TRANSLATION.json')
    figures(r1, r2)
    arrays = {}
    for p in sorted(DATA.rglob('*.npz')):
        arrays[str(p)] = dict(sha256=sha(p), bytes=p.stat().st_size)
    env = dict(python=sys.version, numpy=np.__version__, scipy=scipy.__version__, shapely=shapely.__version__, platform=platform.platform(), threads=1, GPU=False)
    dump(ROOT / 'raw/ENVIRONMENT.json', env)
    cost = dict(acquisition=read(ROOT / 'raw/ACQUISITION.json'), construction_seconds=read(ROOT / 'raw/CONSTRUCTIONS_R1.json')['seconds'], R1_validation_seconds=r1['summary']['seconds'], R2_seconds=r2['summary']['seconds'], R3_seconds=r3['summary']['seconds'] if r3 else None, upstream_scanning_and_generation='UNKNOWN', human_equivalent_discovery_and_tokens='UNKNOWN', fit_calibration=0, physical_measurement='NOT_RUN', intermediate_bytes=sum((x['bytes'] for x in arrays.values())), timing_resolution='PHENOMENOLOGICAL', all_cpu_step_costs_reported=True)
    result = dict(lane='X95-spatial-contact-repair', claim_type='capability', status='SOURCE_SPATIAL_REPAIR_EXECUTABLE_WITH_REMOVAL_FAMILY_REFUTATION_PHYSICAL_UNKNOWN', review_state='PENDING_INDEPENDENT_REVIEW', external_referent=pr['external_referent'], cohort=pr['cohort'], rounds=dict(R1=r1, R2=r2, R3=r3), decision=dict(inherited_50pct_3of10_gate='FAIL', R1_nominal_50pct_cases=0, R2_all_family_possible_cases=1, R2_3of10='IMPOSSIBLE_IN_FROZEN_REMOVAL_FAMILY', full_mask_new_manufacture_required_cases=6, physical_passes=0, commercial_superiority='NOT_ESTABLISHED', algorithm_gain='NOT_ESTABLISHED'), dropout=dict(requested=10, retained=10, rejected=0, rejected_fraction=0.0, reasons={}, failed_case_gate=10, full_shape_and_process_eligibility='UNKNOWN', upstream_selection='Exactly X85R5 frozen10/18 crown sites; remaining8 outside this brief, not replaced or silently rejected', data_records_read='Geometry arrays only; no patient text/images'), summary_sufficiency=sw, controls=dict(R1=ctrl1, R2=ctrl2, R3=None if r3 is None else [dict(key=r['key'], local_order=[o['wall_and_removal']['local_vertical_interval_order_pass'] for o in r['outputs'].values()], downward_crossing_fault=r['downward_crossing_fault']) for r in r3['rows']]), resolution=dict(contact_areas='PER_SURFACE_REGION', gap_displacement_addition='PER_POINT', J_area_integral='PER_TOOTH', case_decisions='PER_TOOTH', counts='POPULATION', cost_estimates='PHENOMENOLOGICAL'), timescale='SIMULTANEOUS', edges=read(ROOT / 'GRAPH_COVERAGE_PROPOSAL.json')['edges'], uncertainty=dict(source_geometry_enclosure='MISSING', rigorous_machine_enclosure='MISSING', pose_boxes='Complete conditional real-arithmetic continuum set bounds for +/-10 and +/-50um, not measured probabilities; R1 maximum discrepancy tolerance0.01mm2', baseline_robust_optimum='R1 gives a nominal scalar minimum bracket. Robust comparison is to this fixed scalar; full boxwise optimal scalar is not certified. R1 fails nominally. R4 uses the shared local constraints and its three nominal gains all fail in explicitly observed frozen pose scenarios.', original_solid_validity='UNKNOWN_STRICT_SELF_INTERSECTION_PROOF', basal_annulus='Excluded from inherited wall bound', clinical_force_adjustment_time='UNKNOWN'), full_cost=cost, environment=env, array_artifacts=arrays, frozen=dict(R1_prereg_sha256=sha(ROOT / 'PREREG_R1.json'), R2_prereg_sha256=sha(ROOT / 'PREREG_R2.json'), R3_prereg_sha256=sha(ROOT / 'PREREG_R3.json'), predictions_sha256=sha(ROOT / 'FROZEN_PREDICTIONS.json'), R3_predictions_sha256=sha(ROOT / 'FROZEN_PREDICTIONS_R3.json') if r3 else None), graph_scope='Definition/missing spatial coverage receipt only. Experiment dispatch at material crown target rejected; source graphs unchanged.', limitations=['Minimum is for a sampled-band lexicographic LP family, not global minimum integrated binary-mask repair', 'External source triangles are used for design and retrospective scoring, not an independent loaded-contact prediction', 'Required added heights are source-grid lower bounds, no addition performed', 'No original-anatomy fidelity guarantee, tool holder clearance, manufacturing transfer or force model', 'No10x time/adjustment/clinical result measured'])
    r4 = read(ROOT / 'raw/ROUND4.json') if (ROOT / 'raw/ROUND4.json').exists() else None
    if r4:
        result['rounds']['R4'] = r4
        result['decision']['R4_nominal_50pct_vs_local_constrained_proxy'] = r4['summary']['nominal_50pct_vs_new_local_proxy']
        result['decision']['R4_robust_cases_not_felled'] = r4['summary']['cases_not_felled_by_observed_poses']
        result['decision']['R4_inherited_gate'] = r4['summary']['inherited50pct3of10']
        result['decision']['R2_proxy_scope'] = 'Limited wall-cap/margin family only; R4 repairs its local-boundary comparator eligibility gap'
        result['full_cost']['R4_seconds'] = r4['summary']['seconds']
        result['frozen']['R4_prereg_sha256'] = sha(ROOT / 'PREREG_R4.json')
        result['frozen']['R4_case_prediction_hashes'] = {p.name: sha(p) for p in ROOT.glob('FROZEN_PREDICTIONS_R4_*.json')}
        result['controls']['R4'] = [dict(key=r['key'], controls=r['controls']) for r in r4['rows']]
    if (ROOT / 'raw/COPIED_REPLAY.json').exists():
        result['copied_replay'] = read(ROOT / 'raw/COPIED_REPLAY.json')
    dump(ROOT / 'results.json', result)
    if not r4:
        raise RuntimeError('Reader report requires all four completed constructions; raw earlier results remain intact')
    latest = ['| Case / tooth type | Uniform R4 mm² | Regional R3 mm² | Vanlig mesh-LP R3 mm² | Missing patch even without depth limit mm² | Nominal ≥ 50% / pose rejects |', '|---|---:|---:|---:|---:|---|']
    previous = ['| Case | Uniform R1 optimumintervall mm² | Regional R1 mm² | R2 floor in older family mm² | Required positive height , max µm |', '|---|---:|---:|---:|---:|']
    for (i, (a, b, c, d)) in enumerate(zip(r1['rows'], r2['rows'], r3['rows'], r4['rows'])):
        yes = lambda v: 'yes' if v else 'no'
        latest.append(f"|{i + 1} {a['case_key'][:8]} / {a['family'].split('_')[0]}|{d['uniform_contact']['symdiff_mm2']:.9g}|{d['regional_R3_contact_mm2']:.9g}|{c['outputs']['ordinary_mesh']['contact']['symdiff_mm2']:.9g}|{d['cap_independent_missing_source_area_mm2']:.9g}|{yes(d['nominal50pct_vs_local_proxy'])} / {yes(d['observed_pose_ratio_failure'])}|")
        ub = a['outputs']['uniform']['optimization']['symdiff_optimum_bracket_mm2']
        previous.append(f"|{i + 1}|[{ub[0]:.9g}; {ub[1]:.9g}]|{a['outputs']['regional']['contact']['symdiff_mm2']:.9g}|{b['nominal']['unavoidable_symdiff_floor_mm2']:.9g}|{1000 * b['nominal']['minimum_added_height_max_mm']:.9g}|")
    latest_table = '\n'.join(latest)
    previous_table = '\n'.join(previous)
    (fig, ax) = plt.subplots(figsize=(11, 5))
    xx = np.arange(10)
    uu = [r['uniform_contact']['symdiff_mm2'] for r in r4['rows']]
    reg = [r['regional_R3_contact_mm2'] for r in r4['rows']]
    ordinary = [r['outputs']['ordinary_mesh']['contact']['symdiff_mm2'] for r in r3['rows']]
    ax.bar(xx - 0.25, uu, 0.25, color='#8796a5', label='Uniform with same local facet constraints')
    ax.bar(xx, reg, 0.25, color='#d48937', label='Regional repair with local facet constraints')
    ax.bar(xx + 0.25, ordinary, 0.25, color='#5e91ad', label='Ordinary mesh LP, same information/constraints')
    ax.plot(xx, [r['cap_independent_missing_source_area_mm2'] for r in r4['rows']], 'ko', label='Missing source contact: even uncapped removal')
    labels = [str(i + 1) + (' *' if r['nominal50pct_vs_local_proxy'] else '') for (i, r) in enumerate(r4['rows'])]
    ax.set_xticks(xx, labels)
    ax.set_xlabel('Case (*: nominal 50% gain; EVERY starred case fails fixed pose scenarios)')
    ax.set_ylabel('Mask mismatch [mm²; PER_SURFACE_REGION]')
    ax.set_title('X95 R4: shared local operation constraints; final robust gate FAIL')
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(ROOT / 'figures/summary_r4.png', dpi=170)
    fig.savefig(ROOT / 'figures/summary_r4.svg')
    plt.close(fig)
    idea = "A executable operation can now suggest local material removal on ten full crowns and point Output contact patches that the operation cannot recover. In 6 / 10 crowns, the source's patches require positive height even if the felling depth would be unlimited. Reference are the registered [Bits2Bites ](https://ditto.ing.unimore.it/bits2bites/) and [Bite2Text ](https://ditto.ing.unimore.it/bite2text/) IOS plots, hash bound per data object. The source geometry is also used in the design; this is a retrospective digital trial, no independent prediction of loaded contact.\n\nThe representation retains the specim, coordinate frame, antagonist and spatial gap mask. The inside and margin are kept exactly. The regional proposal minimizes a declared nodal band error surrogate and then surface sale-integrated change within its field family. It's not a proof of the smallest global repair of the binary contact mask. Contact is defined operationally as 0 ≤ gap ≤ 100 µm ; distance is not converted to Newton power."
    outcome = 'The frozen robust gate is **FAIL**: 3/10 regional proposal reduces nominal mask error at least 50% against the uniform proxy that has the same local facet condition, but all three fails in the frozen pose scenarios. The requirement was at least 3/10 without regression in the inherited gates and with advantage in both ± 10 and ± 50 µm boxes. All ten case is reported; loss 0 / 10 , four original case cluster.\n\nA common, equally informed constrained mesh-LP runs on the same crowns and with the same information and conditions. Its outcome is the same or mixed; no methodological advantage is established. Uniform offset is an explicitly limited proxy. That comparison establishes no superiority against the entire [exocad antagonist-conscious workflow](https://wiki.exocad.com/wiki/index.php/AI_Crown_Design).'
    limits = "What does not holds : no physical contact, adjustment volume, active time, clinical usability or 10 x-saving are validated. The Gap Mask is a geometric definition, not a measured pressure mask. Rigorous floating-point and scanner inclusion is missing; the original soliden's full self-intersection and the entire CAM process is UNKNOWN. The inherited wall boundary excludes the basal annulus. Pose scenarios have no measured probability. Therefore, no physical or manufacturing PASS is given.\n\nThe minimum cost applies to the nodala LP surrogate and specified field family. Surface-integrated change is PER_TOOTH and not actual output volume. Mask areas are PER_SURFACE_REGION , gap/ height PER_POINT and number of cases POPULATION . The time limit to the mask comparison is SIMULTANEOUS . Cost and pose assumptions are PHENOMENOLOGICAL liabilities that need to be replaced by measurement."
    license_text = f"License bound to the data object: {sum((r['dataset'] == 'Bits2Bites' for r in pr['cohort']))}/ 10 crowns references are Bits2Bites , CC BY - NC - SA according to user briefing and inherited local binding ; full archive terms are not shown on the public page. {sum((r['dataset'] == 'Bite2Text' for r in pr['cohort']))}/ 10 is Bite2Text, whose exact archive conditions are UNKNOWN. No redistribution is included. Dataset, locator, file and SHA256 are bound in PREREG_R1.json /cohort, INPUT_LOCK_R1.json and results.json /cohort. Only geometry is read, no patient reports or pictures. Generated crowns are derived digital designs, not measured crowns."
    run_text = _release_expand('```bash\ncd @DENTAL_IMPLEMENTATIONS@/B31_CONTACT\n./run_all.sh\n```\n\nrun all four command designs, facit comparisons, falsifiers and figures. Existing Python environment: `@DENTAL_PYTHON@` with NumPy, SciPy, Shapely and Matplotlib. A CPU wire, no GPU, network is not required. Absolute inputs is hashbound. `X95_DATA_ROOT` specifies the output directory separately; `X95_UPSTREAM_ROOT` can specify the same hashade dependencies on another site . No source tree or full zip extraction is required.')
    measurement = "The next design needs source-owned direct triangle overlap and controlled scan/pose inclusion on a valid solid that the lab owns. Where there is a positive witness to heights, the new proposal needs to be produced. Freeze proposal and prediction before manufacturing ; compare before/after adjustment scan in the same frame and measure active time. [Zhang 2019 , Abstract Methods](https://pubmed.ncbi.nlm.nih.gov/31357225/) describes this separate measuring port; the study's second material and the numbers of individuals do not calibrate this demo. See LAB_PROTOCOL.md and HANDOFF.md . X94 owns any cargo validation."
    replay = read(ROOT / 'raw/COPIED_REPLAY.json') if (ROOT / 'raw/COPIED_REPLAY.json').exists() else None
    replay_text = 'Separately copied rerun : exit 0 and exactly the same mathematical outcome for all four designs, including artifacthashar, gates and interval ends. Se raw/COPIED_REPLAY.json.' if replay and replay['scientific_output_exact_match'] else 'Separately copied rerun still lacks a complete receipt.'
    demo = f"# X95 — spatial contact repair on full crowns\n\n{idea}\n\n{outcome}\n\n{latest_table}\n\nAll areas in the table are mm² on PER_SURFACE_REGION . The missing patch is a lower failure floor for fixed XY - material removal even without depth limit. Zero baselines do not provide relative improvement. Decimals are rounded here; full precision is available in results.json .\n\n![Latest Comparison](figures/summary_r4.png)\n\n{limits}\n\n{run_text}\n\n{replay_text}\n\nResearch STL is available in export/ <key>/R3_regional.stl and R3_ordinary_mesh.stl. Older suggestions are kept in the same catalogue. They're not CAM-validated. Data/ R2 contains spatial polygon witnesses and elevation fields; no elevation is performed by them. All arrays are available on the games disk with hash in results.json.\n\n{license_text}\n\n{measurement}\n\nReview : PENDING_INDEPENDENT_REVIEW . The graph receipt refers to the definition of missing spatial ports; no source graph has been upgraded. Older negative result and contracts are preserved in RESULTS.md and raw/ROUND1–ROUND4.json.\n"
    (ROOT / 'README_DEMO.md').write_text(demo)
    local_caps = [1000 * r['uniform_local_cap_mm'] for r in r4['rows']]
    detailed = f"# Spatial crown repair and its operating limit\n\n{idea}\n\n{outcome}\n\n{latest_table}\n\n![Latest Comparison](figures/summary_r4.png)\n\nMaskareor: PER_SURFACE_REGION. Height and local cap: PER_POINT . Integral of Amendment: PER_TOOTH . Count: POPULATION. SIMULTANEOUS temporal edge. Decimals here are rounded; the results have full precision in results.json .\n\n## What the local trial covers\n\nR3 tests 3 206 503 overlap corners between the relevant upper and lower full crown facets. All ten regional proposal passes the conditional local order lines; four changes compared to R1. Conditions in the vertices of the overlap polygon include the entire affinity polygon under real-time rituals and provided valid original solidity. It does not replace a full self-intersection check.\n\nR4 cuts the scale interval of the uniform change with the same local facet conditions before optimization. Caps Becomes {min(local_caps):.6g}–{max(local_caps):.6g} µm ; these are model values with declared 1e-8 mm tolerance, no measured milling limits. All ten changes have direct mesh/contact control. Inherited minimum wall dimensions follow by old same-hash-boundary minus maximum change; protected vertices have identity error zero . Both wall and order arguments lack rigorous machine closure.\n\nR1 has 60 fully covered continuous scenario boxes with conditional real-ritmetic upper/lower set-bounds. R3 only re-uses them for bitidetic vertices and otherwise recalculates; R4 counts their uniform boxes separately. The tolerance for box width is not more than 0.01 mm². Source/machine closure is MISSING . The nominal scalar optimum of the uniform design is enclosed numerically, but its best scalar optimum over the box is not certified. R4 : s three nominal benefits are dropped by explicitly calculated common poses even before such a robust optimum is needed. 25 f 66 f 15 in R1 can have approximately 7.82 mm² worm error at 50 µm against 0.47344 mm² nominally; thus, the selection of pose can change the answer.\n\n## Preserved previously gates\n\nR1 gave 0/10 nominal 50% enhancements against its previous wall-cap/margin proxy. R2 : s relaxed all-operational floors allowed only one potential case against this particular proxy/family. These outcomes are preserved, but are not a universal or physical impossibility statement. R4 corrects the local operating limit of the uniform control. The stronger positive height witness applies to six cases even without any depth limit, within the fixed XY operation.\n\n{previous_table}\n\nMesh areas are PER_SURFACE_REGION; raising max PER_POINT . Premolar 25 f 66 f 15 missed R1 : s exact 50% - boundary with 0.000184224619 mm² . No threshold was changed to catch this almostoutfall. Increasemax is required height in the source's projected representation, not performed surgery.\n\n![ The source's spatial patches and older suggestions](figures/spatial_masks.png)\n\n## Sufficiency and precipitating controls\n\nTwo state has bitidentical contact area, number, gap p95/mean, change integral and Lipschitz wall boundary. The identity error of the summary is precise {sw['summary_identity_error']}; the applied downstream symmetrical mask difference is changed by: {sw['downstream_difference_mm2']} mm². The minimum extension for this fixed query is overlaparea; a repair operation needs the spatial support and frame of the mask. The witness is a mathematical fixture, not empirically reference.\n\nA truly executed 1 mm mask translation by polygon surgery is folded. All ten thin-intakelio injections fall. All ten local facet crossings fold. Invalid flooring, +3 mm² mask measurement error and modified local-cap is folded by their respective controls. See raw/CONTROLS_R1.json, CONTROLS_R2.json, SUFFICIENCY_AND_TRANSLATION.json and ROUND3/ROUND4.json. The first analytical translation check is replaced as proof of the later performed check; the old one remains.\n\n{limits}\n\n{replay_text}\n\n{license_text}\n\n{measurement}\n\nFreezings, assumption sheets, cost, object-hashar, all cases and failure are available in PREREG_R1 – R4.json , DECOMPOSITION_R1 – R4.json , FROZEN_PREDICTIONS *.json and results.json . Freezing takes place before the subsequent direct digital control; physical measurement is NOT_RUN . See COMMANDS.md for exact commands and preserved failures, CLAIM_CORRECTIONS.md for delimited corrections. PENDING_INDEPENDENT_REVIEW ; graph has only one definition receipt.\n"
    (ROOT / 'RESULTS.md').write_text(detailed)
    return result
if __name__ == '__main__':
    run()
