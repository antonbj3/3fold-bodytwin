"""Reader-facing results, scientific figures and exact handover from an executed run."""
from dental_release.paths import expand as _release_expand
import hashlib, json, math, os
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import trimesh
R = Path(__file__).resolve().parent
DATA = Path(_release_expand('@DENTAL_WORK_ROOT@/X70_end_to_end'))

def read(p):
    return json.loads(Path(p).read_text())

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def write(p, x):
    Path(p).write_text(json.dumps(x, indent=2, sort_keys=True, ensure_ascii=False) + '\n')

def figures(p):
    figdir = R / 'figures'
    figdir.mkdir(exist_ok=True)
    (fig, axs) = plt.subplots(2, 4, figsize=(14, 6), layout='constrained')

    def draw(ax, k):
        if k == 0:
            for (name, color) in [('upper', '#cc6262'), ('lower', '#579abe')]:
                m = trimesh.load_mesh(p / (name + '.stl'), process=False)
                v = m.vertices[::max(1, len(m.vertices) // 3000)]
                ax.scatter(v[:, 0], v[:, 1], s=0.3, c=color, label=name)
            ax.set_aspect('equal')
            ax.set_title('1 Registered scan pair')
            ax.set_xlabel('scan x (mm)')
            ax.set_ylabel('scan y (mm)')
        elif k == 1:
            m = trimesh.load_mesh(p / 'segmented_tooth.stl', process=False)
            v = m.vertices[::max(1, len(m.vertices) // 1500)]
            ax.scatter(v[:, 0], v[:, 2], s=1, c='#5e91a8')
            ax.set_title('2 Predicted tooth 36')
            ax.set_xlabel('local x (mm)')
            ax.set_ylabel('local z (mm)')
        elif k == 2:
            a = np.load(p / 'SDF.npz')
            phi = a['phi_roof_mm']
            h = float(a['pitch_mm'])
            lo = a['origin_mm']
            j = phi.shape[1] // 2
            ax.imshow(phi[:, j, :].T, origin='lower', cmap='coolwarm', extent=[lo[0], lo[0] + h * (phi.shape[0] - 1), lo[2], lo[2] + h * (phi.shape[2] - 1)], vmin=-2, vmax=2, aspect='equal')
            ax.axhline(read(p / 'SDF.json')['preparation_height_mm'], color='black', lw=1, label='virtual plane')
            ax.set_title('3 SDF + virtual preparation')
            ax.set_xlabel('local x (mm)')
            ax.set_ylabel('local z (mm)')
        elif k == 3:
            a = np.load(p / 'DESIGN_FIELDS.npz')
            n = len(a['local_thickness_mm'])
            xyz = a['xyz_mm'][:n]
            q = ax.scatter(xyz[:, 0], xyz[:, 1], c=a['local_thickness_mm'], s=12, cmap='viridis')
            ax.set_aspect('equal')
            ax.set_title('4 X60 roof: local thickness')
            ax.set_xlabel('local x (mm)')
            ax.set_ylabel('local y (mm)')
            plt.colorbar(q, ax=ax, label='mm', shrink=0.65)
        elif k == 4:
            g = read(p / 'GATE.json')
            names = ['mesh', 'scale', 'film max', 'film min', 'contact']
            vals = [g['rules'][n]['status'] for n in ['mesh_health', 'scale', 'film_max', 'film_min', 'occlusal_contact']]
            ax.barh(names, [1] * 5, color=[{'PASS': '#44896b', 'UNKNOWN': '#ca9e43', 'FAIL': '#b54c50'}[s] for s in vals])
            ax.set_xlim(0, 1.6)
            ax.set_xticks([])
            for (i, s) in enumerate(vals):
                ax.text(1.04, i, s, va='center', fontsize=8)
            ax.set_title('5 Design FAIL; transport PASS')
        elif k == 5:
            m = read(p / 'MILLING.json')
            rr = m['rows']
            ax.scatter([r['xyz_mm'][0] for r in rr], [r['xyz_mm'][1] for r in rr], s=20, c=['#44896b' if r['status'] == 'FOUND' else '#b54c50' for r in rr])
            ax.set_aspect('equal')
            ax.set_title('6 Tool poses: 8 points only')
            ax.set_xlabel('local x (mm)')
            ax.set_ylabel('local y (mm)')
        elif k == 6:
            q = read(p / 'PHYSICAL.json')
            rr = q['independent_source_force_control']
            x = np.array([r['observed_N'] for r in rr])
            y = np.array([r['predicted_N'] for r in rr])
            ax.scatter(x, y, c=['#44896b' if r['gate'] else '#b54c50' for r in rr])
            hi = max(x.max(), y.max()) * 1.1
            ax.plot([0, hi], [0, hi], c='#777', lw=0.8)
            ax.set_xlim(0, hi)
            ax.set_ylim(0, hi)
            ax.set_title('7 External force: 4/5')
            ax.set_xlabel('published group mean (N)')
            ax.set_ylabel('source prediction (N)')
        else:
            a = np.load(p / 'metrology/point_fields.npz')
            ax.scatter(a['xyz_mm'][:, 0], a['xyz_mm'][:, 1], c=a['normal_deviation_mm'][0] * 1000, s=10, cmap='coolwarm', vmin=-25, vmax=25)
            ax.set_aspect('equal')
            ax.set_title('8 Digital metrology; lab UNKNOWN')
            ax.set_xlabel('local x (mm)')
            ax.set_ylabel('local y (mm)')
    for (k, ax) in enumerate(axs.flat):
        draw(ax, k)
    fig.savefig(figdir / 'flow.png', dpi=150, metadata={'Software': 'X70 deterministic scientific plot'})
    plt.close(fig)
    for k in range(8):
        (fig, ax) = plt.subplots(figsize=(3.4, 2.6), layout='constrained')
        draw(ax, k)
        fig.savefig(figdir / f'step_{k + 1}.png', dpi=100, metadata={'Software': 'X70 deterministic scientific plot'})
        plt.close(fig)

def main():
    rec = read(R / 'LATEST_RUN.json')
    p = Path(rec['runs'][0]['directory'])
    g = read(p / 'GATE.json')
    gen = read(p / 'GENERATION.json')
    phy = read(p / 'PHYSICAL.json')
    mill = read(p / 'MILLING.json')
    sc = read(p / 'SCAN.json')
    domain = read(p / 'ANTAGONIST_DOMAIN.json')
    meta = read(p / 'METROLOGY.json')
    lab = read(p / 'LAB_ALARM.json')
    figures(p)
    faults = [dict(name='SDF distance', rejected=read(p / 'SDF.json')['wrong_distance_rejected']), dict(name='generator thin wall', rejected=read(p / 'GENERATOR_FAULTS.json')['wall_injection_rejected']), dict(name='generator wrong axial position', rejected=read(p / 'GENERATOR_FAULTS.json')['antagonist_injection_rejected']), dict(name='gate wrong mesh hash', rejected=read(p / 'GATE_FAULTS.json')['wrong_input_hash_rejected']), dict(name='gate wrong units', rejected=read(p / 'GATE_FAULTS.json')['conflicting_units_rejected']), dict(name='holder collision', rejected=mill['holder_fault_rejected']), dict(name='external force/gap wrong value', rejected=all((r['injected_bad_value_rejected'] for r in phy['source_controls']))), dict(name='metrology wrong facet map', rejected=read(p / 'METROLOGY_FAULTS.json')['wrong_region_hash_rejected']), dict(name='X61 wrong patient', rejected=lab['packaged_X61_patient_ID_rejected'])]
    gap = g['rules']['occlusal_contact']['measurement']['minimum_gap_mm']
    summary = phy['summary_sufficiency']
    ext = sc['external_referent']
    outcomes = dict(schema='X70-end-to-end-results-v1', lane='X70-end-to-end', claim_type='capability', review_state='PENDING_INDEPENDENT_REVIEW', outcome='PARTIAL_EXECUTABLE_PATIENT_RESEARCH_CHAIN', full_capability_gate='FAIL_UNKNOWN_FULL_CROWN_AND_PHYSICAL_CALIBRATION', digital_chain=dict(patient_case='27', case_key='14ceebb036f59e99', requested_steps=8, executed_steps=8, scientific_artifacts=len(rec['scientific_artifacts']), bit_identity_gate=rec['byte_identity_gate'], differing_files=len(rec['different_files']), same_generated_mesh_sha256=sha(p / 'crown.stl'), same_patient_sources=sc['sources'], frame=gen['frame'], raw_directory=str(p)), design_gate=dict(verdict=g['verdict'], rules={k: v['status'] for (k, v) in g['rules'].items()}, minimum_antagonist_gap_mm=gap, contact_policy_mm=[0.02, 0.2], policy_resolution='PHENOMENOLOGICAL', gap_resolution='PER_SURFACE_REGION', rule_debt='No empirical justification for scenario0.02–0.20mm contact band; replace with measured lab-specific contact acceptance protocol', export='PASS_TRANSPORT_ONLY'), generator=dict(scope=gen['geometry_scope'], model_ratio=gen['diagnostics']['model_force_quantile_ratio'], ratio_resolution='PER_TOOTH', ratio_status='CONSTITUTIVE_CLOSURE_NOT_CALIBRATED_PATIENT_FORCE', volume_mm3=gen['volume_mm3'], volume_resolution='PER_TOOTH', absolute_force05_N=None), milling=dict(point_tool_queries=mill['tool_point_queries'], found=mill['found_queries'], sample_points=mill['queried_surface_points'], resolution='PER_POINT', area_coverage='NOT_ESTIMATED', physical_status=mill['full_machine_feasibility']), external_force=dict(passed=phy['source_force_pass_count'], compared=phy['source_force_count'], tolerance_log=0.2, resolution='POPULATION', rows=phy['independent_source_force_control']), external_gap=dict(passed=phy['source_gap_pass_count'], compared=phy['source_gap_count'], tolerance_um=20, resolution='PER_SURFACE_REGION source group mean, POPULATION evidence', patient_seated_film='UNKNOWN'), summary_sufficiency=summary, metrology=dict(measurement_kind='OUR_OWN_DIGITAL_CONTROL_NOT_EXTERNAL_FACIT', known_shift_um=25, shift_resolution='PER_POINT', intaglio_signed_bias_um=meta['regions']['intaglio']['signed_bias_um'], bias_resolution='PER_SURFACE_REGION', measured_parts=0, independent_scanner_error='UNKNOWN'), lab_alarm=lab, dropout=dict(cohort_cases=8, selected_before_outputs=1, not_selected=7, selection_explanation='Frozen deterministic first case; not outcome-based rejection', source_scan_reduction=domain, missing_projected_tooth_nodes=dict(count=sc['native_projected_reference_missing'], total=sc['native_projected_reference_nodes'], fraction=sc['missing_fraction'], reason='No native projected triangle at node; nearest measured projected height is a closure'), external_patients_with_matched_physical_measurements=0, force_protocol_support_rejected=1, reason='No matched patient specimen endpoints/setup'), fault_injections=faults, all_faults_rejected=all((x['rejected'] for x in faults)), external_referent=ext, external_referents=[phy['external_referent'], dict(kind='independent_measurement', locator='doi:10.7759/cureus.38688 Table2; doi:10.1155/2023/6698453 Table1', compared_quantity='Dry marginal group gap in um; not patient seated cement film', refutes_us=True)], same_information_control='Unchanged GenCAD all_checks plus independent endpoint log interpolation, naive nearest triangles and published 3MF importer; verification only, no superiority claim.', affine_enclosure='MISSING: no rigorous design sensitivity, continuous SDF, floating-point distance or mesh-convergence enclosure', edges=[dict(producer='registered scan + predicted labels', consumer='source patient site', resolution='PER_POINT', timescale='SIMULTANEOUS'), dict(producer='virtual preparation + X60 shell', consumer='X49/X53', resolution='PER_POINT', timescale='SIMULTANEOUS'), dict(producer='export geometry', consumer='X55 metrology', resolution='PER_POINT', timescale='HANDOVER'), dict(producer='published group force/dry gap', consumer='support gate only', resolution='POPULATION', timescale='SIMULTANEOUS', patient_prediction_transfer='REJECTED'), dict(producer='manufactured matched specimen', consumer='patient-specific calibrated likelihood/alarm', resolution='PER_POINT film and PER_TOOTH force', timescale='HANDOVER', status='UNKNOWN_MEASUREMENT_ABSENT')], cost=dict(wall_seconds_two_runs=rec['wall_seconds'], max_child_RSS_MiB=rec['max_child_RSS_MiB'], own_data_bytes=rec['data_bytes'], threads_max=4, GPU=False, preparation='Directed source read/hash; wall time UNKNOWN', fit='No new fit; historical segmentation/model fit costs UNKNOWN', discovery='One source patient/task; no cohort search for good outcome', validation='Two fresh eight-stage runs, source held controls, exact summary witness and actual mutations', questions=0, fallback='Physical acquisition time/material cost UNKNOWN; lab protocol supplied'), frozen_predictions=dict(path='FROZEN_PREDICTIONS.json', sha256=sha(R / 'FROZEN_PREDICTIONS.json'), physical_outcomes_observed_before_freeze=0), limitations=['Restricted occlusal shell, no complete axial anatomy', 'X11 segment labels and frame not measured PER_POINT ground truth', 'SDF capped projected roof, not enamel/dentin/root', 'X60 support/load/material closure not patient calibrated', 'Actual sinter, holder, fixture, seated film and force unmeasured', 'X61 D1/M1/M2 interface rejects X70 ID; requires distinct measured profile', 'No statistical patient physical benefit inferred'], next_construction='Matched X70 lab acquisition: measured specimen geometry/datum, per-region seated film and failure trace plus independent calibration specimens; or reviewed complete axial generator. Preserve patient/site protocol and freeze new profile before held validation.')
    write(R / 'results.json', outcomes)
    rows = []
    descriptions = [('Registered scan', 'Measured coordinates; pose error UNKNOWN'), ('Segmentation', '17,628 predicted tooth facets; accuracy UNKNOWN'), ('SDF/preparation', '0.25mm grid; virtual plane, no dentin'), ('Generator', 'Closed occlusal shell; full axial crown UNKNOWN'), ('Design/export', 'Gap0.205591mm → FAIL; 3MF transport PASS'), ('Milling', '24/24 point/tool witnesses; fixture UNKNOWN'), ('Force/film', 'Source force4/5, dry gap3/3; patient calibration UNKNOWN'), ('Metrology/lab', 'Digital shift25um recovered; physical data/alarm UNKNOWN')]
    for (i, (label, msg)) in enumerate(descriptions):
        rows.append(f'| {i + 1}. {label} | {msg} | <img src="figures/step_{i + 1}.png" width="130"> |')
    (R / 'README_DEMO.md').write_text('# One patient, one research command\n\nRun `./run_end_to_end.sh`. It extracts only the two local registered scans for **Bits2Bites case27**, invokes reviewed DEMO48_PACKAGE code in read-only mounts, and performs two fresh patient flows. Linux, bubblewrap, the installed scientific Python environment and the local archive are required; no network, GPU or new model fit is used. Output: exact specimen, spatial fields, 3MF, frozen predictions and a lab protocol.\n\n**The executable digital flow is delivered; the complete calibrated patient-crown flow remains PARTIAL.** The existing fallback X60 generates an occlusal shell, with no complete axial anatomy. Published group measurements test source calibration; they do not calibrate this new specimen.\n\n| Step | Result and evidence status | Figure |\n|---|---|---|\n' + '\n'.join(rows) + f"\n\nTwo fresh runs: **{len(rec['scientific_artifacts'])} scientific artifacts byte-identical**, zero differences; {rec['wall_seconds']:.2f}s total, {rec['max_child_RSS_MiB']:.1f}MiB peak child RSS. Scientific coordinates/uncertainties/verdicts are unchanged; runtime clocks are separate telemetry. [Raw paths and exact numbers](results.json), [run receipt](LATEST_RUN.json), [full figure](figures/flow.png).\n\nAll distances/film remain PER_POINT or PER_SURFACE_REGION until the consumer aggregates. The exact-summary test has identity error0 and a45.72% conditional force difference: keep the local thickness, support and load together. R1 normal serialization failure and R2 clock/hash failure remain in [ATTEMPTS.json](ATTEMPTS.json). Nine actual injected faults must reject.\n\nThe external anatomical referent is the [registered Bits2Bites dataset](https://ditto.ing.unimore.it/bits2bites/). Force facit: [Prott Table1](https://doi.org/10.4047/jap.2021.13.5.269), [Chen Table2](https://doi.org/10.3390/ma17020365); dry marginal group gap: [Table2](https://doi.org/10.7759/cureus.38688), [Table1](https://doi.org/10.1155/2023/6698453). One held source force target fails its fixed0.20-log criterion. Controls with the same information verify the quantities; no algorithm superiority is claimed.\n\nWhat does not hold: generated absolute fracture load, physically seated spatial cement film, scanner/segmentation accuracy, whole-machine CAM, clinical eligibility or rigorous sensitivity/mesh-error enclosures. The 0.02–0.20mm contact band and shrink/tool assembly are lab scenarios with measurement debts. Metrology's known25um translation is a software control, not a fabricated part. X61 rejects this patient ID; its old D1/M1/M2 posterior is not updated.\n\n[Frozen predictions](FROZEN_PREDICTIONS.json) precede any physical outcome. Follow [LAB_PROTOCOL.md](LAB_PROTOCOL.md); use `./metrology_lab.sh --scans measured.stl --metadata datum.json` for a real scan, with the frozen specimen/regions. The original export's X49 prediction file describes its gate verification suite, not this patient; the separate X70 freeze is authoritative for patient geometry. No numeric physical prediction has been invented.\n\nLicence: Bits2Bites **CC BY-NC-SA per the supplied local dataset record**; the public landing page does not display the licence, so exact redistribution terms must be checked before distribution. This private demo references local files and distributes no archive. Existing code retains its original research/software licences; article scalar values have per-row attribution. No personal records or images are used. ToothFairy2, Teeth3DS, MMDental, STS and mandibular-defect meshes are not used as this patient's geometry. Review: PENDING_INDEPENDENT_REVIEW.\n")
    (R / 'RESULTS.md').write_text(f"# Registered patient geometry through a research workflow\n\nBits2Bites case27 now runs through eight digital stages with one command. The measured upper/lower scan pair is the anatomical referent, with selected member hashes and source facet correspondence. The output is a restricted occlusal specimen; calibrated complete crown capability remains unsupported.\n\n| Frozen gate | Executed outcome |\n|---|---|\n| Exact scientific rerun | {len(rec['scientific_artifacts'])} files, {rec['byte_identity_gate']}, zero differing files |\n| Patient/source identity | Same case/frame; 17,628 model-predicted tooth facets, accuracy UNKNOWN |\n| SDF vs naive triangle oracle | 0mm difference at32 surface-centroid queries; sampled grid0.25mm, no global error bound |\n| X60 original geometry predicates | PASS for this restricted source task; local rigid-model force ratio {gen['diagnostics']['model_force_quantile_ratio']:.12f}, not calibrated force |\n| X49 design | FAIL: projected minimum gap {gap:.12f}mm lies above fixed0.20mm; wall/film/insertion unknown scopes remain |\n| X53 point/tool | {mill['found_queries']}/{mill['tool_point_queries']} positive poses over8 points; no area fraction or complete machine certificate |\n| External held source force | 4/5 under0.20-log error; every target retained in results.json |\n| External dry group marginal gap | 3/3 under20um; no conversion to this patient's seated film |\n| Generated absolute force / actual film | UNKNOWN; no matched endpoints, manufactured specimen or wet-seating measurement |\n| 3MF | Published importer/metadata transport PASS; failed upstream design retained |\n| X55 metrology | Known25um digital translation: intaglio signed bias {meta['regions']['intaglio']['signed_bias_um']:.9f}um; physical trueness UNKNOWN |\n| X61 | Patient ID rejects; no old profile/posterior borrowed |\n| Adversarial controls | {sum((x['rejected'] for x in faults))}/{len(faults)} rejected |\n\nthe adequacy test: the identity-thickness multiset, volume/min/max/mean has identity error **{summary['identity_error']}**; conditional force ratioB/A **{summary['conditional_quantile_ratio_B_over_A']:.12f}**, downstream difference **{100 * summary['downstream_relative_difference']:.9f}%**. Minimal extension for this pair: spatial thickness together with support and load. This is an existing X60 closure witness, not a physical measurement; no rigorous affine design-box enclosure is available.\n\nR1 retained the entire arch and failed8 stored-normal comparisons. R2 applies the existing v4 XY crop with original facet IDs, excluding {domain['excluded_outside_query']}/{domain['source_total_facets']} facets outside this particular projected query ({100 * domain['excluded_fraction']:.5f}%); zero projected-degenerate exclusions. Binary-coordinate rounding precedes stored-normal evaluation. A separate R2 failure found query timing fields inside hashed reports; its bytes and unmodified scientific outcomes remain in the data directory.\n\nThe independent source endpoint interpolation agrees with the original source predictions; neither provides a patient force law. Every source row has quantity/locator. The generated overall crown is not clinically or physically validated. The fixed gap policy is PHENOMENOLOGICAL and its replacement measurement is a lab acceptance protocol on matched, seated parts.\n\nData and commands: LATEST_RUN.json hashes all declared outputs under {p}; COMMANDS.md preserves process arguments and failures. Physical outcomes0; costs include two fresh runs{rec['wall_seconds']:.3f}s, peak child RSS{rec['max_child_RSS_MiB']:.3f}MiB, existing fitting and physical acquisition costs UNKNOWN. Source modules are invoked unchanged and read-only. No graph scientific admission, no publication and no new physics model.\n")
    (R / 'HANDOFF.md').write_text("X70-end-to-end: eight stages execute for the one frozen Bits2Bites case27. Restricted occlusal shell, virtual preparation, local spatial ports and 3MF are usable research artifacts. Full axial anatomy, calibrated patient force, seated film and actual lab data remain UNKNOWN. Read README_DEMO.md, results.json and LATEST_RUN.json before reuse.\n\nR1 mesh failure (8 antagonist stored-normal disagreements) is preserved; R2 used existing projected crop plus source-facet map and normals from serialized coordinates. Clock fields caused a separate failed byte gate; source raw reports remain. The final two fresh runs have exactly matching scientific files. Actual X49 gap0.2055910485557242mm fails unchanged scenario upper0.20mm; do not tune it to pass. X60 original geometry passes only its restricted roof contract. Source force4/5 and dry gap3/3 are external-source checks, not patient calibration. Exact summary identity0 hides45.7183336% conditional force difference. Keep local geometry/load/support fields. All physical/continuous numerical error enclosures remain scoped or missing.\n\nX61 rejects this patient's ID rather than updating D1/M1/M2. Next physical construction: LAB_PROTOCOL.md, same specimen material/process/die/contact protocol, calibrated datum/triple scan, assembled spatial film and fracture trace, independent calibration specimens and a separately frozen X70 profile before held validation. A new reviewed complete axial generator is the parallel geometric prerequisite; PROOF_LANE_FULL_CROWN was empty at selection. Do not invent an absolute-force posterior from a group mean or a thickness summary.\n\nGraph dispatch to DENT-MFG-ASBUILT-DEVIATION failed while reading a historical result binding with Result hash mismatch. DENT-GOAL packet was unavailable/ambiguous; --query is unsupported by installed rank CLI. Preserve the existing binding. GRAPH_FEEDBACK.json stays PENDING_INDEPENDENT_REVIEW with a coverage proposal and exact outcome/result hash. No source graphs changed. Lab readiness is a reviewable private digital specimen/protocol, not release.\n")
    coverage = dict(proposed_id='DENT-CAP-PATIENT-RESEARCH-HANDOVER', consumer='DENT-MFG-ASBUILT-DEVIATION', scope='Common patient/mesh/frame/protocol ports from registered IOS through generation to matched fabrication/metrology/force/film', existing_target='DENT-MFG-ASBUILT-DEVIATION', status='PENDING_INDEPENDENT_REVIEW', prerequisites=['Reviewed complete axial generator', 'Matched X70 laboratory observation and likelihood profile'], resolution='PER_POINT geometry/film; PER_TOOTH specimen force', timescale='HANDOVER', new_source_information='Common scan-facet and specimen hash context survives module composition; unsupported source population transfer rejected', canonical_graph_mutation=False)
    write(R / 'GRAPH_COVERAGE_PROPOSAL.json', coverage)
    feed = dict(target_id='DENT-MFG-ASBUILT-DEVIATION', claim_type='capability', status='PENDING_INDEPENDENT_REVIEW', review_state='PENDING_INDEPENDENT_REVIEW', result_file=str(R / 'results.json'), sha256=sha(R / 'results.json'), measured_quantity='Same-patient scientific artifact byte identity; projected gap; source protocol support; digital specimen handover', units='bytes, mm, N, um', resolution='PER_POINT / PER_SURFACE_REGION -> PER_TOOTH; external group responses POPULATION', time_scale='SIMULTANEOUS geometry; HANDOVER lab', uncertainty='Measured scanner, segmentation/frame, physical transfer and rigorous numerical sensitivity enclosures UNKNOWN; do not fuse', population_regime='Single Bits2Bites patient27; virtual horizontal preparation and X60 restricted occlusal roof; retrospective published group controls', preregistered_gate=dict(files=['PREREG_X70_R1.json', 'PREREG_X70_R2.json', 'PREREG_RUNTIME_SEPARATION.json'], exact_byte_tolerance=0, full_capability='FAIL_UNKNOWN', source_force_log=0.2, source_gap_um=20), baseline=outcomes['same_information_control'], outcome='PARTIAL_CAPABILITY_WITH_FAILED_DESIGN_AND_UNKNOWN_PHYSICAL_PORTS', negative_result=True, external_referent=ext, graph_registration='Dispatch blocked by previous Result hash mismatch; no old binding repaired', missing_coverage_proposal='GRAPH_COVERAGE_PROPOSAL.json', next_operation=outcomes['next_construction'])
    write(R / 'GRAPH_FEEDBACK.json', feed)
    write(R / 'ARTIFACT_MANIFEST.json', dict(files={str(p / k): v for (k, v) in rec['scientific_artifacts'].items()}, figures={str(q.relative_to(R)): dict(sha256=sha(q), bytes=q.stat().st_size) for q in sorted((R / 'figures').glob('*.png'))}, external_data_root=str(DATA), large_array_policy='All arrays under X70 data directory; no >50MB arrays written to root'))
if __name__ == '__main__':
    main()
