"""Reader report: outcome and applicability, with a figure for each chain step."""
from dental_release.paths import expand as _release_expand
import datetime, hashlib, json, shutil, sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import trimesh
R = Path(__file__).resolve().parents[1]

def read(p):
    return json.loads(Path(p).read_text())

def sha(p):
    h = hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda : f.read(1048576), b''):
            h.update(b)
    return h.hexdigest()

def write(p, obj):
    Path(p).write_text(json.dumps(obj, indent=2, sort_keys=True, ensure_ascii=False) + '\n')
run = read(R / 'LATEST_RUN.json')
work = Path(run['runs'][0]['path'])
sel = read(R / 'PATIENT_SELECTION.json')
scan = read(work / 'SCAN.json')
reg = read(work / 'REGISTRATION.json')
gen = read(work / 'GENERATION.json')
gate = read(work / 'GATE.json')
mill = read(work / 'MILLING.json')
phy = read(work / 'PHYSICAL.json')
unc = read(work / 'UNCERTAINTY.json')
ana = read(work / 'ANATOMICAL_DECISIONS.json')
exp = read(work / 'EXPORT.json')
fields = np.load(work / 'DESIGN_FIELDS.npz')
scene = np.load(work / 'scene.npz')
sdf = np.load(work / 'SDF.npz')
xy = scene['xy']
faces = scene['faces']
(R / 'figures').mkdir(exist_ok=True)
plt.rcParams.update({'font.size': 10, 'axes.spines.right': False, 'axes.spines.top': False})

def save(fig, n, title):
    fig.suptitle(title, fontsize=12)
    fig.tight_layout()
    fig.savefig(R / 'figures' / f'step_{n}.png', dpi=160)
    plt.close(fig)

def roof(ax, z, cmap='viridis'):
    a = ax.tripcolor(xy[:, 0], xy[:, 1], faces, z, shading='gouraud', cmap=cmap)
    ax.set_aspect('equal')
    ax.set_xlabel('local x (conditional mm)')
    ax.set_ylabel('local y (conditional mm)')
    return a
(fig, ax) = plt.subplots(figsize=(10, 3.6))
rows = reg['regions']
vals = [r['p95_mm'] for r in rows]
ax.bar(np.arange(len(rows)), vals, color=['#bc4b39' if not r['pass_gate'] else '#437991' for r in rows])
ax.axhline(0.5, color='black', ls='--', label='frozen 0.5 mm consistency gate')
ax.set_xticks(np.arange(len(rows)))
ax.set_xticklabels([str(r['fdi']) + '/' + str(r['half']) for r in rows], rotation=90, fontsize=7)
ax.set_ylabel('sampled surface p95 (conditional mm)')
ax.legend(fontsize=8)
save(fig, 1, '1. Same-patient regional surface consistency; independent TRE UNKNOWN')
(fig, ax) = plt.subplots(1, 2, figsize=(8.5, 3.6))
a = roof(ax[0], scene['reference_z_mm'])
fig.colorbar(a, ax=ax[0], label='observed projected z')
a = roof(ax[1], fields['virtual_prep_z_mm'])
fig.colorbar(a, ax=ax[1], label='virtual preparation plane z')
save(fig, 2, '2. SDF / virtual preparation; pulp and dentin geometry absent')
(fig, ax) = plt.subplots(1, 2, figsize=(8.5, 3.6))
a = roof(ax[0], fields['outer_z_mm'])
fig.colorbar(a, ax=ax[0], label='generated outer z')
a = roof(ax[1], fields['local_thickness_mm'])
fig.colorbar(a, ax=ax[1], label='local model thickness (mm)')
save(fig, 3, '3. Packaged X60 generated restricted roof; complete axial crown UNKNOWN')
(fig, ax) = plt.subplots(figsize=(8.5, 3.6))
rules = gate['rules']
keys = list(rules)
statuses = [rules[k].get('status', 'UNKNOWN') for k in keys]
mapping = {'PASS': 1, 'FAIL': -1, 'UNKNOWN': 0}
ax.scatter([mapping.get(s, 0) for s in statuses], np.arange(len(keys)), c=[{'PASS': '#26776f', 'FAIL': '#b44c39'}.get(s, '#777777') for s in statuses], s=60)
ax.set_yticks(np.arange(len(keys)))
ax.set_yticklabels(keys)
ax.set_xticks([-1, 0, 1])
ax.set_xticklabels(['FAIL', 'UNKNOWN', 'PASS'])
ax.set_xlim(-1.3, 1.3)
save(fig, 4, '4. Design gate preserves each failure and unknown input')
(fig, ax) = plt.subplots(figsize=(6.5, 3.6))
points = np.array([r['xyz_mm'] for r in mill['rows']])
ax.scatter(points[:, 0], points[:, 1], s=45, alpha=0.5, color='#437991')
ax.set_aspect('equal')
ax.set_xlabel('local x (conditional mm)')
ax.set_ylabel('local y (conditional mm)')
ax.text(0.03, 0.92, f"{mill['found_queries']}/{mill['tool_point_queries']} tool-point witnesses\nFixture / full toolpath UNKNOWN", transform=ax.transAxes, va='top')
save(fig, 5, '5. Nominal tool access on this specimen; physical milling NOT RUN')
(fig, ax) = plt.subplots(1, 2, figsize=(8.5, 3.6))
rr = [r for r in phy['source_controls'] if r['kind'] == 'crown']
ax[0].bar(range(len(rr)), [r['error'] for r in rr], color=['#26776f' if r['gate'] else '#b44c39' for r in rr])
ax[0].axhline(0.2, color='black', ls='--')
ax[0].set_ylabel('source absolute log force error')
a = roof(ax[1], np.asarray(phy['nominal_film_um']))
fig.colorbar(a, ax=ax[1], label='nominal digital gap (um)')
save(fig, 6, '6. Published force controls; patient absolute force and wet cement film UNKNOWN')
fig = plt.figure(figsize=(6.5, 3.6))
ax = fig.add_subplot(111, projection='3d')
m = trimesh.load_mesh(work / 'crown.stl', process=False)
p = m.triangles_center
ax.scatter(p[:, 0], p[:, 1], p[:, 2], s=2, color='#437991')
ax.set_xlabel('x mm')
ax.set_ylabel('y mm')
ax.set_zlabel('z mm')
save(fig, 7, '7. Same-patient STL / 3MF export; transport PASS, design verdict retained')
(fig, ax) = plt.subplots(figsize=(7.5, 3.6))
a = roof(ax, fields['virtual_prep_z_mm'])
fig.colorbar(a, ax=ax, label='prepared point z (mm)')
ax.text(0.02, 0.95, 'PULP DISTANCE: UNKNOWN\nDEJ / PULP NOT IN LOCAL SOURCE\nNo crown exterior substituted as pulp', transform=ax.transAxes, va='top', bbox=dict(facecolor='white', alpha=0.9))
save(fig, 8, '8. Prepared surface retained; no supported anatomical pulp margin')
(fig, ax) = plt.subplots(figsize=(7.5, 3.6))
rp = np.array([r['CBCT_point_mm'] for r in unc['records'] if r['role'].startswith('CBCT')])
ax.scatter(rp[:, 0], rp[:, 2], c='#437991')
ax.set_xlabel('CBCT x (conditional mm)')
ax.set_ylabel('CBCT z (conditional mm)')
ax.text(0.02, 0.95, 'TOOTH ROOT WITNESSES ONLY\nBONE / CANAL / IMPLANT SITE: UNKNOWN', transform=ax.transAxes, va='top', bbox=dict(facecolor='white', alpha=0.9))
save(fig, 9, '9. CBCT root frame; bone and nerve geometry are missing')
(fig, ax) = plt.subplots(figsize=(7.5, 3.6))
qs = unc['records']
ax.plot([r['absolute_error_amplification_upper'] for r in qs], color='#437991')
ax.axvline(len(xy) - 0.5, color='black', ls='--')
ax.set_xlabel('prepared point, then CBCT root witness')
ax.set_ylabel('conditional bound coefficient L')
ax.text(0.02, 0.93, 'Error <= epsilon L + tissue/prep error\nepsilon is UNKNOWN; no physical error inferred', transform=ax.transAxes, va='top')
save(fig, 10, '10. Rigorous conditional pose enclosure, shared registration retained')
frozen = R / 'FROZEN_PREDICTIONS.json'
payload = dict(patient_id=sel['patient_id'], frame=scan['frame'], units=sel['units'], files={n: dict(sha256=sha(work / n), bytes=(work / n).stat().st_size) for n in ['crown.stl', 'prep.stl', 'antagonist.stl', 'DESIGN_FIELDS.npz', 'CONTRACT.json', 'REGIONS.json', 'export/model.3mf', 'export/model.3mf.json', 'EXPORT_FRAME.json']}, conditional_model_force_ratio=gen['diagnostics']['model_force_quantile_ratio'], ratio_resolution='PER_TOOTH conditional closure', absolute_force05_N=None, pulp_margin_mm=None, nerve_clearance_mm=None, bone_width_height_mm=None, seated_film_um=None, independent_registration_bound_mm=None, physical_measurement_status='NOT_RUN', freeze_rule='Before any matched physical measurement; inherited source controls are retrospective and not patient validation')
if frozen.exists():
    if read(frozen)['payload'] != payload or sha(frozen) != (R / 'FROZEN_PREDICTIONS.json.sha256').read_text().strip():
        raise ValueError('Frozen patient prediction changed')
else:
    write(frozen, dict(frozen_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), payload=payload))
    (R / 'FROZEN_PREDICTIONS.json.sha256').write_text(sha(frozen) + '\n')
for name in ['crown.stl', 'prep.stl', 'antagonist.stl', 'EXPORT_FRAME.json']:
    shutil.copyfile(work / name, R / 'exports' / name)
for name in ['model.3mf', 'model.3mf.json', 'DESIGN_GATE.json']:
    shutil.copyfile(work / 'export' / name, R / 'exports' / name)
shutil.copyfile(frozen, R / 'exports/X77_FROZEN_PREDICTIONS.json')
repair = read(work / 'GENERATOR_FAULT_REPAIR.json') if (work / 'GENERATOR_FAULT_REPAIR.json').exists() else None
old_faults = read(work / 'GENERATOR_FAULTS.json')
faults = dict(registration_numeric_fault=reg['injected_distance_plus_1mm_rejected'], **read(work / 'GATE_FAULTS.json'), wall_injection_rejected=old_faults['wall_injection_rejected'], actual_antagonist_collision_rejected=repair['actual_collision_rejected'] if repair else old_faults['antagonist_injection_rejected'], holder_fault=mill['holder_fault_rejected'], source_value_faults=all((r['injected_wrong_value_rejected'] for r in phy['source_controls'])), **unc['controls'])
full = []
for (n, label, folder) in [('FULL_CROWN.json', 'R3', 'full_crown'), ('FULL_CROWN_R4.json', 'R4', 'full_crown_R4')]:
    if (work / n).exists():
        rec = read(work / n)
        full.append(rec)
        for (k, v) in rec.get('controls', {}).items():
            if k.endswith('_rejected') or k.endswith('_pass'):
                faults[label + '_' + k] = v
        (fig, ax) = plt.subplots(figsize=(8, 3.6))
        if (work / folder / 'observed_fields.npz').exists():
            ff = np.load(work / folder / 'observed_fields.npz')
            pp = ff['exterior_points']
            dd = ff['exterior_distances']
            im = ax.scatter(pp[:, 0], pp[:, 2], c=dd, s=9, cmap='magma')
            fig.colorbar(im, ax=ax, label='distance to observed source (mm)')
            ax.set_xlabel('local x (conditional mm)')
            ax.set_ylabel('local z (conditional mm)')
            ax.text(0.02, 0.95, f"Whole exterior p95={rec['observed_exterior_p95_mm']:.4f}mm\n0.35mm gate: {('PASS' if rec['surface_preservation_gate'] else 'FAIL')}", transform=ax.transAxes, va='top', bbox=dict(facecolor='white', alpha=0.9))
        else:
            ax.axis('off')
            ax.text(0.03, 0.75, rec.get('reason', 'No specimen'), transform=ax.transAxes, wrap=True)
        save(fig, 11 if label == 'R3' else 12, label + '. Complete axial construction; missing source support remains explicit')
        if (work / folder / 'crown.stl').exists():
            target = R / 'exports' / label
            target.mkdir(exist_ok=True)
            for pp in [work / folder / 'crown.stl', work / folder / 'prep.stl', work / folder / 'DESIGN_GATE.json', work / folder / 'CONTRACT.json', work / folder / 'export/model.3mf', work / folder / 'export/model.3mf.json']:
                if pp.exists():
                    shutil.copyfile(pp, target / pp.name)
        fp = R / ('FROZEN_PREDICTIONS_' + label + '_V2.json')
        payload_full = dict(round=label, patient_id=sel['patient_id'], outcome=rec['status'], specimen_files={str(pp.relative_to(work / folder)): dict(sha256=sha(pp), bytes=pp.stat().st_size) for pp in sorted((work / folder).rglob('*')) if pp.is_file() and pp.suffix in ('.stl', '.npz', '.3mf', '.json')}, observed_surface_p95_mm=rec.get('observed_exterior_p95_mm'), physical_force05_N=None, pulp_margin_mm=None, canal_clearance_mm=None, seated_film_um=None, measurement_status='NOT_RUN', same_frame='Demo1_local_translated_CBCT', specimen_scope='Failed or unvalidated co-designed research crown/die; actual prepared tooth not measured')
        old_path = R / ('FROZEN_PREDICTIONS_' + label + '.json')
        if old_path.exists():
            old = read(old_path)['payload']
            old_files = old['specimen_files']
            new_files = payload_full['specimen_files']
            changes = [n for n in sorted(set(old_files) | set(new_files)) if old_files.get(n) != new_files.get(n)]
            old_rest = {k: v for (k, v) in old.items() if k != 'specimen_files'}
            new_rest = {k: v for (k, v) in payload_full.items() if k != 'specimen_files'}
            if old_rest != new_rest or any((n != 'export/EXPORT_RECEIPT.json' for n in changes)):
                raise ValueError('R5 repair changed a frozen scientific quantity/asset ' + label + ': ' + str(changes))
            rec['R5_binding_comparison'] = dict(changed_files=changes, scientific_quantities_exactly_unchanged=True, geometry_and_native_asset_hashes_exactly_unchanged=True)
        if fp.exists():
            if read(fp)['payload'] != payload_full or sha(fp) != (R / (fp.name + '.sha256')).read_text().strip():
                raise ValueError('Full specimen freeze drift ' + label)
        else:
            write(fp, dict(frozen_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), payload=payload_full))
            (R / (fp.name + '.sha256')).write_text(sha(fp) + '\n')
allpass = all((v is True for v in faults.values()))
steps = [dict(step=1, name='paired geometry / regional registration', status='SOURCE MEASURED + MODELLED', result=f"8/8 heldout halves pass; {reg['regions_failed']}/{reg['regions_count']} all halves fail", resolution='PER_SURFACE_REGION', figure='step_1.png'), dict(step=2, name='SDF / preparation', status='MODELLED', result='.25mm sampled capped roof; virtual 1.5mm below median plane; no pulp/dentin', resolution='PER_POINT', figure='step_2.png'), dict(step=3, name='generation', status='MODELLED', result=f"Closed roof, {gen['volume_mm3']:.6f} mm3; conditional model force ratio {gen['diagnostics']['model_force_quantile_ratio']:.6f}", resolution='PER_TOOTH conditional model', figure='step_3.png'), dict(step=4, name='design checks', status='NUMERISKT VERIFIERAT', result=gate['verdict'], resolution='PER_POINT / PER_SURFACE_REGION', figure='step_4.png'), dict(step=5, name='milling', status='MODELLED', result=f"{mill['found_queries']}/{mill['tool_point_queries']} witnesses; full machine UNKNOWN", resolution='PER_POINT', figure='step_5.png'), dict(step=6, name='fracture / cement', status='SOURCE CALIBRATION PREVENTED; PATIENT UNKNOWN', result=f"source force {phy['source_force_pass_count']}/{phy['source_force_count']}; dry gap {phy['source_gap_pass_count']}/{phy['source_gap_count']}; absolute patient force / wet film UNKNOWN", resolution='POPULATION sources; PER_POINT nominal film', figure='step_6.png'), dict(step=7, name='export', status='NUMERISKT VERIFIERAT', result=exp['status'] + ' transport; upstream ' + gate['verdict'], resolution='PER_POINT', figure='step_7.png'), dict(step=8, name='pulp / residual dentin', status='UNKNOWN', result='No same-patient pulp / DEJ. Existing X31 event operation abstains.', resolution='PER_POINT inputs, no tooth aggregation', figure='step_8.png'), dict(step=9, name='implant bone / nerve site', status='UNKNOWN', result='No same-patient bone / canal / annotated site', resolution='PER_POINT inputs', figure='step_9.png'), dict(step=10, name='pose-to-decision error', status='DERIVED CONDITIONALLY', result=f"{len(unc['records'])} exact rational point enclosures; measured epsilon UNKNOWN", resolution='PER_POINT', figure='step_10.png')]
for (i, rec) in enumerate(full):
    steps.append(dict(step=11 + i, name=rec['round'] + ' full axial crown/die', status='MODELLED; SOURCE GATE ' + ('PASS' if rec.get('surface_preservation_gate') else 'FAIL'), result=rec['status'] + '; ' + (f"source p95={rec['observed_exterior_p95_mm']:.6f}mm" if 'observed_exterior_p95_mm' in rec else rec.get('reason', 'UNKNOWN')), resolution='PER_POINT exterior; PER_TOOTH specimen', figure='step_' + str(11 + i) + '.png'))
result = dict(schema='X77-same-patient-flow-v1', claim_type='capability', review_state='PENDING_INDEPENDENT_REVIEW', status='SAME_PATIENT_DIGITAL_CROWN_DELIVERED_COMPLETE_FLOW_NOT_ESTABLISHED', full_requested_capability=False, patient=sel, external_referent=read(R / 'PREREG_R1_SAME_PATIENT.json')['external_referent'], other_external_referents=[phy['external_referent'], unc['external_referent']], steps=steps, registration=reg, design_rule_statuses={k: rules[k].get('status', 'UNKNOWN') for k in rules}, generation=gen, physical=phy, anatomical_decisions=ana, uncertainty=dict(exact_points=len(unc['records']), independent_scipy_max_error=unc['max_scipy_weight_error'], measured_registration_epsilon_mm=None, conditional_prepared_coeff_min=min((r['absolute_error_amplification_upper'] for r in unc['records'][:len(xy)])), conditional_prepared_coeff_max=max((r['absolute_error_amplification_upper'] for r in unc['records'][:len(xy)])), shared_relative_pairs=unc['shared_pairs'], sufficiency=unc['sufficiency']), controls=faults, all_faults_rejected=allpass, exclusions=dict(regional_consistency_failed=reg['regions_failed'], regional_consistency_total=reg['regions_count'], regional_failure_fraction=reg['regions_failed'] / reg['regions_count'], anatomical_quantity_inputs_missing=4, anatomical_quantity_inputs_requested=4, anatomical_missing_fraction=1.0, borrowed_other_patient_tissues_rejected=True, missing_FDI24='Inherited source release absence; not synthesized', source_rejected_roof_nodes=scan['roof_missing_nodes'], source_roof_nodes=scan['roof_nodes'], roof_missing_fraction=scan['roof_missing_fraction']), full_cost=dict(execution=run, original_fit_discovery_calibration_cost='UNKNOWN', questions=0, physical_measurements=0, intermediates_bytes=run['data_bytes'], no_speedup_claim=True), frozen_predictions_sha256=sha(frozen), raw_data=dict(directory=str(work), source_lock_sha256=sha(R / 'SOURCE_LOCK.json'), scientific_artifacts=run['scientific_artifacts']), assumption_debts=['STL mm-scale calibration missing', 'Independent directed registration TRE/epsilon absent', 'Interjaw acquisition bite/load state unknown', 'Full axial anatomy missing in packaged X60', 'No clinical finish line/preparation', 'No matched force/wet film/manufacture measurement', 'No rigorous FE/discretization/design-sensitivity enclosure'], next_construction='Attempt complete axial same-source specimen with the available full-prescan operator; independently locate same-patient CBCT tissue before any surgical clearance result.')
result['full_axial_rounds'] = full
result['per_decision_registration_error_map'] = read(work / 'DECISION_ERROR_MAP.json') if (work / 'DECISION_ERROR_MAP.json').exists() else None
result['nearest_surface_summary_counterexample'] = unc.get('nearest_surface_summary_counterexample')
result['initial_fault_instrument_failure_preserved'] = dict(round='R1', inherited_shift_mm=1.0, claimed_collision_rejection=old_faults['antagonist_injection_rejected'], diagnosis='Inherited shift1mm did not reach this patient upper obstacle. This was an inadequate instrument, preserved; R1B injects an actual collision with original checks/tolerance.', repaired=repair)
history = []
for pp in sorted(work.parent.parent.glob('*')):
    failure = pp / 'FAILURE.json'
    receipt = pp / 'RECEIPT.json'
    if failure.exists():
        history.append(dict(run_id=pp.name, status='FAILED_PRESERVED', failure=read(failure), failure_path=str(failure), failure_sha256=sha(failure)))
    elif receipt.exists():
        history.append(dict(run_id=pp.name, status=read(receipt)['byte_identity_gate'], receipt_path=str(receipt), receipt_sha256=sha(receipt)))
failed_binding = [h['failure_path'] for h in history if 'Scientific byte drift' in str(h.get('failure', {}))]
result['artifact_binding_repair'] = dict(round='R5', original_failure='Only full_crown/export/EXPORT_RECEIPT.json differed; stale pre-canonical asset and sidecar SHA', original_frozen_bindings=['FROZEN_PREDICTIONS_R3.json', 'FROZEN_PREDICTIONS_R4.json'], corrected_bindings=['FROZEN_PREDICTIONS_R3_V2.json', 'FROZEN_PREDICTIONS_R4_V2.json'], physical_or_geometry_predictions_changed=False, original_failed_run_artifacts=failed_binding)
write(R / 'ATTEMPTS.json', dict(constructions=[dict(round='R1', outcome='Restricted crown digital chain; complete capability remains PARTIAL', obstacle='Pulp/bone/canal missing, full axial absent'), dict(round='R1B', outcome='Adverse collision instrument repaired; old noncollision remains preserved'), dict(round='R2', outcome='Exact conditional shared-pose point enclosures; physical epsilon UNKNOWN'), *[dict(round=rec['round'], outcome=rec['status'], obstacle=rec.get('reason', 'Unsupported completion exterior, source .35mm gate failed'), raw_result=str(work / ('FULL_CROWN.json' if rec['round'] == 'R3' else 'FULL_CROWN_R4.json'))) for rec in full], dict(round='R5', outcome=run['byte_identity_gate'], obstacle='Old stale pre-canonical receipt repaired with unchanged science')], runs=history))
result['next_construction'] = 'Acquire same-patient full IOS with verified axial/margin ownership and calibrated CBCT pulp/DEJ/bone/canal, plus directed independent registration error. R3/R4 completion hypotheses cannot supply missing measured anatomy.'
write(R / 'results.json', result)
table = '\n'.join(('| ' + str(s['step']) + ' ' + s['name'] + ' | ' + s['result'] + ' | ' + s['status'] + ' | <img src="figures/' + s['figure'] + '" width="160"> |' for s in steps))
readme = _release_expand(f"""# Same patient through the crown chain and anatomical quiz ports\n\nRun `./run_same_patient.sh` . It only uses local Demo 1 data from [Liu/Hao 2023 ](https://doi.org/10.5281/zenodo.8027553), calls unchanged parts in DEMO48_PACKAGE in read-only environment and makes two new runs . No GPU, network access or physical measurement are used. Requires existing Linux/bubblewrap and the Python environment `@DENTAL_EXTERNAL_ROOT@/miniconda3` .\n\n** FDI36 : s digital crown chain is now running in the same CBCT reference as the patient's reconstructed teeth. The entire flow with pulp, bone and nerve is still not established.** The registered release available local contains teeth; the actual CBCT volume and the requested tissue masks are missing for this case. The crown is a limited occlusal shell. All mm - values requires the source's STL scale; physical scale calibration is unknown.\n\n| Step | Results | Knowledge status | Figure |\n|---|---|---|---|\n{table}\n\nRegistration facite is **separate IOS against CBCT -only-tendreconstruction**, never the published fusion result that has already replaced crowns with IOS . All {reg['regions_count']} regional p95 is available in `REGISTRATION_PER_REGION.csv` ; {reg['regions_failed']} fails , including the preserved regional errors in X69 . FDI36 : s two halves provide {reg['selected_tooth_regions'][0]['p95_mm']:.6f} and {reg['selected_tooth_regions'][1]['p95_mm']:.6f} mm. These are the surface consistency of the sample and are **not independent registration error or rigorous bounds**. A single agent/ RMS does not provide a safe pulp or nerve distance.\n\nEach prepared point has an exact rational weight vector. If the four independently targeted anchor measurements later give bounds epsilon_i , "||D(x)|| <= sum |lambda_i(x)| epsilon_i`. Towards a fixed CBCT - tissue, the distance interval is increased by this and the geometry errors of the tissue/preparation. Two objects sharing the same registration use `lambda(x)-lambda(y)`; a common rigid pose does not change the crown's own thickness, gap and tool spacing. Contact with second jaw requires both poses and actually bet state. Epsilon is unknown; no p95 has been inserted in its place. The [external affine reference](https://doc.cgal.org/latest/Barycentric_coordinates_3/namespaceCGAL_1_1Barycentric__coordinates.html) refers to mathematics, not patient accuracy.\n\nThe adequacy test has:  identity error  **0  exact** for  four  ankarmagnituder,  RMS  and  p95 . Two targeted rigid translations give **0,5 mm** difference in distance to a fixed tissue point. The minimum extension for the query is the signed displacement in the direction of the tissue; general queries needs the common targeted pose. This is its own mathematical fixture, not anatomically reference . `UNCERTAINTY.json` contains the entire sample and the rigorous conditional enclosures; the whole design box of the mechanics still lacks rigorous enclosure .\n\nControl outcomes : imported 3MF preserves geometry, the external surface distance sensor gives error {reg['published_triangle_control_error_mm']:.1g} mm , all injected errors rejected={allpass}. Published breaking load groups are tested against [Prott Table 1 ](https://doi.org/10.4047/jap.2021.13.5.269) and [Chen Table 2 ](https://doi.org/10.3390/ma17020365); a group fails the unchanged 0,20 logo. These group data do not calibrate the absolute load of the new shell. Nominal digital gap is not cement film after seating. Tool witnesses cover points and nominal tools, not a complete machine or fixture.\n\n`exports/` contains the same specifics STL, 3MF, retained design report, CBCT transformation and [frozen predictions](FROZEN_PREDICTIONS.json). Manufacturing release is UNKNOWN . The patient-dependent physical predictions are null; this freezes ignorance before future measurement. [ results.json ](results.json) binds raw data /hashes; `MINIMUM_MEASUREMENT.json` indicates what the lab needs to measure.\n\nLicense: use Liu/Hao Demo 1 files are **CC BY 4.0 according to the frozen published Zenodo metadatan in X69**. Figshare001, ToothFairy2, Pulpy3D, MMDental, STS, Bits2Bites, Teeth3DS and mandible defects are not used as this patient's geometry. The previous second datasets appear in inherited code/source control, not as new patient measurements. No personal data or journal fields are used and nothing is published. Package parts retain their original licenses. Status **PENDING_INDEPENDENT_REVIEW**.\n""")
if full:
    block = '\n## Full axial shape: additional constructions\n\n'
    for rec in full:
        block += f"**{rec['round']}: {rec['status']}**. " + (f"Observerad yta → kandidat p95={rec['source_to_ext_p95_mm']:.6f} mm, kandidat → observerad yta p95={rec['ext_to_source_p95_mm']:.6f} mm. Bidirectional gate 0,35 mm : {('PASS' if rec['surface_preservation_gate'] else 'FAIL')}. Slutet forskningsspecimen={rec['complete_axial_specimen']}.\n\n" if 'source_to_ext_p95_mm' in rec else rec.get('reason', 'UNKNOWN') + '\n\n')
    block += "R3 calls existing PROOF_LANE_FULL_CROWN in read-only mode as a complement to the package. R4 only changes the unsaturated source closure from a low base point to the limit local 3D medium; all's original design and design gates are unchanged. A virtual closure is a hypothesis. The next crucial information is an actual full intraoral surface area with verified tooth/margin ownership; raw CBCT/labels must belong to exactly the same patient. No second patients' masks have been inserted.\n\n" + f"R1 : s old +1 mm collision instruments were insufficient for this patient's large jaw gap and failure is preserved. R1B : s actual crossing on {repair['new_shift_mm']:.6f} mm is precipitated by unchanged antagonist control.\n"
    readme += block
(R / 'README_DEMO.md').write_text(readme)
(R / 'README.md').write_text(readme)
(R / 'RESULTS.md').write_text('# The same-patient-chronic chain is dead; tissue issues are unsubstantiated\n\n' + readme.split('**FDI36', 1)[1].split('**', 1)[0] + '\n\n' + f"Registrering: {reg['regions_count']} halvor, {reg['regions_failed']} Failed, 8/8 original holdout halves passed. Full digital design: {gate['verdict']}. Export: {exp['status']} transport. Brottlastfacit: {phy['source_force_pass_count']}/{phy['source_force_count']}; dry-gap-reference {phy['source_gap_pass_count']}/{phy['source_gap_count']}. Pulp / DEJ / bone /channel: 4 / 4 UNKNOWN . No clinical recommendation.\n\nR2 provides {len(unc['records'])} exact conditional query weights, maximum SciPy deviation {unc['max_scipy_weight_error']:.3g}. Satisfaction test: identity error 0 , downstream difference 0,5 mm. See README_DEMO, results.json and raw UNCERTAINTY.\n")
handoff = f"""# X77 R1/R2 handoff\n\n`./run_same_patient.sh` run same-patient Demo1 FDI36 on local data. Read README_DEMO.md , results.json and LATEST_RUN.json . Code/source-lock and prereg is frozen before X77 number. Raw data : {work}. All mm -tal conditional on the STL scale. PENDING_INDEPENDENT_REVIEW.\n\nR1 : the crown chain is executed via unchanged DEMO48_PACKAGE ; 101 roof nodes, closed limited roof, not full axial crown . Gate {gate['verdict']}, export transport {exp['status']}, {mill['found_queries']}/{mill['tool_point_queries']} verktygspunkter. Modellens force-ratio={gen['diagnostics']['model_force_quantile_ratio']:.9g} under uncalibrated Closure; no absolute load. Gruppfacit force {phy['source_force_pass_count']}/{phy['source_force_count']}; dry gap {phy['source_gap_pass_count']}/{phy['source_gap_count']}. The source loss and gate are retained.\n\n54 regional p95 reconstructed exactly from X69 : s immutable vectors; 3 fails , all 8 holdout halves pass. The external resolution is teeth , NOT pulp / bone / channel. No data from Figshare001 / ToothFairy /Pulpy was borrowed. Independent TRE , epsilon, bite state and scale UNKNOWN . 4 / 4 surgical/tissue issues UNKNOWN .\n\nR2: {len(unc['records'])} rational affee queries, {len(unc['shared_pairs'])} gemensamma par. All exact reconstruction and half-band-fault injection tested; SciPy maximum error {unc['max_scipy_weight_error']:.3g}. Common rigid pose leaves its own chronometry invariable; fixed CBCT -tissue requires absolute closure. Sammanfattningsidentitet0 with rigid translation+/-.25 provides fixed-tissue distance difference. 5mm . It's your own fix, not patient facit.\n\nPreserved build errors: freeze.py KeyError on already immutable witness path; first numerical import missing nibabel; physical metadata doubled resolution. All failed run directories exist and COMMANDS binds them. Clean computing environment/installation on new computer is not tested.\n\nNext design: full scan of the axial specific from the same FDI36 with the existing full-prescan operator; do not pretend that X60 roof is whole crown. For surgery, the smallest change is new same-patient CBCT tissues + independent targeted registration anchors with uniform errors. The source's Demo1 tooth-STL cannot identify pulp, DEJ, bone, or channel. Registered teeth alone does not define implant location.\n"""
(R / 'HANDOFF.md').write_text(handoff)
(R / 'HANDOFF_R1_R2.md').write_text(handoff)
if full:
    handoff += '\nR3 / R4 exchange-bound result :\n\n' + '\n'.join((rec['round'] + ': ' + rec['status'] + '; ' + ('p95=' + str(rec['observed_exterior_p95_mm']) + 'mm, source gate=' + str(rec['surface_preservation_gate']) if 'observed_exterior_p95_mm' in rec else rec.get('reason', 'UNKNOWN')) for rec in full)) + '\n\nR3 closed source boundary spans z[-2.99674,1.85524]mm. Full specimen can be closed with a co-designed die but not inferred anatomically valid. R4 changes the unmeasured cap closure and keeps original gates. Do not promote inward-only preservation to full-shell anatomy. Next: same-patient original full IOS and independent margin/region labels, plus CBCT volume/pulp/bone/canal and directed registration anchor errors.\n'
    (R / 'HANDOFF.md').write_text(handoff)
    for rec in full:
        (R / ('HANDOFF_' + rec['round'] + '.md')).write_text(handoff)
    with (R / 'RESULTS.md').open('a') as f:
        f.write('\n\n' + '\n'.join((rec['round'] + ': ' + rec['status'] + '; ' + (f"bidirectional p95={rec['observed_exterior_p95_mm']:.6f}mm versus frozen.35mm." if 'observed_exterior_p95_mm' in rec else rec.get('reason', 'UNKNOWN')) for rec in full)) + '\n')
write(R / 'MINIMUM_MEASUREMENT.json', dict(patient_id=sel['patient_id'], physical_measurements_run=0, required=[dict(quantity='CBCT same-patient pulp/DEJ,bone/canal', resolution='PER_POINT', measurement='Original calibrated-volume coordinates and expert region masks; independent error/scale evidence'), dict(quantity='directed registration error', resolution='PER_POINT', measurement='Four affinely independent same-patient calibration anchors with independent uniform directed error bounds; ideal rigid pose requires >=3 noncollinear landmarks but that alone is not a uniform error certificate'), dict(quantity='bite state', resolution='PER_ARCH plus PER_POINT contact', measurement='Same-session bite registration and source jaw state before contact-to-force use'), dict(quantity='wet cement film/force', resolution='PER_POINT film + PER_TOOTH force', measurement='Matched fabricated X77 specimen with seated-film field, geometry/QC, setup/material/batch/support/contact, force trace and independent noise model')], freeze_sha256=sha(frozen)))
write(R / 'GRAPH_FEEDBACK.json', dict(target_id='DENT-UNC-PROPAGATION', result_file=str(R / 'results.json'), sha256=sha(R / 'results.json'), measured_quantity='Same-source Demo1 regional IOS-CBCT tooth-surface consistency and conditional shared-pose distance enclosure; missing surgical quantities remain unknown', units='conditional mm; dimensionless conditional amplification', uncertainty='Independent TRE/scale and anchor epsilon UNKNOWN; exactly rational algebra does not calibrate patient truth', population_regime='One publisher-derived pair Demo1, sampled regions, virtual restricted FDI36 roof, no clinical inference', preregistered_gate='PREREG_R1_SAME_PATIENT.json + PREREG_R2_CORRELATED_DECISIONS.json; p95 .5mm / exact identity0 / SciPy1e-8 / source force .2log', baseline='Unchanged package geometry checks/3MF importer, published triangle oracle and source endpoint control; capability claim, no algorithm contest', outcome=result['status'], negative_result=True, review_state='PENDING_INDEPENDENT_REVIEW', source_generation=read(R / 'raw/GRAPH_PACKET.json').get('generation'), coverage_proposal='Same-patient semantic tissue ports and directed shared-pose uncertainty, at PER_POINT before aggregation; source tooth exterior not interchangeable with pulp/bone/canal'))
print('Report generated. Full requested physical flow UNKNOWN; digital stages executed.', flush=True)
if not allpass:
    raise RuntimeError('An injected control accepted a false input')
