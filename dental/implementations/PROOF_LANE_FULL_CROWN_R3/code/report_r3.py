from common_r3 import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import time
FAMILIES = ['molar_crown', 'premolar_crown', 'anterior_crown']
NAMES = {'molar_crown': 'Molar', 'premolar_crown': 'Premolar', 'anterior_crown': 'Anterior'}

def run():
    from diagnose_clipping import run as diagnose
    clipping = diagnose()
    from figure_contact import run as plot_contact
    plot_contact()
    rounds = {tag: read(ROOT / 'rounds' / f'{tag}.json') for tag in ['A', 'B', 'C']}
    old = read(OLD / 'rounds/R2.json')
    native = read(ROOT / 'raw/NATIVE_SELF_CONTACT.json')
    labs = read(ROOT / 'LAB_EXPORTS.json')
    source_files = {}
    for r in read(ROOT / 'FROZEN_COHORT.json')['rows']:
        path = V4 / 'payload/whole_private' / r['key'] / 'reference.npz'
        source_files[r['key']] = dict(path=str(path), sha256=sha(path), bytes=path.stat().st_size, resolution='PER_POINT mesh vertices/facets', label_status='X11_INFERRED_NOT_INDEPENDENT_TOOTH_ANNOTATION')
    table = []
    for (label, rr) in [('R2 prior', old), ('A deformation', rounds['A']), ('B explicit boundary', rounds['B']), ('C crossing information', rounds['C'])]:
        for (method, fams) in rr['summary'].items():
            for (fam, s) in fams.items():
                table.append(dict(round=label, method=method, family=fam, **s))
    source_meta = {r['key']: r for r in read(ROOT / 'FROZEN_COHORT.json')['rows']}
    for rr in rounds.values():
        for r in rr['rows']:
            r['split'] = source_meta[r['key']]['split']
    split_tables = []
    for (tag, rr) in rounds.items():
        for method in rr['summary']:
            for fam in FAMILIES:
                for split in ['test', 'auxiliary']:
                    s = [r for r in rr['rows'] if r['participant'] == method and r['family'] == fam and (r['split'] == split)]
                    split_tables.append(dict(round=tag, method=method, family=fam, split=split, requested=len(s), scored=sum((r['status'] == 'SCORED' for r in s)), anatomy=sum((r.get('anatomy_pass', False) for r in s)), function=sum((r.get('functional_pass_sampled', False) for r in s))))
    result = dict(lane='PROOF_LANE-full-crown-r3b', claim_type='capability', review_state='PENDING_INDEPENDENT_REVIEW', scientific_admission=False, external_referent=REFERENT, cohort_sha256=sha(ROOT / 'FROZEN_COHORT.json'), rounds=rounds, table=table, split_tables=split_tables, source_manifest=source_files, native_self_contact=native, premolar_predecessor_regions=read(ROOT / 'raw/PREMOLAR_REGIONS_R2.json'), clipping_sufficiency=read(ROOT / 'raw/CLIPPING_SUFFICIENCY.json'), crossing_acquisition=read(ROOT / 'FROZEN_ACQUISITION_C.json'), clipping_data_diagnosis=clipping, lab_exports=labs, summary_sufficiency=read(ROOT / 'raw/SUFFICIENCY.json'), inherited_summary_witnesses=read(ROOT / 'raw/REUSED_SUMMARY_WITNESSES.json'), negative_controls=read(ROOT / 'raw/NEGATIVE_CONTROLS.json'), scorer_controls=read(ROOT / 'raw/SCORER_CONTROLS.json'), memory_equivalence=read(ROOT / 'raw/MEMORY_EQUIVALENCE.json'), precision_correction=read(ROOT / 'PRECISION_CORRECTION.json'), information_boundary=read(ROOT / 'INFORMATION_BOUNDARY.json'), cost=dict(A_resume=read(ROOT / 'raw/RESUME_A_COST.json'), B_generate=read(ROOT / 'raw/GENERATE_B_COST.json'), C_generate=read(ROOT / 'raw/GENERATE_C_COST.json'), acquisition_A_seconds=read(ROOT / 'FROZEN_ACQUISITION_A.json')['seconds'], acquisition_C_seconds=read(ROOT / 'FROZEN_ACQUISITION_C.json')['seconds'], original_interrupted_peak_rss='UNKNOWN; user reports6GB OOM', completed_original_generation_seconds=19.70149417481989, inherited_scan_label_training_discovery_cost='UNKNOWN', records_evaluated_initial=72, successful_full_mesh_scores_initial=sum((r['status'] == 'SCORED' for rr in rounds.values() for r in rr['rows'])), original_failed_attempt_seconds='UNKNOWN beyond19.701494 seconds for two surviving exports', human_questions=0, new_dataset_downloads=0, hardware_metrics_scope='computational measurements, not dental estimands'), edge=dict(from_quantity='retrospective registered gap sheet', to_quantity='outer boundary nodes', resolution='PER_POINT', time_scale='SIMULTANEOUS', downstream_aggregation='region contact, then tooth conjunction'), debt=[dict(quantity='static gap band as physical contact', resolution='PHENOMENOLOGICAL', replacement_measurement='registered contact footprint and force at specified load'), dict(quantity='virtual co-designed cavity as actual fit', resolution='PHENOMENOLOGICAL', replacement_measurement='independent actual preparation boundary and surface'), dict(quantity='source tooth labels and bite pose', resolution='PER_POINT', replacement_measurement='expert segmentation and calibrated registration')], uncertainty='Geometry/pose +/-0.05mm scenarios and conditional set enclosure only; no source or rigorous floating-point enclosure; sample p95 not continuous Hausdorff certificate', data_files={tag: read(ROOT / f'FROZEN_PREDICTIONS_{tag}.json') for tag in rounds})
    dump(ROOT / 'results.json', result)
    lines = ['# Whole crown R3b: measured contact and full anatomy', '', 'The experiment asks whether a measured local gap prescription can guide a closed crown while preserving the full anatomical surface. The external geometric referent is the local [Bite2Text](https://ditto.ing.unimore.it/bite2text/) / [Bits2Bites](https://ditto.ing.unimore.it/bits2bites/) scan geometry. The prescribed contact data comes from those same scans: reproducing it is a construction test, not independent contact prediction.', '', '| Construction | Tooth type | Scored | Anatomy ≤0.35 mm | Function, sampled wall | Function, facet cover | Positive band restored | Median anatomical p95 mm |', '|---|---|---:|---:|---:|---:|---:|---:|']
    for row in table:
        if row['method'] == 'scalar_control':
            continue
        lines.append(f"| {row['round']} | {NAMES[row['family']]} | {row['scored']}/6 | {row['anatomy_pass']}/6 | {row['functional_pass_sampled']}/6 | {row['functional_pass_cover']}/6 | {row['positive_native_restored']}/6 | {row['median_p95_mm']:.3f} |" if row['median_p95_mm'] is not None else f"| {row['round']} | {NAMES[row['family']]} | 0/6 | 0/6 | 0/6 | 0/6 | 0/6 | — |")
    lines += ['', 'A positive native-band pass only means the fixed absolute spatial tolerances passed; it does not require equality of the whole native band or establish loaded occlusion. Each decision is PER_TOOTH; its contact inputs are PER_POINT and contact metrics PER_SURFACE_REGION. Each type contains four frozen test sites and two auxiliary sites in the same six case clusters. These are retrospective descriptive counts, not population estimates. All rejected constructions stay in the denominator. Medians use scored outputs only. The facet-cover column uses the original0.12mm refinement consistently; R2’s later0.05mm supplement is not substituted into this column. Split tables and every gate are in results.json.', '', '![Anatomy and functional gates](figures/summary.png)', '', '![Premolar band regions](figures/premolar_contact.png)', '', 'The equally informed scalar diagnostic used the same acquired gap sheet and the same shell rebuilding. Its nominal functional counts were2/6 molars,0/6 premolars and4/6 anteriors, with0/6 anatomical passes per type. Complete per-tooth outcomes are retained in results.json. No algorithm superiority is claimed.', '', 'R2 premolars failed contact location at all six sites: three native bands were missed, one expanded, one displaced, and the one native-empty case received spurious contact. Its sampled wall, topology, planar margin and two unsigned neighbour patch gates passed. Three of six native premolar gap maps themselves fail the unchanged interference threshold (native-self comparison in raw/NATIVE_SELF_CONTACT.json). That limits copying the measured template; it does not prove the full design problem impossible.', '', 'A changes column heights locally, then rebuilds both surfaces. B retains explicit contact-bearing triangles and builds only the cavity through an SDF. Its vertical axial skirt is a closure, so native axial anatomy remains independently rejectable. C adds unbounded native gap values only at band-crossing faces, using the identical B geometry operator. All methods co-design a virtual preparation: neither establishes fit to an actual prepared tooth.', '', 'What does not hold: no physical force/pressure, clinical suitability, real preparation fit, scanner accuracy, independently verified tooth labels, global self-intersection proof or patient-disjoint generalization. The 0.35 mm anatomy criterion uses the inherited bidirectional area-stratified p95 probes. Margin planarity does not establish a finish line. Wall facet covers are conservative in exact arithmetic; rigorous rounding enclosure is missing. Contact has ±0.05 mm scenarios and conditional set bounds, not a confidence interval.', '', 'The inherited actual triangle-distance witness has p95=0.25mm in both states (identity error0) but a3mm contact-centroid difference; inherited equal planar margins differ2mm in contour distance. These R2 witnesses are reused with source hashes, not new evidence. The new exact-summary witness has identity error 0 for contact area, count, negative area, mean and minimum gap, yet a 3 mm centroid shift and 2 mm² symmetric difference. One location coordinate is the smallest extension distinguishing this constructed pair; this does not prove that every consumer needs the entire field. Wrong height, location, support, thresholds, anatomy, wall, units and facet data are deliberately rejected by executed controls.', '', 'Run `./run_all.sh` to actually regenerate the three selected crowns in isolated temporary outputs, verify exact vertex/face/role parity, recheck source/frozen hashes, recompute all scored meshes with fresh single-tooth processes, exercise controls, export selected specimens and rebuild this figure. The existing local Python/runtime and source datasets are required; clean-machine installation is not tested. The selected-generator replay is recorded in raw/GENERATOR_REPLAY.json; all original frozen outputs, including the two surviving OOM outputs, stay untouched. Full-cohort generation is the original recorded experiment, not claimed as a second execution. Every completed generator/scorer has a /usr/bin/time -v file under raw/.', '', 'Selected STL and 3MF paths, hashes and rejected gates are in LAB_EXPORTS.json and FROZEN_LAB_PREDICTIONS.json. These are research specimens; physical manufacture and measurement are NOT_RUN. Selection was frozen before scores.', '', 'Data licences: Bits2Bites CC BY-NC-SA according to the task brief. Exact licence for the existing Bite2Text copy remains UNKNOWN. Raw unbounded gap arrays existed in the mounted NPZs but were excluded by the executed generator code; whole private target meshes were not mounted. C changes the gap values actually consumed, not the time of physical acquisition. No new dataset or identifying report/photo is used and no files are published. The prior X11 labelling/training cost is inherited and UNKNOWN. Source pages were checked after the own-data framing and first calculation; they confirmed dataset scope and did not change the construction.', '', 'Graph state: PENDING_INDEPENDENT_REVIEW. This delivery does not admit scientific evidence or modify canonical graphs.']
    joint = sum((r.get('anatomy_pass', False) and r.get('functional_pass_sampled', False) for rr in rounds.values() for r in rr['rows']))
    positive = [rounds[t]['summary'][m]['premolar_crown']['positive_native_restored'] for (t, m) in [('A', 'spatial_sheet'), ('B', 'explicit_sheet'), ('C', 'explicit_sheet')]]
    lines.insert(-2, f"C retained unbounded values at {sum((r['additional_nodes'] for r in clipping['rows']))} of {sum((r['finite_native_nodes'] for r in clipping['rows']))} finite measured nodes. Acquisition-level PL-band error exceeded the frozen1e-7mm² tolerance in {clipping['clipped_fail_exact_count']}/18 clipped sheets and {clipping['augmented_fail_exact_count']}/18 augmented sheets. This is a representation check on source-derived information, not independent contact validation.")
    lines[2] = f'Simultaneous anatomy/function passes: {joint} across72 attempted outputs. Premolar nominal function with a positive native band: A {positive[0]}/6, B {positive[1]}/6, C {positive[2]}/6. ' + lines[2]
    (ROOT / 'README_DEMO.md').write_text('\n'.join(lines) + '\n')
    (ROOT / 'RESULTS.md').write_text('\n'.join(lines[:len(table) + 10]) + '\n\nFull limitations, controls and reproduction: README_DEMO.md. Complete per-tooth raw results: rounds/A.json, rounds/B.json and rounds/C.json.\n')
    (fig, axes) = plt.subplots(1, 3, figsize=(13, 4), constrained_layout=True)
    methods = [('A', 'spatial_sheet', 'A: deformation'), ('A', 'scalar_control', 'Scalar diagnostic'), ('B', 'explicit_sheet', 'B: explicit boundary'), ('C', 'explicit_sheet', 'C: crossing information')]
    for (ax, fam) in zip(axes, FAMILIES):
        for (i, (tag, method, label)) in enumerate(methods):
            s = [r for r in rounds[tag]['rows'] if r['family'] == fam and r['participant'] == method]
            for (j, r) in enumerate(s):
                if r['status'] == 'SCORED':
                    ax.scatter(i + (j - 2.5) * 0.025, r['reconstruction_p95_mm'], marker='o' if r['functional_pass_sampled'] else 'x', c='#187f62' if r['positive_native_band_restored'] else '#45689e', s=45)
        ax.axhline(0.35, color='#aa3333', linestyle='--')
        ax.set_xticks(range(4), ['A', 'Scalar', 'B', 'C'])
        ax.set_title(NAMES[fam])
        ax.set_ylabel('Whole-surface p95 (mm)')
        ax.grid(alpha=0.2)
    fig.suptitle('PER_TOOTH; circle = nominal function pass, green = positive native band restored')
    fig.savefig(ROOT / 'figures/summary.png', dpi=170)
    fig.savefig(ROOT / 'figures/summary.svg')
    plt.close(fig)
    return result
if __name__ == '__main__':
    run()
