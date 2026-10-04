from dental_release.paths import expand as _release_expand
import os
for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[k] = '2'
import csv, json, re, hashlib, time
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from freeze import ROOT, OLD, DATA, write, sha, state, now
from planner import load, dist

def fmt(x):
    return 'UNKNOWN' if x is None else f'{x:.3f}'

def metadata_costs():
    costs = []
    for p in sorted((ROOT / 'logs').glob('*.log')):
        t = p.read_text(errors='replace')
        wall = re.search('Elapsed \\(wall clock\\) time[^\\n]*: ([^\\n]+)', t)
        cpu = re.search('User time \\(seconds\\): ([\\d.]+)', t)
        system = re.search('System time \\(seconds\\): ([\\d.]+)', t)
        rss = re.search('Maximum resident set size \\(kbytes\\): (\\d+)', t)
        if cpu or rss:
            costs.append(dict(log=str(p), wall_clock=wall.group(1) if wall else None, user_seconds=float(cpu.group(1)) if cpu else None, system_seconds=float(system.group(1)) if system else None, peak_rss_KiB=int(rss.group(1)) if rss else None, sha256=sha(p)))
    return costs

def run():
    rounds = {r: json.loads((ROOT / f'RESULTS_{r}.json').read_text()) for r in ('R1', 'R2', 'R3', 'R4')}
    mirror = json.loads((ROOT / 'RESULTS_MIRROR_ALL.json').read_text())
    falsify = json.loads((ROOT / 'CONTROL_FALSIFICATION.json').read_text())
    field = json.loads((ROOT / 'FIELD_DEMO.json').read_text())
    allrows = []
    for (rid, r) in rounds.items():
        for c in r['cases']:
            for (name, m) in c['methods'].items():
                allrows.append(dict(round=rid, case=c['case'], HCL=c['HCL'], quality=c['quality'], method=name, primary_p95_mm=m['primary_p95_mm'], length_mm=m.get('length_mm'), arc_length_error_mm=m.get('arc_length_error_mm'), segment_count=m.get('count'), density_probe_change_mm=m.get('density_probe_change_mm'), non_crossing_miters=m.get('non_crossing_miters'), numerically_available=c['numerically_available']))
    with (ROOT / 'table_cases.csv').open('w') as f:
        w = csv.DictWriter(f, fieldnames=list(allrows[0]))
        w.writeheader()
        w.writerows(allrows)
    classes = sorted(set((c['HCL'] for c in rounds['R1']['cases'])))
    byclass = []
    for label in classes:
        row = dict(HCL=label, n=sum((c['HCL'] == label for c in rounds['R1']['cases'])))
        for (rid, r) in rounds.items():
            vals = [c['methods']['dp']['primary_p95_mm'] for c in r['cases'] if c['HCL'] == label and c['numerically_available']]
            row[rid + '_available'] = len(vals)
            row[rid + '_diagnostic_median_p95_mm'] = float(np.median(vals)) if vals else None
        byclass.append(row)
    write(ROOT / 'TABLE_BY_HCL.json', byclass)
    source = json.loads((OLD / 'CASE_MANIFEST.json').read_text())
    coverage = []
    evaluated_ids = {c['case'] for c in rounds['R4']['cases']}
    for cid in sorted(source['records']):
        roles = source['roles'].get(cid, {})
        present = {k: bool(k in roles and roles[k]['exists']) for k in ('Pre', 'Post', 'Original')}
        coverage.append(dict(case=cid, HCL=source['records'][cid]['HCL classification'], **present, benchmark_status='EVALUATED_FROZEN_COMPLETE_TRIPLET' if cid in evaluated_ids else 'NOT_RUN_INCOMPLETE_TRIPLET', missing_roles=','.join((k for (k, v) in present.items() if not v))))
    with (ROOT / 'CASE_COVERAGE_ALL147.csv').open('w') as f:
        writer = csv.DictWriter(f, fieldnames=list(coverage[0]))
        writer.writeheader()
        writer.writerows(coverage)
    demo_plan = json.loads((DATA / 'R4/001/plan.json').read_text())
    write(ROOT / 'demo_plan.json', demo_plan['grafts']['dp'])
    export_dir = ROOT / 'demo_export'
    export_dir.mkdir(exist_ok=True)
    for p in (DATA / 'R4/001').glob('dp_segment_*.stl'):
        (export_dir / p.name).write_bytes(p.read_bytes())
    maps = {r: {c['case']: c for c in rr['cases']} for (r, rr) in rounds.items()}
    shared = [c for c in maps['R1'] if all((maps[r][c]['numerically_available'] for r in maps))]
    (fig, ax) = plt.subplots(1, 3, figsize=(16, 5), constrained_layout=True)
    for (rid, color) in [('R1', '#6b7280'), ('R2', '#087f8c'), ('R3', '#e07a27'), ('R4', '#7337a5')]:
        vals = sorted((maps[rid][c]['methods']['dp']['primary_p95_mm'] for c in shared))
        label = rid + ' legacy INVALID' if rid in ('R1', 'R2') else rid
        ax[0].plot(vals, np.arange(1, len(vals) + 1) / max(len(vals), 1), label=label, color=color, lw=2)
    mirror_map = {c['case']: c for c in mirror['cases']}
    mv = sorted((mirror_map[c]['primary_p95_mm'] for c in shared))
    ax[0].plot(mv, np.arange(1, len(mv) + 1) / len(mv), label='Raw Pre mirror proxy', color='#308537', lw=2, ls=':')
    ax[0].axvline(3, color='black', ls='--', lw=1, label='Frozen 3 mm gate')
    ax[0].set(xlabel='Defect-region sampled surface p95 (mm)', ylabel='Paired cumulative fraction', title=f'External virtual Post; {len(shared)}/118 paired cases')
    ax[0].legend(fontsize=8)
    ax[0].grid(alpha=0.2)
    cid = '001'
    raw = np.load(DATA / 'R4' / cid / 'evaluation_raw.npz')
    p = raw['post_aligned']
    (pre, _, _) = load(cid, 'Pre')
    plan = json.loads((DATA / 'R4' / cid / 'plan.json').read_text())
    nodes = np.array(plan['grafts']['dp']['nodes_mm'])
    ax[1].scatter(p[::25, 0], p[::25, 1], s=2, alpha=0.25, color='#087f8c', label='expert virtual Post')
    ax[1].scatter(pre[::25, 0], pre[::25, 1], s=2, alpha=0.35, color='gray', label='retained Pre proxy')
    ax[1].plot(nodes[:, 0], nodes[:, 1], 'o-', color='#7337a5', lw=3, label='R4 Pre-only straight graft')
    ax[1].set(xlabel='CT x (assumed mm)', ylabel='CT y (assumed mm)', title='Fixed illustrative case 001: planar view')
    ax[1].axis('equal')
    ax[1].legend(fontsize=8)
    radii = []
    violations = []
    contained = []
    for rr in json.loads((ROOT / 'FROZEN_PREDICTIONS_R3.json').read_text())['predictions']:
        pp = json.loads(Path(rr['path']).read_text())
        if pp['status'] == 'PREDICTED':
            radii.append(pp['selected_radius_mm'] * 2)
            violations.append(max((e['containment_upper_bound_mm'] for e in pp['envelopes'])))
            contained.append(pp['all_envelopes_contained'])
    ax[2].scatter(radii, violations, c=['#087f8c' if b else '#d14b4b' for b in contained], alpha=0.6, s=20)
    ax[2].axhline(0, color='black', ls='--', lw=1)
    ax[2].set(xlabel='Selected donor surrogate diameter (mm)', ylabel='Implant-envelope violation upper bound (mm)', title='R3 scenario plane; occlusion not measured')
    ax[2].grid(alpha=0.2)
    fig.savefig(ROOT / 'FIGURE.png', dpi=170)
    fig.savefig(ROOT / 'FIGURE.pdf')
    plt.close(fig)
    field_path = next((e['path'] for e in field['exports'] if Path(e['path']).name == 'fields.npz'))
    with np.load(field_path) as z:
        ph = z['union']
        prephi = z['pre']
        k = ph.shape[2] // 2
        (fig, ax) = plt.subplots(figsize=(6, 5))
        ax.imshow(ph[:, :, k].T, origin='lower', cmap='coolwarm', vmin=-8, vmax=8)
        ax.contour(prephi[:, :, k].T, levels=[0], colors='black', linewidths=0.8)
        ax.contour(ph[:, :, k].T, levels=[0], colors='#e07a27', linewidths=1)
        ax.set(xlabel='grid x index', ylabel='grid y index', title='Native Pre SDF + analytic graft CSG; research export')
        fig.savefig(ROOT / 'FIELD_FIGURE.png', dpi=170)
        plt.close(fig)
    table = '| Construction | Available / 118 | Diagnostic median p95, mm | Complete-cohort median | Frozen outcome |\n|---|---:|---:|---:|---|\n'
    labels = {'R1': 'Pre mirror, INVALID legacy mesh', 'R2': 'Expert rail, INVALID legacy mesh', 'R3': 'Corrected implant radius/height scenario', 'R4': 'Pre-only 3D centroid rail, corrected mesh'}
    for (rid, r) in rounds.items():
        m = r['methods']['dp']
        table += f"| {labels[rid]} | {r['n_available']} | {fmt(m['diagnostic_median_p95_mm'])} | {fmt(m['complete_cohort_median_p95_mm'])} | {r['outcome']} |\n"
    shape = rounds['R3']['methods']['shape_only']
    table += f"| R3 corrected shape-only graft control | {shape['n']} | {fmt(shape['diagnostic_median_p95_mm'])} | {fmt(shape['complete_cohort_median_p95_mm'])} | matched geometric control |\n"
    table += f"| Raw Pre mirror (surface, all-case baseline) | {mirror['n_available']} | {fmt(mirror['diagnostic_median_p95_mm'])} | {fmt(mirror['complete_cohort_median_p95_mm'])} | comparison only |\n"
    competitor = rounds['R2']['methods']['osteoopt']
    r2 = rounds['R2']
    r3 = rounds['R3']
    demo = _release_expand(f"# Straight graft planning against a virtual mandibular reference\n\nThe idea is to turn a mandibular target shape into 1–3 straight circular graft surrogates with lengths, bend angles, shared osteotomy planes and STL exports. A second operation places finite graft thickness inside segment search; a third searches donor radius and height for a future implant envelope. These are executable geometric research capabilities. The frozen 3 mm surface-accuracy target has not been achieved.\n\nRun from this directory with `./run_all.sh`. Python requires NumPy, SciPy, scikit-image and Matplotlib; the comparison requires Java/JDK. No downloads are made by this command. The installed local dataset and predecessor component caches are read with hashes. All arrays and STLs are on `@DENTAL_WORK_ROOT@/X4_fibula_planner/`. The same command verifies/reuses frozen predictions and completed measurements; move to a new lane/prereg for a different physical input.\n\n`demo_plan.json` and `demo_export/` provide the fixed case001 Pre-only plan and its small STL segments. Reported lengths are centreline lengths; cut normals are oriented along the chain. For segment i, retain `dot(x-node_i,n_i)>=0` and `dot(x-node_(i+1),n_(i+1))<=0`. A donor blank additionally needs the miter extensions `radius*tan(bend/2)` at its ends. These are geometric surrogates, not donor-specific cutting guides. Numerical and figure environments are pinned in `ENVIRONMENT.json`; set `X4_RENDER_PYTHON` to a NumPy-compatible Matplotlib Python when running elsewhere.\n\n{table}\nDiagnostic medians include available cases only and are not complete-cohort successes. Missing cases fail the frozen completeness gate. Every row is in `table_cases.csv`; HCL strata and missing counts are in `TABLE_BY_HCL.json`.\n\n![External comparison and geometric plans](FIGURE.png)\n\nThe facit is [Figshare v2, DOI 10.6084/m9.figshare.28052240.v2](https://doi.org/10.6084/m9.figshare.28052240.v2), interpreted with the [data descriptor](https://doi.org/10.1038/s41597-025-06048-8). Post is expert virtual completion. Original is defective anatomy, useful only for preparation/preservation comparisons. There is no intact same-patient anatomy or postoperative outcome. The local snapshot has 410 STL and 118 complete triplets out of 147 designated cases. Dataset license: **CC BY 4.0**, per local `_metadata.json` and pinned source metadata. Attribution: Wu, Jiang, Shao and colleagues / Mandibular Defect Dataset; Figshare snapshot credited to its deposited creators.\n\nR1 and R4 predict from Pre alone. R4 follows 3D centroid stations of the mirrored missing surface instead of the lower border, and exports corrected straight cylinders. R2 explicitly obtains an inferior contour from the target expert Post: it is assisted CAD planning and cannot support a blind-reconstruction claim. The full added-region surface, including the acquired lower band, is then used for evaluation. R3 uses a superior-anterior surface percentile as a plane scenario, not a labelled occlusal measurement. Condyle landmarks and actual occlusal errors remain UNKNOWN.\n\nAll four rounds use the same118-case cohort. R4 is an adapted repeat-cohort experiment after earlier failures, not independent confirmation on untouched cases; its predictor still opens only Pre. The frozen gates and missing cases are retained without threshold adjustment.\n\nThe equally informed exhaustive chord solver matches the candidate at numerical precision. The published [OsteoOpt++ geometric simplifier](https://github.com/hamidreza-aftabi/OsteoOpt) is actually executed as its pinned, unmodified Java class with a minimal Euclidean-vector adapter. Its R2 diagnostic surface p95 is {fmt(competitor['diagnostic_median_p95_mm'])} mm on {competitor['n']} cases. This is only its geometric simplifier, not the full union/chewing/Bayesian optimizer. That system was not run: MATLAB, ArtiSynth runtime, matched donor anatomy and case-specific mechanics are unavailable. Upstream code license: **PolyForm Noncommercial 1.0.0**; LICENSE and NOTICE are saved in `sources/`. No superiority over the full system is claimed.\n\nWhat fails: R1/R2 also contain a preserved mesh construction bug: independent ring phases twisted side faces, so their surface scores describe the failed legacy mesh rather than a valid straight cylinder. R3 repairs this with axis transport of common miter rings; `GEOMETRY_STRAIGHTNESS_AUDIT.json` and falsifying controls check it. The same-case straightness debugging before the R3 freeze is disclosed as exploratory in `PROTOCOL_DEVIATIONS.json`; it is not prospective confirmation. Shape completion is not the same as a feasible straight donor plan. Circular donor radius is a scenario, not a measured fibula; masks and CT directions are anatomical proxies; sampled distances are not exact surface Hausdorff distances. Primary distances use retained Pre plus graft surface clouds, not an exact filtered union boundary. Poor rails, cut constraints, and very large errors remain in the table. R3 contains all segment implant envelopes in {r3['n_all_envelopes_contained']}/118 cases only under the declared geometry and plane assumptions. This scenario asks for an implant at each segment midpoint; it does not distinguish ramus from implant-bearing body. No cortical support, nerve clearance, vascular pedicle, fixation, union propensity, load or patient suitability is validated.\n\n`FIELD_DEMO.json` and `FIELD_FIGURE.png` show a real Pre mesh → native signed-distance grid → min-CSG graft addition → `Pre_plus_graft_surrogate.stl`. The graft cut field has exact sign but is not Euclidean SDF at edges; the Pre SDF is a finite-grid EDT approximation. Pre interior voxels lost: {field['pre_inside_voxels_lost']}. Manufacture a polymer research model only after reviewing the exported proxy support and seams.\n\nAll predictions are hashed and timestamped before external surface comparison. `CONTROL_FALSIFICATION.json` shows that translated reference geometry, wrong plane/length/objective/radius, a folded short graft and loss of a Pre voxel are rejected. These controls test sensitivity and code contracts; they do not provide physical validation.\n\nThe next decisive inputs are independently labelled recipient rims and endpoint roles, a prosthetic plane with intended implant-bearing body sites, and a measured donor cross-section. `NEXT_CONSTRUCTION_R5.json` puts recipient contact and the full surface objective inside donor-profile sweep search. `implant_variant.py --plane-json` supports a per-case three-point plane format, but frozen runs require a new prereg to incorporate new measurements. See `LAB_MEASUREMENT_SPEC.json` for a prototype comparison with predictions frozen before manufacture/measurement. Results remain PENDING_INDEPENDENT_REVIEW.\n")
    (ROOT / 'README_DEMO.md').write_text(demo)
    costs = metadata_costs()
    write(ROOT / 'COST_ACCOUNTING.json', dict(own_commands=costs, predecessor_prefix=dict(locator=str(OLD / 'COST_ACCOUNTING.json'), sha256=sha(OLD / 'COST_ACCOUNTING.json'), known_preparation_seconds=dict(R1_surface_cache=148.63, R2_components_fields=364.94), charge='inherited complete preparation retained; these include other tasks and are not marginal current costs'), discovery='implementation, source inspection, all failed/incomplete commands included in COMMANDS; model token cost not measured', physical_acquisition='NOT_PERFORMED; no invented zero cost', fallback='UNKNOWN output retained; no human correction', worker_cpu='Log list includes inclusive pipeline-parent CPU and child-stage CPU: DO NOT sum all entries. Use disjoint top-level logs only; per-case rusage cumulative is not added.', disjoint_top_level_logs=['run_all_R1.log', 'R1_evaluate_repaired.log', 'R1_evaluate_repaired_v2.log', 'R1_shard0.log', 'R1_shard1.log', 'R1_aggregate.log', 'R2_predict.log', 'mirror_all.log', 'R2_shard0.log', 'R2_shard1.log', 'run_all_complete.log', 'R4_predict.log', 'R4_shard0.log', 'R4_shard1.log', 'run_all_final.log', 'run_all_verified.log'], resource_limit='two numerical shards, BLAS/OMP each2, total configured<=4; no GPU; own intermediates<=3GB'))
    outcome = dict(status='FOUR_CONSTRUCTIONS_DECIDED', updated_utc=now(), rounds={rid: {k: v for (k, v) in r.items() if k != 'cases'} for (rid, r) in rounds.items()}, raw_mirror_baseline={k: v for (k, v) in mirror.items() if k != 'cases'}, external_referent=dict(kind='published_dataset', locator='https://doi.org/10.6084/m9.figshare.28052240.v2', compared_quantity='expert virtual Post defect-region surface; Original defective preparation/preservation only', refutes_us=True), controls_passed=falsify['passed'], field_export=field, scientific_admission=False, review_state='PENDING_INDEPENDENT_REVIEW', cost_accounting=dict(path=str(ROOT / 'COST_ACCOUNTING.json'), sha256=sha(ROOT / 'COST_ACCOUNTING.json')), reproducible_command='./run_all.sh', figure=dict(path=str(ROOT / 'FIGURE.png'), sha256=sha(ROOT / 'FIGURE.png')), limitations=['no clinical/surgical outcome', 'no donor CT', 'unreviewed component proxies', 'no labelled occlusal/condyle landmarks', 'no measured mechanics'], next_construction='freeze marked target curve, actual donor cross-section and three-point plane, then manufacture polymer plan and test acquired geometry against immutable predictions')
    raw_files = [p for p in DATA.rglob('*') if p.is_file()]
    size = sum((p.stat().st_size for p in raw_files))
    manifest = dict(bytes=size, limit_bytes=3000000000, files=[dict(path=str(p), bytes=p.stat().st_size, sha256=sha(p)) for p in raw_files])
    write(ROOT / 'RAW_MANIFEST.json', manifest)
    outcome['raw_manifest'] = dict(path=str(ROOT / 'RAW_MANIFEST.json'), sha256=sha(ROOT / 'RAW_MANIFEST.json'), bytes=size)
    if size > 3000000000:
        raise RuntimeError('3GB intermediate budget exceeded')
    write(ROOT / 'results.json', outcome)
    (ROOT / 'RESULTS.md').write_text('# X4 — four geometric constructions, external accuracy not achieved\n\n' + table + '\nAll exact gates remain in RESULTS_R1/R2/R3/R4.json; no failed threshold was changed. Clinical validity and anatomical support remain UNKNOWN.\n\nActual changes: Pre-only reflection rail -> supplied expert inferior rail with finite-thickness cuts inside search -> corrected straight cylinders and conditional implant radius/height search -> Pre-only 3D centroid rail. Strong equal-information geometry and scenario-selection controls tie; no algorithmic novelty or 10x superiority is claimed. The full expert surface refutes the declared accuracy target.\n\nR1/R2 mesh surfaces are invalid straight-cylinder representations, preserved with failed gates. R3/R4 use correct cross-section transport. Every primary sample-distance vector is stored under the lane data directory and hashed. FROZEN_PREDICTIONS*.json precede evaluation. table_cases.csv contains lengths, arc deviations, counts, cut-angle information and sampling probes. Condyle and occlusal landmark errors are UNKNOWN; Original is defective.\n\nPublished OsteoOpt geometric class ran; full union optimizer NOT_RUN with enumerated prerequisites in COMPETITOR_BUILD.json. The descriptor/local Figshare release distinction is preserved. 29 designated cases lack complete triplets; no download or release merge was made.\n\nNative field mesh/SDF/CSG/export and falsifying controls ran. This is a research geometric pipeline, not clinical reconstruction or validated manufacture. COST_ACCOUNTING.json includes failed commands and inherited preparation; disk inventory is RAW_MANIFEST.json.\n\nAll147 source records are listed in CASE_COVERAGE_ALL147.csv; 29 lack the preregistered complete triplet, including8 Pre/Post pairs without Original. Next changes information and physical representation: independently labelled recipient rims/condylar roles, target/occlusal curve and donor cross-section; put contact and full surface objective inside the sweep search, then validate a polymer assembly before mechanics.\n')
    feedback = dict(review_state='PENDING_INDEPENDENT_REVIEW', scientific_admission=False, target_id=None, outcome='FAIL_FROZEN_GATES; RUNNABLE_GEOMETRIC_PLANNING_CAPABILITY; SOURCE AND LANDMARK LIMITATIONS', result_file='results/LANE_X4_FIBULA_PLANNER/results.json', sha256=sha(ROOT / 'results.json'), measured_quantity='sampled expert-virtual-surface p95; graft lengths/angles; conditional envelope containment; CSG preservation', units='mm, degrees, counts', uncertainty='component anatomy, CT metric calibration, donor dimensions and occlusal/condyle landmarks UNKNOWN; density probes per case', population_regime='118 complete triplets /147 designated; virtual CAD annotations, circular-donor research scenarios', preregistered_gate='PREREG_R1/R2/R3/R4.json immutable hashes', baseline='Pre mirror, exhaustive equal-info solver, pinned published OsteoOpt geometric class', negative_result=True, graph_binding_status='MISSING_WORKING_COVERAGE_PROPOSAL; no suitable existing fibula planning goal; identity nodes are not used to carry planning claims', proposed_coverage=dict(proposed_id='DENT-RECON-FIBULA-SEGMENT-PLAN', not_registered=True, input_ports=['Pre geometry + independently reviewed jaw support', 'labelled reconstruction/occlusal curve', 'donor cross-section'], output_ports=['1-3 straight segment endpoint/length/cut-plane plan', 'future-implant envelope and uncertainty', 'research manufacture export'], consumers=['CAD/CAM prototype', 'geometric metrology', 'future fixation/mechanics after physical calibration']))
    write(ROOT / 'GRAPH_FEEDBACK.json', feedback)
    handoff = f'# X4 handoff — parent capability open\n\nRun `./run_all.sh`; start with README_DEMO.md and CURRENT_WORK_STATE.json. Four constructions decided; no claim admitted.\n\n{table}\nR1/R4 are Pre-only; R2/R3 explicitly use expert Post geometry as supplied planning information. R1/R2 phased mesh sides are invalid; R3/R4 are corrected. Post is virtual, Original defective. These are not surgical outcomes or intact-anatomy comparisons. Keep every UNKNOWN and all four frozen preregs. Strong exhaustive/scenario-selection controls TIE. Actual published OsteoOpt geometric simplifier ran; full BO not available, prerequisites in COMPETITOR_BUILD.json.\n\nThe mesh/SDF/CSG/export example and all injected-fault controls ran. Figure and table are in README_DEMO. Full output hashes, raw data and complete costs are in results.json/RAW_MANIFEST.json/COST_ACCOUNTING.json. Source cache reads are immutable; all own large arrays stay on games disk ({size / 1000000000.0:.3f} GB). Never overwrite frozen runs with a newly acquired plane.\n\nBinding obstacle: inferior/centroid surface proxies and circular donor shape do not recover a physically constrained graft/implant plan. Recipient contact is not enforced, and the chord objective is a surrogate for the full surface objective. Missing labelled recipient rims, jaw/condyle endpoint roles, occlusion and donor cortex/pedicle dominate. NEXT_CONSTRUCTION_R5.json specifies a donor-profile sweep with signed Pre volume and recipient cap contacts inside search; LAB_MEASUREMENT_SPEC.json describes independent inputs and a frozen polymer-model comparison. Physical measurements have not occurred.\n\nGraph coverage is absent: proposed DENT-RECON-FIBULA-SEGMENT-PLAN only in GRAPH_FEEDBACK.json, not a registered node. Existing dataset/identity nodes have no suitable planning goal path; no broad planning claim is dispatched to them. Coordinator should review native goal/ports before graph admission. PENDING_INDEPENDENT_REVIEW.\n'
    (ROOT / 'HANDOFF.md').write_text(handoff)
    (ROOT / 'HANDOFF_R4.md').write_text(handoff)
    state('four_constructions_decided', 'independent anatomy review and marked-curve/donor measurement construction; see LAB_MEASUREMENT_SPEC.json', '; '.join((r + ':' + v['outcome'] for (r, v) in rounds.items())))
    final = _release_expand("Delivered: `./run_all.sh` , planner for 1 – 3 straight graph segment, table, figure and STL export.\nFour designs were tested against 118 complete cases of dataset 147 .\nSenaste Pre-baserade planeraren: median yt-p95 24,82 mm; reflectionsproxyn: 13,96 mm.\nAll overall frozen gates failed ; equally informed controls gave TIE .\nSegment geometry corrected; previously incorrect mesh and result are preserved. The implant scenario only passed the 3/118 case.\nThe entire execution and raw data verification passed; reference is virtual expert completion, not operation outcome.\nThe biggest obstacle is unmarked contact edges and missing measured donorfibula/protetical plane.\nNext design inserts these measurements and contact conditions before frozen polymer prototype and independent measurement.\nTable, figure and limitations are found in [README_DEMO.md ](@DENTAL_IMPLEMENTATIONS@/X4/README_DEMO.md); independent review remains.\n")
    (ROOT / 'SOL_NIGHT_FINAL_X4-fibula-planner.md').write_text(final)
    print(table)
    print('RAW_BYTES', size, 'CONTROL_FALSIFICATION', falsify['passed'])
if __name__ == '__main__':
    run()
