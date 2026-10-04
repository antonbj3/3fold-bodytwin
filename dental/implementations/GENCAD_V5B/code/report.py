from common import *
import collections, csv, gzip, itertools
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

def pair_analysis(rows):
    valid = [r for r in rows if r.get('geometry_pass') and r.get('force', {}).get('source_force05', {}).get('force05_N') is not None]
    out = []
    for (a, b) in itertools.combinations(valid, 2):
        if a['key'] != b['key'] or a['participant'] == b['participant']:
            continue
        (fa, fb) = (r['force']['source_force05']['force05_N'] for r in [a, b])
        (ia, ib) = (r['force']['source_force05']['parameter_rectangle_N'] for r in [a, b])
        da = a['adapted_shape_p95_mm'] - b['adapted_shape_p95_mm']
        do = a['original_shape_p95_mm'] - b['original_shape_p95_mm']
        df = fa - fb
        out.append(dict(key=a['key'], a=a['participant'], b=b['participant'], a_uid=a['uid'], b_uid=b['uid'], adapted_shape_difference_mm=da, original_shape_difference_mm=do, conditional_force05_difference_N=df, eligible_adapted=abs(da) > 0.01 and abs(df) > 1, eligible_original=abs(do) > 0.01 and abs(df) > 1, adapted_reversal=abs(da) > 0.01 and abs(df) > 1 and (da * df > 0), original_reversal=abs(do) > 0.01 and abs(df) > 1 and (do * df > 0), disjoint_parameter_rectangles=ia[0] > ib[1] or ib[0] > ia[1], physical_robust_ranking='UNKNOWN', resolution='PER_TOOTH'))
    method = []
    names = sorted({r['participant'] for r in valid})
    for (a, b) in itertools.combinations(names, 2):
        qa = {r['key']: r for r in valid if r['participant'] == a}
        qb = {r['key']: r for r in valid if r['participant'] == b}
        keys = sorted(set(qa) & set(qb))
        if not keys:
            continue
        da = np.median([qa[k]['adapted_shape_p95_mm'] - qb[k]['adapted_shape_p95_mm'] for k in keys])
        df = np.median([qa[k]['force']['source_force05']['force05_N'] - qb[k]['force']['source_force05']['force05_N'] for k in keys])
        method.append(dict(a=a, b=b, paired_sites=len(keys), case_clusters=len({k.rsplit('_', 2)[0] for k in keys}), median_shape_difference_mm=float(da), median_conditional_force05_difference_N=float(df), rank_reversal=bool(abs(da) > 0.01 and abs(df) > 1 and (da * df > 0)), scope='Only this pair common accepted in-source support; not global rank', resolution='POPULATION'))
    return dict(rows=out, participant_pairs=method, in_source_geometry_accepted=len(valid), eligible_adapted_pairs=sum((r['eligible_adapted'] for r in out)), adapted_reversals=sum((r['adapted_reversal'] for r in out)), eligible_original_pairs=sum((r['eligible_original'] for r in out)), original_reversals=sum((r['original_reversal'] for r in out)), disjoint_rectangle_reversals=sum((r['adapted_reversal'] and r['disjoint_parameter_rectangles'] for r in out)), physical_rank_change='UNKNOWN', interpretation='Conditional source curve can reorder numerical candidates; exact-summary test refutes median thickness sufficiency for physical strength')

def summarize(rows):
    summary = []
    for name in sorted({r['participant'] for r in rows}):
        rr = [r for r in rows if r['participant'] == name]
        ok = [r for r in rr if r['status'] == 'ADAPTED']
        gg = [r for r in ok if r['geometry_pass']]
        qq = [r for r in gg if r['force'].get('source_force05', {}).get('force05_N') is not None]
        hh = [r for r in ok if r['cement']['status'] == 'CONDITIONAL_REYNOLDS_CONDUCTANCE']
        q50 = [r['force']['source_force05']['force05_N'] for r in qq]
        summary.append(dict(participant=name, requested=len(rr), adapted=len(ok), geometry_accepted=len(gg), force05_in_domain_accepted=len(qq), hydraulic_fields=len(hh), median_shape_before_mm=float(np.median([r['original_shape_p95_mm'] for r in ok])) if ok else None, median_shape_after_mm=float(np.median([r['adapted_shape_p95_mm'] for r in ok])) if ok else None, median_gap_error_p95_um=float(np.median([r['gap_error_p95_um'] for r in ok])) if ok else None, median_conditional_force05_N=float(np.median(q50)) if q50 else None, median_conductance_ratio=float(np.median([r['cement']['spatial_over_mean_ratio'] for r in hh])) if hh else None, force05_calibrated=0, film_physically_measured=0, resolution='POPULATION', force_cohort='Only accepted in-domain subset; participant medians have unequal support and are not global rank'))
    return summary

def figs(rows, summary, rounds, pairs):
    (fig, ax) = plt.subplots(1, 3, figsize=(17, 5), layout='constrained')
    counts = [sum((r['status'] == 'ADAPTED' for r in rows)), sum((r.get('geometry_pass', False) for r in rows)), sum((r.get('cement', {}).get('status') == 'CONDITIONAL_REYNOLDS_CONDUCTANCE' for r in rows)), sum((r.get('force', {}).get('source_force05', {}).get('force05_N') is not None for r in rows)), 0]
    labels = ['STL export', 'Geometry probes', 'Nominal-film solve', 'Source force05', 'Physical calibration']
    ax[0].barh(labels[::-1], counts[::-1], color=['#bb5555', '#517fa4', '#517fa4', '#157f6e', '#157f6e'])
    ax[0].set(xlabel='Complete-crown designs', title='A. Coverage of the183 submitted crowns', xlim=(0, 205))
    for (y, c) in enumerate(counts[::-1]):
        ax[0].text(c + 2, y, str(c), va='center')
    ax[0].text(0.01, 0.02, 'Different query subsets; six original generation failures', transform=ax[0].transAxes, fontsize=7)
    colors = plt.get_cmap('tab20')
    names = sorted({r['participant'] for r in rows})
    legend = []
    for (i, name) in enumerate(names):
        good = [r for r in rows if r['participant'] == name and r.get('geometry_pass') and (r.get('force', {}).get('source_force05', {}).get('force05_N') is not None)]
        if not good:
            continue
        x = np.array([r['adapted_shape_p95_mm'] for r in good])
        y = np.array([r['force']['source_force05']['force05_N'] for r in good])
        band = np.array([r['force']['source_force05']['parameter_rectangle_N'] for r in good])
        ax[1].errorbar(x, y, yerr=np.vstack([y - band[:, 0], band[:, 1] - y]), fmt='o', ms=4, lw=0.6, color=colors(i), label=name)
    ax[1].set(xlabel='Adapted reconstruction p95 (mm)', ylabel='Source-conditioned force05 (N)', title='B. Source-force scenarios')
    ax[1].legend(fontsize=6, loc='best')
    ax[1].text(0.02, 0.97, 'Parameter rectangles only; physical transfer UNKNOWN', transform=ax[1].transAxes, va='top', fontsize=8)
    hs = [r for r in rows if r.get('cement', {}).get('status') == 'CONDITIONAL_REYNOLDS_CONDUCTANCE']
    ax[2].scatter([r['cement']['area_weighted_gap_mean_um'] for r in hs], [r['cement']['spatial_over_mean_ratio'] for r in hs], s=10, alpha=0.5, color='#175c91')
    ax[2].axhline(1, c='black', ls='--', lw=0.8)
    ax[2].set(xlabel='Area-mean nominal gap (µm)', ylabel='Spatial / uniform-mean conductance', title='C. Actual digital film topology matters')
    ax[2].text(0.02, 0.97, '100 Pa s scenario; no observed seated film', transform=ax[2].transAxes, va='top', fontsize=8)
    fig.suptitle('GenCAD v5b — geometry connected to conditional physical queries', fontsize=14)
    fig.savefig(P / 'figures/physical_axes.png', dpi=180)
    fig.savefig(P / 'figures/physical_axes.pdf')
    plt.close(fig)
    r = next((r for r in rows if r.get('geometry_pass')))
    rec = next((x for x in cohort() if x['uid'] == r['uid']))
    old = load_npz(rec['mesh_path'])
    new = load_npz(r['mesh_path'])
    fld = load_npz(r['field_path'])
    fig = plt.figure(figsize=(13, 5))
    aa = fig.add_subplot(121, projection='3d')
    bb = fig.add_subplot(122, projection='3d')
    for (axis, z, title) in [(aa, old, 'Original participant shell'), (bb, new, 'Adapted shell and nominal gap probes')]:
        ff = z['faces']
        tri = z['vertices'][ff[::max(1, len(ff) // 6000)]]
        coll = Poly3DCollection(tri, facecolor='#b7c7cd', edgecolor='none', alpha=0.18)
        axis.add_collection3d(coll)
        v = z['vertices']
        axis.set(xlim=(v[:, 0].min(), v[:, 0].max()), ylim=(v[:, 1].min(), v[:, 1].max()), zlim=(v[:, 2].min(), v[:, 2].max()), xlabel='x (mm)', ylabel='y (mm)', zlabel='z (mm)', title=title)
        axis.set_box_aspect(np.ptp(v, axis=0))
        axis.view_init(25, -65)
    pts = fld['sites_mm']
    sc = bb.scatter(*pts.T, c=fld['gap_mm'] * 1000, s=5, cmap='viridis', vmin=0, vmax=100)
    fig.colorbar(sc, ax=bb, label='Nominal digital gap (µm)', shrink=0.65)
    fig.suptitle(r['participant'] + ' / ' + r['key'] + ' — virtual preparation')
    fig.savefig(P / 'figures/crown_field.png', dpi=180, bbox_inches='tight')
    plt.close(fig)
    return r

def run(rows, tests, accounting, artifact):
    rounds = {}
    for k in range(1, 8):
        rr = read(P / f'raw/R{k}_GEOMETRY.json')['rows']
        ok = [r for r in rr if r['status'] == 'ADAPTED']
        rounds[f'R{k}'] = dict(requested=len(rr), adapted=len(ok), geometry_accepted=sum((r['geometry_pass'] for r in ok)), failed_gates=dict(collections.Counter((g for r in ok for (g, v) in r['gates'].items() if not v))), not_run=sum((r['status'] == 'NOT_RUN_AFTER_PILOT_FAILURE' for r in rr)), adaptation_failures=[dict(uid=r['uid'], reason=r.get('reason')) for r in rr if r['status'] == 'ADAPTATION_FAILED'], resolution='POPULATION')
    summary = summarize(rows)
    pairs = pair_analysis(rows)
    dump(P / 'raw/RANK_COMPARISONS.json', pairs)
    dump(P / 'raw/LEADERBOARD.json', summary)
    source = read(P / 'raw/SOURCE_VALIDATION.json')
    suff = read(P / 'raw/SUFFICIENCY.json')
    ok = [r for r in rows if r['status'] == 'ADAPTED']
    hs = [r for r in ok if r['cement']['status'] == 'CONDITIONAL_REYNOLDS_CONDUCTANCE']
    qs = [r for r in ok if r['force'].get('source_force05', {}).get('force05_N') is not None]
    inherited = json.loads(gzip.decompress((V5 / 'raw/PHYSICAL_ROWS.json.gz').read_bytes()))
    by = {(r.get('track'), r['participant'], r['key']): r for r in rows}
    ports = []
    matched = 0
    for r in inherited:
        target = by.get((r['track'], r['participant'], r['task_id']))
        if target:
            ports.append(dict(r, v5b_uid=target['uid'], v5b_status=target['status'], v5b_geometry_accepted=target.get('geometry_pass'), v5b_physical_record='raw/PHYSICAL_ROWS.json', force05_N=None, cementfilm_um=None))
            matched += 1
        else:
            ports.append(dict(r, v5b_status='NOT_COMPLETE_CROWN_PORT', v5b_reason='Roof/partial design lacks this benchmark complete-shell preparation/intaglio mapping; no physical value fabricated'))
    (P / 'raw/ALL_V5_PORTS.json.gz').write_bytes(gzip.compress(json.dumps(clean(ports), separators=(',', ':')).encode(), mtime=0))
    example = figs(rows, summary, rounds, pairs)
    import trimesh
    from geometry import preparation
    rec = next((x for x in cohort() if x['uid'] == example['uid']))
    (p, pm) = preparation(rec)
    mesh = load_npz(example['mesh_path'])
    trimesh.Trimesh(mesh['vertices'] @ p['source_R'].T + p['source_base'], mesh['faces'], process=False).export(P / 'exports/example_crown.stl')
    trimesh.Trimesh(pm.vertices @ p['source_R'].T + p['source_base'], pm.faces, process=False).export(P / 'exports/example_preparation.stl')
    dump(P / 'exports/EXAMPLE.json', dict(uid=example['uid'], units='mm', frame='original registered source frame', virtual_preparation=True, not_manufacturing_validated=True, crown_sha256=sha(P / 'exports/example_crown.stl'), preparation_sha256=sha(P / 'exports/example_preparation.stl')))
    result = dict(schema='GenCAD-v5b-capability-v1', claim_type='capability', review_state='PENDING_INDEPENDENT_REVIEW', scientific_admission=False, outcome='EXPLICIT_MARGIN_AND_SPATIAL_PHYSICS_PORTS; ABSOLUTE_PHYSICAL_CALIBRATION_UNRESOLVED', rounds=rounds, complete_panel=dict(requested=len(rows), source_scored=183, adapted=len(ok), geometry_accepted=sum((r['geometry_pass'] for r in ok)), case_clusters=len({r['key'].split('_')[0] for r in rows}), sites=len({r['key'] for r in rows}), participants=len(summary), resolution='POPULATION'), physical=dict(force05_source_curve_in_domain=len(qs), force05_source_curve_in_domain_and_geometry_accepted=pairs['in_source_geometry_accepted'], X1B_mean_queries=sum(('X1B_mean_force_diagnostic' in r['force'] for r in ok)), X1B_strict_setup_thickness_supported=sum((r['force'].get('X1B_thickness_in_support', False) for r in ok)), calibrated_generated_force05=0, measured_seated_film=0, spatial_hydraulic_queries=len(hs), hydraulic_mass_max=max((r['cement']['mass_relative_error'] for r in hs)) if hs else None, median_conductance_spatial_over_mean=float(np.median([r['cement']['spatial_over_mean_ratio'] for r in hs])) if hs else None, resolution='POPULATION summaries of PER_POINT/PER_SURFACE_REGION fields and PER_TOOTH diagnostics'), ranking={k: v for (k, v) in pairs.items() if k not in ['rows', 'participant_pairs']}, leaderboard=summary, source_validation=dict(gap_pass=sum((r['gate'] for r in source['rows'] if r['kind'] == 'cement')), gap_total=sum((r['kind'] == 'cement' for r in source['rows'])), force_pass=sum((r['gate'] for r in source['rows'] if r['kind'] == 'crown')), force_total=sum((r['kind'] == 'crown' for r in source['rows'])), X1B_gates=source['X1B_original_validation_gates']), sufficiency=suff, dropouts=dict(source_unscored=sum((r['status'] == 'SOURCE_UNSCORED' for r in rows)), source_unscored_fraction=sum((r['status'] == 'SOURCE_UNSCORED' for r in rows)) / len(rows), adaptation_failures=sum((r['status'] == 'ADAPTATION_FAILED' for r in rows)), adapted_geometry_rejected=sum((not r['geometry_pass'] for r in ok)), force05_outside_source_or_insufficient_rays=len(ok) - len(qs), hydraulics_rejected=len(ok) - len(hs), all_v5_ports=len(ports), matched_complete_ports=matched, roof_or_partial_ports=len(ports) - matched, reason_counts=dict(collections.Counter((r['cement'].get('reason', 'accepted') for r in ok))), no_exclusion_hidden=True), controls=tests, external_referent=read(P / 'PREREG_R7.json')['external_referent'], external_measurement_on_generated_crown=False, uncertainty_debts=[dict(quantity='Absolute crown force', resolution='PHENOMENOLOGICAL', replacement_measurement='Same generated geometry, known batch/material/surface process, support/cement/contact/angle and held-out measured fracture force/origin'), dict(quantity='Seated cement film', resolution='PHENOMENOLOGICAL', replacement_measurement='Independent registered same-part normal gap field after manufacture/seating with pointwise error bounds, rheology/cure and seating load/time'), dict(quantity='Actual anatomical margin', resolution='PHENOMENOLOGICAL', replacement_measurement='Independent marked preparation cervical boundary and datum registration'), dict(quantity='Continuous geometric enclosure', resolution='PER_SURFACE_REGION', replacement_measurement='Certified continuous surface/ray sampling or independent high-resolution same-part metrology')], edges=read(P / 'PREREG_R1.json')['edges'], resources=accounting, artifact_manifest=artifact, frozen_predictions_sha256=sha(P / 'FROZEN_PREDICTIONS.json'), input_lock_sha256=sha(P / 'INPUT_LOCK.json'), linear_sensitivity_enclosure='No rigorous empirical/continuous-geometry enclosure; source parameter rectangles only. No affine sensitivity used as a physical guarantee.', full_cost=dict(preparation='Executed R1-R7 seconds in raw files; input reading/history not fully timed', fit='Frozen X1B fit reused; historical fit cost UNKNOWN', discovery='Reasoning and prior lane cost UNKNOWN', validation='Actual source/numerical/fault checks timed in run_all', queries='per-row geometry and physics seconds retained', fallback='Lab manufacturing/metrology/fracture acquisition cost UNKNOWN'), example_uid=example['uid'])
    result['complete_panel']['absolute_shape_035_pass'] = sum((r['adapted_shape_p95_mm'] <= 0.35 for r in ok))
    result['dropouts']['force_reasons'] = dict(collections.Counter((r['force'].get('source_force05', {}).get('status', r['force']['status']) for r in ok)))
    result['physical']['all_three_spatial_force_scenarios_in_source'] = sum((all((q['force05_N'] is not None for q in r['force'].get('spatial_scenario_force05', []))) and len(r['force'].get('spatial_scenario_force05', [])) == 3 for r in ok))
    source_vectors = {x['uid']: x for x in cohort()}
    full_vectors = [dict(uid=x['uid'], source_v5_vector=source_vectors[x['uid']], adapted_v5b_vector=x, physical_calibration='UNKNOWN', resolution='PER_TOOTH with PER_POINT field locators') for x in rows]
    (P / 'raw/QUALITY_VECTOR_FULL.json.gz').write_bytes(gzip.compress(json.dumps(clean(full_vectors), separators=(',', ':')).encode(), mtime=0))
    dump(P / 'results.json', result)
    fields = ['uid', 'participant', 'key', 'geometry_pass', 'adapted_shape_p95_mm', 'gap_mean_um', 'gap_error_p95_um', 'ray_coverage', 'conditional_force05_N', 'source_rectangle_low_N', 'source_rectangle_high_N', 'X1B_mean_diagnostic_N', 'X1B_interval_low_N', 'X1B_interval_high_N', 'conductance_m3_Pa_s', 'absolute_physical_status']
    with (P / 'raw/QUALITY_VECTOR.csv').open('w') as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for r in ok:
            q = r['force'].get('source_force05', {})
            band = q.get('parameter_rectangle_N') or [None, None]
            x = r['force'].get('X1B_mean_force_diagnostic', {})
            ci = x.get('interval_N', [None, None])
            writer.writerow(dict({k: r.get(k) for k in fields[:8]}, conditional_force05_N=q.get('force05_N'), source_rectangle_low_N=band[0], source_rectangle_high_N=band[1], X1B_mean_diagnostic_N=x.get('point_N'), X1B_interval_low_N=ci[0], X1B_interval_high_N=ci[1], conductance_m3_Pa_s=r['cement'].get('conductance_m3_Pa_s'), absolute_physical_status='UNKNOWN_UNVALIDATED_TRANSFER'))
    with (P / 'FACIT.csv').open('w') as f:
        w = csv.DictWriter(f, fieldnames=['study', 'doi', 'locator', 'quantity', 'predicted', 'observed', 'units', 'gate'])
        w.writeheader()
        w.writerows(({k: r.get(k) for k in w.fieldnames} for r in source['rows']))
    write_text(result, pairs, example)
    return result

def fmt(x, n=2):
    return '—' if x is None else f'{x:.{n}f}'

def write_text(r, pairs, example):
    c = r['complete_panel']
    p = r['physical']
    s = r['source_validation']
    rank = r['ranking']
    drop = r['dropouts']
    n = r['controls']['count']
    gates = c['geometry_accepted']
    table = '| Participant | Requested | Geometry accepted | Force05 in source + geometry | Shape before → after mm | Conditional force05 N | Spatial/mean flow |\n|---|---:|---:|---:|---:|---:|---:|\n'
    for x in r['leaderboard']:
        table += f"| {x['participant']} | {x['requested']} | {x['geometry_accepted']} | {x['force05_in_domain_accepted']} | {fmt(x['median_shape_before_mm'])} → {fmt(x['median_shape_after_mm'])} | {fmt(x['median_conditional_force05_N'], 0)} | {fmt(x['median_conductance_ratio'])} |\n"
    lead = f"GenCAD v5b aligns {gates}/{c['requested']} requested complete crowns to a virtual preparation while retaining local gap and thickness fields. The source-conditioned force axis reverses {rank['adapted_reversals']}/{rank['eligible_adapted_pairs']} eligible within-site adapted-shape orderings; calibrated physical ranking remains UNKNOWN. All final participant counts remain in the denominator. Claim: capability; PENDING_INDEPENDENT_REVIEW."
    (P / 'LEADERBOARD.md').write_text('# GenCAD v5b — physical query axis\n\n' + lead + '\n\n' + table + '\nAll table numbers aggregate PER_TOOTH results to POPULATION; medians use unequal accepted subsets and are not a global ranking. Pair-matched comparisons are in raw/RANK_COMPARISONS.json. ToothCraft remains a one-case adapted pilot. Geometry accepted is a numerical interface result, not milling/insertion/clinical acceptance.\n\n' + f"Original-shape comparisons: {rank['original_reversals']}/{rank['eligible_original_pairs']} eligible pairs reverse. {rank['disjoint_rectangle_reversals']} adapted reversals have disjoint source-parameter rectangles; these omit geometry/material/support transfer error and do not prove a physical rank.\n\n" + 'The full original+updated vector is in raw/QUALITY_VECTOR_FULL.json.gz; the compact vector is in raw/QUALITY_VECTOR.csv and raw/PHYSICAL_ROWS.json, including recalculated antagonist/proximal function after remeshing. Nominal spatial gap and Reynolds conductance are computed; observed seated film and calibrated absolute force are null. Original v5 is read-only. All roof/partial proposal ports remain explicit in raw/ALL_V5_PORTS.json.gz.\n')
    resulttable = f"| Question | Result | Resolution |\n|---|---|---|\n| Complete shells processed | {c['adapted']}/{c['requested']} requested | POPULATION |\n| Final geometry gates | {gates} accepted | PER_TOOTH, aggregated |\n| Source-conditioned force05 | {p['force05_source_curve_in_domain']} in source interval; {p['force05_source_curve_in_domain_and_geometry_accepted']} geometry accepted | PER_TOOTH |\n| X1B mean-force chain | {p['X1B_mean_queries']} queries; {p['X1B_strict_setup_thickness_supported']} in exact material/cement/angle thickness support | PER_TOOTH |\n| Spatial hydraulic queries | {p['spatial_hydraulic_queries']}; median spatial/mean conductance {fmt(p['median_conductance_spatial_over_mean'], 3)} | PER_SURFACE_REGION, aggregated |\n| Adapted shape / source force reversals | {rank['adapted_reversals']}/{rank['eligible_adapted_pairs']} eligible pairs | PER_TOOTH |\n| Published held targets | Dry gap {s['gap_pass']}/{s['gap_total']}; force mean {s['force_pass']}/{s['force_total']} | POPULATION |\n| Fault injections | {n}/{n} rejected | Implementation |\n| Matched manufactured force / seated-film measurements | 0 / 0 | PER_TOOTH / PER_POINT |\n"
    what = 'The useful operation is an explicit finish-plane intersection and annular cap after Boolean cavity construction. R1 point projection collapsed faces; R2 implicit Boolean corners displaced the margin; R3 grid alignment still beveled corners. R4 clips triangles against the exact plane and triangulates the annulus. Domain padding restores complete external-model envelopes. A bounded1µm contour recovery and preclip endpoint snapping within1e-6mm remove STL slivers; every decoded export must pass X55 and closure/volume gates. R6 stopped after3 failed pilots; its remaining186 were not run. Every failed gate and version is retained. X55 proper rigid coordinates/closest-surface operations, X1B preparation and fitted force code, X23 lower-tail convention and X13/BTE1 spatial conductance are reused. No new FE is used as an external reference.'
    limits = f"Absolute generated-crown fracture force and actual seated cement film remain **UNKNOWN**. X1B out-of-study transfer failed before this work; adaptation does not repair material/batch/support/contact/cement identity. The 0.5–1.5mm published force05 curve applies to its own specimen protocol. A median ray thickness cannot make a new anatomy equivalent to that specimen. The exact-summary counterexample has identity error0 and changes predicted force by a factor {r['sufficiency']['thickness']['ratio']:.4f}; the minimum fixed-query extension is local load-weighted hazard. The film counterexample has identity error0 and conductance ratio {r['sufficiency']['film']['ratio']:.4f}; boundary-connected resistance or the full spatial film is needed.\n\nThe absolute0.35mm reconstruction gate remains failed on all final crowns; accepted here means the stated preparation/export interface gates, not adequate anatomy. Exact margin identity refers to the float64 mesh; binary STL adds coordinate rounding, whose margin error is not enclosed. Reloaded topology, face areas and volume are tested. No rigorous continuous geometry or physical-transfer enclosure is available. Source parameter rectangles are mathematical corner ranges, not joint confidence intervals; X1B 90% intervals are provisional diagnostics and not validated crown uncertainty. The nominal gap comes from a coarse, virtual preparation, not a manufactured seating measurement. Marked anatomical margins, self-intersection/insertion and milling feasibility are not certified. A new physical ranking cannot yet be claimed."
    sources = 'External references are pinned local primary data: [Weibull crown loads, Table1](https://doi.org/10.4047/jap.2021.13.5.269), [held crown group loads, Table2](https://doi.org/10.3390/ma17020365), [dry marginal gap, Table2](https://doi.org/10.7759/cureus.38688), [second dry marginal study, Table1](https://doi.org/10.1155/2023/6698453). FACIT.csv names every compared quantity, locator and failed target. Measured tooth surfaces are the frozen local v4 reference.npz files listed in INPUT_LOCK.json; they are an anatomical referent, not optimal crown or physical fracture truth.'
    cost = f"The one-command replay took {r['resources'].get('wall_seconds_before_reporting', 0):.2f}s before reporting, peak RSS {r['resources'].get('peak_rss_MiB', 0):.1f}MiB, at most4 threads, no GPU. Data arrays/meshes occupy {r['artifact_manifest']['total_bytes'] / 1000000.0:.1f}MB under the lane data directory; historical research/lab/fit costs remain UNKNOWN. Full fresh geometry is available with --rebuild; default verifies cached geometry hashes and freshly executes its first accepted case, all force/hydraulic/function queries and fault controls."
    dropout = f"Dropout: {drop['source_unscored']}/{c['requested']} source failures ({drop['source_unscored_fraction']:.2%}), {drop['adaptation_failures']} construction failures, {drop['adapted_geometry_rejected']} rejected adapted geometries, {drop['force05_outside_source_or_insufficient_rays']} outside-source/insufficient force queries and {drop['hydraulics_rejected']} rejected hydraulic inputs. The {drop['roof_or_partial_ports']} roof/partial ports are not complete crown specimens and remain typed as unmapped; no failure is converted to a good score. Reasons and exact cohort support are machine readable."
    readme = '# GenCAD v5b: preparation → physical query vector\n\n' + lead + '\n\nRun `./run_all.sh` in this directory. Use `./run_all.sh --rebuild` to regenerate every R7 shell; all inputs are local. Runtime versions are pinned in RUNTIME_LOCK.json/requirements.lock; DENTAL_PYTHON can select an equivalent local interpreter. `[exports/example_crown.stl]` is a concrete accepted research shell and the paired virtual preparation is supplied. All other meshes/fields are hash-linked in ARTIFACT_MANIFEST.json; `./export_design.sh <uid> <output.stl>` exports a chosen design in the registered frame.\n\n' + resulttable + '\n' + what + '\n\n![Physical query axes](figures/physical_axes.png)\n\n![Actual adapted geometry and field](figures/crown_field.png)\n\n' + limits + '\n\n' + sources + '\n\n' + dropout + '\n\nControls: actual source-value tests, independent X1B algebra, Weibull inverse CDF, X55 datum transform and closed-form Reynolds conductance agree where applicable. Each has a rejected numerical fault; this is no algorithm-superiority claim.\n\n' + cost + '\n\nLicences: local Bits2Bites CC BY-NC-SA per supplied dataset record; exact Bite2Text redistribution terms UNKNOWN, private local research only. ToothCraft checkpoint licence unspecified, model outputs used from v5 without inference/download. X11/Teeth3DS training lineage follows its source terms; no new Teeth3DS meshes copied. STS-derived X1B is reused as code/calibration context, not redistributed source anatomy. Primary article permissions stay with the source (PMC10246932 CC BY; other article permissions per original records). No ToothFairy2/MMDental/mandible-defect dataset used. No personal identifiers or clinical text are exported.\n'
    (P / 'README_DEMO.md').write_text(readme.replace('`[exports/example_crown.stl]`', '[exports/example_crown.stl](exports/example_crown.stl)'))
    (P / 'RESULTS.md').write_text('# GenCAD v5b result\n\n' + lead + '\n\n' + resulttable + '\n' + what + '\n\n' + limits + '\n\n' + sources + '\n\n' + dropout + '\n\n' + cost + '\n')
    (P / 'HANDOFF.md').write_text('# PROOF_LANE-gencad-v5b handoff\n\n' + lead + '\n\n' + resulttable + '\n' + what + '\n\n' + limits + '\n\n' + dropout + '\n\nNext construction: manufacture one accepted exported geometry with an independently marked margin and external datum, acquire signed same-part dry and seated film plus a local uncertainty reference using X55; measure batch/material/flaw origin, contact/support/angle and failure force for matched calibration and held specimens. Freeze a NEW specimen model before measurement. Do not fit another thickness power law or use own FE as truth. LAB_PROTOCOL.md gives the required fields; present numeric predictions remain a retrospective conditional diagnostic.\n\nGraph: existing DENT-VAL-BENCH-CROWN-FIT was dispatched as review because numerical experiment admission is false and inherited hash drift is present. New coverage proposal: DENT-GENCAD-PREP-FILM-FRACTURE, PER_POINT geometry/gap SIMULTANEOUS; PER_SURFACE_REGION seated-state to loaded support HANDOVER; material/batch force calibration POPULATION with explicit per-specimen lineage. No native/source graph changed. GRAPH_FEEDBACK.json is PENDING_INDEPENDENT_REVIEW.\n\nRead raw/R7_GEOMETRY.json, raw/PHYSICAL_ROWS.json, raw/RANK_COMPARISONS.json and FROZEN_PREDICTIONS.json before continuing. All failure versions and exact commands are retained. ' + cost + '\n')
