import csv, xml.etree.ElementTree as ET, copy
from codesign import *
from sliver_scene import SliverScene
from calibration import calibration_controls, displacement_bound
from collision import segment_triangle, segment_segment, point_triangle

def external():
    root = ET.parse(ROOT / 'sources/PMC8391771.xml').getroot()
    table = root.find(".//table-wrap[@id='healthcare-09-00983-t001']")
    rows = table.findall('.//tbody/tr')
    region = None
    out = []
    candidates = 0
    for tr in rows:
        vals = [''.join(td.itertext()).strip() for td in tr]
        if vals and vals[0] in ('Whole\n region', 'Marginal\n region', 'Axial\n region', 'Occlusal\n region'):
            region = ' '.join(vals.pop(0).split())
        elif vals and any((w in vals[0] for w in ('Whole', 'Marginal', 'Axial', 'Occlusal'))):
            region = ' '.join(vals.pop(0).split())
        if not vals or vals[0] not in ('SLA', 'DLP', 'Milling'):
            continue
        candidates += 1
        if vals[0] != 'Milling':
            continue
        mean = float(vals[1].split()[0])
        sd = float(vals[2])
        ci = [float(vals[3]), float(vals[4])]
        out.append(dict(region=region, mean_RMS_um=mean, sd_um=sd, CI95_um=ci, min_RMS_um=float(vals[5]), max_RMS_um=float(vals[6]), units='um', resolution='PER_SURFACE_REGION', population='Published group mean, n=15 crowns per method', locator='https://doi.org/10.3390/healthcare9080983', table_id=table.attrib['id'], quantity='CAD-to-milled scanned intaglio surface RMS trueness; best-fit registration', compared_to_20um_scenario=mean > 20, not_pure_tool_deflection=True))
    if len(out) != 4:
        raise ValueError('Expected4 published milling regional cells, got ' + str(out))
    dump(ROOT / 'raw/EXTERNAL_FACIT.json', dict(rows=out, source_sha256=sha(ROOT / 'sources/PMC8391771.xml'), source_selection=dict(candidates=candidates, accepted=4, rejected=candidates - 4, rejection='SLA/DLP fabrication regime, not milling', dropout_fraction=(candidates - 4) / candidates), interpretation='Refutes universal <=20um full-chain machining error from nominal0.6mm burr. Does NOT refute conditional beam-only inverse algebra; material is PMMA, machine CORITEC250i; no paired local mesh or force/holder data.', paired_beam_facit=dict(accepted=0, candidates=2, rejected=2, dropout_fraction=1.0, reasons=['PMC8391771: full-chain RMS, no paired force/flute/holder measurements', 'PMC6187571: fully sintered zirconia micro-end-milling, published force figures but no paired installed compliance/crown error'])))
    with open(ROOT / 'FACIT.csv', 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(out[0]))
        w.writeheader()
        w.writerows(out)
    return out

def sliver_controls(rows):
    rr = []
    rng = np.random.default_rng(8902)
    for row in rows:
        a = np.load(row['file'], allow_pickle=False)
        s = SliverScene(a['vertices'], a['faces'])
        if not s.sliver.any():
            continue
        tri = s.tri[s.sliver]
        max_upper = -np.inf
        errs = []
        for (j, t) in enumerate(tri):
            for k in range(16):
                p = t.mean(0) + rng.uniform(-0.002, 0.002, 3)
                exact = point_triangle(p, t[None])[0]
                q = segment_segment(p, p, s.a[j:j + 1], s.b[j:j + 1])[0] - s.width[j]
                errs.append(float(q - exact))
                max_upper = max(max_upper, q - exact)
        p = tri[0].mean(0) + np.array([0.0, 0.0, 0.0001])
        exact = point_triangle(p, tri[0:1])[0]
        bound = segment_segment(p, p, s.a[0:1], s.b[0:1])[0] - s.width[0]
        check = lambda proposed: proposed <= exact + 1e-10
        rr.append(dict(id=row['id'], sliver_faces=int(s.sliver.sum()), sliver_face_ids=np.flatnonzero(s.sliver), width_max_mm=float(s.width.max()), max_cover_minus_independent_triangle_mm=max_upper, baseline_pass=max_upper <= 1e-10 and check(bound), injected_deleted_obstacle_rejected=not check(0.05), injection=dict(point_mm=p, independent_distance_mm=exact, baseline_cover_mm=bound, corrupt_cover_mm=0.05), exact_surface_area_mm2=float(s.area[s.sliver].sum()), coverage_status='all slivers retained in collision cover; sampling excludes tiny facets'))
    assert rr and all((x['baseline_pass'] for x in rr))
    return rr

def run():
    pre = lock_prereg('R2')
    lock_prereg()
    start = time.perf_counter()
    rows = intake()
    m = pre['metrics']
    card = read(PARENT / 'NOMINAL_TIP_REFERENCE_CARD.json')['tools']
    r1 = read(ROOT / 'RESULTS_R1.json')
    old = {x['id']: x for x in r1['designs']}
    summary = []
    hashes = {}
    for (ix, row) in enumerate(rows):
        state('R2_RUNNING', f'{len(summary)}/8 mixed-dimensional queries', 'Retain all slivers; coupled local requirements')
        if row['family'] == 'full_crown_R4':
            row['dataset'] = 'Teeth3DS'
            row['license'] = 'UNKNOWN local mesh licence, per PROOF_LANE_FULL_CROWN_R4 README_DEMO.md'
            row['locator'] = 'https://osf.io/xctdy/'
        a = np.load(row['file'], allow_pickle=False)
        scene = SliverScene(a['vertices'], a['faces'])
        goodregions = {k: np.asarray(v)[scene.valid[np.asarray(v)]] for (k, v) in row['regions'].items()}
        pts = samples(scene, goodregions, m['seed'] + ix, m['n_per_region'])
        rr = []
        for p in pts:
            tools = [robust_search(scene, p['point'], p['normal'], t, m) for t in card]
            vi = [t for t in tools if t['status'] == 'CONDITIONAL_LOCAL_WITNESS']
            rr.append(dict(**p, tools=tools, selected_tool=min(vi, key=lambda x: x['nominal_tip_diameter_mm'])['tool_id'] if vi else None, status='CONDITIONAL_LOCAL_WITNESS' if vi else 'NO_LIBRARY_WITNESS_AT_DECLARED_LOAD', resolution='PER_POINT', time_scale='HANDOVER'))
        out = DATA / 'R2' / f"{row['id']}.json"
        dump(out, dict(input={k: v for (k, v) in row.items() if k != 'regions'}, prereg_sha256=sha(ROOT / 'PREREG_R2.json'), rows=rr, sliver_faces=int(scene.sliver.sum())))
        hashes[str(out)] = sha(out)
        ss = dict(id=row['id'], family=row['family'], points=len(rr), found=sum((x['selected_tool'] is not None for x in rr)), ball_only=sum((any((t['ball_only'] for t in x['tools'])) for x in rr)), raw_file=out, sliver_faces=int(scene.sliver.sum()), physical_status='UNKNOWN_MISSING_ASSEMBLY_FORCE_CAM_PROCESS_AND_WHOLE_SURFACE_BOUND', regions={reg: dict(n=sum((x['region'] == reg for x in rr)), found=sum((x['region'] == reg and x['selected_tool'] is not None for x in rr))) for reg in row['regions']})
        summary.append(ss)
        print(ss['id'], ss['found'], '/', ss['points'], flush=True)
        if row['family'] == 'crown_loop':
            before = read(old[row['id']]['raw_file'])['rows']
            assert clean(rr) == before, 'Original well-conditioned crown-loop results changed'
    t1 = time.perf_counter()
    checks = sliver_controls(rows)
    cal = calibration_controls()
    facit = external()
    val_s = time.perf_counter() - t1
    assert all((cal[k] for k in ['injected_missing_curvature_rejected', 'injected_extrapolation_rejected', 'injected_pose_rejected', 'injected50um_rejected', 'injected_understated_enclosure_rejected']))
    result = dict(claim_type='capability', round='R2', status='ALL8_CONDITIONAL_LOCAL_REQUIREMENTS; PHYSICAL_UNKNOWN', designs=summary, tool_mechanics=r1['tool_mechanics'], metrics=m, prereg_sha256=sha(ROOT / 'PREREG_R2.json'), controls=dict(slivers=checks, calibration=cal, beam=r1['controls']), external_referent=dict(kind='independent_measurement', locator='https://doi.org/10.3390/healthcare9080983; PMC8391771 Table1', compared_quantity='Regional CAD-to-milled intaglio RMS trueness [um], published group means', refutes_us=False), external_refutation=dict(assertion='Nominal smallest0.6mm tip assures full-chain local error<=20um', refuted=True, evidence=facit, matched_beam_validation='UNKNOWN_NO_PAIRED_ASSEMBLY'), artifact_sha256=hashes, resolution='PER_POINT', time_scale='HANDOVER', dropout=dict(design_rejected=0, requested=8, physical_qualifications_missing=8, sliver_sampling_excluded=sum((x['sliver_faces'] for x in summary)), sliver_collision_deleted=0), cost=dict(total_wall_s=time.perf_counter() - start, query_s=time.perf_counter() - start - val_s, validation_s=val_s, fit_s=0, questions_s=0, prior_search_preparation_s=None, fallback_s=None, peak_RSS_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, gpu=0), rigorous_floating_enclosure=False, rigorous_calibration_arithmetic=True, physical_uncertainty='UNBOUNDED missing force/neck/holder/CAM/sinter/runout data')
    dump(ROOT / 'RESULTS_R2.json', result)
    dump(ROOT / 'results.json', result)
    dump(ROOT / 'raw/CONTROLS_R2.json', result['controls'])
    state('R2_COMPLETE', '8 designs queried; same5 crown-loop rows; calibration/sliver injections passed', 'Integrate read-only X49 rule; report external mismatch and physical debt; CLI failure probes')
    return result
if __name__ == '__main__':
    run()
