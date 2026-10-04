"""Complete requested per-design surface reporting; no new superiority claim.

R2's twenty frozen intaglio maps are reused unchanged. Other regions use the
identical R2 operation/thresholds. Each region is a separate sparse ball set;
neither this set nor the forced diagnostic is a CAM toolpath stock model.
"""
from dental_release.paths import expand as _release_expand
import time, resource
from common import *
from collision import Scene, search, first_union_ball_entry
from run_r1 import lib
from run_r2 import normal_gaps, closed_cap
import trimesh

def run():
    prepath = ROOT / 'PREREG_COVERAGE_R2.json'
    if not prepath.exists():
        pre = dict(claim_type='capability', kind='Requested coverage completion, not a new construction', capability='Per-design sampled stock/overcut fractions for all sixty meshes and both surface regions; same-point film transport only on intaglio', obstacle='R2 maps only twenty locked lab cases', changed_operation='None: extend identical frozen R2 operator with independently retained region scope', consumer=['DENT-MFG-ASBUILT-DEVIATION', _release_expand('X13')], resolution='PER_POINT', area_resolution='PER_SURFACE_REGION', timescale='HANDOVER', metrics=read(ROOT / 'PREREG_R2.json')['metrics'], frozen_R2_sha256=sha(ROOT / 'FROZEN_PREDICTIONS_R2.json'), decision_criteria='All60 designs have intaglio/exterior maps with unknown rates; original20 intaglio hashes remain unchanged; no physical stock/fit claim', strongest_equal_information_control='Identical R2 scalar sphere-ray substitution plus injected50um rejection, already run; no method competition', external_referent=read(ROOT / 'PREREG_R2.json')['external_referent'], falsifier='Original20 maps change or any unexplained unknown disappears', full_cost=dict(preparation='read frozen R1/R2 fields', fit=0, discovery='40 additional intaglio and60 exterior regions', validation='R2 controls and immutable hash checks', questions=0, fallback='Unknown stock/normal gap retained', physical='No new CAM or measurement'), assumptions=['Independent sparse cutting-ball union per region', 'Nominal fixed-normal film is not seating', 'Same declared vhf5 capsule geometry and.8 scale', 'No rigorous floating arithmetic enclosure; inherited limitations'])
        dump(prepath, pre)
        prepath.with_suffix('.sha256').write_text(sha(prepath) + '\n')
    assert sha(prepath) == prepath.with_suffix('.sha256').read_text().strip()
    pre = read(prepath)
    r1 = read(ROOT / 'RESULTS_R1.json')
    r2 = read(ROOT / 'RESULTS_R2.json')
    assert sha(ROOT / 'FROZEN_PREDICTIONS_R2.json') == pre['frozen_R2_sha256']
    originals = {read(f)['id']: (f, h) for (f, h) in r2['raw_files'].items()}
    manifest = {r['id']: r for r in read(ROOT / 'INPUT_MANIFEST.json')['rows']}
    start = time.perf_counter()
    rows = []
    hashes = {}
    budget()
    state('SURFACE_COVERAGE_RUNNING', 'Frozen coverage contract; original twenty maps unchanged', 'Complete remaining regions with existing R2 operation')
    for (f, h) in r1['raw_files'].items():
        assert sha(f) == h
        old = read(f)
        ident = old['id']
        row = manifest[ident]
        dest = DATA / 'surface_ports' / (ident + '.json')
        if dest.exists():
            out = read(dest)
            assert out['prereg_sha256'] == sha(prepath)
            rows.extend(out['regions'])
            hashes[str(dest)] = sha(dest)
            continue
        assert sha(row['file']) == row['sha256']
        a = np.load(row['file'], allow_pickle=False)
        scene = Scene(a['vertices'], a['faces'])
        variant = next((v for v in old['points_and_tools'] if v['library'] == 'vhf' and v['axes'] == 5))
        tools = lib('vhf', 5)
        regions = []
        allrecords = []
        for region in ['intaglio', 'exterior_and_rim']:
            rec = [p for p in variant['records'] if p['region'] == region]
            if region == 'intaglio' and ident in originals:
                (path, digest0) = originals[ident]
                assert sha(path) == digest0
                out = read(path)
                records = out['records']
                summary = dict(out['summary'], region=region, original_R2_file=path, original_R2_sha256=digest0)
            else:
                P = np.array([p['point'] for p in rec])
                N = np.array([p['normal'] for p in rec])
                selected = []
                centres = []
                radii = []
                for (p0, p, n) in zip(rec, P, N):
                    available = [w for w in p0['tools'] if w['status'] == 'FOUND']
                    if available:
                        w = min(available, key=lambda w: w['diameter_design_mm'])
                    else:
                        w = None
                        for off in pre['metrics']['offset_grid_mm']:
                            for t in tools:
                                v = search(scene, p, n, t, 5, offsets=(off,))
                                if v['status'] == 'FOUND':
                                    w = dict(**v, tool_id=t['id'], diameter_design_mm=t['diameter_mm'])
                                    break
                            if w is not None:
                                break
                    if w is not None:
                        centres.append(w['centre'])
                        radii.append(w['diameter_design_mm'] / 2)
                    selected.append(w)
                stock = first_union_ball_entry(P, N, centres, radii)
                forced = first_union_ball_entry(P, N, P + tools[0]['diameter_mm'] / 2 * N, np.full(len(P), tools[0]['diameter_mm'] / 2))
                if region == 'intaglio':
                    if row.get('die'):
                        assert sha(row['die']) == row['die_sha256']
                        die = trimesh.load_mesh(row['die'], process=True)
                    else:
                        die = closed_cap(a['prep_vertices'], a['prep_faces'])
                    h0 = normal_gaps(die, P, N)
                else:
                    h0 = np.full(len(P), np.nan)
                records = [dict(point=P[i], normal=N[i], face_id=rec[i]['face_id'], normal_gap_nominal_um=h0[i] * 1000, stock_safe_ball_union_um=stock[i] * 1000, forced_signed_stock_um=forced[i] * 1000, normal_gap_after_safe_union_um=(h0[i] - stock[i]) * 1000, normal_gap_after_forced_union_um=(h0[i] - forced[i]) * 1000, safe_pose_witness=selected[i], stock_status='FINITE_BALL_WITNESS_SET' if np.isfinite(stock[i]) else 'UNKNOWN_NO_RAY_CUT', fit_status='CONDITIONAL_FIXED_NORMAL_RAY' if np.isfinite(h0[i]) else 'UNKNOWN_OR_NOT_APPLICABLE_EXTERIOR', geometry_resolution='PER_POINT') for i in range(len(P))]
                summary = dict(id=ident, family=row['family'], region=region, n=len(P), stock_gt25um_fraction=float(np.mean(stock > 0.025)), safe_union_overcut_gt25um_fraction=float(np.mean(stock < -0.025)), forced_overcut_gt25um_fraction=float(np.mean(forced < -0.025)), unknown_stock_fraction=float(np.mean(~np.isfinite(stock))), unknown_nominal_gap_fraction=float(np.mean(~np.isfinite(h0))), empirical_fit='UNKNOWN')
            summary.update(found_zero_offset_fraction=sum((p['smallest_found_green_mm'] is not None for p in rec)) / len(rec), not_found_zero_offset_fraction=sum((p['smallest_found_green_mm'] is None for p in rec)) / len(rec), area_mm2=sum((p['area_weight_mm2'] for p in rec)), fraction_resolution='PER_SURFACE_REGION', point_resolution='PER_POINT', hoeffding_95_halfwidth=float(np.sqrt(np.log(40) / (2 * len(rec)))), stock_population_scope='All sampled points, including unknowns in denominator; safe/forced sparse ball unions independently per region; not CAM', physical_gate='UNKNOWN')
            for p in records:
                p = dict(p, region=region)
                allrecords.append(p)
            regions.append(summary)
        dump(dest, dict(id=ident, prereg_sha256=sha(prepath), source_sha256=row['sha256'], regions=regions, records=allrecords))
        rows.extend(regions)
        hashes[str(dest)] = sha(dest)
        state('SURFACE_COVERAGE_RUNNING', str(len(rows)) + '/120 regions saved', 'Continue without changing frozen twenty intaglios')
    assert len(rows) == 120
    for (f, h) in r2['raw_files'].items():
        assert sha(f) == h
    result = dict(claim_type='capability', coverage_only=True, requested_designs=60, regions=120, points=7680, rows=rows, raw_files=hashes, seconds=time.perf_counter() - start, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, original20_unchanged=True, physical_gate='UNKNOWN', external_referent=pre['external_referent'])
    dump(ROOT / 'RESULTS_SURFACE_COVERAGE.json', result)
    freeze(ROOT / 'FROZEN_PREDICTIONS_SURFACE_COVERAGE.json', dict(prereg_sha256=sha(prepath), result_sha256=sha(ROOT / 'RESULTS_SURFACE_COVERAGE.json'), point_files=hashes, code_sha256=sha(Path(__file__)), physical_measurement_performed=False))
    print('surface coverage', len(rows), result['seconds'])
if __name__ == '__main__':
    run()
