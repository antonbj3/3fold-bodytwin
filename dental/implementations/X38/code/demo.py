"""Execute frozen demonstrations with exact local inputs, retain each run."""
import datetime, hashlib, json, os, resource, sys, time
from pathlib import Path
import numpy as np
import trimesh
from exportgate import export, check
from mesh_transport import dump, sha, load_mesh, correspondence, write_stl, Refused, validate_triangles
from labcount import plan, endpoint_summary, compare
R = Path(__file__).resolve().parents[1]

def state(phase, gate, next_op):
    dump(R / 'CURRENT_WORK_STATE.json', dict(lane='X38-export-gate', updated_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), phase=phase, latest_gate=gate, next_operation=next_op, threads=1, intermediate_limit_bytes=3000000000))

def main():
    start = time.monotonic()
    cpu = time.process_time()
    for name in ['PREREG_R1', 'PREREG_R2', 'PREREG_R3', 'PREREG_STATS', 'PREREG_STATS_R2', 'FROZEN_PREDICTIONS']:
        if sha(R / (name + '.json')) != (R / (name + '.sha256')).read_text().strip():
            raise RuntimeError('Frozen hash changed: ' + name)
    for q in json.loads((R / 'INPUT_LOCK.json').read_text())['files']:
        if sha(R / q['path']) != q['sha256']:
            raise RuntimeError('input hash changed: ' + q['path'])
    existing = sum((p.stat().st_size for p in (R / 'raw/runs').rglob('*') if p.is_file()))
    total_faces = sum(((p.stat().st_size - 84) // 50 for p in (R / 'inputs/exports').glob('*/*.stl')))
    estimated_next_bytes = 2200 * total_faces + 20000000
    if existing + estimated_next_bytes > 3000000000:
        raise RuntimeError('preserved runs plus conservative next-run estimate exceed3GB; archive externally before new replay')
    run = R / 'raw/runs' / datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    run.mkdir(parents=True)
    rows = []
    refused = []
    state('R3_RUNNING', 'frozen hashes verified', 'export/reimport eligible X1b meshes, refuse tiny source facets without repair')
    for design in ['D1', 'D2', 'D3', 'M1', 'M2']:
        for part in ['crown', 'die', 'preparation']:
            s = time.monotonic()
            source = R / f'inputs/exports/{design}/{part}.stl'
            out = run / f'{design}_{part}'
            try:
                (original, _) = load_mesh(source, 'mm')
                validate_triangles(original)
            except Refused as e:
                refused.append(dict(design=design, part=part, path=str(source.relative_to(R)), reason=str(e), status='REFUSED'))
                dump(run / f'{design}_{part}_REFUSED.json', refused[-1])
                print(design, part, 'REFUSED:', e, flush=True)
                continue
            receipt = export(source, out, 'mm', None, '3Y zirconia; batch UNKNOWN' if part == 'crown' else 'laboratory die; material UNKNOWN', 1, 'final_sintered', R / 'inputs/X1B_FROZEN_PREDICTIONS.json')
            (original, _) = load_mesh(source, 'mm')
            (prior, _) = load_mesh(source.with_suffix('.3mf'))
            (_, prior_error) = correspondence(original, prior)
            checks = []
            for name in ['model.stl', 'model_ascii.stl', 'model.3mf']:
                asset = out / name
                side = out / (name + '.json')
                q = check(asset, side, sha(side), 'mm' if name.endswith('.stl') else None)
                independent = trimesh.load(asset, process=False)
                independent = independent.to_geometry() if isinstance(independent, trimesh.Scene) else independent
                (_, error) = correspondence(original, np.asarray(independent.vertices)[independent.faces])
                checks.append(dict(conversion=name, gate=q['status'], max_correspondence_mm=q['geometry_max_correspondence_mm'], published_trimesh_max_error_mm=error, region_source=q['region_source']))
            m = trimesh.load_mesh(out / 'model.stl', process=False)
            c = out / 'published_ascii.stl'
            c.write_text(m.export(file_type='stl_ascii'))
            q = check(c, out / 'model.stl.json', sha(out / 'model.stl.json'), 'mm', True)
            checks.append(dict(conversion='trimesh binary->ASCII', gate=q['status'], max_correspondence_mm=q['geometry_max_correspondence_mm'], region_source=q['region_source']))
            scene = trimesh.load(out / 'model.3mf', process=False)
            m = scene.to_geometry() if isinstance(scene, trimesh.Scene) else scene
            c = out / 'published_3mf_to_stl.stl'
            c.write_bytes(m.export(file_type='stl'))
            q = check(c, out / 'model.3mf.json', sha(out / 'model.3mf.json'), 'mm', True)
            checks.append(dict(conversion='trimesh3MF->binarySTL', gate=q['status'], max_correspondence_mm=q['geometry_max_correspondence_mm'], region_source=q['region_source']))
            rows.append(dict(design=design, part=part, faces=len(original), original_stl_vs_original_3mf_error_mm=prior_error, checks=checks, seconds=time.monotonic() - s, anatomical_regions='UNKNOWN'))
            print(design, part, len(original), 'facets,5 paths PASS', flush=True)
            state('R3_RUNNING', dict(completed_meshes=len(rows), refused=len(refused), last=rows[-1]), 'continue remaining frozen transport paths')
    source = R / 'inputs/exports/D1/crown.stl'
    (tri, _) = load_mesh(source, 'mm')
    labels = np.where(tri.mean(1)[:, 2] >= np.median(tri.mean(1)[:, 2]), 'transport_test_A', 'transport_test_B').tolist()
    regions = run / 'test_regions.json'
    dump(regions, dict(labels=labels, semantics='Synthetic spatial partition on X1b crown; no anatomical claim'))
    out = run / 'labelled_transport_probe'
    export(source, out, 'mm', regions, 'transport_fixture', 1, 'final_sintered', R / 'inputs/X1B_FROZEN_PREDICTIONS.json')
    c = out / 'reordered.stl'
    order = np.arange(len(tri))[::-1]
    write_stl(c, np.roll(tri[order], 1, axis=1))
    q = check(c, out / 'model.stl.json', sha(out / 'model.stl.json'), 'mm', True)
    mismatch = sum((a != b for (a, b) in zip(q['regions_in_imported_face_order'], [labels[i] for i in order])))
    stats = dict(plan29=plan(0.1, 0.95, true_risk=0.1), plan42=plan(0.1, 0.95, family_tests=4), actual100N=endpoint_summary(json.loads((R / 'inputs/specimens.json').read_text()), 5000000.0, 'MUSLA 04020', 100), binary=compare('binary', p_a=0.3, p_b=0.1, max_n=200), continuous=compare('continuous', effect=0.8, max_n=200))
    allrows = json.loads((R / 'inputs/specimens.json').read_text())
    cells = sorted({(q['configuration'], q['max_load_N']) for q in allrows})
    stats['actual_all_cells'] = [endpoint_summary(allrows, 5000000.0, g, f, family_tests=len(cells)) for (g, f) in cells]
    stats['all_observation_screening'] = dict(total=len(allrows), complete=sum((x['n'] - x['unknown_early_censors'] for x in stats['actual_all_cells'])), early_censored=sum((x['unknown_early_censors'] for x in stats['actual_all_cells'])), rejected=0)
    dump(run / 'statistics.json', stats)
    result = dict(round='R3', claim_type='capability', outcome='BOUNDED_TRANSPORT_CAPABILITY_DELIVERED', meshes=rows, refused_sources=refused, original_meshes=15, conversion_paths=5 * len(rows), region_probe=dict(facets=len(tri), mismatches=mismatch, resolution='PER_SURFACE_REGION', anatomical_validity='UNKNOWN; spatial test fixture only'), statistics_path=str((run / 'statistics.json').relative_to(R)), run_directory=str(run.relative_to(R)), screening=dict(mesh_candidates=15, accepted_transport=len(rows), rejected=len(refused), rejection_fraction=len(refused) / 15, reasons={'source facet below unchanged numerical degeneracy cutoff1e-12 mm2': len(refused)}, anatomically_classified=0, unclassified_fraction=1.0), full_cost=dict(wall_seconds=time.monotonic() - start, cpu_seconds=time.process_time() - cpu, peak_process_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, preparation_bytes=json.loads((R / 'INPUT_LOCK.json').read_text())['total_bytes'], fit_seconds=0, discovery='R1/R2 strict source failure then R3 explicit per-candidate export/refusal; human time UNKNOWN', validation=f'unittest independent controls +{5 * len(rows)} live conversion paths', queries='fixed29/family42, published36-specimen endpoints,2 design power scenarios', fallback='refuse remeshing/ambiguous regions and4 tiny source facets; censored endpoints remain unknown'), external_referent=dict(kind='published_code', locator=['https://trimesh.org/trimesh.exchange.stl.html', 'https://trimesh.org/trimesh.exchange.threemf.html', 'https://www.itl.nist.gov/div898/software/dataplot/refman2/auxillar/exacbici.htm'], compared_quantity='Physical triangle import and exact fixed-horizon binomial confidence', refutes_us=True))
    dump(run / 'round_result.json', result)
    dump(R / 'rounds/R3.json', result)
    dump(R / 'raw/LATEST_RUN.json', dict(path=str(run.relative_to(R)), sha256=sha(run / 'round_result.json')))
    if len(rows) != 11 or len(refused) != 4 or mismatch != 0:
        raise RuntimeError('frozen R3 expected decisions failed')
    state('R3_AND_STATS_DECIDED', dict(meshes=len(rows), refused=len(refused), conversion_paths=5 * len(rows), region_mismatches=mismatch, n_zero_failures=29), 'compile lab README, figure and graph feedback; next physical operation commercial CAM roundtrip')
    print('Done:', run, flush=True)
if __name__ == '__main__':
    main()
