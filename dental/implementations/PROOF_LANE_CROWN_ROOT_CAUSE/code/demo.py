"""One-command reproduction in an isolated data directory, without source edits."""
from dental_release.paths import expand as _release_expand
from joint_solid import *
import joint_solid as js
import oriented_crown as oc
import score_ray as sc
import score_joint as sj
import validate as val
import figure as fig
import shutil

def main():
    start = time.perf_counter()
    root = R
    data = D
    stamp = datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    work = data / ('demo_' + stamp)
    (work / 'raw').mkdir(parents=True)
    (work / 'data').mkdir()
    (work / 'figures').mkdir()
    if sum((p.stat().st_size for p in data.rglob('*') if p.is_file() and (not p.is_symlink()))) > 2500000000.0:
        raise RuntimeError('Lane budget reserve below 0.5 GB; no large writes made')
    for name in ['PREREG_C.json', 'PREREG_D.json', 'RESULTS_A2.json', 'RESULTS_A.json']:
        shutil.copyfile(root / name, work / name)
    for name in ['close_native', 'boolean']:
        (work / 'data' / name).symlink_to(data / name)
    (work / _release_expand('B31_FULL_CROWN_R6')).symlink_to(data.parent / _release_expand('B31_FULL_CROWN_R6'), target_is_directory=True)
    for mod in [js, oc, sc, sj, val, fig]:
        mod.R = work
    for mod in [js, oc, sc, sj]:
        mod.D = work / 'data'
    checks = val.run()
    frozen = {r['key']: r for r in read(root / 'RESULTS_D.json')['rows']}
    inputs = read(B / 'PROOF_LANE_FULL_CROWN_R4/FROZEN_INPUTS_B.json')['records'][:3]
    rows = []
    for inp in inputs:
        rec = oc.generate_d(inp)
        pred = work / ('PREDICTION_' + inp['key'] + '.json')
        dump(pred, dict(frozen_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), record=rec, scope='Digital geometry; frozen before independent scorer'))
        pred.with_suffix('.sha256').write_text(sha(pred) + '  ' + pred.name + '\n')
        got = sc.evaluate_c(rec)
        old = frozen[inp['key']]
        assert got['geometric_conjunction']
        assert got['chosen'] == old['chosen']
        old_arrays = load(old['mesh_path'])
        new_arrays = load(got['mesh_path'])
        identity = all((np.array_equal(old_arrays[k], new_arrays[k]) for k in old_arrays))
        assert identity, 'Regenerated mesh arrays changed'
        got['mesh_arrays_identical_to_frozen'] = identity
        rows.append(got)
        print('REPLAY', inp['key'], 'all mesh arrays identical; geometric gates pass', flush=True)
    dump(work / 'RESULTS_D.json', dict(rows=rows))
    fig.run()
    result = dict(status='PASS', scope='Three prespecified native-region geometric witnesses, not clinical crowns', work_directory=str(work), regenerated=len(rows), exact_mesh_array_replays=sum((r['mesh_arrays_identical_to_frozen'] for r in rows)), controls_pass=checks['all_pass'], wall_seconds=time.perf_counter() - start, peak_RSS_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, source_prereg_sha256=sha(root / 'PREREG_D.json'), rows=[dict(key=r['key'], gates=r['gates'], film=r['film']) for r in rows])
    dump(work / 'DEMO_RESULT.json', result)
    dump(root / 'raw/DEMO_LATEST.json', result)
    print(json.dumps({k: v for (k, v) in result.items() if k != 'rows'}, indent=2), flush=True)
if __name__ == '__main__':
    main()
