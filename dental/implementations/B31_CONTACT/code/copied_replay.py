"""Only lane-local code/contracts copied. Upstream trees stay read-only in place."""
from dental_release.paths import expand as _release_expand
from common import *
import shutil, subprocess, time, resource

def scientific(j):
    out = dict(decision=j['decision'], sufficiency=dict(identity=j['summary_sufficiency']['summary_identity_error'], downstream=j['summary_sufficiency']['downstream_difference_mm2']), rows={})
    for phase in ['R1', 'R2', 'R3', 'R4']:
        rr = j['rounds'][phase]
        arr = []
        for r in rr['rows']:
            if phase == 'R1':
                arr.append(dict(key=r['key'], gate=r['digital_case_gate'], values={n: dict(hash=x['sha256'], contact=x['contact']['symdiff_mm2'], wall=x['wall_and_removal']['wall_lower_mm'], boxes=x['boxes']) for (n, x) in r['outputs'].items()}))
            elif phase == 'R2':
                arr.append(dict(key=r['key'], hash=r['artifact_sha256'], nominal=r['nominal'], impossible=r['all_removal_50pct_impossible'], new_manufacture=r['full_source_mask_answer']))
            elif phase == 'R3':
                arr.append(dict(key=r['key'], values={n: dict(hash=x['sha256'], contact=x['contact']['symdiff_mm2'], wall=x['wall_and_removal']['wall_lower_mm'], order=x['wall_and_removal']['local_vertical_interval_order_pass']) for (n, x) in r['outputs'].items()}))
            else:
                arr.append(dict(key=r['key'], hash=r['artifact_sha256'], local_cap=r['uniform_local_cap_mm'], contact=r['uniform_contact']['symdiff_mm2'], missing_without_cap=r['cap_independent_missing_source_area_mm2'], nominal=r['nominal50pct_vs_local_proxy'], pose_failure=r['observed_pose_ratio_failure'], gate=r['case_gate'], controls=r['controls']))
        if phase == 'R1':
            for r in arr:
                for x in r['values'].values():
                    x['boxes'] = {b: {k: v for (k, v) in z.items() if k != 'seconds'} for (b, z) in x['boxes'].items()}
        out['rows'][phase] = arr
    return out

def run():
    start = time.perf_counter()
    base = Path(_release_expand('@DENTAL_EXTERNAL_ROOT@/storage/tmp_dental_audit/X95-spatial-contact-repair'))
    base.mkdir(parents=True, exist_ok=True)
    dest = base / ('replay_' + datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%SZ'))
    dest.mkdir()
    for pat in ['*.json', '*.sha256', '*.md', 'run_all.sh']:
        for p in ROOT.glob(pat):
            shutil.copy2(p, dest / p.name)
    shutil.copytree(ROOT / 'code', dest / 'code')
    (dest / 'raw').mkdir()
    env = os.environ.copy()
    env.update(X95_DATA_ROOT=str(dest / 'data'), X95_UPSTREAM_ROOT=str(UP), OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', NUMEXPR_NUM_THREADS='1', PYTHONDONTWRITEBYTECODE='1', MPLBACKEND='Agg')
    print('copied replay', dest, flush=True)
    with (dest / 'RUN.log').open('w') as log:
        cp = subprocess.run(['bash', str(dest / 'run_all.sh')], cwd=dest, env=env, stdout=log, stderr=subprocess.STDOUT)
    original = read(ROOT / 'results.json')
    copied = read(dest / 'results.json') if cp.returncode == 0 else None
    same = copied is not None and scientific(original) == scientific(copied)
    receipt = dict(audit_root=str(dest), command='X95_DATA_ROOT=' + str(dest / 'data') + ' X95_UPSTREAM_ROOT=' + str(UP) + ' bash ' + str(dest / 'run_all.sh'), exit_code=cp.returncode, scientific_output_exact_match=same, comparison='All gates, R1/R2/R3/R4 artifact hashes, contact/wall/cap/missing-region values, continuum bracket endpoints and sufficiency; runtime/path excluded', seconds=time.perf_counter() - start, threads=1, source_trees_copied=False, source_input_lock_sha256=sha(ROOT / 'INPUT_LOCK_R1.json'), copied_results_sha256=sha(dest / 'results.json') if copied else None, run_log_sha256=sha(dest / 'RUN.log'), physical_validation='NOT_RUN')
    dump(ROOT / 'raw/COPIED_REPLAY.json', receipt)
    dump(dest / 'REPLAY_COMPARISON.json', receipt)
    if not same:
        dump(ROOT / 'raw/REPLAY_FAILURE.json', dict(receipt=receipt, original_scientific=scientific(original), copied_scientific=scientific(copied) if copied else None))
        raise RuntimeError('Copied replay failed; see audit log and REPLAY_FAILURE.json')
    print(json.dumps(receipt), flush=True)
if __name__ == '__main__':
    run()
