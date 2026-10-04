"""One-command execution/resume, with isolated replay artifacts and finite disk budget."""
from dental_release.paths import expand as _release_expand
import os, pathlib, json, subprocess, sys, argparse, datetime, hashlib, shutil
SRC = pathlib.Path(__file__).resolve().parent
BASE = pathlib.Path(_release_expand('@DENTAL_WORK_ROOT@/X3_crown_benchmark'))

def run():
    ap = argparse.ArgumentParser()
    ap.add_argument('--replay', action='store_true', help='new complete execution, original frozen outputs retained')
    ap.add_argument('--dataset-root', type=pathlib.Path, help='local Teeth3DS root with identical frozen source bytes')
    ap.add_argument('--data-dir', type=pathlib.Path, help='large-output base directory on a data disk')
    ap.add_argument('--output-dir', type=pathlib.Path, help='new report directory; must not exist')
    a = ap.parse_args()
    env = dict(os.environ)
    replay = a.replay or any([a.dataset_root, a.data_dir, a.output_dir])
    base = (a.data_dir or BASE).resolve()
    for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
        env[k] = '1'
    if not replay and all(((SRC / n).exists() for n in ['RESULTS_R1.json', 'RESULTS_R2.json', 'RESULTS_R3.json', 'README_DEMO.md', 'benchmark_figure.png'])):
        env.update(X3_RUN_DIR=str(SRC), X3_DATA_DIR=str(base))
        code = "import crownbench,pathlib,os; crownbench.P=pathlib.Path(os.environ['X3_RUN_DIR']); from verify_delivery import verify; verify(write_receipt=False)"
        subprocess.run([sys.executable, '-c', code], cwd=SRC, env=env, check=True)
        print('READY', SRC / 'README_DEMO.md', flush=True)
        return
    if replay:
        stamp = datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%SZ')
        out = a.output_dir.resolve() if a.output_dir else SRC / 'replays' / stamp
        data = base / 'replays' / stamp
        out.mkdir(parents=True)
        data.mkdir(parents=True)
        for n in ['PREREG_R1.json', 'SPLIT.json', 'DECOMPOSITION.json', 'GRAPH_FEEDBACK.json', 'IMPLEMENTATION_FROZEN_R1.json', 'SOURCE_LITERATURE.json', 'PREREG_NUMERICAL_KERNEL.json', 'PREREG_NUMERICAL_KERNEL_K2.json', 'NUMERICAL_KERNEL_FAILED_V1.json', 'NUMERICAL_KERNEL_VERIFICATION_K2.json', 'ABORTED_RTREE_EVALUATION.json']:
            shutil.copy2(SRC / n, out / n)
            if n == 'NUMERICAL_KERNEL_VERIFICATION_K2.json':
                (out / n).rename(out / 'INHERITED_NUMERICAL_KERNEL_VERIFICATION_K2.json')
        os.symlink(SRC / 'sources', out / 'sources')
        initial = {'lane': 'X3-crown-benchmark-replay', 'phase': 'INITIALIZED', 'updated_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'data_directory': str(data), 'max_intermediate_bytes': 3000000000}
        (out / 'CURRENT_WORK_STATE.json').write_text(json.dumps(initial, indent=2) + '\n')
        (out / 'COMMANDS.md').write_text('# Replay exact commands\n\nOriginal protocol and split reused. This is a reproducibility run, not fresh independent data.\n')
    else:
        out = SRC
        data = base
    if a.dataset_root:
        newroot = a.dataset_root.resolve()
        pr = json.load(open(out / 'PREREG_R1.json'))
        oldroot = pathlib.Path(pr['data']['root'])
        split = json.load(open(out / 'SPLIT.json'))
        for row in split['scans']:
            for key in ['mesh', 'labels']:
                row[key] = str(newroot / pathlib.Path(row[key]).relative_to(oldroot))
                h = hashlib.sha256()
                with open(row[key], 'rb') as src:
                    for block in iter(lambda : src.read(1048576), b''):
                        h.update(block)
                if h.hexdigest() != row[key + '_sha256']:
                    raise RuntimeError('Frozen source bytes differ: ' + row[key])
        split['split_sources'] = {str(newroot / pathlib.Path(k).relative_to(oldroot)): v for (k, v) in split['split_sources'].items()}
        split['relocation'] = {'original_split_sha256': hashlib.sha256((SRC / 'SPLIT.json').read_bytes()).hexdigest(), 'new_root': str(newroot), 'source_bytes_verified': True}
        (out / 'SPLIT.json').write_text(json.dumps(split, indent=2) + '\n')
        pr['relocation'] = {'original_prereg_sha256': hashlib.sha256((SRC / 'PREREG_R1.json').read_bytes()).hexdigest(), 'operators_metrics_thresholds_unchanged': True}
        pr['data']['root'] = str(newroot)
        pr['split_sha256'] = hashlib.sha256((out / 'SPLIT.json').read_bytes()).hexdigest()
        pr['frozen_utc'] = datetime.datetime.now(datetime.timezone.utc).isoformat()
        (out / 'PREREG_R1.json').write_text(json.dumps(pr, indent=2) + '\n')
    env.update(X3_RUN_DIR=str(out), X3_DATA_DIR=str(data), X3_CODE_DIR=str(SRC))
    manifest = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in SRC.glob('*.py')}
    (out / 'RUN_CODE_MANIFEST.json').write_text(json.dumps({'recorded_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'code_sha256': manifest, 'replay': replay}, indent=2) + '\n')
    prefix = "import os,pathlib; import crownbench; crownbench.P=pathlib.Path(os.environ['X3_RUN_DIR']); "

    def stage(name, code, done):
        if (out / done).exists():
            print('REUSE', name, flush=True)
            return
        total = sum((p.stat().st_size for p in base.rglob('*') if p.is_file()))
        if total > 2850000000:
            raise RuntimeError('Lane total disk near3GB. Replay needs budget; no automatic deletion of evidence.')
        cmd = [sys.executable, '-c', prefix + code]
        log = out / (name + '.log')
        print('RUN', name, flush=True)
        with open(log, 'w') as f:
            r = subprocess.run(cmd, cwd=SRC, env=env, stdout=f, stderr=subprocess.STDOUT)
        with open(out / 'COMMANDS.md', 'a') as f:
            f.write('\n' + json.dumps({'argv': cmd, 'env': {'X3_RUN_DIR': str(out), 'X3_DATA_DIR': str(data), 'threads': 1}, 'exit_code': r.returncode, 'log': str(log)}) + '\n')
        if r.returncode:
            raise RuntimeError(f'{name} failed: see {log}')
        if not (out / done).exists():
            raise RuntimeError(name + ' missing completion artifact')
    generate_extra = "; from crownbench import P,put,sha; import json; f=P/('FROZEN_PREDICTIONS_'+rn+'.json'); j=json.load(open(f)); j['implementation_manifest']=[{'path':str(q),'sha256':sha(q)} for q in pathlib.Path(os.environ['X3_CODE_DIR']).glob('*.py')];put(f,j); put(P/'FROZEN_PREDICTIONS.json',j) if rn=='R1' else None"
    stage('TRAIN_REPLAY', 'from prepare import train; train()', 'TRAIN_FIT.json')
    for (rn, group) in [('R1', 'test_R1'), ('R2', 'reserve_R2')]:
        if rn == 'R2':
            stage('FREEZE_R2_REPLAY', 'from freeze_r2 import freeze; freeze()', 'PREREG_R2.json')
        stage('PREP_' + rn + '_REPLAY', f'from prepare import contexts; contexts({group!r})', f'PREPARATION_{group}.json')
        stage('GENERATE_' + rn + '_REPLAY', f'rn={rn!r};from generate import generate; generate({group!r},rn)' + generate_extra, f'FROZEN_PREDICTIONS_{rn}.json')
        if rn == 'R1':
            stage('VERIFY_KERNEL_REPLAY', 'import verify_kernel_k2', 'NUMERICAL_KERNEL_VERIFICATION_K2.json')
        stage('EVALUATE_' + rn + '_REPLAY', f'from evaluate_fast import run_eval; run_eval({rn!r})', f'RESULTS_{rn}.json')
        subprocess.run([sys.executable, '-c', prefix + 'from make_report import report; report()'], cwd=SRC, env=env, check=True)
    stage('R3_COLLAR_REPLAY', 'from round3_collar import run; run()', 'RESULTS_R3.json')
    subprocess.run([sys.executable, '-c', prefix + 'from make_report import report; report(); from finish_delivery import finish; finish()'], cwd=SRC, env=env, check=True)
    from_script = [sys.executable, '-s', str(SRC / 'plot_figure.py')]
    subprocess.run(from_script, cwd=SRC, env=env, check=True)
    verify = [sys.executable, '-c', prefix + 'from verify_delivery import verify; verify()']
    subprocess.run(verify, cwd=SRC, env=env, check=True)
    print('READY', out / 'README_DEMO.md', flush=True)
if __name__ == '__main__':
    run()
