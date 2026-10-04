"""Regenerate every operative participant and compare scientific output bytes."""
from common import *
from integrity import verify
from generate import command
import subprocess, shutil

def run():
    verify()
    start = time.perf_counter()
    root = DATA / 'recomputation'
    root.mkdir(exist_ok=True)
    receipts = []
    for (name, runner, script, extra, expected, extensions) in [('roofs', 'runner', 'worker.py', [], DATA / 'predictions', {'.npz'}), ('whole', 'whole_runner', 'whole_worker.py', [('whole_inputs', '/whole_inputs'), ('predictions', '/roofs')], DATA / 'whole_predictions', {'.npz', '.stl'}), ('prepared', 'prep_runner', 'prep_worker.py', [('whole_inputs', '/whole_inputs')], DATA / 'prep_predictions', {'.npz', '.stl'})]:
        out = root / name
        out.mkdir(exist_ok=True)
        cmd = command(out, '/runner/' + script)
        cmd[cmd.index(str(DATA / 'runner'))] = str(DATA / runner)
        i = cmd.index('/runtime/bin/python3')
        opts = []
        for (a, b) in extra:
            opts += ['--ro-bind', str(DATA / a), b]
        cmd[i:i] = opts
        st = time.perf_counter()
        with (ROOT / 'raw' / ('recompute_' + name + '.log')).open('w') as f:
            cp = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT)
        if cp.returncode:
            raise RuntimeError('recompute failed ' + name)
        old = {str(p.relative_to(expected)): sha(p) for p in expected.rglob('*') if p.is_file() and p.suffix in extensions}
        new = {str(p.relative_to(out)): sha(p) for p in out.rglob('*') if p.is_file() and p.suffix in extensions}
        diff = sorted((k for k in set(old) | set(new) if old.get(k) != new.get(k)))
        receipts.append(dict(track=name, files=len(old), differences=diff, all_bytes_equal=not diff, seconds=time.perf_counter() - st, expected_hashes=old, actual_hashes=new))
        dump(ROOT / 'raw/RECOMPUTATION_RECEIPT.json', dict(tracks=receipts, total_seconds=time.perf_counter() - start))
        if diff:
            raise RuntimeError('Recomputed bytes differ ' + name)
    verify()
    dump(ROOT / 'raw/RECOMPUTATION_RECEIPT.json', dict(passed=True, tracks=receipts, total_seconds=time.perf_counter() - start))
    shutil.rmtree(root)
    print('all participant recomputation exact')
if __name__ == '__main__':
    run()
