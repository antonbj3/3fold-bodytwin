"""Evaluate a third-party natural-tooth generator on the frozen public contexts.
Example: python3 benchmark_method.py --name my_method --command 'python3 /path/method.py {input} {output}'
The subprocess receives only masked-context and predicted-output paths. Python methods
may import crownbench primitives, but must not read original/hidden data. This is an
information contract, not an OS-level sandbox against malicious disk access.
"""
from crownbench import *
from evaluate import metric, balanced, METRICS
import argparse, subprocess, shlex, collections, re

def run(name, command, round_name):
    if not re.fullmatch('[a-zA-Z0-9_-]+', name):
        raise ValueError('method name must be simple')
    group = 'test_R1' if round_name == 'R1' else 'reserve_R2'
    public = json.load(open(DATA / 'public' / f'{group}_INDEX.json'))
    base = P / 'submissions' / name / round_name
    if (base / 'FROZEN.json').exists():
        raise RuntimeError('Frozen submission already exists; choose new method version name.')
    base.mkdir(parents=True, exist_ok=True)
    generated = []
    tokens = shlex.split(command)
    for r in public:
        out = DATA / 'submissions' / name / round_name / f"{r['key']}.npz"
        row = {'key': r['key'], 'case': r['case'], 'jaw': r['jaw'], 'fdi': r['fdi'], 'eligible': r['eligible']}
        if r['eligible']:
            out.parent.mkdir(parents=True, exist_ok=True)
            args = [a.format(input=r['context'], output=str(out)) for a in tokens]
            ts = time.perf_counter()
            proc = subprocess.run(args, capture_output=True, text=True, timeout=600)
            log = base / (r['key'] + '.log')
            log.write_text(proc.stdout + '\nSTDERR\n' + proc.stderr)
            row.update(exit_code=proc.returncode, seconds=time.perf_counter() - ts, argv=args)
            if proc.returncode == 0 and out.exists():
                (v, f) = read_surface(out)
                if len(v) and len(f) and np.isfinite(v).all():
                    row.update(path=str(out), sha256=sha(out))
            if not row.get('path'):
                row['failure'] = 'NO_VALID_OUTPUT'
        else:
            row['failure'] = r['reason']
        generated.append(row)
    put(base / 'FROZEN.json', {'frozen_utc': now(), 'split_sha256': sha(P / 'SPLIT.json'), 'command': command, 'predictions': generated, 'before_hidden_evaluation': True})
    hidden = {r['key']: r for r in json.load(open(DATA / 'hidden' / f'{group}_INDEX.json'))}
    pose = {r['key']: r for r in json.load(open(DATA / 'predictions' / round_name / 'INDEX.json'))}
    rows = []
    for r in generated:
        row = {k: r[k] for k in ('key', 'case', 'jaw', 'fdi', 'eligible')}
        row['metrics'] = {}
        if r.get('path'):
            h = hidden[r['key']]
            assert sha(r['path']) == r['sha256']
            assert sha(h['reference']) == h['reference_sha256']
            with np.load(h['context']) as d:
                cv = d['vertices']
                cf = d['faces']
                cl = d['labels']
            k = r['fdi']
            assert k not in np.unique(cl)
            neighbors = [select(cv, cf, cl, k - 1), select(cv, cf, cl, k + 1)]
            n = np.array(pose[r['key']]['frame'])[:, 2]
            row['metrics'][name] = metric(read_surface(r['path']), read_surface(h['reference']), neighbors, n, r['key'])
        else:
            row['failure'] = r['failure']
        rows.append(row)
    put(base / 'RAW.json', rows)
    summary = {'method': name, 'scheduled': len(public), 'completed': sum((name in r['metrics'] for r in rows)), 'coverage': sum((name in r['metrics'] for r in rows)) / len(public), 'metrics_patient_balanced': {q: float(np.mean(list(balanced(rows, name, q).values()))) if balanced(rows, name, q) else None for q in METRICS}, 'by_type': {}, 'external_referent': json.load(open(P / 'PREREG_R1.json'))['external_referent'], 'frozen_sha256': sha(base / 'FROZEN.json'), 'antagonist_contact': 'UNKNOWN_NO_VERIFIED_BITE_TRANSFORM', 'review_state': 'PENDING_INDEPENDENT_REVIEW'}
    for (typ, idx) in [('incisor', 2), ('canine', 3), ('premolar', 5), ('molar', 6)]:
        sub = [r for r in rows if r['fdi'] % 10 == idx]
        summary['by_type'][typ] = {q: float(np.mean(list(balanced(sub, name, q).values()))) if balanced(sub, name, q) else None for q in METRICS}
    put(base / 'RESULTS.json', summary)
    print(json.dumps(summary, indent=2))
if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--name', required=True)
    ap.add_argument('--command', required=True)
    ap.add_argument('--round', default='R1', choices=['R1', 'R2'])
    a = ap.parse_args()
    run(a.name, a.command, a.round)
