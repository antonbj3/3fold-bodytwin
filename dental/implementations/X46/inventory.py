"""Direct, one-level metadata inventory; never traverses source/input trees."""
from dental_release.paths import expand as _release_expand
import collections, datetime, hashlib, json, pathlib
ROOT = pathlib.Path(_release_expand('@DENTAL_INPUT_ROOT@/workspace'))
HERE = pathlib.Path(__file__).resolve().parent
CUTOFF = datetime.datetime.fromisoformat('2026-10-02T12:00:00+02:00').timestamp()

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    state_path = ROOT / 'tasks/swarm48/STATE.json'
    state = json.loads(state_path.read_text())
    (rows, rejected) = ([], collections.Counter())
    for (jid, job) in state['jobs'].items():
        p = ROOT / 'results' / jid
        (md, js) = (p / 'RESULTS.md', p / 'results.json')
        if job.get('planner') or '-REVIEW-' in jid:
            rejected['planner_or_prior_review'] += 1
            continue
        if not md.is_file() or not js.is_file():
            rejected['missing_result_pair'] += 1
            continue
        if max(md.stat().st_mtime, js.stat().st_mtime) < CUTOFF:
            rejected['before_cutoff_mtime'] += 1
            continue
        try:
            result = json.loads(js.read_text())
        except Exception:
            rejected['unparseable_results_json'] += 1
            continue
        text = md.read_text(errors='replace')
        outcome = result.get('outcome', {}) if isinstance(result, dict) else {}
        universal = 'for_all_in_declared_domain' in json.dumps(result)
        keys = job.get('source_keys', [])
        cov = any((k.startswith('COV') for k in keys))
        family = 'LIT_CROWN' if 'LIT_CROWN' in keys else 'COV' if cov else keys[0] if keys else 'UNSPECIFIED'
        topfiles = sorted((f.name for f in p.iterdir() if f.is_file()))
        runtime = {}
        for n in ['JOB.json', 'job.json', 'RUN.json', 'runtime.json', 'CLOUD_COLLECTION.json']:
            q = p / n
            if q.is_file() and q.stat().st_size < 200000:
                try:
                    runtime[n] = json.loads(q.read_text())
                except Exception:
                    pass
        rows.append(dict(job_id=jid, source_path=str(p.resolve()), source_alias=str(p), family=family, target_id=job.get('target_id'), category=job.get('category'), source_keys=keys, decision=job.get('decision'), universal=universal, priority=int(universal) * 4 + int('LIT_CROWN' in keys) * 3 + int(cov) * 2, report_sha256=sha(md), result_sha256=sha(js), result_mtime=js.stat().st_mtime, result_keys=list(result) if isinstance(result, dict) else [], outcome=outcome, title=text.splitlines()[0] if text else '', opening=text[:1100], topfiles=topfiles, runtime=runtime))
    rows.sort(key=lambda r: (-r['priority'], -r['result_mtime'], r['job_id']))
    payload = dict(cutoff='2026-10-02T12:00:00+02:00', time_selection='result-pair latest mtime, completion proxy; completion metadata checked separately', state_sha256=sha(state_path), state_jobs=len(state['jobs']), eligible=len(rows), rejected=dict(rejected), rows=rows)
    (HERE / 'INVENTORY.json').write_text(json.dumps(payload, ensure_ascii=False, indent=2) + '\n')
    print(json.dumps({k: v for (k, v) in payload.items() if k != 'rows'}, ensure_ascii=False))
    for r in rows[:90]:
        print(r['job_id'], r['family'], r['category'], 'ALL' if r['universal'] else '', r['title'][:145], r['source_keys'])
if __name__ == '__main__':
    main()
