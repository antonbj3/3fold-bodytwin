from common_r3 import *
import subprocess, time, collections, sys

def run(tag):
    start = time.perf_counter()
    records = read(DATA / tag / 'RECORDS.json')
    rows = []
    assert (ROOT / f'FROZEN_PREDICTIONS_{tag}.json').exists()
    if not (ROOT / f'FROZEN_SCORER_{tag}.json').exists():
        freeze(ROOT / f'FROZEN_SCORER_{tag}.json', dict(files={str(p): sha(p) for p in [ROOT / 'code/score_worker.py', ROOT / 'code/score_round.py', ROOT / 'code/memory_height.py', OLD / 'code/score.py']}, source='Unchanged R2 full surface scorer and v6 spatial thresholds; only bounded candidate-list allocation'))
    for rec in records:
        key = rec['key']
        method = rec['participant']
        ident = key + '__' + method
        dest = ROOT / 'raw' / f'{tag}_SCORE_{ident}.json'
        if not dest.exists():
            cmd = ['/usr/bin/time', '-v', '-o', str(ROOT / 'raw' / f'{tag}_SCORE_MEMORY_{ident}.txt'), PYTHON, str(ROOT / 'code/score_worker.py'), tag, key, method]
            dump(ROOT / 'raw' / f'{tag}_SCORE_CMD_{ident}.json', dict(argv=cmd))
            with (ROOT / 'raw' / f'{tag}_score.log').open('a') as f:
                p = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT)
            if p.returncode:
                raise RuntimeError('scorer process failed ' + ident)
        rows.append(read(dest))
        state(tag + '_SCORING', str(len(rows)) + '/' + str(len(records)) + ' records scored', 'Complete full frozen denominator')
    summary = {}
    for method in sorted({r['participant'] for r in records}):
        summary[method] = {}
        for fam in ['molar_crown', 'premolar_crown', 'anterior_crown']:
            rr = [r for r in rows if r['participant'] == method and r['family'] == fam]
            ss = [r for r in rr if r['status'] == 'SCORED']
            summary[method][fam] = dict(requested=len(rr), scored=len(ss), anatomy_pass=sum((r['anatomy_pass'] for r in ss)), functional_pass_sampled=sum((r['functional_pass_sampled'] for r in ss)), functional_pass_cover=sum((r['functional_pass_cover'] for r in ss)), joint_pass=sum((r['anatomy_pass'] and r['functional_pass_sampled'] for r in ss)), positive_native_restored=sum((r['positive_native_band_restored'] for r in ss)), median_p95_mm=float(np.median([r['reconstruction_p95_mm'] for r in ss])) if ss else None, gate_counts={k: sum((r['gates'][k] for r in ss)) for k in ss[0]['gates']} if ss else {})
    result = dict(round=tag, claim_type='capability', external_referent=REFERENT, rows=rows, summary=summary, dropout=dict(requested=len(rows), rejected=sum((r['status'] != 'SCORED' for r in rows)), fraction=sum((r['status'] != 'SCORED' for r in rows)) / len(rows), reasons=dict(collections.Counter((r.get('reason', '') for r in rows if r['status'] != 'SCORED')))), seconds=time.perf_counter() - start, peak_worker_rss_MiB=max((r.get('score_peak_rss_MiB', 0) for r in rows)), resolution='PER_POINT distances -> PER_SURFACE_REGION contact -> PER_TOOTH conjunction; across-case counts retain family and split', rigorous_floating_enclosure='MISSING')
    dump(ROOT / 'rounds' / f'{tag}.json', result)
    text = json.dumps(clean(summary), indent=2)
    (ROOT / 'history' / f'HANDOFF_{tag}.md').write_text(text + '\nFrozen full anatomy and functional gates retained. No physical contact validation.\n')
    (ROOT / 'HANDOFF.md').write_text('Round ' + tag + ' scored.\n' + text + '\nNext: diagnose residuals and change load-bearing operation.\n')
    state(tag + '_DECIDED', summary, 'Diagnose residuals and preregister changed construction')
    print(text)
    return result
if __name__ == '__main__':
    run(sys.argv[1])
