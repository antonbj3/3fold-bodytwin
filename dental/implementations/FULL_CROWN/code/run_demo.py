from fc_common import *
import argparse, subprocess, resource

def run():
    parser = argparse.ArgumentParser()
    parser.add_argument('--regenerate', action='store_true')
    args = parser.parse_args()
    st = time.perf_counter()
    if (ROOT / 'RELEASE_CODE.json').exists():
        for (rel, meta) in read(ROOT / 'RELEASE_CODE.json')['files'].items():
            if sha(ROOT / rel) != meta['sha256']:
                raise ValueError('Code/protocol drift: ' + rel)
    if args.regenerate:
        import regenerate
        regenerate.run()
    from score_predictions import run as score
    for tag in ['R2', 'R3', 'R4', 'R5', 'R6', 'R7']:
        score(tag)
    from validate_demo import run as validate
    val = validate()
    unexpected = [r for r in val['checks'] if not r['pass_gate'] and (not r['name'].startswith('STL_closed_'))]
    if unexpected:
        raise RuntimeError('Unexpected validation failures:' + str(unexpected))
    from export_lab import run as export
    if not export()['all_pass']:
        raise RuntimeError('Selected3MF export failed')
    from report_demo import run as report
    result = report()
    dump(ROOT / 'raw/DEMO_RUN.json', dict(seconds=time.perf_counter() - st, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, regenerated=args.regenerate, anatomical_capability_passed=False, retained_STL_failures=[r['name'] for r in val['checks'] if not r['pass_gate']], status='REPRODUCED_NEGATIVE_ANATOMY_WITH_REVIEWABLE_SPECIMENS'))
    print('Demo complete: anatomical gate remains FAILED; see README_DEMO.md and results.json')
if __name__ == '__main__':
    run()
