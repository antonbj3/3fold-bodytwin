import sys, os, time, argparse, platform
from pathlib import Path
from util import *
from integrity import Bundle

def main():
    start = time.perf_counter()
    ap = argparse.ArgumentParser()
    ap.add_argument('--phase', choices=['all', 'verify', 'controls', 'score', 'report'], default='all')
    ap.add_argument('--recompute', action='store_true')
    a = ap.parse_args()
    if sys.version_info[:3] != (3, 10, 12):
        raise RuntimeError('CPython3.10.12 required; no fallback')
    bundle = Bundle()
    phase = 'all' if a.phase == 'report' else a.phase
    if phase == 'verify':
        print('VERIFIED', len(bundle.files), 'locked files')
        return
    if a.recompute:
        from sandbox import generate
        dest = PAYLOAD / 'recomputed'
        generate(dest, receipt_name='RECOMPUTATION_PROCESS.json', log_name='recomputation.log')
        fr = read(PAYLOAD / 'predictions/FROZEN_PREDICTIONS.json')
        differences = []
        for (rel, info) in fr['files'].items():
            if not (dest / rel).exists() or sha(dest / rel) != info['sha256']:
                differences.append(rel)
        dump(ROOT / 'raw/RECOMPUTATION_RECEIPT.json', dict(identical=not differences, differences=differences, files=len(fr['files'])))
        if differences:
            raise RuntimeError('regenerated predictions differ; frozen originals retained')
    if phase in ['all', 'controls']:
        from controls import geometry_controls, integrity_controls, sandbox_control, lp_controls
        from cohort_audit import run as cohort_audit
        controls = dict(geometry=geometry_controls(), integrity=integrity_controls(bundle), sandbox=sandbox_control(bundle), cohort=cohort_audit(bundle), lp=lp_controls(bundle))
        dump(ROOT / 'raw/CONTROL_REPORT.json', controls)
    if phase in ['all', 'score']:
        from scorer import score
        from source_audit import run as source_audit
        from preparation import evaluate as prep2
        from preparation_replay import run as prep3
        from crown_track import evaluate as crown
        from source_solid_replay import run as solid
        scores = score()
        source = source_audit(bundle)
        dump(ROOT / 'raw/SOURCE_REFERENCE_AUDIT.json', source)
        if not source['all_pass']:
            raise RuntimeError('independent source reference mismatch')
        r2 = prep2(bundle)
        r3 = prep3(bundle)
        r4 = crown(bundle)
        r5 = solid(bundle, 5)
        r6 = solid(bundle, 6)
        r7 = solid(bundle, 7)
    if phase in ['all', 'report']:
        from report import generate, figure
        out = generate(scores, controls, source, r2, r3, r4, r5, r6, r7, bundle)
        figure(scores, r3, r4, r6, bundle)
        bundle.verify()
    dump(ROOT / 'raw/LAST_RUN.json', dict(argv=sys.argv, seconds=time.perf_counter() - start, python=sys.version, platform=platform.platform(), network='No network calls made by replay; isolation independently tested', threads=2))
    print('DONE', phase, 'seconds', round(time.perf_counter() - start, 2))
if __name__ == '__main__':
    try:
        main()
    except Exception as e:
        dump(ROOT / 'raw/LAST_RUN.json', dict(status='REJECTED', exception=type(e).__name__, reason=str(e), utc=now(), argv=sys.argv))
        raise
