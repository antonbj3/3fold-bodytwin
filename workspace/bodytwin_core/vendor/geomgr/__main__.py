"""CLI: python3 -m geomgr <command> ...

  demo        --subject z001 [--obs S6|S2|S4|S1|S3] [--out DIR] [--data DIR] [--lex-only]
  instantiate --features L_mech=431,CCD=128 [--landmarks lm.json] [--bone femur_r|tibia_r] --out DIR
  batch       --stage numeric|local --n 100 --data DIR --out DIR
  adversarial [--subject z001]
  evaluate    ... (see geomgr.evaluate)   summarize --in DIR... --out DIR   data --out DIR
  outlier     --data DIR --out results.json
  pipeline    --subject z001 [--obs S6] [--stage numeric|local|all] [--draws 50] [--register] [--out DIR]   (N7c)
  pipeline-batch --subjects z009,z013,... --stage numeric --out DIR                                     (N7c, cloud)
"""
import argparse
import json
import sys
from pathlib import Path

DATA = 'external_media'


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv:
        print(__doc__)
        return 2
    cmd, rest = argv[0], argv[1:]
    if cmd == 'evaluate':
        from . import evaluate
        return evaluate.main(rest)
    if cmd == 'summarize':
        from . import summarize
        return summarize.main(rest)
    if cmd == 'data':
        from . import data
        return data.main(rest)
    if cmd in ('pipeline', 'pipeline-batch'):
        from . import pipeline
        return pipeline.main(rest)
    if cmd in ('outlier', 'landmark-outlier'):
        from . import outlier
        return outlier.main(rest)
    ap = argparse.ArgumentParser(prog=f'geomgr {cmd}')
    ap.add_argument('--data', default=DATA)
    ap.add_argument('--out', default=None)
    if cmd == 'demo':
        ap.add_argument('--subject', default='z001')
        ap.add_argument('--obs', default='S6')
        ap.add_argument('--lex-only', action='store_true')
        a = ap.parse_args(rest)
        from . import demo
        s = demo.run(a.subject, a.obs, a.out, a.data, validate_full=not a.lex_only)
        print(json.dumps(dict(accepted=s['accepted'], error_vs_truth=s['error_vs_truth'],
                              opensim_ok=s['exports']['opensim']['load']['ok'], sdf_ok=s['sdf'].get('ok'),
                              timings=s['timings']), indent=1))
        return 0
    if cmd == 'adversarial':
        ap.add_argument('--subject', default='z001')
        a = ap.parse_args(rest)
        from . import adversarial
        r = adversarial.run(a.subject, a.out, a.data)
        print(json.dumps({k: v for k, v in r.items() if k not in ('cases', 'fold_scan')}, indent=1))
        return 0
    if cmd == 'batch':
        ap.add_argument('--stage', choices=('numeric', 'local'), required=True)
        ap.add_argument('--n', type=int, default=100)
        ap.add_argument('--kappa', default=None, help='calibration.json')
        a = ap.parse_args(rest)
        from . import batch
        if a.stage == 'numeric':
            k = None
            if a.kappa:
                kk = json.loads(Path(a.kappa).read_text())['kappa']['femur_r']
                k = {sc: v.get('cov3_placed', v.get('cov3')) for sc, v in kk.items()}
            r = batch.stage_numeric(a.data, a.out, a.n, kappa=k)
        else:
            r = batch.stage_local(a.data, a.out, a.n)
        print(json.dumps({k: v for k, v in r.items() if k != 'rows'}, indent=1, default=float))
        return 0
    if cmd == 'instantiate':
        ap.add_argument('--bone', default='femur_r')
        ap.add_argument('--features', default='')
        ap.add_argument('--landmarks', default=None, help='JSON {name: [x,y,z]} in mm')
        ap.add_argument('--subject', default='cli')
        a = ap.parse_args(rest)
        from . import core as C
        from .api import GeometryManager
        from . import export as X
        gm = GeometryManager(a.data, bone=a.bone, cache_dir='external_media')
        f = {kv.split('=')[0]: float(kv.split('=')[1]) for kv in a.features.split(',') if kv}
        lm = json.loads(Path(a.landmarks).read_text()) if a.landmarks else None
        r = gm.instantiate(features=f, feature_units={k: C.feature_units(k) for k in f}, landmarks=lm, subject=a.subject)
        out = Path(a.out or '.')
        out.mkdir(parents=True, exist_ok=True)
        if r.instance is not None:
            X.to_json(r, out / f'{a.subject}_instance.json')
        print(json.dumps(r.summary(), indent=1, default=float))
        return 0 if r.accepted else 1
    print(__doc__)
    return 2


if __name__ == '__main__':
    sys.exit(main())
