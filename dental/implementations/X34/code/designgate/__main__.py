import argparse
import json
import os
from pathlib import Path
for k in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS']:
    os.environ[k] = '1'
from .gate import check
from .report import html_report

def main():
    p = argparse.ArgumentParser(description='Independent STL geometric design checks with scoped PASS/FAIL/UNKNOWN.')
    sub = p.add_subparsers(dest='command', required=True)
    c = sub.add_parser('check')
    for k in ['prep', 'crown', 'antagonist']:
        c.add_argument('--' + k, required=True, type=Path)
    c.add_argument('--neighbor', type=Path)
    c.add_argument('--material', default='3Y')
    c.add_argument('--contract', type=Path)
    c.add_argument('--units', choices=['mm', 'um'])
    c.add_argument('--contact-axis', nargs=3, type=float)
    c.add_argument('--common-frame', action='store_true')
    c.add_argument('--json', type=Path, default=Path('designgate-report.json'))
    c.add_argument('--html', type=Path, default=Path('designgate-report.html'))
    args = p.parse_args()
    try:
        out = check({k: getattr(args, k) for k in ('prep', 'crown', 'antagonist', 'neighbor')}, args.material, args.contract, args.units, args.contact_axis, args.common_frame, round_version=os.environ.get('DESIGNGATE_ROUND', 'R2'))
        code = 0 if out['verdict'] == 'PASS' else 2 if out['verdict'] == 'FAIL' else 3
    except (ValueError, KeyError, OSError, TypeError, OverflowError) as e:
        out = dict(schema='designgate-0.1', claim_type='capability', verdict='UNKNOWN', input_error=str(e), rules={})
        code = 4
    for f in (args.json, args.html):
        f.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(json.dumps(out, indent=2, allow_nan=False) + '\n')
    args.html.write_text(html_report(out))
    print(json.dumps(dict(verdict=out['verdict'], rules={k: v['status'] for (k, v) in out['rules'].items()}, json=str(args.json), html=str(args.html))))
    return code
if __name__ == '__main__':
    raise SystemExit(main())
