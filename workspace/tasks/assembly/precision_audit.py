"""How many significant figures does each number we carry actually have?

FOUND BY THE SWARM, verified here. A creative job on edge H-E29 took the number
2.220261437908497 from our own net and reconstructed it as exactly 3397/1530 -- the nearest double
to a ratio of two four-significant-figure readings, with zero error across denominators from 1e4 to
1e7. So the number carries FOUR significant figures and the remaining eleven digits are floating
point representation, not measurement. Quoting it with fifteen is false precision.

That is a class, not an instance. A chain that divides measured quantities produces long decimals
from short measurements, and once a long decimal is written down it reads as precise. The same
mistake appeared in my own work twice today: a residual reported as -0.000000 when the inputs
supported 3e-5, and a bracket quoted to ten digits that was 1.81 float32 ULPs wide.

WHAT THIS DOES. For every number in the net and in the assembly outputs, find the smallest
denominator whose fraction reproduces the number EXACTLY as a double. A small denominator means the
number is a ratio of short readings and its honest precision is the shorter of the two. A number
that needs a huge denominator is either a genuine long computation or an irrational, and is left
alone.

WHAT IT DOES NOT DO. It cannot tell a four-figure measurement from a four-figure assumption, and a
number with few significant figures is not wrong -- it is just not precise. The output is a precision
label per number, not a defect list. PENDING_INDEPENDENT_REVIEW.
"""
from __future__ import annotations

import json
import math
import pathlib
from fractions import Fraction

W = pathlib.Path('')


def honest_digits(x: float) -> dict | None:
    """Smallest exact rational, and the significant figures that implies."""
    if x == 0 or not math.isfinite(x) or abs(x) > 1e12:
        return None
    for limit in (10, 100, 1000, 10**4, 10**5, 10**6):
        f = Fraction(x).limit_denominator(limit)
        if float(f) == x:
            # significant figures = the longer of numerator and denominator digit counts
            digits = max(len(str(abs(f.numerator))), len(str(abs(f.denominator))))
            return {'exact_rational': f'{f.numerator}/{f.denominator}',
                    'denominator_limit': limit, 'significant_figures': digits,
                    'printed_digits': len(repr(x).replace('-', '').replace('.', '').lstrip('0')),
                    'overstated_by': len(repr(x).replace('-', '').replace('.', '').lstrip('0')) - digits}
    return None


def walk(obj, path=''):
    if isinstance(obj, dict):
        for k, v in obj.items():
            yield from walk(v, f'{path}/{k}')
    elif isinstance(obj, list):
        for i, v in enumerate(obj):
            yield from walk(v, f'{path}/{i}')
    elif isinstance(obj, float):
        yield path, obj


def main() -> None:
    rows = []
    sources = [W / 'CONSTRAINT_NETS.json'] + sorted(W.glob('results/ASSEMBLY_*/*.json'))
    for f in sources:
        try:
            doc = json.loads(f.read_text())
        except Exception:
            continue
        tag = f.name if f.name == 'CONSTRAINT_NETS.json' else f.parent.name
        for path, val in walk(doc):
            h = honest_digits(val)
            # The first filter caught 24 026 numbers and they were almost all 1/6, 4/3 and 8/3
            # printed as long floats. A simple fraction printed long is not false precision; it is
            # an exact value in a format that looks imprecise, which is the opposite problem. The
            # dangerous signature is the one the swarm found: TWO SHORT READINGS DIVIDED, where both
            # numerator and denominator are real measurements of three or four figures. That shows
            # up as a denominator in the hundreds or thousands, not a single digit.
            if not h:
                continue
            den = int(h['exact_rational'].split('/')[1])
            if den >= 50 and 3 <= h['significant_figures'] <= 6 and h['overstated_by'] >= 8:
                rows.append(dict(h, source=tag, field=path[-70:], value=val))
    rows.sort(key=lambda r: -r['overstated_by'])
    out = W / 'results/ASSEMBLY_PRECISION_AUDIT'
    out.mkdir(parents=True, exist_ok=True)
    (out / 'AUDIT.json').write_text(json.dumps(
        {'found_by': 'swarm job on edge H-E29, verified here', 'n': len(rows), 'rows': rows,
         'review_state': 'PENDING_INDEPENDENT_REVIEW'}, indent=1, ensure_ascii=False))
    print(f'numbers with more printed digits than they carry: {len(rows)}')
    print(f'{"source":<30} {"digits":>8} {"printed":>10}  exact ratio')
    for r in rows[:14]:
        print(f'  {r["source"][:28]:<30} {r["significant_figures"]:>6} {r["printed_digits"]:>10}  '
              f'{r["exact_rational"][:24]:<26} {r["field"][-34:]}')
    print('wrote', out / 'AUDIT.json')


if __name__ == '__main__':
    main()
