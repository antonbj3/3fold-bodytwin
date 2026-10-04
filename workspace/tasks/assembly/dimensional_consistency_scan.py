"""Check every stated record against its own arithmetic: does A times B equal the C it also states?

WHERE THIS CAME FROM. Four constraint-net edges said "At load=1000 N" while stating a contact area
and a peak pressure whose product is 132 N, 168 N, 115.5 N and 133 N. A peak pressure acting over
the contact area cannot carry less than the total load, so the label was unsupported -- and the
arithmetic that showed it needed no source, no network and no judgement. Anton's point is that this
cannot be the only one.

It is not a special case. A record that states several quantities WITH UNITS constrains itself:
area times pressure is a force, stiffness times length is a force, concentration times volume is an
amount, flux times area times time is an amount. Where a record states all three legs of such a
triple, the third is checkable against the other two. This file does that check over every record we
hold, and it needs nothing the bunnies have to be given -- which is why it belongs here and not in a
brief. They run against the web with text excerpts; the repository is local.

WHAT IT DOES NOT DO. It does not decide which leg is wrong. A disagreement means the record cannot
be read as stated; whether the load label, the area or the pressure is the error needs the source.
It also does not flag a record whose legs are consistent but wrong. Dimensional self-consistency is
necessary, not sufficient. PENDING_INDEPENDENT_REVIEW.
"""
from __future__ import annotations

import json
import pathlib
import re

W = pathlib.Path('')
NET = W / 'CONSTRAINT_NETS.json'

# Units as (dimension, factor to the SI-ish base used below).
UNITS = {
    'N': ('force', 1.0), 'mN': ('force', 1e-3), 'uN': ('force', 1e-6), 'µN': ('force', 1e-6),
    'nN': ('force', 1e-9), 'kN': ('force', 1e3),
    'Pa': ('pressure', 1.0), 'kPa': ('pressure', 1e3), 'MPa': ('pressure', 1e6),
    'GPa': ('pressure', 1e9), 'mmHg': ('pressure', 133.322),
    'mm2': ('area', 1e-6), 'mm²': ('area', 1e-6), 'cm2': ('area', 1e-4), 'cm²': ('area', 1e-4),
    'um2': ('area', 1e-12), 'µm2': ('area', 1e-12), 'um²': ('area', 1e-12), 'µm²': ('area', 1e-12),
    'm2': ('area', 1.0), 'm²': ('area', 1.0),
    'mm': ('length', 1e-3), 'cm': ('length', 1e-2), 'um': ('length', 1e-6), 'µm': ('length', 1e-6),
    'nm': ('length', 1e-9), 'm': ('length', 1.0),
    'N/mm': ('stiffness', 1e3), 'N/m': ('stiffness', 1.0),
    'uM': ('concentration', 1e-3), 'µM': ('concentration', 1e-3), 'mM': ('concentration', 1.0),
    'nM': ('concentration', 1e-6), 'M': ('concentration', 1e3),
    'L': ('volume', 1.0), 'mL': ('volume', 1e-3), 'uL': ('volume', 1e-6), 'µL': ('volume', 1e-6),
}
# (left, right) -> product dimension
PRODUCTS = {
    ('pressure', 'area'): 'force',
    ('stiffness', 'length'): 'force',
    ('concentration', 'volume'): 'amount',
}
NUM = re.compile(r'(?<![\w.])(-?\d+(?:[.,]\d+)?(?:[eE][+-]?\d+)?)\s*([A-Za-zµ°/²³]+[0-9²³]?)')
TOL = 0.05          # 5 percent: looser than rounding, tighter than a unit error


def quantities(text: str) -> list[tuple[float, str, str]]:
    """Every number with a unit we know, as (value in base units, dimension, literal)."""
    out = []
    for val, unit in NUM.findall(text):
        u = UNITS.get(unit)
        if not u:
            continue
        try:
            v = float(val.replace(',', '.'))
        except ValueError:
            continue
        out.append((v * u[1], u[0], f'{val} {unit}'))
    return out


def check(text: str) -> list[dict]:
    """Every product triple the text states all three legs of, with its disagreement."""
    qs = quantities(text)
    findings = []
    for i, (a, da, la) in enumerate(qs):
        for j, (b, db, lb) in enumerate(qs):
            if i == j:
                continue
            target = PRODUCTS.get((da, db))
            if not target:
                continue
            implied = a * b
            for k, (c, dc, lc) in enumerate(qs):
                if dc != target or k in (i, j):
                    continue
                if implied == 0 or c == 0:
                    continue
                ratio = implied / c
                if abs(ratio - 1.0) > TOL:
                    findings.append({'left': la, 'right': lb, 'stated': lc,
                                     'implied_over_stated': round(ratio, 6),
                                     'relation': f'{da} x {db} -> {target}'})
    # keep the worst per stated leg so one record does not emit the same thing many times
    best = {}
    for f in findings:
        key = (f['stated'], f['relation'])
        if key not in best or abs(f['implied_over_stated'] - 1) > abs(best[key]['implied_over_stated'] - 1):
            best[key] = f
    return list(best.values())


def main() -> None:
    net = json.loads(NET.read_text())['bodytwin']['tissue_constraint_net']
    rows = []
    for e in net['edges']:
        text = ' '.join(str(e.get(k, '')) for k in ('constraint', 'evidence'))
        for f in check(text):
            rows.append(dict(f, record=e['id'], status=e.get('status'), kind='net_edge'))
    docs = sorted(pathlib.Path('source_repository/docs').glob('MECHANISM_*.md'))
    for p in docs:
        try:
            lines = p.read_text(errors='replace').splitlines()
        except OSError:
            continue
        for n, line in enumerate(lines, 1):
            if len(line) > 400:
                continue
            for f in check(line):
                rows.append(dict(f, record=f'{p.name}:L{n}', status=None, kind='detail_layer_line'))
    rows.sort(key=lambda r: -abs(r['implied_over_stated'] - 1))
    out = W / 'results/ASSEMBLY_DIMENSIONAL_SCAN'
    out.mkdir(parents=True, exist_ok=True)
    (out / 'FINDINGS.json').write_text(json.dumps(
        {'tolerance': TOL, 'n_findings': len(rows), 'documents_scanned': len(docs),
         'rows': rows, 'review_state': 'PENDING_INDEPENDENT_REVIEW'}, indent=1, ensure_ascii=False))
    edges = [r for r in rows if r['kind'] == 'net_edge']
    print(f'documents scanned: {len(docs)}   findings: {len(rows)}   '
          f'of which in net edges: {len(edges)}')
    for r in rows[:14]:
        print(f"  {str(r['record'])[:46]:<48} {r['left']:>12} x {r['right']:>12} "
              f"-> {r['stated']:>12}  factor {r['implied_over_stated']:.4g}")
    print('wrote', out / 'FINDINGS.json')


if __name__ == '__main__':
    main()
