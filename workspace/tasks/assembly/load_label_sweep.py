#!/usr/bin/env python3
"""Every place in the nets where a pressure is stated next to a load, and whether the two agree.

Why this exists. Three separate findings on 2026-10-04 turned out to be the same fault, each found by
hand: the meniscus rows multiply area by peak pressure to 3450 N, 3120 N and 132 N against a stated
1000 N (factor 26 between them, and the 132 N row is impossible because it needs a peak below the
mean); the disc comparison puts our multiplier 2.09x above the in-vivo facit, but the model's own
reference load is sigma*A = 0.9045*1800 = 1628.1 N while our multiplier is consistent with the measured
pressure at 636-900 N; and the IOL edge named a field that counts something else entirely. A fault
found three times by hand is a family, not three accidents, so it gets enumerated.

What it does. For each edge in each net it pulls out the pressures, the areas and the forces that the
edge's own text states, forms every area*pressure product, and reports the ratio to each stated force.
It decides nothing: a ratio above 1 is expected, because peak pressure times area must exceed the load
that the mean pressure carries. The two readings worth looking at are a ratio BELOW 1, which is
impossible for a peak, and two ratios on the same edge that differ by more than the peak-to-mean factor
a contact can plausibly have.

Units are read from the text, not assumed: mm^2 with MPa gives newtons directly, and anything else is
converted or skipped with the reason named.

Two traps the first version fell into, both fixed here and both worth naming. It scanned the fields
this coordinator had itself written onto the edges (`load_label_unsupported`, which carries an
`implied_load_N`), so every edge reported a ratio of exactly 1.000 against my own note -- a tool
reading its own output, the same fault `precision_audit.py` had. And it paired the dispersions with
each other: "+/- 8 mm^2" times "+/- 0.2 MPa" is 1.6 N, which against a stated 1000 N looks like a
ratio of 0.0016 and means nothing. Dispersions are now dropped, and only the edge's own stated
quantities are read.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

W = Path(__file__).resolve().parents[2]
NETS = ('CONSTRAINT_NETS.json', 'data/CONSTRAINT_NET_TISSUE.json')

# A number immediately followed by a unit. The lookbehind keeps exponents and ids out.
NUM = r'(?<![\d.A-Za-z])(\d+(?:[.,]\d+)?)'
PRESSURE = re.compile(NUM + r'\s*(MPa|kPa|Pa)\b')
AREA = re.compile(NUM + r'\s*(mm\^?2|mm²|cm\^?2|cm²|m\^?2|m²)\b')
FORCE = re.compile(NUM + r'\s*(kN|N)\b(?![a-z])')

TO_MPA = {'MPa': 1.0, 'kPa': 1e-3, 'Pa': 1e-6}
TO_MM2 = {'mm^2': 1.0, 'mm2': 1.0, 'mm²': 1.0,
          'cm^2': 100.0, 'cm2': 100.0, 'cm²': 100.0,
          'm^2': 1e6, 'm2': 1e6, 'm²': 1e6}
TO_N = {'N': 1.0, 'kN': 1e3}

# Peak over mean for a contact pressure distribution. Below 1 is impossible; the upper end is the
# loosest value the meniscus rows support (3.45), kept as a flag threshold rather than a constant.
PEAK_OVER_MEAN_MAX = 3.5


def nums(pat, text, scale):
    out = []
    for m in pat.finditer(text):
        try:
            v = float(m.group(1).replace(',', '.'))
        except ValueError:
            continue
        u = m.group(2)
        f = scale.get(u)
        if f is None:
            continue
        out.append((v * f, m.group(0)))
    return out


def iter_edges(obj):
    if isinstance(obj, dict):
        if isinstance(obj.get('edges'), list):
            for e in obj['edges']:
                if isinstance(e, dict):
                    yield e
        for v in obj.values():
            yield from iter_edges(v)
    elif isinstance(obj, list):
        for v in obj:
            yield from iter_edges(v)


def main() -> int:
    rows = []
    for net in NETS:
        p = W / net
        if not p.exists():
            continue
        for e in iter_edges(json.loads(p.read_text())):
            # Only the edge's own claim and its locator. The annotation fields are this
            # coordinator's writing, and reading them back is a feedback loop, not evidence.
            text = ' '.join(str(e.get(k, '')) for k in ('constraint', 'evidence'))
            # A value written as a dispersion is not a value. "110 +/- 8 mm^2" states one area.
            text = re.sub(r'(?:\+/-|±|\+-)\s*\d+(?:[.,]\d+)?\s*'
                          r'(?:MPa|kPa|Pa|mm\^?2|mm²|cm\^?2|cm²|m\^?2|m²|kN|N)\b', ' ', text)
            pres = nums(PRESSURE, text, TO_MPA)
            area = nums(AREA, text, TO_MM2)
            force = nums(FORCE, text, TO_N)
            if not (pres and force):
                continue
            # MPa * mm^2 = N exactly, so no further factor enters here.
            # The same number may be written several times in one edge; each distinct value counts
            # once, or the ratio list is dominated by repetitions.
            pres = sorted({v: t for v, t in pres}.items())
            area = sorted({v: t for v, t in area}.items())
            force = sorted({v: t for v, t in force}.items())
            products = [(a * q, f'{at} x {pt}') for q, pt in pres for a, at in area]
            ratios = [(prod / fo, label, fot) for prod, label in products for fo, fot in force if fo]
            rows.append({'net': net, 'id': e.get('id'), 'status': e.get('status'),
                         'pressures': [t for _, t in pres], 'areas': [t for _, t in area],
                         'forces': [t for _, t in force],
                         'ratios': [{'ratio': r, 'product': l, 'against': f} for r, l, f in ratios]})

    # One edge per id. Both nets carry the harvested edges, and the tissue net is nested inside the
    # main one in places, so an id can appear four times.
    uniq = {}
    for r in rows:
        uniq.setdefault(r['id'], r)
    rows = list(uniq.values())

    impossible, spread, bare = [], [], []
    for r in rows:
        rs = [x['ratio'] for x in r['ratios']]
        if not rs:
            bare.append(r)
            continue
        if min(rs) < 1.0:
            impossible.append(r)
        if max(rs) / min(rs) > PEAK_OVER_MEAN_MAX:
            spread.append(r)

    out = {'edges_with_pressure_and_force': len(rows),
           'edges_with_an_area_so_a_product_exists': len(rows) - len(bare),
           'edges_with_a_product_below_its_stated_load': len(impossible),
           'edges_whose_products_spread_more_than_peak_over_mean': len(spread),
           'peak_over_mean_flag_threshold': PEAK_OVER_MEAN_MAX,
           'reading': ('a ratio below 1 is impossible for a peak pressure; a spread above the '
                       'threshold means the rows cannot all be at the same load'),
           'review_state': 'PENDING_INDEPENDENT_REVIEW',
           'impossible': impossible, 'spread': spread, 'no_area_stated': bare}
    d = W / 'results/ASSEMBLY_LOAD_LABEL_SWEEP'
    d.mkdir(parents=True, exist_ok=True)
    (d / 'SWEEP_V1.json').write_text(json.dumps(out, indent=1, ensure_ascii=False))

    print(f'Edges with both pressure and load: {len(rows)}')
    print(f'  of which with an area, so a product is available: {len(rows) - len(bare)}')
    print(f'  produkt UNDER sin angivna last (impossible for a top): {len(impossible)}')
    print(f'  products that span more than {PEAK_OVER_MEAN_MAX}x: {len(spread)}')
    seen = set()
    for r in impossible + spread:
        if r['id'] in seen:
            continue
        seen.add(r['id'])
        rs = sorted(x['ratio'] for x in r['ratios'])
        print(f"   {r['id']}  {r['status']}  kvoter {rs[0]:.4g} .. {rs[-1]:.4g}")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
