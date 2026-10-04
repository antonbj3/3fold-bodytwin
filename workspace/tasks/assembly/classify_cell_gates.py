"""Classify every cell's frozen acceptance criterion: scalar, band, or distributional.

RESULTS_2026-10-03.md section 5 reports that of 43 cells, 21 accept on a scalar point-or-band test,
only 7 carry any distributional statistic, and of eight tested against a distribution 3 flipped. The
instruction is to run the distribution gate on all 43. The first column of that table is which gate
each cell actually has, and that is readable from each cell's own frozen criterion rather than
assumed -- so it is done here, mechanically, before any resampling.

Classification is by what the criterion TEXT commits to, not by what the cell computes:

  distributional  the criterion names draws, a confidence or credible interval, a percentile, a
                  quantile, a standard deviation or an explicit distribution
  band            the criterion accepts inside an interval or within a tolerance, with no dispersion
  scalar          the criterion compares a single number to a threshold
  none_found      no frozen criterion text located in the cell directory

A cell can carry several; the strongest present is reported, and every hit is listed so a band that
sits inside a distributional criterion is visible rather than hidden. PENDING_INDEPENDENT_REVIEW.
"""
from __future__ import annotations

import json
import pathlib
import re

W = pathlib.Path('')
CELLS = W / 'tasks/free48/sources'

DIST = re.compile(r'\b(\d[\d\s,]*draws?|seeded draws|confidence interval|credible interval|'
                  r'\d+\s*%\s*(CI|interval)|percentile|quantile|standard deviation|\bSD\b|\bSEM\b|'
                  r'bootstrap|posterior|distribution)\b', re.I)
BAND = re.compile(r'\b(within\s+[\d.]|inside the|tolerance|between\s+[\d.]|band|range of)\b', re.I)
SCALAR = re.compile(r'\b(at least\s+[\d.]|exceeds?\s+[\d.]|greater than\s+[\d.]|less than\s+[\d.]|'
                    r'below\s+[\d.]|above\s+[\d.]|threshold)\b', re.I)
CRIT_FILE = re.compile(r'(PREREG|prereg|CRITERI|criteri|ACCEPT|accept)', re.I)


def criterion_text(cell_dir: pathlib.Path) -> tuple[str, str]:
    """The frozen-criterion section if there is one, else the whole prereg, else empty."""
    for f in sorted(cell_dir.rglob('*')):
        if not f.is_file() or f.suffix not in ('.md', '.txt') or not CRIT_FILE.search(f.name):
            continue
        txt = f.read_text(errors='replace')
        # The first version matched only "frozen acceptance", "acceptance criterion" and "frozen
        # criterion", and fell back to the first 2500 characters. That reported 28 of 45 cells as
        # having no criterion language, which was my parser and not the cells: Q019 has "## Frozen
        # acceptance criteria" (plural) and Q049 has "## Frozen comparisons and criteria". Match the
        # heading on the WORDS rather than on a phrase, and fall back to the WHOLE file rather than
        # its opening, because a criterion section usually sits two thirds of the way down.
        heads = [m for m in re.finditer(r'^#+[^\n]*\b(criteri|acceptance|gate|falsif)[^\n]*$',
                                        txt, re.I | re.M)]
        if heads:
            start = heads[0].start()
            nxt = re.search(r'^#+\s', txt[heads[0].end():], re.M)
            end = heads[0].end() + (nxt.start() if nxt else len(txt))
            return txt[start:end], f.name
        return txt, f.name
    return '', ''


def main() -> None:
    rows = []
    for d in sorted(p for p in CELLS.iterdir() if p.is_dir()):
        text, src = criterion_text(d)
        hits = {'distributional': DIST.findall(text), 'band': BAND.findall(text),
                'scalar': SCALAR.findall(text)}
        if not text:
            kind = 'none_found'
        elif hits['distributional']:
            kind = 'distributional'
        elif hits['band']:
            kind = 'band'
        elif hits['scalar']:
            kind = 'scalar'
        else:
            kind = 'no_criterion_language'
        rows.append({'cell': d.name, 'gate_kind': kind, 'criterion_source': src,
                     'distributional_hits': len(hits['distributional']),
                     'band_hits': len(hits['band']), 'scalar_hits': len(hits['scalar']),
                     'needs_distribution_gate': kind in ('band', 'scalar', 'no_criterion_language',
                                                         'none_found')})
    counts = {}
    for r in rows:
        counts[r['gate_kind']] = counts.get(r['gate_kind'], 0) + 1
    out = W / 'results/ASSEMBLY_CELL_GATES'
    out.mkdir(parents=True, exist_ok=True)
    (out / 'CLASSIFICATION.json').write_text(json.dumps(
        {'n_cells': len(rows), 'counts': counts, 'rows': rows,
         'review_state': 'PENDING_INDEPENDENT_REVIEW'}, indent=1, ensure_ascii=False))
    print(f'celler: {len(rows)}   {counts}')
    print(f"need distribution gate: {sum((r['needs_distribution_gate'] for r in rows))}")
    print()
    for r in rows:
        if r['needs_distribution_gate']:
            print(f"  {r['cell']:<26} {r['gate_kind']:<22} band={r['band_hits']} skalar={r['scalar_hits']}  {r['criterion_source'][:28]}")
    print('wrote', out / 'CLASSIFICATION.json')


if __name__ == '__main__':
    main()
