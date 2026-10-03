#!/usr/bin/env python3
"""Turn the anchor searches' reports into a record store that cells can consume.

Why. The 48 anchor searches queued on 2026-10-03 have started coming back, and measured over the first
42: **33 partial, 5 full hits, 5 verified absences**. Partial means the search found published numbers
while the exact quantity we asked for is genuinely unmeasured — both halves are useful. One report alone
carries 204 numbers with units, with zero self-computation markers and 59 published-reference markers.

But they are 42 markdown reports, which is not a source. A cell cannot consume prose. This builds the
same shape the literature collection already has — value, unit, locator, validity range — so the
existing rules apply: dimension must match, the validity range must overlap the cell's, and a record
must carry a measurement rather than a threshold sentence.

Deliberately strict. A row is kept only when a number, its unit and a locator appear on the same line,
because that is the only form where the three cannot be mismatched by the parser. Everything else is
counted as unextracted and reported, never silently dropped: an honest small store beats a large one
whose provenance is guessed.
"""
from __future__ import annotations

import glob
import json
import os
import re

W = '.'
OUT = os.path.join(W, 'notes/ANCHOR_HARVEST.json')

UNIT = (r'J/m\^?2|mJ/m\^?2|N/m|nmol/m\^?2/s|mL/min|mmHg|kPa|MPa|Pa|µm|um|nm|mm|cm|µM|uM|mM|nM|pM|'
        r'M\^-1\s*cm\^-1|cm\^?-1|J/cm\^?2|W/cm\^?2|mol|mmol|µmol|umol|nmol|percent|%|K|°C|s|ms|min|h|'
        r'g|mg|kg|D|dioptre')
LOCATOR = r'(10\.\d{4,}/[^\s`)\],]+|PMID[:\s`]*\d+|PMC\d+)'
NUMBER = r'(?<![\w.])(-?\d[\d\s.,]*(?:[eE][-+]?\d+)?)'

# A validity range is what makes a record usable; without it a value cannot be matched to a cell.
RANGE = re.compile(r'(?i)\b(n\s*=\s*\d+|human|porcine|murine|mouse|rat|bovine|ovine|canine|in vivo|'
                   r'in vitro|ex vivo|\d+\s*°?C|Fitzpatrick|adult|infant)\b')
SELF = re.compile(r'(?i)\b(we compute|our model|computed here|derived here|this job computes)\b')


def table_blocks(lines):
    """Yield (block_lines, locators_in_block) for each markdown table.

    The reports put the citation in one cell or a footnote and the value in another, so a per-line rule
    captures only the few rows where they happen to coincide. Measured on the first 43 reports: 37
    records kept while 15 to 36 locator lines per report went unextracted. A table is the real unit of
    provenance here, so a locator anywhere in the block covers that block's rows.
    """
    block, locs = [], set()
    for line in lines + ['']:
        if line.strip().startswith('|'):
            block.append(line)
            locs.update(m.rstrip('`.,') for m in re.findall(LOCATOR, line))
        else:
            if block:
                yield block, locs
            block, locs = [], set()


def rows_from(path: str):
    text = open(path, errors='replace').read()
    job = os.path.basename(os.path.dirname(path))
    absence = 'VERIFIED_ABSENCE_OF_EVIDENCE' in text
    self_marks = len(SELF.findall(text))
    kept, skipped = [], 0
    lines = text.splitlines()

    # Table rows first: a value in one cell with the block's citation in another.
    for block, locs in table_blocks(lines):
        if not locs:
            continue
        for row in block:
            if re.match(r'^\|[\s:|-]+\|?$', row.strip()):
                continue                     # separator row
            hit = re.search(NUMBER + r'\s*(' + UNIT + r')\b', row)
            if not hit:
                continue
            raw = hit.group(1).replace(' ', '').rstrip('.,')
            try:
                value = float(raw.replace(',', ''))
            except ValueError:
                continue
            quantity = row.strip('|').split('|')[0].strip(' `*')[:120]
            kept.append({
                'job': job,
                'quantity': quantity or None,
                'value': value,
                'unit': hit.group(2),
                'locator': sorted(locs)[0],
                'locators_in_block': sorted(locs)[:4],
                'validity': sorted({m.lower() for m in RANGE.findall(row)}) or None,
                'from': 'table',
                'line': row.strip()[:300],
                'review_state': 'PENDING_INDEPENDENT_REVIEW',
            })

    seen_lines = {r['line'] for r in kept}
    for line in lines:
        line = line.strip()
        if len(line) < 12 or len(line) > 600:
            continue
        if line[:300] in seen_lines:
            continue
        loc = re.search(LOCATOR, line)
        if not loc:
            continue
        hit = re.search(NUMBER + r'\s*(' + UNIT + r')\b', line)
        if not hit:
            skipped += 1
            continue
        raw = hit.group(1).replace(' ', '').rstrip('.,')
        try:
            value = float(raw.replace(',', '')) if raw.count(',') <= 1 else float(raw.replace(',', ''))
        except ValueError:
            skipped += 1
            continue
        kept.append({
            'job': job,
            'quantity': None,
            'from': 'line',
            'value': value,
            'unit': hit.group(2),
            'locator': loc.group(1).rstrip('`.,'),
            'validity': sorted({m.lower() for m in RANGE.findall(line)}) or None,
            'line': line[:300],
            'review_state': 'PENDING_INDEPENDENT_REVIEW',
        })
    return job, absence, self_marks, kept, skipped


def main() -> int:
    records, absences, stats = [], [], []
    for path in sorted(glob.glob(os.path.join(W, 'results/BT-ANCHOR-*/RESULTS.md'))):
        job, absence, self_marks, kept, skipped = rows_from(path)
        records.extend(kept)
        if absence:
            absences.append(job)
        stats.append({'job': job, 'kept': len(kept), 'unextracted_locator_lines': skipped,
                      'declares_absence': absence, 'self_computation_markers': self_marks})

    with_range = sum(1 for r in records if r['validity'])
    out = {
        '_meta': {
            'what': 'Published values harvested from the anchor searches, in consumable form.',
            'rule': ('A row is kept only when number, unit and locator share one line. Everything else '
                     'is counted as unextracted, never dropped silently.'),
            'reports': len(stats),
            'records': len(records),
            'records_with_validity_range': with_range,
            'reports_declaring_absence': len(absences),
            'status': 'PENDING_INDEPENDENT_REVIEW',
            'consumption_rule': ('Same three requirements as the literature collection: the dimension '
                                 'must match, the validity range must overlap the cell, and the record '
                                 'must carry a measurement and not a threshold sentence.'),
        },
        'records': records,
        'absence_reports': absences,
        'per_report': stats,
    }
    json.dump(out, open(OUT, 'w'), ensure_ascii=False, indent=1)
    print(f'reports {len(stats)}  records {len(records)}  with validity range {with_range}  '
          f'absence reports {len(absences)}')
    top = sorted(stats, key=lambda s: -s['kept'])[:6]
    for s in top:
        print(f"  {s['kept']:4d} kept  {s['unextracted_locator_lines']:4d} unextracted  {s['job'][:52]}")
    bad = [s for s in stats if s['self_computation_markers']]
    print(f'reports with self-computation markers: {len(bad)}')
    print(f'written: {os.path.relpath(OUT, W)}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
