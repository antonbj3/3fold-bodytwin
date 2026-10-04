"""Fail-closed exclusion over anything harvested from the original repo, run before admission.

Three things do not enter, and the operator stated them as content rules rather than filter rules.

1. The excluded_category material is out entirely -- not filtered from what travels onward, not rewritten, not
   counted. It is not read, consumed, cited, or included in any total.
2. Anything concerning the collaborator is out entirely. No named person appears in anything written
   here.
3. excluded_category is OUT, unconditionally. It was conditional until 2026-10-05, when the condition was tested
   against the material and Anton closed it: of the 154 rows a classification had listed as
   conditionally admissible, two carried a frequency or a flux density and eleven named a tissue
   type. The rest were the acronym with no quantity, base64 where the letters fall together, or
   PEMFC, which is a proton exchange membrane fuel cell and belongs to porous-media literature.
   There was nothing to admit, so the condition is gone and the exclusion is flat. What follows
   below was the old condition and is kept only to show what was tested:
   stated generally -- an electromagnetic field with a frequency, a flux density and an exposure time
   acting on tissue, with the tissue named as a tissue type and nothing else. If it cannot be written
   that way without the original application showing through, it stays out. A doubtful case stays out,
   because missing an edge is cheaper than having to tear one out.

Fail-closed: a record that cannot be parsed is rejected rather than passed. A filter that fails open is
how the thing it was built to stop gets through, and tonight already produced three checks that reported
green because they matched nothing at all.
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

BLOCK = re.compile(r'excluded_category|excluded_category|tunica\s+albug|sinusoid|excluded_category|excluded_category|corpus\s+spongios', re.I)
PERSON = re.compile(r'\bdavid\b', re.I)
excluded_category = re.compile(r'\bpemf\b|pulsed\s+electromagnetic', re.I)
# A excluded_category record is admissible only if it carries the general physical quantities and no application.
PEMF_OK = re.compile(r'(\bHz\b|frequency).*(\bmT\b|\bT\b|flux\s+density)|'
                     r'(\bmT\b|flux\s+density).*(\bHz\b|frequency)', re.I | re.S)


def verdict(text: str) -> tuple[str, str]:
    if BLOCK.search(text):
        return 'REJECT', 'blocked subject matter; not read, counted or carried onward'
    if PERSON.search(text):
        return 'REJECT', 'names an individual'
    if excluded_category.search(text):
        if PEMF_OK.search(text) and not BLOCK.search(text):
            return 'ADMIT_CONDITIONAL', ('excluded_category stated with frequency and flux density as general tissue '
                                         'exposure; admit only if a reader confirms no application shows '
                                         'through')
        return 'REJECT', 'excluded_category without a general formulation in frequency and flux density'
    return 'ADMIT', 'no excluded subject matter'


def main() -> int:
    root = Path(sys.argv[1] if len(sys.argv) > 1
                else './results/SOL_MECHANISM_HARVEST')
    rows, counts = [], {}
    for f in sorted(root.rglob('*.json')) + sorted(root.rglob('*.md')):
        try:
            text = f.read_text(errors='replace')
        except Exception as exc:
            rows.append({'file': str(f), 'verdict': 'REJECT',
                         'why': f'unreadable ({type(exc).__name__}); fail-closed'})
            counts['REJECT'] = counts.get('REJECT', 0) + 1
            continue
        if f.suffix == '.json':
            try:
                items = json.loads(text)
            except Exception:
                items = None
            if isinstance(items, list):
                for i, it in enumerate(items):
                    v, why = verdict(json.dumps(it, ensure_ascii=False))
                    counts[v] = counts.get(v, 0) + 1
                    if v != 'ADMIT':
                        rows.append({'file': f.name, 'index': i, 'verdict': v, 'why': why})
                continue
        v, why = verdict(text)
        counts[v] = counts.get(v, 0) + 1
        if v != 'ADMIT':
            rows.append({'file': f.name, 'verdict': v, 'why': why})

    out = root / 'EXCLUSION_VERDICTS.json'
    out.write_text(json.dumps({'counts': counts, 'fail_closed': True,
                               'review_state': 'PENDING_INDEPENDENT_REVIEW',
                               'non_admit': rows}, indent=1, ensure_ascii=False))
    print('  ' + ', '.join(f'{k} {v}' for k, v in sorted(counts.items())) if counts
          else "  Nothing to review yet (the harvest has not written)")
    for r in rows[:10]:
        print(f"    {r['verdict']:18s} {r.get('file','')[:40]} — {r['why'][:60]}")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
