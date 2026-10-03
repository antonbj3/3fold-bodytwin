#!/usr/bin/env python3
"""Index the completed work cells in the older BodyTwin project so current work can consume them.

Why this exists. Anton asked repeatedly tonight whether we are supposed to have masses of Q-families
and whether work was lost on the road from cad-to-simulation to bodytwin to this workspace. Measured: the
older project at ~/projects/bodytwin holds 1309 distinct work cells with a decisive number, of which
696 are MEASURED and 431 REFUTED. The current workspace has 43 cells and cites none of them. Nothing
was lost; it was disconnected. An index is what reconnects it, and it is cheaper than rebuilding.

The 431 REFUTED cells matter as much as the measured ones: a preserved negative stops a later lane
from re-running a question that already failed, and five lanes closed tonight on questions that were
ill-posed rather than hard.

Read-only on the source project. This writes one file into the workspace and touches nothing else.

Usage:
  python3 tasks/build_night/index_old_cells.py            # writes notes/OLD_CELL_INDEX.json
  python3 tasks/build_night/index_old_cells.py --quantity flow   # grep the index by quantity text
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

OLD = Path('source_repository/data')
OUT = Path(__file__).resolve().parents[2] / 'notes' / 'OLD_CELL_INDEX.json'
MAX_BYTES = 2_000_000          # a few domain files are multi-megabyte dumps, not cell records
FIELDS = ('cell', 'decisive_number', 'decisive_quantity', 'status', 'verdict',
          'registered_prediction_was', 'written_utc', 'gate', 'consumes')


def harvest() -> list[dict]:
    rows = []
    for domain in sorted(p for p in OLD.iterdir() if p.is_dir()):
        for f in sorted(domain.glob('*.json')):
            try:
                if f.stat().st_size > MAX_BYTES:
                    continue
                d = json.loads(f.read_text())
            except (OSError, ValueError, UnicodeDecodeError):
                continue
            if not isinstance(d, dict):
                continue
            # A cell record is one that committed to a number or a verdict. Anything else is data.
            if not ('decisive_number' in d or 'decisive_quantity' in d
                    or ('cell' in d and 'verdict' in d)):
                continue
            row = {k: d[k] for k in FIELDS if k in d}
            row['domain'] = domain.name
            row['path'] = str(f.relative_to(OLD.parent))
            # Keep the verdict readable but bounded; the path is there for the full text.
            for k in ('verdict', 'decisive_quantity', 'gate', 'consumes'):
                if isinstance(row.get(k), str) and len(row[k]) > 400:
                    row[k] = row[k][:400] + ' …'
            rows.append(row)
    return rows


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument('--quantity', help='print index rows whose decisive_quantity matches this text')
    args = ap.parse_args()

    if args.quantity and OUT.exists():
        rows = json.loads(OUT.read_text())['cells']
        q = args.quantity.lower()
        hits = [r for r in rows if q in str(r.get('decisive_quantity', '')).lower()]
        print(f'{len(hits)} cells whose crucial quantity matches {args.quantity!r}:')
        for r in hits[:40]:
            print(f"  [{r.get('status')}] {r.get('cell')} ({r['domain']}): "
                  f"{r.get('decisive_number')} — {str(r.get('decisive_quantity'))[:90]}")
        return 0

    if not OLD.is_dir():
        print(f'source project missing: {OLD}', file=sys.stderr)
        return 1

    rows = harvest()
    from collections import Counter
    status = Counter(str(r.get('status')) for r in rows)
    domains = Counter(r['domain'] for r in rows)
    numeric = sum(1 for r in rows if isinstance(r.get('decisive_number'), (int, float)))
    out = {
        '_meta': {
            'what': 'Completed work cells in the older BodyTwin project, indexed so current work can '
                    'consume them instead of re-deriving. Read-only harvest; source unchanged.',
            'source': str(OLD),
            'cells': len(rows),
            'distinct_cell_ids': len({str(r.get('cell')) for r in rows if r.get('cell')}),
            'with_numeric_decisive_number': numeric,
            'status_counts': dict(status.most_common()),
            'domains': len(domains),
            'caveat': 'Statuses are the source cells own and are PENDING_INDEPENDENT_REVIEW here. '
                      'A REFUTED cell is a preserved negative and is as useful as a MEASURED one: it '
                      'stops a later lane re-running a question that already failed.',
        },
        'top_domains': dict(domains.most_common(25)),
        'cells': rows,
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, ensure_ascii=False, indent=1))
    print(f'skrev {OUT}')
    print(f"  celler {len(rows)}, distinkta ID {out['_meta']['distinct_cell_ids']}, "
          f"med numeriskt tal {numeric}")
    print(f"  status {dict(status.most_common(6))}")
    print(f'  domains {len(domains)}; largest: {list(domains.most_common(5))}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
