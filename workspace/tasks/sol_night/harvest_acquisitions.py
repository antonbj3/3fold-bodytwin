#!/usr/bin/env python3
"""Harvest acquisition requirements from all lane ports to EN list.

Why. On the night of 2–3 October, five independent lanes converged on the same conclusion: what binds is
measurements we do not have, not resolution and not method. The resolution atlas gave 96,08 % of the undecidables
the magnitudes without nameable other magnitude; expert source rendered 175 of 188 records unprovable; the eyeball landed
on a single missing greatness; the anchor cells were closed on a missing receipt. But the requirements remain in each
lane's own files — mentioned in twelve lanes — while the combined list has nine entries. This crop makes them
to a list.

Rule that makes the output useful: a record WITHOUT entity is not searchable, so entity is the only hard thing
the requirement for an item to be considered orderable. Entries without units are listed separately as incomplete,
aldrig tysta.
"""
from __future__ import annotations

import json
import os
import re
import sys

W = ''
SLOTS = os.path.join(W, 'tasks/build_night/slots.txt')

# Keys that in practice carry an acquisition requirement in the lanes' gates. The list is intentionally broad;
# the precision comes from the unit requirement below, not from the key name.
KEY_HINTS = re.compile(
    r'acquisition|acquire|required_measurement|measurement_required|reopening_receipt|'
    r'needs_measurement|unmeasured|missing_measurement|forvarv|data_law',
    re.I)
UNIT_PAT = re.compile(
    "\\b(m|mm|cm|um|µm|nm|s|ms|min|h|K|Pa|kPa|MPa|mmHg|N|J|W|mol|mmol|nmol|umol|µmol|M|mM|uM|µM|nM|pM|mL|L|g|mg|kg|D|dioptr\\w*|cm\\^-1|per_min|1/min|J/m\\^2|J/m2|mL/min|m/s|nmol/m\\^2/s|percent|%|dimensionless|dimensionless)\\b")


def texts_of(obj, path=''):
    """Flatten to (path, text) for each string and each dict with a claim key."""
    out = []
    if isinstance(obj, dict):
        for k, v in obj.items():
            p = f'{path}/{k}'
            if KEY_HINTS.search(k) and not isinstance(v, (dict, list)):
                out.append((p, str(v)))
            elif KEY_HINTS.search(k) and isinstance(v, (dict, list)):
                out.append((p, json.dumps(v, ensure_ascii=False)))
            out.extend(texts_of(v, p))
    elif isinstance(obj, list):
        for i, v in enumerate(obj):
            out.extend(texts_of(v, f'{path}[{i}]'))
    elif isinstance(obj, str) and KEY_HINTS.search(obj):
        out.append((path, obj))
    return out


def main() -> int:
    lanes = [l.strip() for l in open(SLOTS) if l.strip()]
    found, incomplete, wellformed = [], [], []
    seen = set()
    # Lane-local ACQUISITION_TARGETS files are already in the correct schema and should go in unchanged.
    for lane in lanes:
        base = os.path.join(W, 'results', lane)
        for root, _dirs, fns in os.walk(base) if os.path.isdir(base) else []:
            for fn in fns:
                if fn.startswith('ACQUISITION_TARGETS') and fn.endswith('.json'):
                    try:
                        items = json.load(open(os.path.join(root, fn)))
                    except Exception:
                        continue
                    if isinstance(items, dict):
                        items = items.get('targets') or []
                    for it in items if isinstance(items, list) else []:
                        if isinstance(it, dict) and it.get('quantity'):
                            it = dict(it)
                            it['_lane'] = lane
                            it['_file'] = os.path.relpath(os.path.join(root, fn), W)
                            wellformed.append(it)
    for lane in lanes:
        base = os.path.join(W, 'results', lane)
        if not os.path.isdir(base):
            continue
        files = []
        for pat in ('PORT_R', 'night_rounds/r'):
            d = os.path.join(base, os.path.dirname(pat)) if '/' in pat else base
            if not os.path.isdir(d):
                continue
            for fn in os.listdir(d):
                if fn.startswith(os.path.basename(pat)) and fn.endswith('.json'):
                    files.append(os.path.join(d, fn))
        # newest first, and take no more than six files per lane: requirements are repeated between rounds
        files.sort(key=os.path.getmtime, reverse=True)
        for f in files[:6]:
            try:
                d = json.load(open(f))
            except Exception:
                continue
            for p, t in texts_of(d):
                t = ' '.join(t.split())
                if len(t) < 12 or len(t) > 600:
                    continue
                norm = re.sub(r'[^a-z0-9]+', '', t.lower())[:110]
                if norm in seen:
                    continue
                seen.add(norm)
                unit = UNIT_PAT.search(t)
                rec = {'lane': lane, 'where': os.path.relpath(f, W) + p,
                       'requirement': t, 'unit_found': unit.group(0) if unit else None}
                (found if unit else incomplete).append(rec)

    out = {
        '_meta': {
            'what': "Acquisition requirements harvested from the lane ports to a list.",
            'rule': "A post without a device cannot be searched and ends up incomplete, never silent.",
            'lanes_scanned': len(lanes),
            'wellformed_orderable': len(wellformed),
            'with_unit': len(found),
            'without_unit': len(incomplete),
            'status': 'PENDING_INDEPENDENT_REVIEW',
        },
        'wellformed_orderable': wellformed,
        'with_unit': found,
        'without_unit': incomplete,
    }
    dst = os.path.join(W, 'notes/ACQUISITION_HARVEST.json')
    json.dump(out, open(dst, 'w'), ensure_ascii=False, indent=1)
    print(f'lanes: {len(lanes)}  well-shaped: {len(wellformed)}  med enhet: {len(found)}  without unit: {len(incomplete)}')
    per = {}
    for r in found:
        per[r['lane']] = per.get(r['lane'], 0) + 1
    for lane, n in sorted(per.items(), key=lambda x: -x[1]):
        print(f'  {n:3d}  {lane}')
    print(f'skrivet: {os.path.relpath(dst, W)}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
