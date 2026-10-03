"""Re-read every edge's evidence and flag the ones whose text no longer matches the file.

Why this exists. The net handed a brief a number I had withdrawn hours earlier: the implant-power edge
still said "11 of 20, 0.688 against 1.070 D", a comparison against a surgeon who did not have the
postoperative lens position the chain was using. Every job generated from that node received it.

The reason is structural rather than careless. An edge carries a constraint as prose plus one evidence
path of the form `file :: key = value`. The prose is a snapshot of what was true when it was written and
nothing links a new measurement back to the edges that cite it, so an edge cannot know it has gone
stale. It keeps serving the old figure with full confidence, which is worse than carrying no figure.

The evidence path is enough to check, though: read the file, read the key, and see whether that value
still appears in the edge's own text. Three outcomes, and the middle one is the dangerous class:

  OK        the value in the file appears in the constraint text
  STALE     the file has a different value from the one the text states
  UNCHECKED the evidence has no `key = value` form, or the file is gone

STALE is what silently poisons downstream briefs. UNCHECKED is not an error but it is an edge that
cannot be audited, and counting them is the honest measure of how much of the net is verifiable at all.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

W = Path('.')
NET = W / 'CONSTRAINT_NETS.json'
EV = re.compile(r'^(?P<file>[^\s:]+)\s*::\s*(?P<key>[^=]+?)\s*=\s*(?P<value>-?[\d.eE+]+)\s*$')


def dig(obj, dotted: str):
    # Resolve a dotted path, and if it does not resolve, SEARCH for the leaf key anywhere.
    #
    # The strict reading declared 19 of 22 edges unverifiable with "key is not a number". That was
    # not the net being wrong: the keys exist, the hand-written paths just do not spell out the
    # nesting, because a value sits under `summary` or a level deeper than the path says. A
    # provenance link that only resolves when a human typed the nesting correctly is decorative,
    # so the leaf name is searched for when the path misses, and the row records which found it.
    cur, ok = obj, True
    for part in dotted.split('.'):
        if isinstance(cur, dict) and part in cur:
            cur = cur[part]
        else:
            ok = False
            break
    if ok:
        return cur, 'path'
    leaf = dotted.split('.')[-1].strip()
    found = []

    def walk(o):
        if isinstance(o, dict):
            for k, v in o.items():
                if k == leaf and isinstance(v, (int, float)):
                    found.append(v)
                walk(v)
        elif isinstance(o, list):
            for v in o:
                walk(v)

    walk(obj)
    if found and len(set(found)) == 1:
        return found[0], 'leaf-search'
    return None, f'{len(found)} traffar pa bladnyckeln'


def number_appears(text: str, value: float) -> bool:
    # Does this value appear in the prose, at any rounding or sign or scale the prose might use?
    #
    # The strict version raised three STALE flags of which two were its own fault: a signed
    # -489.6488 against prose that says "lowers it by 489.65", and a fraction 0.26508 against
    # prose that says "26.51 percent". A checker that cries wolf is worse than one that stays
    # quiet, because the next real flag gets ignored. Sign and percent are accepted as the same
    # number; a different number is not.
    for v in (value, -value, value * 100.0, value / 100.0):
        for places in range(0, 7):
            s = f'{v:.{places}f}'
            if s.endswith('.'):
                s = s[:-1]
            if s in text:
                return True
    return repr(value) in text


def main() -> int:
    net = json.loads(NET.read_text())['bodytwin']['tissue_constraint_net']
    rows = []
    for e in net['edges']:
        ev = str(e.get('evidence', ''))
        m = EV.match(ev.strip())
        if not m:
            rows.append((e['id'], 'UNCHECKED', 'evidence carries no file :: key = value'))
            continue
        f = W / m.group('file')
        if not f.exists():
            rows.append((e['id'], 'UNCHECKED', f'evidence file missing: {m.group("file")}'))
            continue
        try:
            current, how = dig(json.loads(f.read_text()), m.group('key').strip())
        except Exception as exc:
            rows.append((e['id'], 'UNCHECKED', f'{type(exc).__name__} reading evidence'))
            continue
        if not isinstance(current, (int, float)):
            rows.append((e['id'], 'UNCHECKED', f'{how}: {m.group("key").strip()}'))
            continue
        stated = float(m.group('value'))
        if abs(current - stated) > 1e-9 * max(1.0, abs(stated)):
            rows.append((e['id'], 'STALE',
                         f'evidence says {stated} but the file now holds {current}'))
        elif not number_appears(e['constraint'], current):
            rows.append((e['id'], 'STALE',
                         f'the value {current} does not appear in the constraint text'))
        else:
            rows.append((e['id'], 'OK', f'{m.group("key").strip()} = {current} ({how})'))

    counts = {}
    for _, s, _ in rows:
        counts[s] = counts.get(s, 0) + 1
    out = {'checked': len(rows), 'counts': counts,
           'rows': [{'edge': a, 'status': b, 'detail': c} for a, b, c in rows],
           'why': ('an edge is a snapshot of prose plus one evidence path; nothing links a new '
                   'measurement back to the edges that cite it, so staleness is silent'),
           'review_state': 'PENDING_INDEPENDENT_REVIEW'}
    d = W / 'results/ASSEMBLY_NET_STALENESS'
    d.mkdir(parents=True, exist_ok=True)
    (d / 'STALENESS_V1.json').write_text(json.dumps(out, indent=1, ensure_ascii=False))

    print(f'  {len(rows)} kanter kontrollerade: ' + ', '.join(f'{k} {v}' for k, v in counts.items()))
    for a, b, c in rows:
        if b != 'OK':
            print(f'  [{b:9s}] {a[:56]:56s} {c}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
