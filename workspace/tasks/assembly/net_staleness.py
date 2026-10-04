"""Re-read every edge's evidence and flag the ones whose text no longer matches the file.

READ THIS FIRST: sign and percent-versus-fraction are the same number here. The first version of
this script raised three STALE alarms and two were its own fault -- a signed -489.6488 against prose
that says "lowers it by 489.65", and a fraction 0.26508 against prose that says "26.51 percent". An
alarm that cries wolf twice means the third one gets ignored, so any new comparison added below must
be checked against that before it is allowed to flag.

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
import pathlib
import re
from pathlib import Path

import importlib.util as _ilu

_spec = _ilu.spec_from_file_location(
    'relocate_stale_evidence',
    str(Path(__file__).resolve().parent / 'relocate_stale_evidence.py'))
_rel = _ilu.module_from_spec(_spec)
_spec.loader.exec_module(_rel)
# carries() is the relocator's test: literal match, or agreement to the last digit the COARSER
# side states. Keeping a second, literal-only copy of that rule here flagged HARVEST-E0065 as
# stale for citing 4.072 where its source says 4.07 -- a correct edge reported as poison. One
# rule, one file.
carries = _rel.carries

W = Path('.')
NET = W / 'CONSTRAINT_NETS.json'
# The value class used to be [\d.eE+], which accepts 1.25e+09 and rejects 1.25e-09. Four edges
# whose provenance resolved perfectly were therefore reported as "evidence carries no
# file :: key = value": the exponent sign, not the evidence, was the defect. Exponents of either
# sign are now accepted, and the grammar is spelled out rather than being a character bag.
EV = re.compile(r'^(?P<file>[^\s:]+)\s*::\s*(?P<key>[^=]+?)\s*=\s*'
                r'(?P<value>[+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?)\s*$')


def dig(obj, dotted: str):
    # Resolve a dotted path, and if it does not resolve, SEARCH for the leaf key anywhere.
    #
    # The strict reading declared 19 of 22 edges unverifiable with "key is not a number". That was
    # not the net being wrong: the keys exist, the hand-written paths just do not spell out the
    # nesting, because a value sits under `summary` or a level deeper than the path says. A
    # provenance link that only resolves when a human typed the nesting correctly is decorative,
    # so the leaf name is searched for when the path misses, and the row records which found it.
    # List indices count as nesting too. The path walker only stepped through dicts, so an edge
    # whose number sits in a list -- one device row of several -- could not be addressed at all,
    # and the leaf search then found one value per row and called the leaf ambiguous. A path
    # segment of the form name[i] now steps into the list as well.
    cur, ok = obj, True
    for part in dotted.split('.'):
        part = part.strip()
        idx = re.findall(r'\[(\d+)\]', part)
        name = part.split('[')[0]
        if isinstance(cur, dict) and name in cur:
            cur = cur[name]
        else:
            ok = False
            break
        for i in idx:
            if isinstance(cur, list) and int(i) < len(cur):
                cur = cur[int(i)]
            else:
                ok = False
                break
        if not ok:
            break
    if ok:
        return cur, 'path'
    leaf = dotted.split('.')[-1].strip().split('[')[0]
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
    # A rendering only counts if it still IS the value. At zero decimal places every quantity
    # below 0.5 renders as "0", and "0" occurs in almost any prose, so the loop used to clear such
    # an edge on a digit that carries none of its information -- a silent pass, the mirror of the
    # false alarm this module was written to stop. Measured on the 56-edge net: all 53 verifiable
    # edges match a rendering within 1 percent, so no current verdict rests on the loose reading
    # and this bound changes no count today; it stops a future edge from being cleared by a zero.
    for v in (value, -value, value * 100.0, value / 100.0):
        for places in range(0, 7):
            s = f'{v:.{places}f}'
            if s.endswith('.'):
                s = s[:-1]
            if s in text and abs(float(s) - v) <= 0.01 * max(abs(v), 1e-300):
                return True
    # Prose writes a small number in scientific notation and it writes it short: an evidence value
    # of 1.255019203345805e-09 appears in the text as "1.255e-9". Fixed-point rendering of that
    # value is "0.000000" at every width the loop above tries, so a correct edge was one rounding
    # convention away from being called STALE -- the cry-wolf failure this module exists to avoid.
    # Only renderings OF THIS VALUE are accepted, so a different mantissa or exponent still flags.
    for digits in range(1, 8):
        sci = f'{value:.{digits}e}'
        mant, _, exp = sci.partition('e')
        mant = mant.rstrip('0').rstrip('.') if '.' in mant else mant
        sign = '-' if exp[0] == '-' else ''
        n = exp[1:].lstrip('0') or '0'
        for e_txt in (f'e{sign}{n}', f'e{sign}{int(n):02d}'):
            if f'{mant}{e_txt}' in text:
                return True
    return repr(value) in text


def main() -> int:
    net = json.loads(NET.read_text())['bodytwin']['tissue_constraint_net']
    rows = []
    for e in net['edges']:
        # An edge whose evidence was SEARCHED FOR and not found is a different state from one nobody
        # examined, and lumping them under UNCHECKED hides the difference. T-E6's count exists only as
        # the string "98/98", T-E9's numbers appear nowhere in its lane, and T-E18's index file is not
        # under this workspace -- all three carry evidence_unresolved with a note saying what was
        # searched. UNCHECKED should mean unexamined, and that number is the one that ought to be zero.
        if e.get('evidence_unresolved'):
            rows.append((e['id'], 'UNRESOLVED_DECLARED',
                         str(e.get('evidence_unresolved_note', 'searched, not found'))[:120]))
            continue
        ev = str(e.get('evidence', ''))
        # The harvested detail-layer edges cite a LINE in a markdown document rather than a key in
        # a JSON file: `path :: L145 = <the text>`. My checker read only JSON keys, so admitting 99
        # of them turned 99 verifiable edges into 99 unchecked ones. Read the line and compare the
        # text, which is the same audit in a different file format.
        mline = re.match(r'^(?P<file>[^\s:]+)\s*::\s*L(?P<line>\d+)\s*=\s*(?P<text>.+)$',
                         ev.strip(), re.S)
        if mline:
            f = pathlib.Path(mline.group('file'))
            if not f.is_absolute():
                f = W / mline.group('file')
            if not f.exists():
                rows.append((e['id'], 'UNCHECKED', f"evidence file missing: {mline.group('file')}"))
                continue
            try:
                lines = f.read_text(errors='replace').splitlines()
            except Exception as exc:
                rows.append((e['id'], 'UNCHECKED', f'{type(exc).__name__} reading evidence'))
                continue
            idx = int(mline.group('line')) - 1
            # Normalise the dash family before comparing. The cited text writes -90 with an ASCII
            # hyphen while the source table uses the Unicode minus U+2212, and 25 of 30 stale
            # flags came from that one difference in one file. Fourth sign-related bug in this
            # checker tonight: exponent sign, hyphen-as-minus, heading-level citation, and now
            # the dash encoding. Signs are where this kind of code keeps failing.
            DASH = {'\u2212': '-', '\u2013': '-', '\u2014': '-', '\u2010': '-', '\u00ad': '-'}

            def norm(s):
                for a, b in DASH.items():
                    s = s.replace(a, b)
                return s

            claimed = norm(' '.join(mline.group('text').split()))
            # The harvested citations point at the SECTION the claim lives in, not at the exact
            # line: L61 is '## 2. Method', L148 is '## 4. Falsifier A'. A plus-or-minus-three-line
            # window flagged 41 correctly cited edges as stale for that reason alone. When the
            # cited line is a heading, read to the next heading, which is what the citation means.
            if lines[idx].lstrip().startswith('#'):
                end = idx + 1
                while end < len(lines) and not lines[end].lstrip().startswith('#'):
                    end += 1
                window = norm(' '.join(' '.join(lines[idx:end]).split()))
            else:
                window = norm(' '.join(' '.join(lines[max(idx - 3, 0):idx + 4]).split()))
            # the cited text, or every number in it, must appear within a few lines of the cite
            # A hyphen between digits is a RANGE, not a minus. Reading "1.176-1.357" as 1.176 and
            # -1.357 made the second number unfindable and flagged 43 correctly cited edges as
            # stale, including one I had verified by hand. Same family as the exponent-sign bug an
            # audit found in this file earlier: my number extraction keeps mis-reading signs.
            nums = re.findall(r'(?<![\d.])-?\d+[.,]?\d*', claimed)
            if carries(window, claimed)[0]:
                rows.append((e['id'], 'OK', f"line {mline.group('line')} matches"))
            else:
                rows.append((e['id'], 'STALE',
                             f"line {mline.group('line')} does not carry the cited text or numbers"))
            continue
        # A third evidence format, carried by 24 of the 99 admitted detail-layer edges and read by
        # nothing until now: `file :: /a/b/c = <free text with numbers>`. The locator is a JSON
        # pointer and the value is prose, so both the dotted-key walker and the numeric value class
        # reject it, and all 24 counted as unauditable. They are not: resolve the pointer, render
        # the subtree, and run the same text-and-number comparison the markdown branch runs.
        mptr = re.match(r'^(?P<file>[^\s:]+)\s*::\s*(?P<ptr>/[^=]+?)\s*=\s*(?P<text>.+)$',
                        ev.strip(), re.S)
        if mptr and not EV.match(ev.strip()):
            f = pathlib.Path(mptr.group('file'))
            if not f.is_absolute():
                f = W / mptr.group('file')
            if not f.exists():
                rows.append((e['id'], 'UNCHECKED', f"evidence file missing: {mptr.group('file')}"))
                continue
            try:
                cur = json.loads(f.read_text())
            except Exception as exc:
                rows.append((e['id'], 'UNCHECKED', f'{type(exc).__name__} reading evidence'))
                continue
            how, missing = 'pointer', None
            for seg in [x for x in mptr.group('ptr').split('/') if x]:
                if isinstance(cur, dict) and seg in cur:
                    cur = cur[seg]
                elif isinstance(cur, list) and seg.isdigit() and int(seg) < len(cur):
                    cur = cur[int(seg)]
                else:
                    missing = seg
                    break
            if missing is not None:
                # The pointer missing a segment is not the same as the claim being absent: fall
                # back to the whole document, and say in the row that the address did not resolve.
                cur = json.loads(f.read_text())
                how = f'pointer broke at /{missing}, searched whole file'
            rendered = _rel.norm(json.dumps(cur, ensure_ascii=False))
            claimed = _rel.norm(mptr.group('text'))
            if carries(rendered, claimed)[0]:
                rows.append((e['id'], 'OK', f'{how}: numbers present'))
            else:
                rows.append((e['id'], 'STALE', f'{how}: cited numbers not in that subtree'))
            continue
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
        # An evidence pointer quotes the value at the precision a reader can carry, so the comparison
        # has to be made at the precision STATED and not to nine decimals. Measured 2026-10-05: three
        # pointers quoting 0.8048, 0.4425 and 1.6641 were reported STALE against files holding
        # 0.8047619047619047, 0.4425149700598802 and 1.6641005886756877 -- all three correct
        # roundings. The checker was comparing a rounded quotation against a full-precision value and
        # printing a truncation of the latter, which read as a different number.
        raw = m.group('value')
        digits = len(raw.split('.')[1]) if '.' in raw and 'e' not in raw.lower() else None
        if digits is not None:
            agrees = f'{current:.{digits}f}' == f'{stated:.{digits}f}'
        else:
            agrees = abs(current - stated) <= 1e-9 * max(1.0, abs(stated))
        if not agrees:
            rows.append((e['id'], 'STALE',
                         f'evidence says {stated} but the file now holds {current}'))
        elif not number_appears(e['constraint'], current):
            rows.append((e['id'], 'STALE',
                         f'the value {current} does not appear in the constraint text'))
        else:
            rows.append((e['id'], 'OK', f'{m.group("key").strip()} = {current} ({how})'))

    # Write a stamp back onto each verified edge so a LATER run can see that the evidence moved,
    # without anyone having to declare a supersession. That is the half the graph lane's
    # consumers_of_superseded() cannot cover: it needs the supersession declared, and the failure
    # that started this was precisely that nobody declared it. A stamp turns movement into something
    # the file itself reports.
    full = json.loads(NET.read_text())
    net_w = full['bodytwin']['tissue_constraint_net']
    by_id = {e['id']: e for e in net_w['edges']}
    moved = []
    for eid, status, detail in rows:
        e = by_id.get(eid)
        if e is None:
            continue
        m = EV.match(str(e.get('evidence', '')).strip())
        if not m:
            continue
        f = W / m.group('file')
        if not f.exists():
            continue
        stamp = round(f.stat().st_mtime, 3)
        prev = e.get('evidence_mtime_at_last_check')
        if prev is not None and stamp != prev and status == 'OK':
            moved.append((eid, prev, stamp))
        e['evidence_mtime_at_last_check'] = stamp
        e['last_checked_status'] = status
    NET.write_text(json.dumps(full, indent=1, ensure_ascii=False))

    counts = {}
    for _, s, _ in rows:
        counts[s] = counts.get(s, 0) + 1
    if moved:
        counts['MOVED_SINCE_LAST_CHECK'] = len(moved)
    out = {'checked': len(rows), 'counts': counts,
           'evidence_moved_since_last_check': [
               {'edge': a, 'was': b, 'now': c} for a, b, c in moved],
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
    for a, b, c in moved:
        print(f'  [MOVED    ] {a[:56]:56s} beviset andrades efter forra kontrollen')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
