"""Re-resolve a STALE line citation by searching the source file for the cited text.

Measured 2026-10-04: all five distinct STALE flags in the admitted detail-layer edges are
LOCATOR drift, not bad evidence. HARVEST-E0001 cites L61 of MECHANISM_INTERVERTEBRAL_DISC.md;
'A_disc | 1800 mm2' sits on line 143. HARVEST-E0065 cites L147 of the Fahraeus document; the
D=6.84 -> 1.247 pair sits on 119 and 149. The numbers are in the file. The line number is from
a different pass over a document that has since grown, so the checker was right to flag and
wrong about the cause -- and a line number is the wrong address for a claim in a living file.

So: for each flagged edge, search EVERY line of the cited file with the same window rule the
checker uses, and rewrite Lnnn to the line that actually carries the claim. An edge whose text
is nowhere in the file is a different state and gets evidence_unresolved with what was searched,
never a silent repair. Nothing is admitted or upgraded here; only the address changes.
"""
from __future__ import annotations

import json
import pathlib
import re
import sys
from pathlib import Path

W = Path('.')
NET = (W / 'CONSTRAINT_NETS.json').resolve()
DASH = {'−': '-', '–': '-', '—': '-', '‐': '-', '­': '-'}
CITE = re.compile(r'^(?P<file>[^\s:]+)\s*::\s*L(?P<line>\d+)\s*=\s*(?P<text>.+)$', re.S)


def norm(s: str) -> str:
    for a, b in DASH.items():
        s = s.replace(a, b)
    return ' '.join(s.split())


def window_at(lines, idx: int) -> str:
    """The checker's own rule: a heading cites its section, any other line cites +-3 lines."""
    if lines[idx].lstrip().startswith('#'):
        end = idx + 1
        while end < len(lines) and not lines[end].lstrip().startswith('#'):
            end += 1
        return norm(' '.join(lines[idx:end]))
    return norm(' '.join(lines[max(idx - 3, 0):idx + 4]))


# A digit welded to a letter is an exponent in a unit, not a quantity: the 2 in "mm2" and the 3
# in "cm3". Ten meniscus claims could only be matched by a window wide enough to contain some
# unrelated literal 2 -- which the section HEADING "## 2. Falsifier A" supplied -- so every one of
# them cited a whole section when its own table row was the source. Units are excluded from the
# number list; a quantity is not preceded by a letter.
NUM = re.compile(r'(?<![\d.A-Za-z\u00b2\u00b3])-?\d+[.,]?\d*')


def half_unit(lit: str) -> float:
    """Half of the last stated digit: 4.07 means 4.07 +- 0.005, 1800 means 1800 +- 0.5."""
    frac = len(lit.split('.')[1]) if '.' in lit else 0
    return 0.5 * 10 ** (-frac)


def num_match(lit: str, window_nums: list[str]) -> tuple[bool, str | None]:
    """A cited number matches a source number literally, or by rounding at the COARSER precision.

    HARVEST-E0065 cites 4.072; the source document states 4.07 and nothing finer. The claim is
    sound and the extra digit is the harvest's, not the source's, so an exact-literal test calls a
    correct edge unfindable while a loose numeric test would accept a genuinely different value.
    Rule: the two agree only if they agree to the last digit the COARSER one states, and when they
    do, the source's literal is the one that gets written back.
    """
    if lit in window_nums:
        return True, None
    try:
        a = float(lit.replace(',', '.'))
    except ValueError:
        return False, None
    for w in window_nums:
        if '.' not in w or len(w) >= len(lit):
            # The coarser side must still state a decimal. Without that guard "1" matched the
            # cited "1.3" (|1.3-1| <= 0.5) and the repair rewrote a Nachemson factor of 1.3-1.5
            # into 1-1. Sixth loose-matcher failure in this pipeline; the guard is the specific
            # thing that was missing, so it is stated rather than widened.
            continue
        try:
            b = float(w.replace(',', '.'))
        except ValueError:
            continue
        frac = len(w.split('.')[1])
        if f'{a:.{frac}f}' == f'{b:.{frac}f}' and abs(a - b) <= half_unit(w):
            return True, w
    return False, None


def carries(window: str, claimed: str) -> tuple[bool, dict]:
    # A hyphen between digits is a range, not a minus -- the same trap as in net_staleness.py.
    nums = NUM.findall(claimed)
    if claimed[:60] in window:
        return True, {}
    if not nums:
        return False, {}
    wnums = NUM.findall(window)
    fixes = {}
    for lit in nums[:4]:
        ok, src = num_match(lit, wnums)
        if not ok:
            return False, {}
        if src:
            fixes[lit] = src
    return True, fixes


def main() -> int:
    apply = '--apply' in sys.argv
    net = json.loads(NET.read_text())
    edges = net['bodytwin']['tissue_constraint_net']['edges']
    moved, unresolved, already = [], [], 0
    for e in edges:
        if e.get('evidence_unresolved'):
            continue
        m = CITE.match(str(e.get('evidence', '')).strip())
        if not m:
            continue
        f = pathlib.Path(m.group('file'))
        if not f.is_absolute():
            f = W / m.group('file')
        if not f.exists():
            continue
        lines = f.read_text(errors='replace').splitlines()
        claimed = norm(m.group('text'))
        cited = int(m.group('line')) - 1
        if 0 <= cited < len(lines) and carries(window_at(lines, cited), claimed)[0]:
            already += 1
            continue
        hits = [(i, fx) for i in range(len(lines))
                for ok, fx in [carries(window_at(lines, i), claimed)] if ok]
        if hits:
            # nearest hit to the original citation: the document grew, the claim did not move far
            best, fixes = min(hits, key=lambda h: abs(h[0] - cited))
            text = m.group('text')
            for lit, src in fixes.items():
                text = text.replace(lit, src)
            moved.append((e['id'], cited + 1, best + 1, len(hits), fixes))
            if apply:
                e['evidence'] = f"{m.group('file')} :: L{best + 1} = {text}"
                note = f"L{cited + 1} -> L{best + 1} on 2026-10-04; line drift"
                if fixes:
                    note += ('; cited precision the source does not state, corrected to the '
                             'source literal: ' + ', '.join(f'{a}->{b}' for a, b in fixes.items()))
                    for lit, src in fixes.items():
                        e['constraint'] = str(e.get('constraint', '')).replace(lit, src)
                else:
                    note += ', text unchanged'
                e.setdefault('evidence_relocated', []).append(note)
        else:
            unresolved.append((e['id'], m.group('file'), cited + 1))
            if apply:
                e['evidence_unresolved'] = True
                e['evidence_unresolved_note'] = (
                    f"cited text not found anywhere in {m.group('file')} on 2026-10-04 "
                    f"(searched every line with the heading/+-3 window rule); cite was L{cited + 1}")
    print(f"line citations already correct: {already}")
    print(f"relocated: {len(moved)}")
    for i, a, b, n, fx in moved:
        extra = ('  precision corrected ' + ', '.join(f'{k}->{v}' for k, v in fx.items())) if fx else ''
        print(f"  {i:<16} L{a} -> L{b}  ({n} line(s) carry the claim){extra}")
    print(f"text not in file: {len(unresolved)}")
    for i, f, a in unresolved:
        print(f"  {i:<16} {f} L{a}")
    if apply:
        NET.write_text(json.dumps(net, ensure_ascii=False, indent=2) + '\n')
        print("written")
    else:
        print("dry run; pass --apply to write")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())


# ---------------------------------------------------------------------------
# Second locator class: `file :: /json/pointer = free text`.
#
# Measured 2026-10-04 after teaching the staleness checker to read this format: 19 of the 24
# edges carrying it resolve their pointer perfectly and do NOT contain the cited claim. The first
# one examined by hand, HARVEST-E0092, cites /citations_verified_live_this_session/
# graziosi_2007_rca_vs_lad_doppler/numbers for a swine autoregulatory gain of 0.46 +- 0.11; that
# subtree is a 24-human intracoronary Doppler study, and the swine gain is in the SIBLING entry
# berwick_2012_autoregulatory_gain. Right file, right document, wrong citation key -- the harvest
# attached each claim to a neighbouring entry. Flagging all 19 as STALE would be literally true
# and useless, because the evidence exists one key away. So: search the document for the subtree
# that does carry the claim, report which, and move the pointer there.

def _walk_json(obj, path=''):
    yield path or '/', obj
    if isinstance(obj, dict):
        for k, v in obj.items():
            yield from _walk_json(v, f'{path}/{k}')
    elif isinstance(obj, list):
        for i, v in enumerate(obj):
            yield from _walk_json(v, f'{path}/{i}')


def find_pointer(doc, claimed: str):
    """The smallest subtree that carries the claim, so the pointer gains precision, not loses it."""
    best = None
    for path, sub in _walk_json(doc):
        rendered = norm(json.dumps(sub, ensure_ascii=False))
        ok, fixes = carries(rendered, claimed)
        if ok and (best is None or len(rendered) < best[2]):
            best = (path, fixes, len(rendered))
    return best


def companion_md(f: Path) -> Path:
    """`X_evidence.json` and `X.md` are one document in two files; the harvest confused them.

    Eleven pointer-cited edges carry a claim that is in neither the cited subtree nor any sibling:
    HARVEST-E0092's swine gain of 0.46 is on line 63 of MECHANISM_CORONARY_BLOOD_FLOW.md, cited
    there to Berwick 2012 (reference 5), while the pointer addresses the evidence JSON. So the
    second place to look is the other half of the same document, and a claim found there becomes
    an ordinary line citation.
    """
    stem = f.name[:-len('_evidence.json')] if f.name.endswith('_evidence.json') else f.stem
    return f.with_name(stem + '.md')


def relocate_pointers(edges, apply: bool) -> dict:
    """Repair pointer citations. A move must keep or gain precision: a claim whose numbers are
    only found by widening to an ancestor container is NOT repaired, because the ancestor carries
    every sibling's numbers and verifies nothing. Those are downgraded with that reason stated."""
    PTR = re.compile(r'^(?P<file>[^\s:]+)\s*::\s*(?P<ptr>/[^=]+?)\s*=\s*(?P<text>.+)$', re.S)
    out = {'ok': 0, 'repointed': [], 'to_markdown': [], 'downgraded': []}
    for e in edges:
        if e.get('evidence_unresolved'):
            continue
        m = PTR.match(str(e.get('evidence', '')).strip())
        if not m:
            continue
        f = pathlib.Path(m.group('file'))
        if not f.is_absolute():
            f = W / m.group('file')
        if not f.exists():
            continue
        try:
            doc = json.loads(f.read_text())
        except Exception:
            continue
        claimed = norm(m.group('text'))
        segs = [x for x in m.group('ptr').split('/') if x]
        cur, broke = doc, False
        for seg in segs:
            if isinstance(cur, dict) and seg in cur:
                cur = cur[seg]
            elif isinstance(cur, list) and seg.isdigit() and int(seg) < len(cur):
                cur = cur[int(seg)]
            else:
                broke = True
                break
        if not broke and carries(norm(json.dumps(cur, ensure_ascii=False)), claimed)[0]:
            out['ok'] += 1
            continue
        hit = find_pointer(doc, claimed)
        if hit and len([x for x in hit[0].split('/') if x]) >= len(segs):
            out['repointed'].append((e['id'], m.group('ptr'), hit[0]))
            if apply:
                e['evidence'] = f"{m.group('file')} :: {hit[0]} = {m.group('text')}"
                e.setdefault('evidence_relocated', []).append(
                    f"pointer {m.group('ptr')} -> {hit[0]} on 2026-10-04; the harvest attached the "
                    f"claim to a neighbouring citation entry in the same document")
            continue
        md = companion_md(f)
        found_md = None
        if md.exists():
            lines = md.read_text(errors='replace').splitlines()
            # Prefer the SMALLEST window that carries the claim. Ten meniscus claims first landed
            # on one line -- the heading of the section whose table holds all ten rows. True, but a
            # citation that fits every row of a table identifies none of them, so the narrow
            # three-line window around the claim's own row wins when one exists, and the section
            # heading is only the fallback.
            cands = [(len(w), i + 1) for i in range(len(lines))
                     for w in [window_at(lines, i)] if carries(w, claimed)[0]]
            found_md = min(cands)[1] if cands else None
        if found_md:
            out['to_markdown'].append((e['id'], f.name, md.name, found_md))
            if apply:
                e['evidence'] = f"{md} :: L{found_md} = {m.group('text')}"
                e.setdefault('evidence_relocated', []).append(
                    f"moved from {f.name} pointer {m.group('ptr')} to {md.name} L{found_md} on "
                    f"2026-10-04; the claim is in the prose half of the document, not the evidence JSON")
            continue
        why = ('the claim is not in the cited subtree, and its numbers are only found by widening '
               'to an ancestor container that carries every sibling citation'
               if hit else 'the claim is in neither the cited JSON nor its companion markdown')
        out['downgraded'].append((e['id'], why))
        if apply:
            e['evidence_unresolved'] = True
            e['evidence_unresolved_note'] = (
                f"{why}; searched every subtree of {f.name} and all of {md.name} on 2026-10-04 "
                f"(cite was {m.group('ptr')})")
    return out
