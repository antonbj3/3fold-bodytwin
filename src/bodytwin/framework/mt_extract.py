#!/usr/bin/env python3
"""Extract the cell JSON from a contributor's raw output.

extract_cell_json(text) finds an MT-JSON-BEGIN {...} MT-JSON-END envelope, or a bare JSON
object, in a raw transcript or file, handles one level of backslash escaping, and returns
the parsed dict (None when nothing parseable is present).
"""
import json
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mt_id_alias import declared_id, normalize_id_alias

SPAN = re.compile(r"MT-JSON-BEGIN\s*(.*?)\s*MT-JSON-END", re.S)

def _variants(s, max_rounds=4):
    """Every way the payload might need unwrapping, cheapest first.

    ★ UNESCAPING IS ITERATIVE, NOT SINGLE-SHOT. Measured : a span captured
    from a transcript OF a transcript arrives DOUBLE-escaped, so one round of
    json.loads('"'+s+'"') leaves `{\\"proposed_cell\\":` — which fails with exactly the
    same error as the un-unescaped form, so a single-round extractor reports "no
    parseable span" on a payload that is entirely intact. Peel until it parses or
    stops changing.
    """
    cur = s
    for _ in range(max_rounds):
        yield cur
        try:
            nxt = json.loads('"' + cur + '"')
        except Exception:
            break
        if nxt == cur:
            break
        cur = nxt
    try:
        yield s.encode("utf-8", "ignore").decode("unicode_escape", "ignore")
    except Exception:
        pass

def _is_cell(d):
    """A cell payload, wrapped or flat.

    ★ WHY THIS IS NOT JUST `"proposed_cell" in d`. Measured : an agent emitted a
    fully intact cell with id/status/claim/cert_design at the TOP LEVEL, no wrapper. The
    extractor required the wrapper both to accept a parse AND to locate the object, so it
    printed "NO PARSEABLE MT-JSON SPAN FOUND" — the exact output it gives for an agent that
    produced nothing at all. A whole result was one keystroke from being discarded as an
    empty run, and nothing downstream could have told the difference. The fold path already
    lifts a top-level cert_design when the wrapper is absent; this makes the EXTRACTOR agree
    with it, so the two ends of the pipeline accept the same shapes.

    ★ THIRD MEMBER OF THE SAME CLASS, measured . A wrapped, fully intact cell declared
    its id as `proposed_cell.cell_id` instead of `proposed_cell.id`. _is_cell ACCEPTED it (the
    wrapper is there) — and main() still printed `id=None`, while the fold path's four independent
    `.get("id")` lookups each rejected it, so no artifact file and no node were ever written and
    the run was byte-indistinguishable from an agent that produced nothing. The flat branch below
    was the same trap one level down: `"id" in d` is a KEY test, and the key varies. The id lookup
    is now ONE shared definition (scripts/mt_id_alias.py) instead of a literal key spelled out at
    each site. The rule this encodes is the one this file already carries for parse tolerances,
    applied to KEY NAMES: any site that asks "what is this cell's id" must ask the SAME function.
    """
    if not isinstance(d, dict):
        return False
    if isinstance(d.get("proposed_cell"), dict):
        return True
    return "claim" in d and declared_id(d) is not None      # flat cell (id or cell_id)

def _normalize(d):
    """Always hand back the wrapped shape AND a usable `id`, so every caller sees one shape."""
    if isinstance(d, dict) and not isinstance(d.get("proposed_cell"), dict):
        d = {"proposed_cell": d}
    return normalize_id_alias(d)

def _parse(s):
    for v in _variants(s):
        for strict in (True, False):             # strict=False tolerates raw newlines
            try:
                return json.loads(v, strict=strict)
            except Exception:
                continue
    return None

def _by_object_anchor(raw):
    """Fallback: ignore the delimiters, find the OBJECT.

    ★ WHY THE DELIMITERS ARE NOT TRUSTWORTHY. Measured : the literal string
    MT-JSON-BEGIN also appears in the agent's PROMPT and planning text inside the
    transcript, and the real payload can be split across JSONL records. A non-greedy
    BEGIN...END match then captures a span that starts in one record and ends in
    another, so the "extracted" text is transcript metadata with a fragment of the cell
    glued on — and it fails to parse for a reason that looks exactly like a truncated
    agent. Anchoring on `"proposed_cell"` and letting raw_decode find the object's
    end is immune to both problems: trailing junk is simply ignored.

    ★ AND IT MUST WORK PER-LINE, NOT PER-FILE. Measured  on an 860 KB transcript:
    unescaping the WHOLE file fails, because the file is JSONL — real newlines and record
    boundaries make it an invalid JSON string as a whole. So only the raw variant got
    searched, and in the raw variant the marker is the ESCAPED \\"proposed_cell\\", which
    the bare search never matches. The payload was intact and the extractor reported
    nothing found. Each JSONL record is one line, so the line is the right unit to unescape.
    ★ AND A JSONL LINE IS A RECORD, NOT A STRING. json.loads('"' + line + '"') fails on the
    leading brace, so the escaped payload inside it never gets unwrapped. The line must be
    PARSED as JSON and its string VALUES searched -- that is where the agent's text lives,
    already unescaped by the parser. Measured on an 860 KB transcript whose payload was
    entirely intact while the extractor reported nothing found.
    """
    best = None
    dec = json.JSONDecoder(strict=False)
    ANCHORS = ('"proposed_cell"', '"claim"')   # flat cells have no wrapper to anchor on

    def _strings(o):
        if isinstance(o, str):
            yield o
        elif isinstance(o, dict):
            for v in o.values():
                yield from _strings(v)
        elif isinstance(o, list):
            for v in o:
                yield from _strings(v)

    for line in raw.splitlines():
        if not any(a.strip('"') in line for a in ANCHORS):
            continue
        candidates = [line]
        try:
            candidates.extend(_strings(json.loads(line)))   # JSONL record -> its own text
        except Exception:
            pass
        for cand in candidates:
            for text in _variants(cand):
                for anchor in ANCHORS:
                    pos = 0
                    while True:
                        i = text.find(anchor, pos)
                        if i < 0:
                            break
                        pos = i + 1
                        # walk left over candidate opening braces: a flat cell's "claim" may sit
                        # far from its own "{", unlike the wrapper key which is adjacent to it.
                        j = i
                        while True:
                            j = text.rfind("{", 0, j)
                            if j < 0:
                                break
                            try:
                                obj, _ = dec.raw_decode(text, j)
                            except Exception:
                                continue
                            if _is_cell(obj):
                                best = obj                # keep scanning: take the LAST
                            break
    return best

def extract_cell(path):
    """Return the LAST parseable MT-JSON object in the file, or None.

    Last rather than first: an agent that retries emits several spans and the final one
    is the corrected payload. Taking the first silently folds a superseded draft.
    """
    raw = open(path, encoding="utf-8", errors="replace").read()
    for s in SPAN.findall(raw)[::-1]:
        d = _parse(s)
        if _is_cell(d):
            return _normalize(d)
    got = _by_object_anchor(raw)
    return _normalize(got) if got else None

def main():
    if len(sys.argv) < 3:
        print(__doc__)
        return 2
    d = extract_cell(sys.argv[1])
    if not d:
        print("NO PARSEABLE MT-JSON SPAN FOUND")
        return 1
    cell = d.get("proposed_cell", d)
    json.dump(d, open(sys.argv[2], "w"), ensure_ascii=False, indent=1)
    print(f"{cell.get('id')}  claim={len(cell.get('claim',''))} chars -> {sys.argv[2]}")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
