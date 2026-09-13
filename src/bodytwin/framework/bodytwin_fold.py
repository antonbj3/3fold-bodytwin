#!/usr/bin/env python3
"""Fold adapter: turn a contributor's cell JSON into a graph node.

Input: a raw result file containing an MT-JSON-BEGIN {...} MT-JSON-END
envelope (or a bare JSON object) with a proposed_cell block, plus a cell id.
Output: an upserted node in the dependency graph (default $BODYTWIN_GRAPH, else
ANCHOR_GRAPH.json in the working directory); blocked submissions are appended to a
fold queue JSONL instead.

A submission without an executed-result receipt is booked as a HYPOTHESIS, i.e. status
OPEN. A submission carrying {verdict, rerun:{command, artifact, checks}} is booked as a
measured RESULT. Claim text is passed through a dialect sanitizer before it is stored.

CLI:
  python bodytwin_fold.py <raw.json> <CELL-ID> [cluster] [layer]
  python bodytwin_fold.py --dry <raw.json> <CELL-ID>      # gate only, no write

Requires the fold gate modules (fold_gate, fold_gate_v2) on the import path; they are
not part of this repository. Everything except fold_one() imports and
runs without them.
"""
import json, os, re, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
REPO = os.environ.get("BODYTWIN_ROOT", os.getcwd())
BODYTWIN_GRAPH = os.environ.get("BODYTWIN_GRAPH", os.path.join(REPO, "ANCHOR_GRAPH.json"))
FOLD_QUEUE = os.environ.get("BODYTWIN_FOLD_QUEUE", os.path.join(REPO, "fold_queue.jsonl"))
# Guarded graph writer (atomic os.replace + flock + reload-verify).
from bodytwin_graph_io import write_graph, graph_lock

# Shared with any cron driver so the skip predicate cannot drift between call sites.
from bodytwin_fold_skip_rules import skip_reason
from mt_id_alias import ID_KEYS as _ID_KEYS, declared_id as _declared_id, normalize_id_alias

# Remove configured domain-specific wording before storing a public claim.
_JARGON_EXTRA = [p for p in os.environ.get("BODYTWIN_JARGON", "").split(",") if p.strip()]
# Claim text is rewritten before it is stored so that a submission cannot carry a private
# dialect into the public graph. Substrate-specific patterns are supplied via $BODYTWIN_JARGON.
_JARGON = [
    (r"\bx0\s*\+\s*V_d\s*@\s*z\b", "a low-rank latent parametrization"),
    (r"\bV_d\s*@\s*z\b", "a low-rank latent code"),
    (r"\bV_d\b", "the latent basis"),
    (r"\bin-manifold\b", "in-latent-space"),
    (r"\bmanifold\b", "latent space"),
    (r"[σ]_?min\b|\bsigma_min\b", "min-sensitivity"),
    (r"\breduced-reps?\b", "low-rank"),
    (r"\bgradient-projection\b", "gradient method"),
] + [(p, "[redacted]") for p in _JARGON_EXTRA]

# The sanitizer must be a superset of the detector in bodytwin_hardening_check: read the detector's
# own pattern out of that module's source (no import, so the check does not run at import time).
def _detector_pattern():
    try:
        src = open(os.path.join(HERE, "bodytwin_hardening_check.py")).read()
        m = re.search(r'CS_DIALECT\s*=\s*re\.compile\(\s*r"([^"]+)"', src)
        return re.compile(m.group(1), re.I) if m else None
    except Exception:
        return None
_DETECTOR = _detector_pattern()

def _sanitize(claim):
    for pat, rep in _JARGON:
        claim = re.sub(pat, rep, claim, flags=re.I)
    if _DETECTOR is not None:                      # closure guarantee: scrubber ⊇ detector, by construction
        claim = _DETECTOR.sub("the inherited method", claim)
    return re.sub(r"\s{2,}", " ", claim).strip()

def _load_gates():
    """Import the separate fold-gate dependencies only when folding is requested."""
    from fold_gate import fold_gate
    from fold_gate_v2 import fold_gate_v2
    return fold_gate, fold_gate_v2

_REQUIRED_ANY = ("proposed_cell", "decorrelated_anchor", "honest_gaps", "data_sources")
_SENTINEL = re.compile(r"MT-JSON-BEGIN\s*(\{.*?\})\s*MT-JSON-END", re.S)

def _unicode_escape(s):
    """Decode a backslash-escaped envelope; return None when decoding fails."""
    try:
        return s.encode("utf-8", "ignore").decode("unicode_escape")
    except Exception:
        return None

def _largest_balanced(text):
    """largest top-level balanced {...} that parses (string-aware; ignores braces in prose)."""
    best = None
    for i, ch in enumerate(text):
        if ch != '{':
            continue
        depth = 0; ins = False; esc = False
        for j in range(i, len(text)):
            c = text[j]
            if esc: esc = False; continue
            if c == '\\': esc = True; continue
            if c == '"': ins = not ins
            if ins: continue
            if c == '{': depth += 1
            elif c == '}':
                depth -= 1
                if depth == 0:
                    try:
                        o = json.loads(text[i:j + 1])
                        if isinstance(o, dict) and len(o) > len(best or {}): best = o
                    except Exception: pass
                    break
    return best

def extract_cell_json(raw_text):
    """Extract the last valid sentinel envelope, a JSON fence, or a balanced object.

    Require at least four fields and a recognized research field, then normalize
    identifier aliases. Return None when no acceptable record is present."""
    if not isinstance(raw_text, str):
        return None
    cand = None
    # Use the last parseable envelope so quoted earlier records cannot replace the final result.
    for m in reversed(list(_SENTINEL.finditer(raw_text))):
        span = m.group(1)
        for variant in (span, _unicode_escape(span)):     # raw span first, then unicode-escape decode
            if variant is None:
                continue
            try:
                c = json.loads(variant)
            except Exception:
                continue
            if isinstance(c, dict):
                cand = c
                break
        if cand is not None:
            break
    if cand is None:
        m = re.search(r"```json\s*(\{.*?\})\s*```", raw_text, re.S)
        if m:
            try: cand = json.loads(m.group(1))
            except Exception: cand = None
    if cand is None:
        cand = _largest_balanced(raw_text)
    if not isinstance(cand, dict):
        return None
    # The >=4-key thinness test counts keys at the level it happens to be handed. A correctly
    # WRAPPED cell puts its 4+ keys one level down, so `{"proposed_cell": {...}, "result_paths":
    # [...]}` = 2 keys and was refused as "too thin" — measured  on the 5.6 KB cell
    # NODE-0-ENERGY-HEAT-ROUTES-CLOSURE-CLOTHED-EXT, which only survived because autofold's
    # repair fallback happened to catch it. Count where the cell actually IS.
    _inner = cand.get("proposed_cell") if isinstance(cand.get("proposed_cell"), dict) else None
    if len(cand) < 4 and not (_inner and len(_inner) >= 4):
        return None                       # refuse: too thin to be a real cell
    if not any(k in cand for k in _REQUIRED_ANY):
        return None                       # refuse: no research payload -> don't fabricate a claim
    # EXTRACTION BOUNDARY: mirror an aliased id (`cell_id`) onto `id` here, so no consumer
    # downstream has to know the alias exists. See scripts/mt_id_alias.py.
    return normalize_id_alias(cand)

# Root-cause fix (HOLE-CELL-CRITERION-FALSIFIABILITY-CENSUS): the claim field was hard-cut at
# 400 chars, truncating 1776/2382 (74.6%) of all cells MID-WORD. Measured consequence: a cell's pass
# criterion was severed at "FAIL i" — one letter before its condition — so "no criterion visible" became
# indistinguishable from "no criterion written", and every text-based census over this field read a
# systematically incomplete string. Raised well past the observed natural claim length.
# ⚠ , second incident, same shape: raising 400→4000 did NOT fix it — 16 nodes sat at EXACTLY
# 4000, i.e. still severed, and one had its decisive quantified negative (a ceiling-vs-noise-floor
# comparison) cut mid-sentence, so the folded node CONTRADICTED its own raw evidence file and anyone
# reading only the graph missed the result entirely. The bug was never the number; it was the SILENCE.
# A truncated claim must be impossible to mistake for a complete one, so truncation now announces itself.
CLAIM_MAX = 16000
_TRUNC_MARK = " …[CLAIM TRUNCATED AT CLAIM_MAX — the full text is in this cell's evidence file; do not read this claim as complete]"

def _cap_claim(s):
    """Cap a claim, but never silently: a cut claim says so and points at the untruncated source."""
    if len(s) <= CLAIM_MAX:
        return s
    return s[:CLAIM_MAX - len(_TRUNC_MARK)] + _TRUNC_MARK

def _cap_field(s, limit, name):
    """Cap ANY content field the same way a claim is capped — announced, never silent.

    : claim had this protection and no other field did. Measured on the live graph, three
    fold-path fields carried SILENT fixed slices, and roughly a third of every populated instance sat
    exactly at its ceiling: verdict[:600] (71 nodes), regime_note[:800] (227 of 614 populated, 37.0%),
    decorrelated_anchor[:800] (747 of 1964, 38.0%). Union: 898 of 2833 nodes, 31.7% of the whole graph,
    carrying at least one invisibly-severed field.

    decorrelated_anchor is the worst of the three to lose, because it is the field the method relies on
    to tell a genuine external anchor from a tautology gate — a reader cannot judge an anchor they can
    only see the first 800 characters of, and nothing told them it was cut.

    Raising a ceiling alone does not fix this; the same bug returns invisibly the day any field exceeds
    the new one. So the announcement is the fix and the ceiling is only a backstop.
    """
    if not isinstance(s, str) or len(s) <= limit:
        return s
    mark = (" …[%s TRUNCATED AT %d CHARS — full text is in this cell's evidence file; "
            "do not read this field as complete]" % (name.upper(), limit))
    return s[:max(0, limit - len(mark))] + mark

def _claim_from(out):
    """Derive a compact claim from flat or nested result fields."""
    for k in ("proposed_cell", "overall", "decorrelated_anchor"):
        v = out.get(k)
        # : len>40 alone let a bare TOPIC SLUG win. Measured: a record whose proposed_cell was
        # the 42-char string "hippocampus_trisynaptic_and_cognitive_map" returned here immediately, so
        # its genuinely rich sibling fields (trisynaptic_computation, replay_mechanism, honest_gaps, …)
        # were never reached — the node landed carrying a label where a claim belonged. This is the
        # INVERSE of the scalars-only loss fixed in tier 3 below: there a full record rendered as empty,
        # here a slug renders as complete, and a confidently-wrong claim is the more dangerous of the two
        # because nothing about it invites a second look. A real claim is prose: require whitespace, or
        # enough length that a single identifier token is implausible.
        if isinstance(v, str) and len(v) > 40 and (" " in v.strip() or len(v) > 80):
            return _cap_claim(_sanitize(v))    # scrub CS reduced-rep dialect in the LIVE fold path
    # Search nested claim fields before falling back to the complete structured record.
    for src in (out.get("proposed_cell"), out.get("cert_design"), out):
        if not isinstance(src, dict):
            continue
        for k in ("claim", "overall", "summary", "verdict"):
            v = src.get(k)
            if isinstance(v, str) and len(v) > 40:
                # DOMINANCE GUARD (): returning the first named key that clears 40 chars silently
                # discards every sibling field. Measured: a flat record with mechanism=900 chars and
                # verdict=120 returned the 120 and dropped the 900. Only take this shortcut when the field
                # really IS the record; otherwise fall through to the full longest-first stitch below, which
                # still leads with this field. Nested cells (proposed_cell/cert_design) keep the old
                # behaviour — there the claim genuinely is the claim.
                if src is not out:
                    return _cap_claim(_sanitize(v))
                _tot = sum(len(str(x)) for x in out.values() if isinstance(x, (str, int, float)))
                if _tot and len(v) >= 0.6 * _tot:
                    return _cap_claim(_sanitize(v))
                break
    # PREVIEW-COLLAPSE FIX (an audit):
    # the old fallback stitched `k: str(v)[:80]` over only the FIRST 3 scalar fields, producing a ~281-char
    # label-joined PREVIEW while the evidence file held the full record. Measured damage: 103 nodes, one of
    # them showing 700,892 chars on disk vs 11,504 in the graph (1.64% visible). Two silent cuts caused it —
    # 80 chars per field and 3 fields — neither announced. Fix: keep EVERY scalar field, longest first (the
    # real claim is almost always the longest), full text, and let the existing CLAIM_MAX/_cap_claim net do
    # the one truncation that DOES announce itself. Ad hoc schemas (grand_challenge_knee, honest_gaps,
    # hidden_state, mechanism, …) are covered because nothing is keyed on a field-name vocabulary here.
    scalars = [(k, str(v)) for k, v in out.items()
               if isinstance(v, (str, int, float)) and k not in ("topic",)]
    # : scalars-only was a SILENT TOTAL LOSS for one whole record shape. Measured: 9 records
    # whose full_output has ZERO top-level scalars besides the excluded "topic" — every other field is a
    # dict/list of research sections. All three tiers correctly find nothing, so this returned the
    # sentinel while 6–16 kB of PMID-cited research sat one level down, untouched. Those nodes then read
    # "(no claim extractable)" in the graph and were RECOMMENDED FOR DELETION on that basis; reading the
    # evidence file one hop deeper showed 9/9 were real. A field that renders a full record as "empty"
    # does not merely lose data — it invites its destruction.
    # Fix reuses this file's _flatten_text (defined below, already used for decorrelated_anchor):
    # when nothing scalar survives, harvest the structured fields instead. Same longest-first join, same
    # CLAIM_MAX/_cap_claim announced-truncation net. Cannot regress the 121 nodes where tier-3 already
    # works, because this only fires when `scalars` is empty.
    if not scalars:
        scalars = [(k, _flatten_text(v)) for k, v in out.items()
                   if isinstance(v, (list, dict)) and k not in ("topic",) and _flatten_text(v)]
    scalars.sort(key=lambda kv: len(kv[1]), reverse=True)
    if scalars:
        return _cap_claim(_sanitize("\n\n".join(f"{k.upper()}: {v}" for k, v in scalars)))
    return "(no claim extractable)"

def _as_list(v):
    if v is None: return []
    if isinstance(v, list):
        return [x if isinstance(x, str) else (json.dumps(x) if isinstance(x, (dict, list)) else str(x)) for x in v]
    if isinstance(v, dict):
        return [f"{k}: {json.dumps(val) if isinstance(val,(dict,list)) else val}" for k, val in v.items()]
    return [str(v)]

def _as_datapoints(v):
    """Normalize dictionary or list measurements into a datapoint list."""
    if isinstance(v, dict):
        return [{"name": str(k), "value": val} for k, val in v.items()]
    if isinstance(v, list):
        return v
    # second stranding shape (measured , HOLE-STRANDED-NUMBERS-CENSUS): 125 nodes emit
    # derived_values as a prose STRING; dropping it loses the content entirely, so keep it verbatim.
    if isinstance(v, str) and v.strip():
        return [{"name": "summary", "value": v.strip()[:2000]}]
    return []

def _flatten_text(v):
    """Render a field that may arrive as a str, a dict of sub-fields, or a list into one readable string.

    Agents legitimately emit e.g. decorrelated_anchor as {"leg1_method":..., "convergence":..., "refs":[...]}.
    str(dict) would store an unreadable repr, and the old isinstance(v, str) test dropped it entirely."""
    if v is None:
        return ""
    if isinstance(v, str):
        return v
    if isinstance(v, dict):
        return "; ".join(f"{k}: {_flatten_text(x)}" for k, x in v.items() if x not in (None, "", [], {}))
    if isinstance(v, (list, tuple)):
        return "; ".join(_flatten_text(x) for x in v if x not in (None, "", [], {}))
    return str(v)

def _lift_cert_design(out, layer=""):
    """Preserve flat and nested result fields in the canonical cert_design structure."""
    pc = out.get("proposed_cell") if isinstance(out.get("proposed_cell"), dict) else {}
    # Root-cause fix (HOLE-FOLD-PATH-SILENTLY-EMPTIES-CERT-DESIGN-ARRAYS): when the record is
    # NOT wrapped in proposed_cell, pc stays {} and cd stays {}, so pick() searched pc / cd / out and never
    # descended into a TOP-LEVEL cert_design — every array (legs/datapoints/couples_to/conflicts/
    # honest_gaps) was silently dropped and an EMPTY SHELL written, while fold still reported
    # "wrote_node": true and the hardening check still reported PASS. Measured by controlled A/B probe on
    # Flat and wrapped submissions carry the same content at different nesting levels; a flat one
    # must be lifted, not silently accepted as an empty node.
    if not pc and isinstance(out.get("cert_design"), dict):
        pc = out                                        # flat submission: the record IS the cell
    cd = pc.get("cert_design") if isinstance(pc.get("cert_design"), dict) else {}   # cert fields can be nested here
    # Include nested verify dictionaries when resolving verdicts and supporting fields.
    _verifies = [s.get("verify") for s in (pc, cd, out) if isinstance(s.get("verify"), dict)]
    def pick(*keys):
        for src in (pc, cd, out, *_verifies):
            for k in keys:
                v = src.get(k)
                if v not in (None, "", [], {}):
                    return v
        return None
    couples = pick("couples_to") or []
    if isinstance(couples, str): couples = [couples]
    elif isinstance(couples, dict): couples = [str(x) for x in couples.values()]
    elif isinstance(couples, list): couples = [x if isinstance(x, str) else str(x) for x in couples]
    lifted = {
        "hidden_state": _sanitize(str(pick("hidden_state") or "")),
        "occluded_truth": True,
        "legs": [_sanitize(s) for s in _as_list(pick("legs"))],
        "anchor": _sanitize(str(pick("anchor") or "")),
        "composition_mode": str(pick("composition_mode") or ""),
        "couples_to": couples,
        "regime_note": _cap_field(_sanitize(str(pick("regime_note") or "")), 4000, "regime_note"),
        "datapoints": _as_datapoints(pick("derived_values", "datapoints")),
        # Preserve these fields when they are nested under proposed_cell.
        "decorrelated_anchor": _cap_field(_sanitize(_flatten_text(pick("decorrelated_anchor"))), 8000, "decorrelated_anchor"),
        "conflicts": _as_list(pick("conflicts")),
        "honest_gaps": _as_list(pick("honest_gaps")),
        "data_sources": _as_list(pick("data_sources")),
        "couples_to_inherited": "none",
        # Carry the supplied grade and verdict; use explicit defaults only when absent.
        "bodytwin_grade": _grade_from(pick("bodytwin_grade"), pick("verdict")),
        # : the cap was [:600] and was SLICING THE BOTTOM-LINE FIELD MID-WORD. Measured across the
        # graph: 169 nodes sat at exactly 600 chars, and for the 122 whose pre-fold source survived, a total of
        # 59,181 characters had been cut (median 384, max 1722). That is not cosmetic — verdict is the
        # fold_gate-facing field downstream cells read, and a correction appended to a long verdict was being
        # silently discarded, which is exactly how 4 kidney cells kept asserting a superseded elasticity in
        # their verdict while their claim/conflicts carried the fix. Raised to 6000; the cap exists only to
        # stop a runaway blob from becoming the verdict, and 6000 is far above any real one.
        "verify": {"verdict": _cap_field(_sanitize(_flatten_text(pick("verdict"))), 6000, "verdict")
                              or "SEED-DESIGN (fold_gate ALLOW as hypothesis; cited literature, not executed)",
                   "layer": layer or "biology"},
    }
    # Preserve additional cert_design keys that are not in the canonical field list.
    _known = set(lifted) | {"verdict", "derived_values", "verify", "cert_design", "evidence",
                            "claim", "id", "status", "type", "cluster", "risk", "load",
                            "depends_on", "actionable_by", "proposed_cell"}
    # ⚠ The extras scan must cover EVERY level pick() covers, plus a nested proposed_cell — measured the
    # same day: BUILD-FROM-CITED-SOURCE-BATCH13 carried its cert_design one level deeper and lost
    # void_floor_pass_rates while a sibling cell folded the same hour kept all six of its extras. A
    # partial fix to a nesting-level defect is the defect.
    _inner = pc.get("proposed_cell") if isinstance(pc.get("proposed_cell"), dict) else {}
    _inner_cd = _inner.get("cert_design") if isinstance(_inner.get("cert_design"), dict) else {}
    for src in (cd, pc, _inner_cd, _inner, *[v for v in _verifies if isinstance(v, dict)]):
        if not isinstance(src, dict):
            continue
        for k, v in src.items():
            if k in _known or v in (None, "", [], {}):
                continue
            if isinstance(v, str):
                lifted[k] = _cap_field(_sanitize(v), 20000, k)
            elif isinstance(v, (list, dict)):
                lifted[k] = v            # structure is the payload here — keep it, do not flatten
            else:
                lifted[k] = v
            _known.add(k)
    return lifted

# Resolve grade and verdict together so they cannot contradict each other.
_DECISIVE = re.compile(r"\b(PASS|FAIL|REFUTED|CONFIRMED|WEAKENED|MIXED|PARTIAL|MEASURED|EXECUTED|"
                       r"RESOLVED|FALSIFIED|ROBUST|FLIP|CLEAN-NEGATIVE|HYPOTHESIS-awaiting-QC)\b", re.I)
_TEMPLATE = re.compile(r"<\s*(one word|verdict|anchor|claim|number|value)\s*>", re.I)

def _grade_from(explicit, verdict):
    """Grade must be consistent with the verdict sitting next to it, never independently defaulted."""
    v = _flatten_text(verdict) or ""
    if _TEMPLATE.search(v):                       # an unfilled placeholder is NOT a measured result
        return "TEMPLATE-UNFILLED"
    if explicit:
        return str(explicit)
    if v.strip() and not v.lstrip().upper().startswith("SEED-DESIGN") and _DECISIVE.search(v):
        return "EXECUTED-VERDICT-PRESENT"         # derived, not asserted — the verdict carries a decision
    return "SEED-DESIGN"

def _is_shell(node):
    """A node whose cert_design carries neither legs nor an anchor is an un-filled shell.

    Also (): a node carrying the literal '(no claim extractable)' sentinel is a shell REGARDLESS of
    legs/anchor. Those two fields are rescued by the depth-aware pick(), so a degraded node looks non-shell and
    a re-fold with a FIXED extractor would silently no-op via the idempotent-skip path — leaving the degraded
    state frozen while the repairer believes it succeeded. Treat the sentinel as overwritable."""
    cd = node.get("cert_design") or {}
    if "(no claim extractable)" in str(node.get("claim") or ""):
        return True
    return (not (cd.get("legs") or [])) and (not str(cd.get("anchor") or "").strip())

# Accept the declared identifier aliases consistently with the extraction and provenance checks.

def fold_one(raw_path, cell_id, cluster="", layer="", dry=False):
    doc = json.load(open(raw_path))
    out = doc.get("full_output", doc)
    # SKIP GATE (generalised from the autofold-cron-only checks — see the import above):
    # applied BEFORE the provenance guard and BEFORE building a submission, exactly like the cron script
    # applies it before minting an id — a cell we're about to discard doesn't need its id-provenance
    # validated. Fires regardless of `dry`, so a dry-run also reports what would be skipped.
    _reason = skip_reason(out)
    if _reason is not None:
        return {"cell": cell_id, "decision": "SKIP", "skip_reason": _reason, "wrote_node": False,
                "note": ("dup-signal: input text signals this would duplicate an existing node"
                          if _reason == "dup-signal" else
                          "design-stub: imperative next-step prose with no data_sources -- a TODO, not "
                          "a researched cell") +
                         " (GRAPH-STUB-SKIP-EXISTS-ONLY-IN-THE-UNCERTIFIED-FOLD-PATH-)"}
    # Refuse placeholder and instruction-shaped envelopes before an idempotent write can persist them.
    _PLACEHOLDER = re.compile(r"^\s*(\.\.\.|…|<[^>]{1,40}>)\s*$")
    _INSTRUCTION = re.compile(r"<\s*(one word|one sentence|real ids|ids|anchor|numbers|gaps)\b", re.I)
    _pc = out.get("proposed_cell", out) if isinstance(out, dict) else {}
    if isinstance(_pc, dict):
        _cd = _pc.get("cert_design") or {}
        # ⚠ The instruction pattern is only evidence of a stub in a SHORT field. A long claim that
        # QUOTES a placeholder is legitimate — the census cell reporting this very defect was refused
        # by the first version of this guard for quoting the text it was reporting on. A guard that
        # blocks the audit of the thing it guards against is too broad. Whole-field match always
        # counts; instruction-shaped text only counts when the field is too short to be a real claim.
        _ph = [k for k, v in list(_pc.items()) + list(_cd.items() if isinstance(_cd, dict) else [])
               if isinstance(v, str) and (_PLACEHOLDER.match(v)
                                          or (_INSTRUCTION.search(v) and len(v) < 200))]
        if _ph:
            raise ValueError(
                "TEMPLATE-ECHO REFUSED: %s carries literal '...' placeholders in %s — this is the "
                "input-format template, not a result record. It reaches the extractor "
                "when the input contains only a template. Supply the "
                "completed result envelope before folding." % (raw_path, ", ".join(sorted(_ph))))
    #  PROVENANCE GUARD. A prior repair pass re-extracted each node's claim FROM its cited
    # evidence, gated ONLY on "produces >200 more chars than what's stored" — with no check that the
    # evidence belonged to the node being written. Where an evidence link pointed at another cell's
    # artifact, that cell's claim was copied wholesale under the wrong id. Measured: 179 nodes carry
    # that pass's tag, 40 have an evidence file declaring a DIFFERENT id, 7 are confirmed verbatim
    # cross-cell swaps — id and hidden_state intact and correct, claim belonging to someone else, and
    # 3 of the 7 wearing a confident PASS-flavoured verdict. That is the dangerous direction: an empty
    # claim invites a second look, a confident wrong one does not, and nothing in the schema check can
    # see it because the schema is perfectly satisfied.
    # The artifact states its own id one JSON key away. Length gain is not provenance — ask.
    _decl = _declared_id(out)
    if _decl and _decl.upper() != str(cell_id).upper():
        raise ValueError(
            "PROVENANCE MISMATCH: %s declares id %r but is being folded as %r. Refusing — this is the "
            "exact shape that produced 7 cross-cell claim swaps. If the evidence link is genuinely "
            "shared (a comparison reference, which is legitimate and common for adjudication cells), "
            "fold from the artifact that belongs to this cell instead." % (raw_path, _decl, cell_id))
    # a literature-grounded cell is a HYPOTHESIS until it has a forward-model + local re-run artifact
    submission = {
        "hypothesis": True,                              # => fold_gate ALLOW as='hypothesis'
        "cell": cell_id,
        "claim": _claim_from(out),
        "evidence": [os.path.relpath(raw_path, REPO)],
        "note": "SEED-DESIGN: cited literature anchor, not yet executed against data",
    }
    # Apply the forward gate; hypothesis submissions do not require executed-result coverage.
    gate = fold_gate_v2(json.dumps(submission), base_dir=REPO)
    decision = gate.get("decision")
    if decision != "ALLOW":
        rec = {"cell": cell_id, "raw": os.path.relpath(raw_path, REPO), "decision": decision,
               "missing": gate.get("missing")}
        if not dry:
            with open(FOLD_QUEUE, "a") as f: f.write(json.dumps(rec) + "\n")
        return {"cell": cell_id, "decision": decision, "wrote_node": False, "missing": gate.get("missing")}

    # Upsert an allowed hypothesis into the canonical graph without duplicating its identifier.
    if dry:
        return {"cell": cell_id, "decision": "ALLOW", "as": gate.get("as"), "wrote_node": False, "dry": True}
    # ⚠ RACE FIX : this read-decide-write cycle MUST be inside ONE lock. Locking only the
    # write (the old behaviour) serialized the writes but not the decisions: two concurrent folds each
    # read the same snapshot, each appended their own node to their own copy, and the second write
    # replaced the file wholesale — silently erasing the first fold's node, undetectable by
    # reload_verify() since it compares the file against the writer's count. graph_lock() is NOT
    # reentrant, so every write_graph() inside this block passes _locked=True.
    with graph_lock():
        g = json.load(open(BODYTWIN_GRAPH))
        nodes = g if isinstance(g, list) else (g.get("nodes") or g.get("cells"))
        nodes_is_list = isinstance(nodes, list)   # the CONTAINER of nodes, not g itself
        cert = _lift_cert_design(out, layer)      # preserve the complete result structure
        rich_incoming = bool(cert["legs"]) or bool(str(cert["anchor"]).strip())
        by_id = {n.get("id"): n for n in nodes if isinstance(n, dict)}

        def _agent_out_paths(ev):
            """Return the source JSON paths that identify the folded analysis."""
            return {e for e in (ev or []) if isinstance(e, str) and "agent_outputs" in e}

        # Refuse repeated evidence paths with a matching claim prefix even when the identifier differs.
        _inc_out = _agent_out_paths(submission["evidence"])
        if cell_id not in by_id and _inc_out:
            _inc_claim = (submission["claim"] or "")[:200]
            for n in nodes if nodes_is_list else nodes.values():
                if not isinstance(n, dict):
                    continue
                if _agent_out_paths(n.get("evidence")) & _inc_out and (n.get("claim") or "")[:200] == _inc_claim:
                    return {"cell": cell_id, "decision": "ALLOW", "wrote_node": False,
                            "note": f"same analysis already in graph as {n.get('id')} (evidence-path dedup)"}

        if cell_id in by_id:
            ex = by_id[cell_id]
            # UPSERT policy: FILL an empty shell (the fold used to strand this content), but NEVER
            # clobber a node that already carries real cert content. Idempotent once filled -> a later
            # re-fold of the same cell sees a non-shell and skips, so this cannot loop.
            if _is_shell(ex) and rich_incoming:
                ex["cert_design"] = cert
                if cert["regime_note"]: ex["regime_note"] = cert["regime_note"]
                if len(submission["claim"]) > 20: ex["claim"] = submission["claim"]
                for e in submission["evidence"]:
                    if e not in (ex.get("evidence") or []): ex.setdefault("evidence", []).append(e)
                write_graph(g, BODYTWIN_GRAPH, _locked=True)   # atomic + flock-guarded + reload-verify (bodytwin_graph_io)
                return {"cell": cell_id, "decision": "ALLOW", "wrote_node": True, "note": "filled shell (lifted cert_design)"}
            return {"cell": cell_id, "decision": "ALLOW", "wrote_node": False, "note": "already in graph (idempotent skip)"}
        node = {"id": cell_id, "claim": submission["claim"], "type": "EMPIRICAL",
                "status": "OPEN",                                # CS-enum: designed, NOT measured
                "evidence": submission["evidence"], "depends_on": [], "load": 0,
                "risk": "MED", "regime_note": cert["regime_note"], "actionable_by": "bodytwin-reviewer",
                "cluster": cluster, "cert_design": cert}
        if nodes_is_list: nodes.append(node)
        else: nodes[cell_id] = node
        write_graph(g, BODYTWIN_GRAPH, _locked=True)   # atomic + flock-guarded + reload-verify (bodytwin_graph_io)
        return {"cell": cell_id, "decision": "ALLOW", "as": gate.get("as"), "wrote_node": True}

def assert_landed(expected_id, evidence_path=None):
    """★ THE DESTINATION CHECK () — the actual antidote to cell substitution.

    A node count that moves after a fold proves A cell landed, not that YOURS did: the
    substitution bug measured tonight had `fold_one` report `wrote_node: true` for a real,
    well-formed node — it just wasn't the node the caller thought it was folding. The only
    thing that catches that is asking the graph, by the id the CALLER expected, not the id
    the (possibly wrong) extraction produced.

    Raises LOUDLY (never returns a soft False) on either failure mode:
      - expected_id absent from the graph entirely -> the fold silently landed nothing, or
        landed under a different id.
      - present, but its evidence does NOT include the path the caller folded FROM -> the
        node with this id already existed and was NOT updated by this fold (a same-id
        collision hid a no-op), so "landed" would be a false positive of a different shape.
    """
    g = json.load(open(BODYTWIN_GRAPH))
    nodes = g if isinstance(g, list) else (g.get("nodes") or g.get("cells"))
    node = None
    for n in (nodes if isinstance(nodes, list) else nodes.values()):
        if isinstance(n, dict) and n.get("id") == expected_id:
            node = n
            break
    if node is None:
        raise AssertionError(
            "DESTINATION CHECK FAILED: expected id %r is NOT in %s. A count moving elsewhere "
            "is not proof; this id specifically did not land." % (expected_id, BODYTWIN_GRAPH))
    if evidence_path is not None:
        rel = os.path.relpath(evidence_path, REPO) if os.path.isabs(evidence_path) else evidence_path
        if rel not in (node.get("evidence") or []):
            raise AssertionError(
                "DESTINATION CHECK FAILED: %r is in the graph but its evidence %r does NOT "
                "include %r -- this fold did not touch it (stale node, not a landed cell)."
                % (expected_id, node.get("evidence"), rel))
    return node

if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if a != "--dry"]
    dry = "--dry" in sys.argv
    if len(args) < 2:
        print(__doc__); sys.exit(1)
    raw, cell_id = args[0], args[1]
    cluster = args[2] if len(args) > 2 else ""
    layer = args[3] if len(args) > 3 else ""
    raw_abs = raw if os.path.isabs(raw) else os.path.join(REPO, raw)
    r = fold_one(raw_abs, cell_id, cluster, layer, dry=dry)
    print(json.dumps(r, indent=1))
    print(f"[isolation] target graph = {os.path.relpath(BODYTWIN_GRAPH, REPO)} (NEVER data/ANCHOR_GRAPH.json)")
    if not dry and r.get("wrote_node"):
        assert_landed(cell_id, raw_abs)   # ⚠ destination check: a moved count is not proof it's YOURS
        print(f"[destination-check] PASS — {cell_id} confirmed in graph, evidence matches.")
