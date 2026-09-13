#!/usr/bin/env python3
"""The single definition of the DESIGN-STUB and DUP-SIGNAL fold skips.

skip_reason(payload, graph_nodes) returns a string naming why a submission should not be
folded (a template echo with no measured content, or a duplicate of an existing node), or
None when the submission may proceed. Both the fold adapter and any cron driver call this
one function so the predicate cannot drift between call sites.
"""
import re

def flat(v):
    """Coerce proposed_cell / decorrelated_anchor (str, dict, or list) to text for matching."""
    if isinstance(v, str):
        return v
    if isinstance(v, dict):
        for k in ("design", "claim", "cell", "summary", "test"):
            if isinstance(v.get(k), str):
                return v[k]
        return " ".join(flat(x) for x in v.values())
    if isinstance(v, list):
        return " ".join(flat(x) for x in v)
    return "" if v is None else str(v)

# DUP-SIGNAL ( origin): an agent that concluded "this is a duplicate, create NO new node" must
# not be reified into a phantom node (measured: 3 such 'AUTO-NO-NEW-...' / 'AUTO-DO-NOT-CREATE-...'
# phantoms). Pattern and 600-char/lower-case window copied verbatim from the cron script.
_DUP_SIGNAL_RE = re.compile(
    r"(no new (top|node|graph)|no new top-level|do not create|don't create"
    r"|would duplicate|is a duplicate|duplicate-|duplicate --|already covered"
    r"|already established|already in the graph|dup_of|topic is saturated"
    r"|no dedicated node needed|not warrant a (new|separate) node)")

# DESIGN-STUB ( origin): a proposed cell whose CLAIM is imperative next-step prose
# ("Build/Pull/Add X") AND that cites NO data_sources is a TODO, not a researched cell (measured 63 such
# stubs already folded, inflating the node count with intentions). Pattern copied verbatim.
_STUB_CLAIM_RE = re.compile(
    r"^(build|pull|add|locate|test whether|do not|do n|wire|split|extend|make|"
    r"write|use|set|compute|derive)\b", re.I)

def is_dup_signal(out):
    """True iff the agent's text (claim/verdict/proposed_cell) signals this would duplicate an
    existing node. `out` is the extracted cell dict (fold_one's `out = doc.get("full_output", doc)`,
    same shape autofold_pending's `out` already has)."""
    if not isinstance(out, dict):
        return False
    sig = (str(out.get("claim") or "") + " " + str(out.get("verdict") or "")
           + " " + flat(out.get("proposed_cell")))[:600].lower()
    # ⚠ NARROWED , MEASURED. Phrasing alone was refusing real work: a corpus scan of all
    # 3414 agent_outputs found 19 live dup-signal hits, 12 of them (63%) carrying a real structured
    # payload. The concrete casualty was FIX-IRON-DUPLICATE-LEG-DEINFLATION- — an executed
    # meta-audit finding with 13 measured datapoints — refused for the sole reason that a cell ABOUT
    # a duplicate pair must use the word "duplicate". That is the second guard today to block the
    # audit of the thing it guards against (cf. the template-echo guard refusing the template census
    # for quoting a placeholder). A cell that SIGNALS duplication and carries NO measured payload is
    # the true positive; one that signals it and carries legs+datapoints is a finding about
    # duplication. Same discriminator is_design_stub already relies on — content, not phrasing.
    return bool(_DUP_SIGNAL_RE.search(sig)) and _payload_size(out) == 0

def _cert(out):
    pcx = out.get("proposed_cell")
    cd = (pcx.get("cert_design") if isinstance(pcx, dict) else None) or out.get("cert_design")
    return cd if isinstance(cd, dict) else {}

def _payload_size(out):
    """How much MEASURED content the cell carries, independent of how its claim is phrased.

    ★ WHY THIS EXISTS (measured). The data_sources clause below was checking two paths
    that are structurally dead -- 0 of 3064 graph nodes carry a top-level `data_sources`, because the
    house convention nests it at proposed_cell.cert_design.data_sources, exactly where
    bodytwin_fold._lift_cert_design already looks. So the clause never fired and the rule was
    effectively "imperative claim" alone, which is why two attempts to reproduce it failed.
    Adding the correct path fixes that -- but on its own it OVER-flags: three AUTO- cells whose claim
    is a terse `seed: <topic>` marker carry 12-13 REAL datapoints each (named quantities with values
    and units), and skipping them would have discarded that content.
    Measured discriminator: all 58 confirmed stubs have payload EXACTLY 0 (max 0, not merely a low
    median). So gating on payload costs 0/58 sensitivity and spares every cell that carries measured
    values behind a terse claim. Content, not phrasing, is the thing that separates a reified
    intention from a researched cell.
    """
    cd = _cert(out)
    return len(cd.get("datapoints") or []) + len(cd.get("legs") or [])

def is_design_stub(out):
    """True iff the claim is imperative next-step prose ('Build/Pull/Add X'), the cell cites no
    data_sources at ANY of the three real paths, AND it carries no measured payload -- a TODO, not a
    researched cell. See _payload_size for why the third condition is not optional."""
    if not isinstance(out, dict):
        return False
    claim = str(out.get("claim") or "").strip()
    ds = out.get("data_sources")
    if not (isinstance(ds, list) and ds):
        pcx = out.get("proposed_cell")
        if isinstance(pcx, dict) and isinstance(pcx.get("data_sources"), list):
            ds = pcx.get("data_sources")
    if not (isinstance(ds, list) and ds):                 # the path that actually holds citations
        ds = _cert(out).get("data_sources")
    return (bool(_STUB_CLAIM_RE.match(claim))
            and not (isinstance(ds, list) and ds)
            and _payload_size(out) == 0)

def skip_reason(out):
    """Single entry point fold_one() calls: None if the cell should fold, else a short tag naming which
    rule fired ('dup-signal' checked first, matching the cron script's check ORDER)."""
    if is_dup_signal(out):
        return "dup-signal"
    if is_design_stub(out):
        return "design-stub"
    return None
