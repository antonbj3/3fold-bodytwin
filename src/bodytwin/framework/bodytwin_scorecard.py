#!/usr/bin/env python3
"""Scorecard over the graph: node counts by status, type and cluster.

Input: $BODYTWIN_GRAPH. Output: a summary table on stdout.
"""
import json, re
from pathlib import Path
import os
REPO = Path(os.environ.get("BODYTWIN_ROOT", os.getcwd()))
G = json.load(open(os.environ.get("BODYTWIN_GRAPH", REPO / "ANCHOR_GRAPH.json")))

def verdict(n):
    return ((n.get("cert_design") or {}).get("verify") or {}).get("verdict", "") or ""

_REFUTE_RE = re.compile(r"\bREFUTE[SD]?(_H\d+)?\b")

def tier(n):
    st = n.get("status", "OPEN"); v = verdict(n).upper()
    ev = bool(n.get("evidence"))
    cd = n.get("cert_design") or {}
    exec_note = (cd.get("execution_status") or "") + " " + ((cd.get("verify") or {}).get("residual") or "")
    if st == "PROVEN": return "MEASURED-PROVEN"
    if st == "ASSUMED" and ev: return "MEASURED-B"
    # ⚠ , two bugs found by an audit that TESTED this function against hand-verified nodes
    # instead of reading it. Both made a genuinely-refuted cell invisible in this scorecard:
    #
    # (1) ORDERING. The BLOCKED branch matched "gated" anywhere in exec_note — which includes
    #     verify.residual — and sat BEFORE the refutation branch. MITO-CELL's residual mentions an
    #     unrelated FUTURE candidate anchor as "dbGaP-gated", so that substring stole priority and
    #     bucketed a cleanly-refuted cell as BLOCKED (waiting on access) rather than NEEDS-REDESIGN.
    #     A refuted hypothesis does not stop being refuted because some other dataset is gated.
    #     Refutation is now decided FIRST; access-gating only classifies what is still live.
    # (2) LEXICAL. "REFUTED" as an exact substring missed the variant "REFUTE_H1", dropping that
    #     node into the 2773-entry OTHER/un-audited catch-all where nothing looks at it.
    #
    # Why this mattered: the graph's documented boot command is verdict-blind and marks 99.45% of
    # nodes actionable, so THIS scorecard is the only thing that would have flagged them.
    refuted = bool(_REFUTE_RE.search(v)) or "WEAKENED" in v or "RE-SCOPED" in v
    if "FENCED" in v or "HONEST-NEG" in v: return "FENCED"
    if refuted: return "NEEDS-REDESIGN"
    if st == "DEFERRED" or "BLOCK" in exec_note.upper() or "gated" in exec_note.lower(): return "BLOCKED"
    if v.startswith("REAL"): return "READY-TO-EXECUTE"
    if "PROPOSED" in v: return "PROPOSED"
    return "OTHER/un-audited"

ORDER = ["MEASURED-PROVEN", "MEASURED-B", "READY-TO-EXECUTE", "BLOCKED", "NEEDS-REDESIGN", "FENCED", "PROPOSED", "OTHER/un-audited"]
buckets = {t: [] for t in ORDER}
for n in G["nodes"]:
    buckets[tier(n)].append(n)

print("=" * 78)
print("BODYTWIN CERTIFICATION SCORECARD — measurement state  (design->execute->prove)")
print("=" * 78)
for t in ORDER:
    cells = buckets[t]
    if not cells: continue
    print(f"\n### {t}  ({len(cells)})")
    for n in cells:
        ev = "  ⟨evidence⟩" if n.get("evidence") else ""
        note = verdict(n)[:70]
        print(f"  {n['id']:<42} {note}{ev}")
print("\n" + "-" * 78)
print("tally: " + " · ".join(f"{t}={len(buckets[t])}" for t in ORDER if buckets[t]))
print(f"total cells: {len(G['nodes'])}")
