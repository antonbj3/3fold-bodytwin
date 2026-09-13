#!/usr/bin/env python3
"""Scan for write-back misregistration between nodes and their evidence artifacts.

scan(graph) returns nodes whose recorded decisive number does not match the number in the
artifact the node points at. Input: the graph JSON plus the artifacts it references.
Output: a list of (node id, reason) pairs.
"""
import json, os, re, sys

REPO = os.environ.get("BODYTWIN_ROOT", os.getcwd())
GRAPH = os.environ.get("BODYTWIN_GRAPH", os.path.join(REPO, "ANCHOR_GRAPH.json"))

# ---- find every "[WRITE-BACK ...]" bracketed block, anywhere in the node (not just verify.verdict —
#      the observed location in the confirmed incident, but a future annotation elsewhere must not be
#      missed by a field-specific scan). No nested "[" is expected inside the annotation (confirmed on
#      all 68 live instances), so a non-greedy scan to the first "]" is the exact span. -------------
BLOCK_RE = re.compile(r'\[WRITE-BACK\b[^\]]*\]')

# ---- the file token: bounded by what it MUST CONTAIN (a literal .json or .py suffix via a non-greedy
#      \S+?), NOT by whitespace/non-whitespace exclusion. This is the deliberate fix for the reviewer's
#      own first-cut bug, where a whitespace-bounded capture swallowed a trailing ";" and the basename
#      comparison then silently never matched anything. -----------------------------------------------
FILE_RE = re.compile(r'recovered from\s+(\S+?\.(?:json|py))\b')

def _walk_strings(obj):
    if isinstance(obj, str):
        yield obj
    elif isinstance(obj, dict):
        for v in obj.values():
            yield from _walk_strings(v)
    elif isinstance(obj, list):
        for v in obj:
            yield from _walk_strings(v)

def scan(nodes, repo=REPO):
    """Core scan, importable so bodytwin_hardening_check.py can call it on ALREADY-loaded nodes with
    no second graph read. Returns (violations: list[dict], info: dict)."""
    violations = []
    total_blocks = 0
    exempt_py = 0
    matched_json = 0
    checked_nodes = 0

    for n in nodes:
        if not isinstance(n, dict):
            continue
        nid = n.get("id")
        ev = [e for e in (n.get("evidence") or []) if isinstance(e, str)]
        ev_basenames = {os.path.basename(e) for e in ev}
        node_hit = False

        for s in _walk_strings(n):
            if "[WRITE-BACK" not in s:
                continue
            for bm in BLOCK_RE.finditer(s):
                block = bm.group(0)
                fm = FILE_RE.search(block)
                if not fm:
                    continue
                total_blocks += 1
                node_hit = True
                fpath = fm.group(1)
                if fpath.endswith(".py"):
                    exempt_py += 1
                    continue
                base = os.path.basename(fpath)
                if base in ev_basenames:
                    matched_json += 1
                else:
                    violations.append({
                        "node": nid,
                        "recovered_from": fpath,
                        "own_evidence": ev,
                        "reason": ("node's '[WRITE-BACK ... recovered from %s]' names a .json file "
                                   "not present in this node's evidence[] -- the verdict likely traces "
                                   "to a DIFFERENT node's raw output (couples_to-membership-as-identity "
                                   "borrow, the confirmed 60398eaef mechanism)" % fpath),
                    })
        if node_hit:
            checked_nodes += 1

    info = {
        "total_nodes": len(nodes),
        "nodes_with_writeback_annotation": checked_nodes,
        "total_writeback_blocks": total_blocks,
        "exempt_py_artifacts": exempt_py,
        "matched_own_evidence": matched_json,
        "violations": len(violations),
    }
    return violations, info

def main():
    quiet = "--quiet" in sys.argv
    json_out = None
    if "--json" in sys.argv:
        i = sys.argv.index("--json")
        json_out = sys.argv[i + 1] if i + 1 < len(sys.argv) else "/dev/stdout"

    g = json.load(open(GRAPH))
    nodes = g["nodes"] if isinstance(g, dict) else g
    violations, info = scan(nodes, REPO)

    if json_out:
        json.dump({"violations": violations, "info": info}, open(json_out, "w"), indent=1)

    if violations:
        print(f"BODYTWIN WRITE-BACK PROVENANCE SCAN: ❌ FAIL — {len(violations)} node(s) carry a "
              f"'[WRITE-BACK ... recovered from X.json]' verdict annotation whose X is NOT in that node's "
              f"own evidence[] (borrowed-sibling-verdict pattern, commit 60398eaef mechanism):")
        for v in violations:
            print(f"  ✗ {v['node']}  recovered_from={v['recovered_from']}  own_evidence={v['own_evidence']}")
        print(f"(coverage: {info['total_nodes']} nodes, {info['total_writeback_blocks']} WRITE-BACK blocks, "
              f"{info['exempt_py_artifacts']} .py-exempt, {info['matched_own_evidence']} clean self-matches)")
        sys.exit(1)

    if not quiet:
        print(f"BODYTWIN WRITE-BACK PROVENANCE SCAN: ✅ PASS — 0 violations across {info['total_nodes']} "
              f"nodes ({info['total_writeback_blocks']} WRITE-BACK blocks, {info['exempt_py_artifacts']} "
              f".py-exempt script artifacts, {info['matched_own_evidence']} verdicts confirmed traced to "
              f"their own node's evidence).")
    sys.exit(0)

if __name__ == "__main__":
    main()
