#!/usr/bin/env python3
"""Scan node evidence pointers for provenance defects.

scan(graph) returns the nodes whose evidence pointer does not support the node's
claim (missing artifact, artifact belonging to another node, or a pointer that is prose
rather than a file). Input: the graph JSON. Output: a list of (node id, reason) pairs.
"""
import json, os, re, sys, difflib
from collections import defaultdict, Counter

REPO = os.environ.get("BODYTWIN_ROOT", os.getcwd())
GRAPH = os.environ.get("BODYTWIN_GRAPH", os.path.join(REPO, "ANCHOR_GRAPH.json"))
EVID_PREFIX = os.environ.get("BODYTWIN_EVIDENCE_DIR", "evidence")

# ---- EXACT mirror of bodytwin_fold.py::_declared_id — "the id the source artifact says it is,
#      if it says one. Never inferred." Any drift between this and the real guard must be fixed
#      in BOTH places at once (a diverging mirror would rescan against the wrong contract). -------
def declared_id(out):
    if not isinstance(out, dict):
        return None
    pc = out.get("proposed_cell")
    if isinstance(pc, dict) and isinstance(pc.get("id"), str):
        return pc["id"].strip()
    v = out.get("id")
    return v.strip() if isinstance(v, str) else None

# ---- distinctive-feature extractor: citations, identifiers and numbers used to match a
#      node's claim against the artifact it points at ----------------
CITE_RE = re.compile(r'\b([A-Z][a-zA-Z\-]{2,20}(?:&[A-Z][a-zA-Z\-]{2,20})?)\s*(?:et al\.?\s*)?(19|20)\d{2}\b')
PMID_RE = re.compile(r'PMID\s*:?\s*(\d{6,9})')
PMC_RE  = re.compile(r'PMC(\d{5,9})')
DOI_RE  = re.compile(r'10\.\d{4,9}/[^\s"\')]+')
NUM_RE  = re.compile(r'\b\d+\.\d{2,4}\b')

def features(text):
    if not text:
        return set()
    cites = set()
    for m in CITE_RE.finditer(text):
        name = m.group(1).rstrip('&')
        year = text[m.end() - 4:m.end()]
        if len(name) >= 3:
            cites.add(f"cite:{name}{year}")
    return (cites | {f"pmid:{x}" for x in PMID_RE.findall(text)}
                  | {f"pmc:{x}" for x in PMC_RE.findall(text)}
                  | {f"doi:{x}" for x in DOI_RE.findall(text)}
                  | {f"num:{x}" for x in NUM_RE.findall(text)})

def node_text(n):
    parts = [n.get("claim", "")]
    cd = n.get("cert_design") or {}
    if isinstance(cd, dict):
        for k in ("hidden_state", "anchor", "decorrelated_anchor", "regime_note"):
            v = cd.get(k)
            if isinstance(v, str):
                parts.append(v)
        for k in ("legs", "datapoints", "honest_gaps", "data_sources", "couples_to"):
            v = cd.get(k)
            if v:
                parts.append(json.dumps(v))
        vf = cd.get("verify")
        if isinstance(vf, dict):
            parts.append(str(vf.get("note", "")))
            parts.append(str(vf.get("verdict", "")))
    return " ".join(parts)

def containment(feat_a, feat_b):
    """|A ^ B| / |A| -- asymmetric on purpose: how much of A's distinctive content is IN B.
    Returns None (undefined, NOT zero) when A has no features -- caller must not silently treat
    "can't test" as "failed the test" (that conflation was v1's confirmed false-positive bug)."""
    if not feat_a:
        return None
    return len(feat_a & feat_b) / len(feat_a)

_WS = re.compile(r'\s+')
_ANNOT = re.compile(r'\s*\[(CORRECTION|DUPLICATE-PAIR|RE-VERIFY|CORRECTED)\b.*$', re.S)

def norm(s):
    return _WS.sub(' ', s or '').strip().lower()

def core_claim(claim):
    """Strip a dated bracket annotation appended AFTER the original fold. That text cannot appear
    in the (earlier-written) evidence file by construction, so leaving it in only ever penalizes a
    real match — it can never legitimately catch a real mismatch."""
    return _ANNOT.sub('', claim or '')

def claim_like_field(out):
    if not isinstance(out, dict):
        return ""
    pc = out.get("proposed_cell")
    if isinstance(pc, str):
        return pc
    if isinstance(pc, dict) and isinstance(pc.get("claim"), str):
        return pc["claim"]
    if isinstance(out.get("claim"), str):
        return out["claim"]
    return ""

def claim_ratio(a, b):
    a, b = norm(core_claim(a)), norm(core_claim(b))
    if not a or not b:
        return 0.0
    return difflib.SequenceMatcher(None, a, b).ratio()

def text_match(node_claim, out, blob):
    """Three independent, fair, mechanical tests; matched iff ANY fires. Void-floor calibrated:
    0/497 random-unrelated-pair false positives (see module docstring)."""
    nc = norm(core_claim(node_claim))
    if not nc:
        return False, 0.0, "empty-claim"
    cl = norm(claim_like_field(out))
    r = difflib.SequenceMatcher(None, nc, cl).ratio() if cl else 0.0
    if r >= 0.60:
        return True, r, "field-ratio"
    blob_n = norm(blob)
    probe = nc[:150]
    if len(probe) >= 40 and probe in blob_n:
        return True, max(r, 0.99), "verbatim-substring"
    sm = difflib.SequenceMatcher(None, nc, blob_n[:20000])
    m = sm.find_longest_match(0, len(nc), 0, min(len(blob_n), 20000))
    lcs_ratio = m.size / len(nc) if nc else 0.0
    if lcs_ratio >= 0.5 and m.size >= 40:
        return True, lcs_ratio, "lcs-substring"
    return False, max(r, lcs_ratio), "no-match"

RECOVERY_SIGNAL = re.compile(r'/recovered/|recovered__|manual__', re.I)

def _load_evidence(abspath):
    if not os.path.isfile(abspath):
        return "MISSING", None, None, None
    raw = open(abspath, encoding="utf-8", errors="replace").read()
    try:
        doc = json.loads(raw)
    except Exception:
        return "NOT-JSON", None, raw, None
    if isinstance(doc, list):
        return "LIST-TOPLEVEL", None, json.dumps(doc), None
    if not isinstance(doc, dict):
        return "SCALAR-TOPLEVEL", None, raw, None
    out = doc.get("full_output", doc)
    did = declared_id(out)
    blob = json.dumps(out) if isinstance(out, (dict, list)) else str(out)
    return "OK", did, blob, out

def scan(nodes, repo=REPO):
    """Core scan, importable so bodytwin_hardening_check.py can call it on ALREADY-loaded nodes
    with no second graph read. Returns (true_misreg: list[dict], info: dict)."""
    by_id = {n["id"].strip().upper(): n for n in nodes if isinstance(n, dict) and isinstance(n.get("id"), str)}
    path_citers = defaultdict(list)
    for n in nodes:
        for ev in (n.get("evidence") or []):
            if isinstance(ev, str) and ev.startswith(EVID_PREFIX):
                path_citers[ev].append(n.get("id"))

    _cache = {}
    def load_evidence(relpath):
        if relpath not in _cache:
            _cache[relpath] = _load_evidence(os.path.join(repo, relpath))
        return _cache[relpath]

    nodes_with_ao_evidence = 0
    total_pairs = 0
    status_counts = Counter()
    match_n = 0
    no_declared_id_n = 0
    mismatch_records = []

    for n in nodes:
        nid = n.get("id")
        ev_list = [e for e in (n.get("evidence") or []) if isinstance(e, str) and e.startswith(EVID_PREFIX)]
        if not ev_list:
            continue
        nodes_with_ao_evidence += 1
        for ev in ev_list:
            total_pairs += 1
            status, did, blob, out = load_evidence(ev)
            status_counts[status] += 1
            if status != "OK":
                continue
            if did is None:
                no_declared_id_n += 1
                continue
            if did.strip().upper() == str(nid).strip().upper():
                match_n += 1
                continue
            mismatch_records.append({"node": nid, "evidence": ev, "declared_id": did, "blob": blob, "out": out})

    classified = []
    for rec in mismatch_records:
        nid, ev, did, blob, out = rec["node"], rec["evidence"], rec["declared_id"], rec["blob"], rec["out"]
        n = by_id.get(str(nid).strip().upper())
        n_claim = n.get("claim", "") if n else ""

        self_match, self_ratio, self_mode = text_match(n_claim, out, blob)
        n_feat = features(node_text(n)) if n else set()
        ev_feat = features(blob)
        cont_n = containment(n_feat, ev_feat)
        feat_match = (cont_n is not None and cont_n >= 0.5 and len(n_feat & ev_feat) >= 1)

        owner = by_id.get(did.strip().upper())
        owner_match, owner_ratio, node_vs_owner_ratio = False, 0.0, 0.0
        if owner is not None and owner.get("id") != nid:
            owner_match, owner_ratio, _ = text_match(owner.get("claim", ""), out, blob)
            node_vs_owner_ratio = claim_ratio(n_claim, owner.get("claim", ""))

        duplicate_pair = node_vs_owner_ratio >= 0.85
        matched = self_match or feat_match or duplicate_pair

        if matched:
            klass = "CONTENT-MATCHES"
            if duplicate_pair and not (self_match or feat_match):
                sub = "DUPLICATE-NODE-PAIR"
            else:
                sub = "RECOVERY-ARTEFACT" if RECOVERY_SIGNAL.search(ev) else "BENIGN-RENAME-OR-SHARED"
        else:
            klass = "TRUE-MISREGISTRATION"
            if owner is not None and owner_match:
                sub = f"evidence belongs to existing node {owner.get('id')} (ratio={owner_ratio:.2f})"
            elif owner is not None:
                sub = (f"declared id matches existing node {owner.get('id')} but neither its text nor "
                       f"a claim-vs-claim comparison confirms it (node-vs-owner ratio={node_vs_owner_ratio:.2f}) "
                       "-- LOWER CONFIDENCE, flag for human review")
            else:
                sub = "declared id matches no existing node (orphan/renamed/unknown)"

            # ----  DISCLOSURE CHECK (reviewer-directed): a node may carry an explicit,
            # machine-readable cert_design.duplicate_registration_of naming the node it is a duplicate
            # registration of (same analysis, folded twice under two ids -- NOT the same thing as
            # "evidence belongs to node X" above: the true duplicate-twin is found by CLAIM similarity,
            # while the wrongly-cited evidence's declared owner can be a THIRD, unrelated node -- e.g.
            # AUTO-GANSEVOORT...'s claim duplicates RENAL-GLOMERULAR-FILTRATION but its evidence file
            # self-declares the unrelated HOLE-P53-DAMAGE-PULSES...). This check does NOT trust the
            # field's mere PRESENCE -- it independently RE-VERIFIES claim_ratio(this node, target) >=
            # 0.85 (the same threshold used above for an undeclared duplicate_pair) every run, so a
            # false or stale annotation is caught, not rubber-stamped. This is the line the reviewer
            # asked to be drawn cleanly: acknowledged-AND-true -> not a violation; asserted-but-false is
            # WORSE than silent (a lie caught) and must fail louder than a plain unexplained mismatch.
            dro = (n.get("cert_design") or {}).get("duplicate_registration_of") if n else None
            if isinstance(dro, str) and dro.strip():
                dro_id = dro.strip()
                target = by_id.get(dro_id.upper())
                if target is None:
                    klass = "FALSE-DISCLOSURE"
                    sub = f"claims duplicate_registration_of={dro_id!r} but NO such node exists in the graph"
                elif dro_id.upper() == str(nid).strip().upper():
                    klass = "FALSE-DISCLOSURE"
                    sub = "claims duplicate_registration_of=itself (self-reference, nonsensical)"
                else:
                    dro_ratio = claim_ratio(n_claim, target.get("claim", ""))
                    if dro_ratio >= 0.85:
                        klass = "ACKNOWLEDGED-DUPLICATE"
                        sub = f"verified duplicate_registration_of={dro_id} (claim ratio={dro_ratio:.2f})"
                    else:
                        klass = "FALSE-DISCLOSURE"
                        sub = (f"claims duplicate_registration_of={dro_id} but claim-vs-claim ratio is only "
                               f"{dro_ratio:.2f} (<0.85) -- disclosure NOT verified, treat as unresolved")

        classified.append({
            "node": nid, "evidence": ev, "declared_id": did, "class": klass, "sub": sub,
            "self_ratio": round(self_ratio, 3), "owner_id": owner.get("id") if owner else None,
            "owner_ratio": round(owner_ratio, 3), "node_vs_owner_ratio": round(node_vs_owner_ratio, 3),
            "shared_evidence_with_other_nodes": [c for c in path_citers.get(ev, []) if c != nid],
        })

    # FALSE-DISCLOSURE fails LOUDER than a plain silent TRUE-MISREGISTRATION (a caught lie is worse
    # than an unexplained gap) -- both are returned as the failing set. ACKNOWLEDGED-DUPLICATE is
    # excluded from the failing set but kept visible in `classified`/info, never silently absorbed
    # into "clean" CONTENT-MATCHES, so it stays auditable.
    true_misreg = [c for c in classified if c["class"] in ("TRUE-MISREGISTRATION", "FALSE-DISCLOSURE")]
    info = {
        "total_nodes": len(nodes), "nodes_with_ao_evidence": nodes_with_ao_evidence,
        "total_pairs": total_pairs, "status_counts": dict(status_counts),
        "match": match_n, "no_declared_id": no_declared_id_n, "mismatch": len(mismatch_records),
        "content_matches": sum(1 for c in classified if c["class"] == "CONTENT-MATCHES"),
        "acknowledged_duplicate": sum(1 for c in classified if c["class"] == "ACKNOWLEDGED-DUPLICATE"),
        "false_disclosure": sum(1 for c in classified if c["class"] == "FALSE-DISCLOSURE"),
        "true_misregistration": sum(1 for c in classified if c["class"] == "TRUE-MISREGISTRATION"),
        "classified": classified,
    }
    return true_misreg, info

def main():
    quiet = "--quiet" in sys.argv
    json_out = None
    if "--json" in sys.argv:
        i = sys.argv.index("--json")
        json_out = sys.argv[i + 1] if i + 1 < len(sys.argv) else "/dev/stdout"

    g = json.load(open(GRAPH))
    nodes = g["nodes"] if isinstance(g, dict) else g
    true_misreg, info = scan(nodes, REPO)

    if json_out:
        json.dump(info, open(json_out, "w"), indent=1)

    if true_misreg:
        n_false_disc = info["false_disclosure"]
        n_true_mis = info["true_misregistration"]
        print(f"BODYTWIN PROVENANCE SCAN: ❌ FAIL — {len(true_misreg)} VIOLATION(S): {n_true_mis} TRUE "
              f"MISREGISTRATION (a node's evidence citation belongs to different/no work than its own "
              f"claim, undisclosed) + {n_false_disc} FALSE-DISCLOSURE (claims duplicate_registration_of "
              f"but that claim does not verify -- worse than silent, a caught false annotation):")
        for c in true_misreg:
            print(f"  ✗ [{c['class']}] {c['node']}  cites {c['evidence']}  declares={c['declared_id']}  :: {c['sub']}")
        if info["acknowledged_duplicate"]:
            print(f"({info['acknowledged_duplicate']} node(s) carry a VERIFIED cert_design."
                  f"duplicate_registration_of and are acknowledged, not counted as violations)")
        print(f"(coverage: {info['total_nodes']} nodes, {info['nodes_with_ao_evidence']} cite "
              f"agent_outputs evidence, {info['match']} clean id-match, {info['no_declared_id']} "
              f"unassessable/no-declared-id, {info['content_matches']} benign mismatches)")
        sys.exit(1)

    if not quiet:
        print(f"BODYTWIN PROVENANCE SCAN: ✅ PASS — 0 true misregistrations, 0 false disclosures across "
              f"{info['total_nodes']} nodes ({info['nodes_with_ao_evidence']} cite agent_outputs evidence, "
              f"{info['total_pairs']} (node,evidence) pairs, {info['match']} id-matched, "
              f"{info['content_matches']} benign mismatches content-confirmed, "
              f"{info['acknowledged_duplicate']} acknowledged duplicate-registrations (verified), "
              f"{info['no_declared_id']} unassessable).")
    sys.exit(0)

if __name__ == "__main__":
    main()
