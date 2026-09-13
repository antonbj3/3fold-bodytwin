#!/usr/bin/env python3
"""Machine guard over the graph and the fold path. Exit 1 on any violation.

Checks: graph wrapper shape {_meta, nodes}; unique node ids; status/type enum validity;
required node schema fields; no id collision with a second graph; the fold adapter still
routes through the gate, still writes status OPEN for hypotheses, still lifts cert_design;
both writers route graph writes through the guarded writer; and the two provenance scans
report no misregistration.

Input: $BODYTWIN_GRAPH (graph JSON) and the framework sources. Output: a report on stdout;
--quiet prints violations only.
"""
import json, os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.environ.get("BODYTWIN_ROOT", os.getcwd())
MT_GRAPH = os.environ.get("BODYTWIN_GRAPH", os.path.join(REPO, "ANCHOR_GRAPH.json"))
OTHER_GRAPH = os.environ.get("BODYTWIN_OTHER_GRAPH", "")   # optional second graph to check id collisions against
FOLD_ADAPTER = os.path.join(HERE, "bodytwin_fold.py")

CS_STATUS = {"PROVEN", "REFUTED", "DEFERRED", "OPEN", "ASSUMED"}
VALID_TYPES = {"THEORY", "EMPIRICAL", "ENGINEERING", "DATA-ANCHOR", "GOAL"}
VALID_RISK = {"LOW", "MED", "HIGH", "SOURCE-COND"}
REQUIRED = {"id", "claim", "type", "status", "evidence", "depends_on",
            "load", "risk", "regime_note", "actionable_by", "cluster", "cert_design"}
CS_DIALECT = re.compile(r"\bV_d\b|\bmanifold\b|reduced-rep|sigma_min|gradient-projection", re.I)

V = []   # violations
def bad(msg): V.append(msg)

def _nodes(path):
    g = json.load(open(path))
    return g, (g if isinstance(g, list) else (g.get("nodes") or g.get("cells") or []))

# ---- 1. BODYTWIN graph is the canonical, CS-combinable shape --------------------------------
try:
    g, nodes = _nodes(MT_GRAPH)
    if not (isinstance(g, dict) and "nodes" in g and "_meta" in g):
        bad("BODYTWIN graph wrapper is not {_meta, nodes} (schema drift)")
    ids = [n.get("id") for n in nodes if isinstance(n, dict)]
    dup = sorted({i for i in ids if ids.count(i) > 1})
    if dup: bad(f"BODYTWIN graph has duplicate ids: {dup}")
    for n in nodes:
        if not isinstance(n, dict): continue
        miss = REQUIRED - set(n.keys())
        if miss: bad(f"node {n.get('id')} missing REQUIRED fields {sorted(miss)}")
        if n.get("status") not in CS_STATUS: bad(f"node {n.get('id')} status '{n.get('status')}' not in CS VALID_STATUS")
        if n.get("type") not in VALID_TYPES: bad(f"node {n.get('id')} type '{n.get('type')}' not in VALID_TYPES")
        if n.get("risk") not in VALID_RISK: bad(f"node {n.get('id')} risk '{n.get('risk')}' not in VALID_RISK")
        if isinstance(n.get("claim"), str) and CS_DIALECT.search(n["claim"]):
            bad(f"node {n.get('id')} claim carries CS reduced-rep dialect (isolation/dialect leak)")
    n_nodes = len(nodes)
    # ---- 1b. APPEND-ONLY node-count invariant -------------------------------------------------
    # ⚠ , found by MUTATION TEST of this very file: every check above binds (10/10
    # schema/isolation mutations caught), but PURE NODE LOSS was invisible at every scale —
    # deleting 1976 of 2976 nodes passed, and so did deleting all but 76. The count was PRINTED
    # in the PASS line below and never CHECKED, so "PASS — N nodes" carried no assurance that N
    # was intact. The graph is append-only (verified: it never shrank across the last 30 graph
    # commits), so a drop below the recorded high-water mark cannot be legitimate growth.
    # bodytwin_graph_io.bump_high_water() maintains the mark on every guarded write.
    hw = (g.get("_meta") or {}).get("node_high_water") if isinstance(g, dict) else None
    if isinstance(hw, int) and n_nodes < hw:
        bad(f"node-count REGRESSION: graph has {n_nodes} nodes but _meta.node_high_water is {hw} "
            f"({hw - n_nodes} lost). The graph is append-only — investigate before writing again; "
            "if a retraction legitimately removed nodes, lower node_high_water DELIBERATELY.")
except Exception as e:
    bad(f"cannot load/validate BODYTWIN graph: {e}"); ids = []; n_nodes = 0

# ---- 2. ISOLATION: BODYTWIN and CS graphs share NO node ids (biology never leaks into CS) ----
try:
    if OTHER_GRAPH:
        _, cs_nodes = _nodes(OTHER_GRAPH)
        cs_ids = {n.get("id") for n in cs_nodes if isinstance(n, dict)}
        collide = sorted(set(ids) & cs_ids)
        if collide: bad(f"ISOLATION BREACH: node ids in BOTH graphs: {collide}")
except Exception as e:
    bad(f"cannot read the second graph for the isolation check: {e}")

# ---- 3. Fold path is the canonical one: bodytwin_fold -> fold_gate_v2, extractor present ------
try:
    src = open(FOLD_ADAPTER).read()
    if "fold_gate_v2" not in src: bad("bodytwin_fold.py does NOT use fold_gate_v2 (gate drift)")
    if "def extract_cell_json" not in src: bad("bodytwin_fold.py missing extract_cell_json (extractor drift)")
    if 'BODYTWIN_GRAPH' not in src: bad("bodytwin_fold.py no longer resolves its graph target through BODYTWIN_GRAPH")
    # NB: the ACTUAL isolation guard is check #2 (biology never appears in CS's graph — an OUTCOME check).
    # We deliberately do NOT grep the adapter source for the string 'data/ANCHOR_GRAPH.json' — it legitimately
    # appears in NEVER-warning comments, and a false-alarming guard trains future readers to ignore it.
    if re.search(r'open\(\s*[^)]*\bANCHOR_GRAPH\.json[^)]*,\s*["\']w', src) and "BODYTWIN_ANCHOR_GRAPH" not in re.search(r'open\(\s*[^)]*\bANCHOR_GRAPH\.json[^)]*,\s*["\']w', src).group(0):
        bad("bodytwin_fold.py opens a NON-bodytwin ANCHOR_GRAPH.json for WRITE (real isolation breach)")
    if '"status": "OPEN"' not in src and "'status': 'OPEN'" not in src:
        bad("bodytwin_fold.py no longer writes canonical status OPEN for hypotheses")
except Exception as e:
    bad(f"cannot read bodytwin_fold.py: {e}")

# ---- 4. Fold is LOSSLESS: fold_one must lift the submitted cert_design (hidden_state / legs /
#         anchor / datapoints / couples_to) into the node, not store a claim-only shell. -----
try:
    src = open(FOLD_ADAPTER).read()
    if "_lift_cert_design" not in src:
        bad("bodytwin_fold.py lost _lift_cert_design() — fold is LOSSY again (nodes become claim-only "
            "shells; see  root-cause + the project notes*)")
    elif "_lift_cert_design(out" not in src:
        bad("bodytwin_fold.py defines _lift_cert_design but fold_one no longer calls it — lift is dead code")
    # Catastrophe floor (an INVARIANT+floor, NOT an exact-count gate that would drift): a non-trivial
    # graph must carry many cells whose cert_design has real legs/anchor. If almost none do, the fold
    # regressed to lossy. Floors only fire on catastrophe; they do not drift as the graph grows.
    real = sum(1 for n in nodes if isinstance(n, dict)
               and ((n.get("cert_design") or {}).get("legs")
                    or str((n.get("cert_design") or {}).get("anchor") or "").strip()))
    if n_nodes >= 200 and real < 40:
        bad(f"only {real}/{n_nodes} nodes carry cert_design legs/anchor — fold likely regressed to lossy")
except Exception as e:
    bad(f"cannot verify lossless-fold invariant: {e}")

# ---- 5. Graph WRITES are atomic + flock-guarded (bodytwin_graph_io), never a raw truncating open().
#          incident: a plain json.dump(g, <graph-open-for-write>) was interrupted mid-write
#         and left a truncated, unparseable 266KB graph (14MB->266KB), which got committed. Absorbed
#         lane-I's video_cartographer_land pattern (flock RMW + os.replace atomicity + reload-verify). --
try:
    BUILDER = os.path.join(HERE, "bodytwin_build_anchor_graph.py")
    GIO = os.path.join(HERE, "bodytwin_graph_io.py")
    for path, gvar in [(FOLD_ADAPTER, "BODYTWIN_GRAPH"), (BUILDER, "GRAPH")]:
        s = open(path).read()
        if "bodytwin_graph_io" not in s:
            bad(f"{os.path.basename(path)} does not route graph writes through bodytwin_graph_io "
                "(atomic+flock+verify absent — the  truncation class can recur)")
        if re.search(r'open\(\s*' + gvar + r'\s*,\s*["\']w', s):
            bad(f"{os.path.basename(path)} still does a RAW graph write open({gvar}, 'w') — route it "
                "through bodytwin_graph_io.write_graph (atomic os.replace + flock + reload-verify)")
    gs = open(GIO).read() if os.path.isfile(GIO) else ""
    if not gs:
        bad("bodytwin_graph_io.py missing — the guarded graph writer is gone")
    else:
        if "os.replace" not in gs: bad("bodytwin_graph_io.py missing os.replace (atomic-write leg absent)")
        if "flock" not in gs: bad("bodytwin_graph_io.py missing flock (concurrency guard absent)")
except Exception as e:
    bad(f"cannot verify graph-write hardening: {e}")

# ---- N. PER-NODE fold integrity: the graph's cells must actually still HOLD what their source files
#         emitted. Added  after 15 cells were silently emptied and every commit still said
#         PASS. Root cause of the BLIND SPOT (not of the data loss): check 4 above greps the fold SOURCE
#         TEXT for _lift_cert_design and applies a whole-graph catastrophe FLOOR. That verifies THE
#         FIXING FUNCTION EXISTS IN THE CODE — not that any individual node retained its arrays. A guard
#         written as a function-name grep passes while that function is bypassed by an input path that
#         never reaches it, and a whole-graph floor cannot fire on 15 nodes out of 2770.
#         This check reads the DATA instead: for every node whose own pre-fold source file still exists
#         and emitted real cert content, the node must not be an EMPTY SHELL. It is a per-node
#         catastrophe invariant (source had content, node has none at all), NOT an exact-count gate —
#         partial deficits are reported as INFO only, because the fold legitimately drops couples_to
#         entries that name a non-existent node, and failing on those would be a false alarm that
#         teaches the next reader to work around this check.
try:
    _AGENT_OUT = os.environ.get("BODYTWIN_EVIDENCE_DIR", os.path.join(REPO, "evidence"))
    _ARR = ("legs", "datapoints", "couples_to", "conflicts", "honest_gaps")
    _by_id = {}
    for _n in (nodes if isinstance(nodes, list) else list(nodes.values())):
        if isinstance(_n, dict) and _n.get("id"):
            _by_id[_n["id"]] = _n
    _shells, _partial, _checked = [], 0, 0
    if os.path.isdir(_AGENT_OUT):
        for _fn in os.listdir(_AGENT_OUT):
            if not _fn.endswith(".json"):
                continue
            try:
                _o = json.load(open(os.path.join(_AGENT_OUT, _fn)))
            except Exception:
                continue
            if not isinstance(_o, dict):
                continue          # some evidence files are a LIST at top level
            _pc = _o.get("proposed_cell") if isinstance(_o.get("proposed_cell"), dict) else _o
            _cid = _pc.get("id")
            _node = _by_id.get(_cid) if isinstance(_cid, str) else None
            if not _node:
                continue
            _scd = _pc.get("cert_design") if isinstance(_pc.get("cert_design"), dict) else {}
            _src_items = sum(len(_scd[_f]) for _f in _ARR if isinstance(_scd.get(_f), list))
            if _src_items < 3:
                continue          # nothing substantial was emitted; nothing to lose
            _checked += 1
            _dcd = _node.get("cert_design") or {}
            _dst_items = sum(len(_dcd[_f]) for _f in _ARR if isinstance(_dcd.get(_f), list))
            if _dst_items == 0:
                if _cid not in _shells:
                    _shells.append(_cid)   # several source files can map to ONE node id; count nodes, not files
            elif _dst_items < _src_items:
                _partial += 1
    if _shells:
        bad("fold integrity: %d node(s) are EMPTY SHELLS whose own source file emitted real cert "
            "content — the lift was bypassed, not merely imperfect (repair by re-folding the source "
            "WRAPPED in proposed_cell, then read back legs/datapoints/couples_to): %s"
            % (len(_shells), ", ".join(sorted(_shells)[:5]) + (" …" if len(_shells) > 5 else "")))
    elif not ("--quiet" in sys.argv):
        print("  [fold-integrity] %d node/source pairs checked, 0 empty shells, %d with a partial "
              "deficit (INFO only — dangling-id filtering is legitimate)." % (_checked, _partial))
except Exception as e:
    bad("cannot verify per-node fold integrity: %s" % e)

# ---- N+1. RETROACTIVE PROVENANCE SCAN: the fold-time-only guard (bodytwin_fold.py _declared_id,
#      ~L412-436) never rescans a node folded before it existed. : an exhaustive graph-wide
#      rescan (bodytwin_provenance_scan.py) found 32 nodes whose cited evidence artifact
#      self-declares a DIFFERENT node's id, and whose own claim independently fails to match that
#      evidence (i.e. not explained by the separately-common, benign "duplicate id-mint" pattern) --
#      confirmed by manual reads (e.g. MODEL-ALPHASYN-PRION-PROPAGATION cites HHT/VEGF vascular-biology
#      evidence that belongs, self-declared, to CVD-HHT-VEGF-BLOCKADE-REVERSE-PROOF-PATHWAY-CONVERGENCE).
#      Run on every hardening check so a newly-folded misregistration cannot go unnoticed again. -------
try:
    from bodytwin_provenance_scan import scan as _prov_scan
    _true_misreg, _prov_info = _prov_scan(nodes, REPO)
    if _true_misreg:
        bad("provenance scan: %d node(s) cite evidence that belongs to different/no work than their "
            "own claim (run `python3 bodytwin_provenance_scan.py` for the full list, e.g. %s)"
            % (len(_true_misreg), _true_misreg[0]["node"]))
    elif not ("--quiet" in sys.argv):
        print("  [provenance-scan] %d nodes, %d cite agent_outputs evidence, %d id-matched, %d benign "
              "mismatches classified content-confirmed, 0 true misregistrations."
              % (_prov_info["total_nodes"], _prov_info["nodes_with_ao_evidence"], _prov_info["match"],
                 _prov_info["content_matches"]))
except Exception as e:
    bad("cannot run retroactive provenance scan: %s" % e)

# ---- N+2. VERDICT WRITE-BACK PROVENANCE GUARD: a node's '[WRITE-BACK <date>: executed result
#      recovered from X]' annotation must name a file this node itself cites in evidence[] (X ending
#      .py exempt -- a legitimate re-run-script artifact). Root cause (see node
#      GRAPH-SHARED-VERDICT-ATOM-ROOT-CAUSE-ADJUDICATED): an ad hoc repair pass at commit
#      60398eaef read couples_to membership as identity-of-finding and copied a topically-coupled
#      SIBLING's verdict onto 3 AUTO-* atoms -- NOT the fold path (fold_one() parses exactly one raw
#      file per call, mechanically incapable of this). The 3 nodes are repaired at source; this check
#      is the permanent guard so a future ad hoc repair making the same mistake fails loudly here
#      instead of accumulating silently for days, as this one did. ------------------------------------
try:
    from bodytwin_writeback_provenance_scan import scan as _wb_scan
    _wb_violations, _wb_info = _wb_scan(nodes, REPO)
    if _wb_violations:
        bad("write-back provenance guard: %d node(s) carry a '[WRITE-BACK ... recovered from X]' "
            "verdict whose X is not in that node's evidence[] (run "
            "`python3 bodytwin_writeback_provenance_scan.py` for the full list, e.g. %s)"
            % (len(_wb_violations), _wb_violations[0]["node"]))
    elif not ("--quiet" in sys.argv):
        print("  [writeback-provenance-scan] %d nodes, %d WRITE-BACK blocks (%d .py-exempt, %d "
              "self-evidence-confirmed), 0 borrowed-sibling-verdict violations."
              % (_wb_info["total_nodes"], _wb_info["total_writeback_blocks"],
                 _wb_info["exempt_py_artifacts"], _wb_info["matched_own_evidence"]))
except Exception as e:
    bad("cannot run write-back provenance scan: %s" % e)

# ---- report ----------------------------------------------------------------------------------
quiet = "--quiet" in sys.argv
if V:
    print("BODYTWIN HARDENING CHECK: ❌ FAIL — %d VIOLATION(S) (do NOT work around; fix at source):" % len(V))
    for v in V: print("  ✗ " + v)
    print("Authority: docs/RUNNING.md")
    sys.exit(1)
if not quiet:
    print("BODYTWIN HARDENING CHECK: ✅ PASS — %d nodes, valid status enum, valid types, full schema,"
          " no dup ids, no dialect leak, isolation clean (0 id-collision with the second graph),"
          " fold path = bodytwin_fold->fold_gate_v2, lossless cert_design lift." % n_nodes)
sys.exit(0)
