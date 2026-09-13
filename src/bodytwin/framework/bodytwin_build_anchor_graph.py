#!/usr/bin/env python3
"""Build a dependency graph from a cell-design catalog.

Input: a JSON catalog of designed cells ($BODYTWIN_DESIGNS, default
designs/cell_designs.json). Output: the graph JSON ($BODYTWIN_GRAPH) written through
the guarded writer in bodytwin_graph_io, plus a human-readable cell inventory markdown
file ($BODYTWIN_INVENTORY). All nodes are written with status OPEN (designed, not
measured).
"""
import json
import os
import sys
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bodytwin_graph_io import write_graph  # guarded writer (atomic+flock+verify); never json.dump(open(...,"w"))

REPO = Path(os.environ.get("BODYTWIN_ROOT", os.getcwd()))
DESIGNS = Path(os.environ.get("BODYTWIN_DESIGNS", REPO / "designs" / "cell_designs.json"))
GRAPH = Path(os.environ.get("BODYTWIN_GRAPH", REPO / "ANCHOR_GRAPH.json"))
INV = Path(os.environ.get("BODYTWIN_INVENTORY", REPO / "CELL_INVENTORY.md"))

def legs_str(legs):
    out = []
    for leg in legs:
        if isinstance(leg, list) and len(leg) == 2:
            name, ids = leg
            out.append(f"{name} [{', '.join(ids) if ids else 'stub'}]")
        else:
            out.append(str(leg))
    return " · ".join(out)

def main():
    d = json.load(open(DESIGNS))
    nodes = []
    def add(cluster, c):
        n = c["node"]; cert = c.get("cert", {})
        nodes.append({
            "id": n["id"], "claim": n["claim"], "type": n["type"], "status": n["status"],
            "evidence": n.get("evidence", []), "depends_on": n.get("depends_on", []),
            "load": n.get("load", 0), "risk": n.get("risk", "MED"),
            "regime_note": n.get("regime_note", ""), "actionable_by": "bodytwin-reviewer",
            "cluster": cluster,
            "cert_design": {
                "hidden_state": cert.get("hidden", ""),
                "occluded_truth": cert.get("occluded", False),
                "legs": cert.get("legs", []),
                "anchor": cert.get("anchor_src") or cert.get("anchor", ""),
                "composition_mode": cert.get("mode", ""),
                "couples_to_inherited": cert.get("couples_to_inherited", "none"),
                "verify": cert.get("_verify"),
            },
        })
    for c in d.get("foundational", []):
        add("foundational", c)
    for cl, blk in d.get("clusters", {}).items():
        for c in blk.get("cells", []):
            add(cl, c)

    # load = number of nodes that depend on this one (matches anchor_graph_tools' computed load)
    dependents = {n["id"]: 0 for n in nodes}
    for n in nodes:
        for dep in n["depends_on"]:
            if dep in dependents:
                dependents[dep] += 1
    for n in nodes:
        n["load"] = dependents[n["id"]]

    graph = {
        "_meta": {
            "name": "BODYTWIN_ANCHOR_GRAPH — bio build-goal cells (example designs, all HYPOTHESIS/OPEN)",
            "note": "Synthesized from the cell-design catalog. status OPEN = designed, NOT measured. Work a node by measuring its legs in-regime and folding the result against the held-out anchor.",
            "build_order": "roots (no deps): NODE-0-ENERGY, KNEE-CELL, ORG-WAVEFORM-ANCHOR + the standalone DATA-ANCHORs; then the cells that depend on them. Cheapest first-close: MOL-SAMESAMPLE-PATCHSEQ-CALIBRATION (LOW risk, open DANDI, same-sample).",
        },
        "nodes": nodes,
    }
    # BOOTSTRAP RE-RUN GUARD. This script SYNTHESIZES a fresh graph from the
    # seed catalog and hands write_graph a whole new object -- it is a one-time bootstrap, not
    # an append. Run against a populated graph it silently replaces every folded cell with the 51
    # seed cells. That happened: 3323 nodes -> 51, and the loss reached a commit before anyone
    # looked, because the only signal was a node COUNT and nothing compares it to the file's
    # prior value. write_graph's node_high_water is monotonic and so cannot be lowered by a
    # truncating write -- but nothing REFUSED the write, and the hardening check reads the graph
    # it is handed. The append-only invariant has to be enforced at the writer that breaks it.
    if GRAPH.exists():
        try:
            existing = json.load(open(GRAPH))
            prior = len(existing.get("nodes", []) if isinstance(existing, dict) else existing)
        except Exception:
            prior = None
        if prior is not None and prior > len(nodes):
            raise SystemExit(
                f"REFUSING TO BOOTSTRAP OVER A POPULATED GRAPH: {GRAPH} already holds {prior} "
                f"nodes and this bootstrap would write {len(nodes)}. The graph is append-only. "
                f"If you genuinely intend to re-seed, move the existing file aside first -- this "
                f"script must never be the thing that decides to discard {prior - len(nodes)} cells."
            )
    write_graph(graph, GRAPH)   # atomic + flock-guarded + reload-verify (bodytwin_graph_io)

    # inventory markdown
    L = []
    def verd(n):
        v = (n["cert_design"].get("verify") or {}).get("verdict", "un-audited")
        return v
    def bucket(v):
        v = v.upper()
        if v.startswith("REAL"): return "REAL"
        if "FENCED" in v: return "FENCED"
        if "HONEST-NEG" in v: return "HONEST-NEG"
        if "RE-SCOPED" in v: return "RE-SCOPED"
        if "REFUTED" in v: return "REFUTED"
        if "WEAKENED" in v: return "WEAKENED"
        return "un-audited"
    from collections import Counter
    tally = Counter(bucket(verd(n)) for n in nodes)
    real = [n["id"] for n in nodes if bucket(verd(n)) == "REAL"]

    L.append("# BODYTWIN CELL INVENTORY — the bio build-goal map\n")
    L.append(f"Synthesized from the cell-design catalog into the dependency graph. "
             f"**{len(nodes)} cells, adversarially audited (none MEASURED yet — verdicts are on the DESIGN).**\n")
    L.append("## Certifiability (design-audit verdicts)\n")
    L.append("| bucket | count | meaning |")
    L.append("|---|---|---|")
    L.append(f"| REAL | {tally['REAL']} | design survives adversarial verify + external fact-check → ready to EXECUTE |")
    L.append(f"| WEAKENED | {tally['WEAKENED']} | fixable defect, residual noted per cell |")
    L.append(f"| RE-SCOPED | {tally['RE-SCOPED']} | honest weaker claim (archive lacks the ideal instrument) |")
    L.append(f"| REFUTED | {tally['REFUTED']} | anchor/leg defect not yet resolved |")
    L.append(f"| FENCED / HONEST-NEG | {tally['FENCED']+tally['HONEST-NEG']} | no valid anchor in archive (=PASS, non-promotable) |")
    L.append(f"| un-audited | {tally['un-audited']} | (honest-neg no-anchor cells / node-0) |")
    L.append(f"\n**REAL (execute these first):** {', '.join(real)}\n")
    L.append("Each cell = an occluded hidden state cornered by ≥2 decorrelated legs vs a held-out external anchor. "
             "`occluded✓` = needs a consequence leg.\n")
    order = ["foundational", "musculoskeletal", "energy_metabolism", "neural_sensory", "organ_systemic", "molecular_program"]
    by_cluster = {}
    for n in nodes:
        by_cluster.setdefault(n["cluster"], []).append(n)
    gaps = {cl: blk.get("honest_gaps", []) for cl, blk in d.get("clusters", {}).items()}
    for cl in order:
        cells = by_cluster.get(cl, [])
        if not cells:
            continue
        L.append(f"\n## {cl}  ({len(cells)} cells)\n")
        for n in cells:
            cd = n["cert_design"]
            occ = "occluded✓" if cd["occluded_truth"] else "direct"
            vb = bucket(verd(n))
            L.append(f"### {n['id']}  · **[{vb}]** · {n['type']} · risk={n['risk']} · {occ}\n")
            vfy = cd.get("verify") or {}
            if vfy.get("note"):
                L.append(f"- **verdict:** {vfy.get('verdict','')} — {vfy['note']}")
            if vfy.get("residual") and vfy["residual"] != "none":
                L.append(f"- **residual:** {vfy['residual']}")
            L.append(f"- **claim:** {n['claim']}")
            L.append(f"- **hidden state:** {cd['hidden_state']}")
            L.append(f"- **legs:** {legs_str(cd['legs'])}")
            L.append(f"- **anchor (held-out):** {cd['anchor']}  ·  **mode:** {cd['composition_mode']}")
            if n["depends_on"]:
                L.append(f"- **depends_on:** {', '.join(n['depends_on'])}")
            if cd["couples_to_inherited"] and cd["couples_to_inherited"] != "none":
                L.append(f"- **couples into CS:** {cd['couples_to_inherited']}")
            if n["regime_note"]:
                L.append(f"- **regime:** {n['regime_note']}")
            L.append("")
        if gaps.get(cl):
            L.append(f"**honest gaps ({cl}):**")
            for g in gaps[cl]:
                L.append(f"- {g}")
            L.append("")
    INV.write_text("\n".join(L))
    print(f"wrote {GRAPH}  ({len(nodes)} nodes)")
    print(f"wrote {INV}")
    # quick coupling summary
    real = sorted({t for n in nodes for t in [n['cert_design']['couples_to_inherited']]
                   if t and t != 'none'})
    print(f"cells with cross-graph coupling annotation: {sum(1 for n in nodes if n['cert_design']['couples_to_inherited'] not in ('','none'))}")

if __name__ == "__main__":
    main()
