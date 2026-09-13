#!/usr/bin/env python3
"""Cell-id alias handling at the extraction boundary.

ID_KEYS lists the keys a contributor envelope may use for the cell id; declared_id(payload)
returns the id under any of them; normalize_id_alias(payload) rewrites the payload so every
downstream reader sees the same key. Fixing only one reader leaves the others rejecting a
valid cell, so the lookup lives here and is shared.

--selftest runs the built-in assertions (requires bodytwin_fold and mt_extract on the path).
"""

# `id` first — an artifact that declares BOTH is taken at its `id`.
ID_KEYS = ("id", "cell_id")

def declared_id(obj):
    """The id the artifact says it is (wrapped or flat), or None. Never inferred, never minted."""
    if not isinstance(obj, dict):
        return None
    pc = obj.get("proposed_cell")
    if isinstance(pc, dict):
        for k in ID_KEYS:
            v = pc.get(k)
            if isinstance(v, str) and v.strip():
                return v.strip()
    for k in ID_KEYS:
        v = obj.get(k)
        if isinstance(v, str) and v.strip():
            return v.strip()
    return None

def normalize_id_alias(obj):
    """Mirror an aliased id onto `id` IN PLACE, at the level where it was found, and return obj.

    Applied at every EXTRACTION BOUNDARY so no consumer downstream has to know the alias exists.
    Idempotent; a no-op when `id` is already present or when no id is declared at all. The alias
    key is left in place — it is provenance, and removing it would make a re-run of this diagnosis
    impossible.
    """
    if not isinstance(obj, dict):
        return obj
    for level in (obj.get("proposed_cell"), obj):
        if not isinstance(level, dict):
            continue
        cur = level.get("id")
        if isinstance(cur, str) and cur.strip():
            continue
        for k in ID_KEYS[1:]:
            v = level.get(k)
            if isinstance(v, str) and v.strip():
                level["id"] = v.strip()
                break
    return obj

def _selftest():
    """Reproduces the SPECIFIC failing payload shapes, end to end through the real extractors."""
    import json, os, sys, tempfile
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    ok = True

    def chk(name, cond):
        nonlocal ok
        print(("  PASS  " if cond else "  FAIL  ") + name)
        ok = ok and bool(cond)

    CID = "NODE-0-ENERGY-HEAT-ROUTES-CLOSURE-CLOTHED-EXT"
    # The measured shape: wrapper present, id declared as cell_id, NO top-level `claim` key,
    # a sibling `result_paths` key, as one earlier submission returned it.
    payload = {"proposed_cell": {"cell_id": CID, "status": "OPEN", "parent_cells": ["NODE-0-ENERGY"],
                                 "summary": "clothed-extension heat-route closure, 15-31 C sweep",
                                 "cert_design": {"legs": ["a", "b"], "anchor": "x",
                                                 "datapoints": [1], "honest_gaps": ["g"]}},
               "result_paths": ["/tmp/x.json"]}
    flat = {"cell_id": CID + "-FLAT", "claim": "c" * 900,
            "cert_design": {"legs": ["a"], "anchor": "x", "datapoints": [1], "honest_gaps": []}}

    chk("declared_id reads the wrapped alias", declared_id(payload) == CID)
    chk("declared_id reads the flat alias", declared_id(flat) == CID + "-FLAT")
    chk("declared_id: no id declared -> None", declared_id({"proposed_cell": {"status": "OPEN"}}) is None)
    chk("declared_id: `id` wins over `cell_id`",
        declared_id({"proposed_cell": {"id": "A", "cell_id": "B"}}) == "A")
    chk("normalize is idempotent",
        normalize_id_alias(normalize_id_alias({"proposed_cell": {"cell_id": CID}}))["proposed_cell"]["id"] == CID)
    chk("normalize never invents an id",
        "id" not in normalize_id_alias({"proposed_cell": {"status": "OPEN"}})["proposed_cell"])
    chk("normalize keeps the alias as provenance",
        normalize_id_alias({"cell_id": "Z"}).get("cell_id") == "Z")

    d = tempfile.mkdtemp(prefix="mt_id_alias_")
    tpath = os.path.join(d, "t.txt")
    with open(tpath, "w") as f:
        f.write("prose before\nMT-JSON-BEGIN " + json.dumps(payload) + " MT-JSON-END\nprose after\n")

    import mt_extract
    got = mt_extract.extract_cell(tpath)
    chk("mt_extract.extract_cell finds the aliased cell", isinstance(got, dict))
    chk("mt_extract yields a usable `id` (was None -> the drop)",
        isinstance(got, dict) and got.get("proposed_cell", {}).get("id") == CID)

    # A FLAT aliased cell has no "id" and may have no wrapper: _is_cell must still accept it.
    fpath = os.path.join(d, "f.txt")
    with open(fpath, "w") as f:
        f.write("MT-JSON-BEGIN " + json.dumps(flat) + " MT-JSON-END\n")
    gf = mt_extract.extract_cell(fpath)
    chk("mt_extract accepts a FLAT cell whose id is `cell_id`",
        isinstance(gf, dict) and gf.get("proposed_cell", {}).get("id") == CID + "-FLAT")

    import bodytwin_fold as mf
    chk("bodytwin_fold._declared_id agrees", mf._declared_id(payload) == CID)
    chk("bodytwin_fold.extract_cell_json normalizes at the boundary",
        (mf.extract_cell_json("MT-JSON-BEGIN " + json.dumps(payload) + " MT-JSON-END")
         or {}).get("proposed_cell", {}).get("id") == CID)

    print("SELFTEST", "PASS" if ok else "FAIL")
    return 0 if ok else 1

if __name__ == "__main__":
    raise SystemExit(_selftest())
