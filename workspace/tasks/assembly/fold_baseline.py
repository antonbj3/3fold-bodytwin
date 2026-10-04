#!/usr/bin/env python3
"""Pre-rewrite baseline for by-name evidence pointers. Run BEFORE filter-repo; nothing reconstructs it after.

A fold id resolves by NAME, so a rewrite that changes a fold row's content while keeping its id leaves every
citing pointer intact and now meaning something else. No staleness check and no reproducibility re-run sees
that: the pointer resolves and the script reruns against the new content. The only guard is a content hash
taken before the rewrite.

Two design points that matter:

  `commit` is EXCLUDED from the hash. A history rewrite changes commit SHAs legitimately, so hashing a row
  with its commit in it reports every such row as changed for a correct reason -- and a gate that fires on
  the healthy case teaches you to ignore it. The commit value is recorded SEPARATELY so it can be remapped
  through filter-repo's commit-map and checked on its own.

  The hash is over a CANONICALISED row (sorted keys, no whitespace), so a reformat is not a content change
  while any value change is.

Usage:  fold_baseline.py capture <repo> <out.json>
        fold_baseline.py verify  <repo> <out.json>   [after the rewrite]
"""
import hashlib
import json
import pathlib
import sys

HASH_EXCLUDE = ("commit",)          # legitimately changed by a rewrite


def canon(row: dict) -> str:
    keep = {k: v for k, v in row.items() if k not in HASH_EXCLUDE}
    return json.dumps(keep, sort_keys=True, separators=(",", ":"), ensure_ascii=False)


def fold_rows(repo: pathlib.Path):
    for p in sorted(repo.rglob("FOLD_LEDGER*.jsonl")):
        for ln, line in enumerate(p.read_text(errors="ignore").splitlines(), 1):
            if not line.strip():
                continue
            try:
                row = json.loads(line)
            except ValueError:
                yield str(p.relative_to(repo)), ln, None, line, None
                continue
            fid = next((str(row[k]) for k in ("id", "fold_id", "name") if k in row), None)
            yield str(p.relative_to(repo)), ln, fid, line, row


def graphs(repo: pathlib.Path):
    for p in sorted(repo.rglob("*ANCHOR_GRAPH*.json")):
        try:
            d = json.loads(p.read_text())
        except ValueError:
            continue
        n = d["nodes"] if isinstance(d, dict) and "nodes" in d else d
        yield p, [x for x in (n.values() if isinstance(n, dict) else n) if isinstance(x, dict)]


def capture(repo: pathlib.Path) -> dict:
    rows, dup, unparsable = {}, [], []
    for path, ln, fid, line, row in fold_rows(repo):
        if row is None:
            unparsable.append({"file": path, "line": ln})
            continue
        h = hashlib.sha256(canon(row).encode()).hexdigest()
        key = fid or f"{path}:{ln}"
        if key in rows:
            dup.append(key)
        rows[key] = {"file": path, "line": ln, "sha256": h,
                     "commit": row.get("commit"), "has_id": fid is not None}
    cited, byname_only = {}, []
    for gp, nodes in graphs(repo):
        for nd in nodes:
            ev = [e for e in (nd.get("evidence") or []) if isinstance(e, str)]
            if not ev:
                continue
            ids = [e for e in ev if "/" not in e]
            for e in ids:
                cited.setdefault(e, []).append(nd["id"])
            if ids and len(ids) == len(ev):
                byname_only.append({"node": nd["id"], "graph": str(gp.relative_to(repo)), "fold_ids": ids})
    missing = sorted(e for e in cited if e not in rows)
    return {"rows": rows, "cited_fold_ids": {k: sorted(v) for k, v in cited.items()},
            "evidence_is_by_name_only": byname_only, "cited_but_not_in_any_ledger": missing,
            "duplicate_keys": sorted(set(dup)), "unparsable_lines": unparsable,
            "n_rows": len(rows), "n_cited": len(cited), "n_byname_only_nodes": len(byname_only),
            "hash_excludes": list(HASH_EXCLUDE),
            "note": "capture BEFORE the rewrite; `commit` is excluded from the hash and recorded separately"}


def verify(repo: pathlib.Path, base: dict) -> dict:
    now = capture(repo)
    changed, vanished, commit_moved = [], [], []
    for key, old in base["rows"].items():
        new = now["rows"].get(key)
        if new is None:
            vanished.append(key)
            continue
        if new["sha256"] != old["sha256"]:
            changed.append({"fold": key, "file": new["file"], "was": old["sha256"][:12], "now": new["sha256"][:12]})
        if old.get("commit") != new.get("commit"):
            commit_moved.append({"fold": key, "was": old.get("commit"), "now": new.get("commit")})
    cited = set(base["cited_fold_ids"])
    critical = {n["node"]: n["fold_ids"] for n in base["evidence_is_by_name_only"]}
    hit = sorted({nd for c in changed + [{"fold": v} for v in vanished]
                  for nd, ids in critical.items() if c["fold"] in ids})
    return {"ok": not changed and not vanished,
            "content_changed_while_id_kept": changed,
            "cited_fold_vanished": [v for v in vanished if v in cited],
            "vanished_total": len(vanished),
            "commit_field_remapped": len(commit_moved),
            "byname_only_nodes_affected": hit,
            "verdict": ("CLEAN" if not changed and not vanished else
                        "POINTERS NOW MEAN SOMETHING ELSE -- see content_changed_while_id_kept")}


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print(__doc__)
        sys.exit(2)
    cmd, repo, out = sys.argv[1], pathlib.Path(sys.argv[2]), pathlib.Path(sys.argv[3])
    if cmd == "capture":
        d = capture(repo)
        out.write_text(json.dumps(d, indent=1, ensure_ascii=False))
        print(json.dumps({k: d[k] for k in ("n_rows", "n_cited", "n_byname_only_nodes",
                                            "cited_but_not_in_any_ledger", "duplicate_keys")}, indent=1)[:900])
    elif cmd == "verify":
        print(json.dumps(verify(repo, json.loads(out.read_text())), indent=1, ensure_ascii=False)[:1200])
    else:
        print(__doc__)
        sys.exit(2)
