#!/usr/bin/env python3
"""The one guarded writer for the graph JSON.

write_graph(graph, path) writes atomically (temp file + os.replace), under an exclusive
flock, and re-reads the file to verify node count and wrapper shape before returning.
graph_lock(path) exposes the same lock for read-modify-write callers. A raw json.dump into
an open()-for-write truncates the graph if the process is interrupted; this module exists
so that never happens.
"""
import os
import json
import fcntl
import tempfile
from contextlib import contextmanager

REPO = os.environ.get("BODYTWIN_ROOT", os.getcwd())
GRAPH_PATH = os.environ.get("BODYTWIN_GRAPH", os.path.join(REPO, "ANCHOR_GRAPH.json"))
# In-repo lock (shared by cl/cl2 in the same worktree), gitignored, NEUTRAL name — no biology term
# and no bodytwin identity in a filesystem path (isolation §1.3 neutral-dialect).
LOCK_PATH = os.environ.get("BODYTWIN_GRAPH_LOCK", os.path.join(REPO, ".graph.lock"))

@contextmanager
def graph_lock():
    """Exclusive flock held across a whole read-modify-write (open lock -> LOCK_EX -> yield -> unlock).

    ⚠ NOT REENTRANT. Each entry open()s a NEW fd; flock(LOCK_EX) on a second fd of the same file
    blocks forever even within one process. So `with graph_lock(): ... write_graph(G)` SELF-DEADLOCKS,
    because write_graph() already takes this lock. Measured : two runs hung silently to the
    timeout with zero output (stdout buffering ate the progress prints, so it looked like a slow parse
    of the 21 MB graph — it was a deadlock; the same edit ran in 0.4 s once the outer lock was removed).
    Read-modify-write the graph as:  G = json.load(open(GRAPH_PATH)); ...mutate...; write_graph(G)
    Take graph_lock() yourself ONLY around a critical section that does NOT call write_graph().
    """
    with open(LOCK_PATH, "a+") as lf:
        fcntl.flock(lf, fcntl.LOCK_EX)
        try:
            yield lf
        finally:
            fcntl.flock(lf, fcntl.LOCK_UN)

def _count(g):
    nodes = g if isinstance(g, list) else (g.get("nodes") or g.get("cells") or [])
    return len(nodes)

def atomic_write_json(obj, path, indent=2):
    """Serialize obj to `path` ATOMICALLY: write a sibling tmp, fsync, os.replace (atomic rename).
    Does NOT take the lock — the caller holds graph_lock() when concurrency matters. ensure_ascii is
    left at the json default (True) to stay byte-compatible with the existing graph serialization."""
    path = os.fspath(path)
    d = os.path.dirname(path) or "."
    fd, tmp = tempfile.mkstemp(dir=d, prefix=".graph.", suffix=".tmp")
    try:
        with os.fdopen(fd, "w") as f:
            json.dump(obj, f, indent=indent)
            f.flush()
            os.fsync(f.fileno())
        os.replace(tmp, path)          # atomic on POSIX (same filesystem)
    except BaseException:
        try:
            os.unlink(tmp)
        except OSError:
            pass
        raise

def reload_verify(path, expected_nodes=None):
    """Re-read the just-written file: it MUST parse, and (if given) carry expected_nodes nodes.
    Raises (JSONDecodeError / IOError) on a bad write so corruption is caught, never trusted."""
    with open(os.fspath(path)) as f:
        back = json.load(f)
    got = _count(back)
    if expected_nodes is not None and got != expected_nodes:
        raise IOError(f"graph write verify FAILED: expected {expected_nodes} nodes, read back {got}")
    return got

def bump_high_water(obj):
    """Record the graph's node-count high-water mark in _meta, monotonically.

    ★ WHY. Mutation-tested : the hardening check PRINTS the node count in its
    PASS line but never CHECKS it, so deleting 1976 of 2976 nodes passed the gate
    silently — as did deleting all but 76. Ten schema/isolation mutations all bind; pure
    node loss was invisible at every scale. The count was reported, not verified.

    The invariant that fixes it is the one the project established the same day for census
    claims: the graph is APPEND-ONLY, so the count can never legitimately fall. Measured
    over the last 30 graph commits — it never shrank once. Storing the high-water mark
    here, in the one guarded writer every path is already forced through, means any
    write that loses nodes leaves evidence the hardening check can fail on.

    Monotonic on purpose: a caller that hands us a truncated graph cannot lower the bar
    by writing it. Never lowers, only raises.
    """
    if not isinstance(obj, dict):
        return None                                # bare-list graph: schema check owns it
    meta = obj.setdefault("_meta", {})
    prev = meta.get("node_high_water")
    n = _count(obj)
    meta["node_high_water"] = max(n, prev) if isinstance(prev, int) else n
    return meta["node_high_water"]

def write_graph(obj, path=GRAPH_PATH, indent=2, _locked=False):
    """Single guarded write: flock -> atomic tmp+fsync+os.replace -> reload+node-count verify.
    Returns the verified node count. Use this everywhere instead of json.dump(open(path,"w")).

    ⚠ A LOCKED WRITE IS NOT ENOUGH FOR READ-MODIFY-WRITE. Locking only the write serializes the
    writes but not the DECISIONS: two processes can each json.load() the same pre-change snapshot,
    each append their own node to their OWN in-memory copy, and each call this — the second write
    replaces the whole file with its stale object, silently erasing the first writer's node.
    reload_verify() cannot catch it, because it compares the file against the WRITER'S OWN count,
    which matches. Measured  in bodytwin_fold.fold_one, which read the graph unlocked at
    line 355 and wrote here ~45 lines later; the concurrency precondition was observed live (the
    node count changed under a reader that took no write action).

    Callers doing read-modify-write MUST hold graph_lock() across the whole cycle and pass
    _locked=True here — graph_lock() is NOT reentrant (see its docstring), so holding it and calling
    this WITHOUT _locked=True self-deadlocks.
    """
    bump_high_water(obj)                          # append-only invariant, see its docstring
    if _locked:                                   # caller already holds graph_lock()
        atomic_write_json(obj, path, indent=indent)
        return reload_verify(path, _count(obj))
    with graph_lock():
        atomic_write_json(obj, path, indent=indent)
        return reload_verify(path, _count(obj))

if __name__ == "__main__":
    # self-test: round-trip the live graph through the guarded writer on a scratch COPY (never the real
    # file), proving atomicity + verify without touching production. Prints the verified node count.
    import shutil
    src = GRAPH_PATH
    tmpcopy = os.path.join(tempfile.gettempdir(), "graph_io_selftest.json")
    shutil.copy(src, tmpcopy)
    g = json.load(open(tmpcopy))
    n = write_graph(g, tmpcopy)
    print(f"graph_io self-test: wrote+verified {n} nodes to a scratch copy (atomic+flock+reload OK)")
    os.remove(tmpcopy)
