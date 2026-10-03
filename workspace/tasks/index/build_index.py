#!/usr/bin/env python3
"Indexes all results and datasets so no data point is lost.\nWriting tasks/index/index.json (each file: path, size, sha256, mtime; PREREG before/after run)\nand tasks/index/DATASETS.json. Runs at every wake-up. Does not modify result files."
import hashlib, json, os, time, glob
W = "."
DS = "/mnt/shared_data/datasets/bodytwin"
def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""): h.update(b)
    return h.hexdigest()
out = {"built": time.strftime("%Y-%m-%dT%H:%M:%S%z"), "experiments": {}}
for d in sorted(glob.glob(f"{W}/results/*/")):
    eid = os.path.basename(d.rstrip("/"))
    files = []
    for root, dirs, fs in os.walk(d):
        dirs[:] = [x for x in dirs if x != "__pycache__"]
        for f in fs:
            p = os.path.join(root, f); st = os.stat(p)
            files.append({"path": os.path.relpath(p, W), "bytes": st.st_size,
                          "mtime": time.strftime("%H:%M:%S", time.localtime(st.st_mtime)),
                          "sha256": sha(p) if st.st_size < 200_000_000 else None})
    pre = os.path.join(d, "PREREG.md")
    res = [x for x in files if x["path"].endswith(".json")]
    prereg_first = None
    if os.path.exists(pre) and res:
        pm = os.stat(pre).st_mtime
        prereg_first = all(os.stat(os.path.join(W, x["path"])).st_mtime >= pm for x in res)
    out["experiments"][eid] = {"n_files": len(files), "has_readme": os.path.exists(os.path.join(d, "README.md")),
                               "has_prereg": os.path.exists(pre), "prereg_before_all_json": prereg_first, "files": files}
json.dump(out, open(f"{W}/tasks/index/index.json", "w"), indent=1)
ds = {}
for m in sorted(glob.glob(f"{DS}/*/MANIFEST.json")):
    try: ds[os.path.basename(os.path.dirname(m))] = json.load(open(m))
    except Exception as e: ds[m] = {"error": str(e)}
json.dump(ds, open(f"{W}/tasks/index/DATASETS.json", "w"), indent=1)
for k, v in out["experiments"].items():
    print(f"{k:4s} filer={v['n_files']:4d} README={v['has_readme']!s:5s} PREREG={v['has_prereg']!s:5s} prereg_ffirst={v['prereg_before_all_json']}")
print("dataset:", ", ".join(ds))
