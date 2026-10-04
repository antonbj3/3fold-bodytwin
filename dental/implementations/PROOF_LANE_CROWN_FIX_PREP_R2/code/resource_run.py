import sys, time, json, subprocess, psutil, pathlib, datetime
label = sys.argv[1]
cmd = sys.argv[2:]
start = time.monotonic()
p = subprocess.Popen(cmd)
peak = 0
abort = False
while p.poll() is None:
    try:
        ps = [psutil.Process(p.pid)] + psutil.Process(p.pid).children(recursive=True)
        rss = sum((q.memory_info().rss for q in ps if q.is_running()))
        peak = max(peak, rss)
        if rss > 3500 * 1024 ** 2:
            abort = True
            for q in reversed(ps):
                q.kill()
    except (psutil.NoSuchProcess, psutil.AccessDenied):
        pass
    time.sleep(0.05)
rc = p.wait()
dest = pathlib.Path(__file__).resolve().parents[1] / 'raw' / (label + '_RESOURCES.json')
dest.write_text(json.dumps(dict(command=cmd, exit_code=rc, seconds=time.monotonic() - start, peak_tree_RSS_MiB=peak / 1024 ** 2, limit_MiB=3500, killed_for_memory=abort, timestamp=datetime.datetime.now(datetime.timezone.utc).isoformat()), indent=2) + '\n')
sys.exit(rc)
