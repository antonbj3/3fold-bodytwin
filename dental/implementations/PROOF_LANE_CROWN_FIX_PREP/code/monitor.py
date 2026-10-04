"""Run one worker, measure its process tree, abort own worker before3.8GiB."""
import os, sys, time, subprocess, json, pathlib, signal
import psutil
R = pathlib.Path(__file__).resolve().parents[1]
tag = sys.argv[1]
cmd = sys.argv[2:]
t = time.monotonic()
p = subprocess.Popen(cmd, start_new_session=True)
root = psutil.Process(p.pid)
peak = 0
maxthreads = 0
count = 0
aborted = False
while p.poll() is None:
    try:
        procs = [root] + root.children(recursive=True)
        rss = sum((q.memory_info().rss for q in procs if q.is_running()))
        peak = max(peak, rss)
        count += 1
        if rss > 3800 * 1024 ** 2:
            os.killpg(p.pid, signal.SIGTERM)
            aborted = True
            break
    except psutil.NoSuchProcess:
        pass
    time.sleep(0.2)
rc = p.wait()
(R / 'raw' / f'{tag}_RESOURCE.json').write_text(json.dumps(dict(command=cmd, returncode=rc, seconds=time.monotonic() - t, peak_sampled_tree_RSS_MiB=peak / 1024 ** 2, sampling_interval_s=0.2, samples=count, aborted=aborted, hard_worker_RLIMIT_AS_MiB=3500), indent=2) + '\n')
sys.exit(rc)
