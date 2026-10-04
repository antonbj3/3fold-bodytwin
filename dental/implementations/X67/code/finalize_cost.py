from common import *
import re
p = ROOT / 'raw/PIPELINE_TIME.txt'
t = p.read_text()

def number(label):
    return float(re.search(re.escape(label) + ':\\s*([0-9.]+)', t).group(1))
elapsed = next((s.rsplit(': ', 1)[-1] for s in t.splitlines() if 'Elapsed (wall clock)' in s))
parts = list(map(float, elapsed.split(':')))
wall = sum((v * 60 ** i for (i, v) in enumerate(parts[::-1])))
r = load('results.json')
r['full_cost']['metered_one_command_pipeline'] = dict(wall_seconds=wall, user_CPU_seconds=number('User time (seconds)'), system_CPU_seconds=number('System time (seconds)'), max_RSS_KiB=int(number('Maximum resident set size (kbytes)')), timing_file='raw/PIPELINE_TIME.txt', timing_sha256=sha(p), mode='source response basis reused iff exact source/code/prereg/measurement contract matches; all bisection, fits, certificates, plot and verification executed', tail_cost='This cost-binding write and final binding verification are outside the measured pipeline; unmetered bookkeeping, no comparison gain claimed')
r['full_cost']['lane_bytes'] = sum((p.stat().st_size for p in ROOT.rglob('*') if p.is_file()))
dump('results.json', r)
f = load('GRAPH_FEEDBACK.json')
f['sha256'] = sha(ROOT / 'results.json')
dump('GRAPH_FEEDBACK.json', f)
print('Metered wall', wall, 's; max RSS', r['full_cost']['metered_one_command_pipeline']['max_RSS_KiB'], 'KiB')
