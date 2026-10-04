"""Count BodyTwin work eligible for the cloud driver, including active jobs."""
from pathlib import Path
from collections import Counter
import re

ROOT = Path(__file__).resolve().parents[1]


def supply(root=ROOT):
    lanes = root / 'tasks/lanes'
    remote = Counter(re.findall(r'ovhstart (\S+) ', (lanes / 'ovh_agents/queue_ovh.log').read_text(errors='replace')))
    local = Counter(re.findall(r'\] start (\S+) \((\S+) (\S+)\)', (lanes / 'bt_queue.log').read_text(errors='replace')))
    eligible = set()
    for line in (lanes / 'bt_queue.txt').read_text().splitlines():
        fields = line.split()
        if len(fields) != 3 or fields[0].startswith('#'):
            continue
        profile, model, job = fields
        if job.startswith('BT-DW48-'):
            continue
        directory = root / 'results' / job
        report = directory / 'RESULTS.md'
        if not directory.is_dir() or (report.is_file() and report.stat().st_size):
            continue
        # A running second attempt still supplies work; an exhausted idle job does not.
        if (directory / '.ovh_claim').exists():
            eligible.add(job)
        elif remote[job] < 2 and local[(job, profile, model)] < 3:
            eligible.add(job)
    return len(eligible)


if __name__ == '__main__':
    print(supply())
