from common_r3 import *
import shlex
base = read(ROOT / 'FROZEN_GENERATOR_A.json')['argv']
commands = []
for row in read(DATA / 'A/RECORDS.json'):
    if 'peak_rss_MiB' not in row:
        continue
    ident = row['key'] + '__' + row['participant']
    cmd = ['/usr/bin/time', '-v', '-o', str(ROOT / 'raw' / ('MEMORY_' + ident + '.txt'))] + base[:-1] + ['/runner/resume_worker.py', row['key'], row['participant']]
    commands.append(dict(stage='resume_A', argv=cmd))
for p in sorted((ROOT / 'raw').glob('*_CMD_*.json')):
    commands.append(dict(source_file=p.name, **read(p)))
dump(ROOT / 'raw/EXACT_COMMANDS.json', commands)
(ROOT / 'raw/EXACT_COMMANDS.sh').write_text('#!/bin/bash\n# Recorded commands, for inspection; do not overwrite frozen outputs.\n' + '\n'.join((shlex.join(c['argv']) for c in commands)) + '\n')
