"""Apply shared research-value advice without mutating or dropping queue rows."""
from pathlib import Path
import re,sys,signal
signal.signal(signal.SIGPIPE,signal.SIG_DFL)
import os
# The covardation is located in the shared dashboard outside the tradet; the sokvagen star in
# ~/.bodytwin/runners.env so repot didn't carry a local sokvag.
sys.path.insert(0, os.environ.get('RESEARCH_VALUE_DIR',
                                  os.path.expanduser('~/research/AGENT_DASHBOARD_20260930')))
from research_value import order
lines=Path(sys.argv[1]).read_text().splitlines(keepends=True)
for line in order(lines)[0]:
    sys.stdout.write(line)
