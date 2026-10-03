"""Apply shared research-value advice without mutating or dropping queue rows."""
from pathlib import Path
import re,sys,signal
signal.signal(signal.SIGPIPE,signal.SIG_DFL)
sys.path.insert(0,'~/research/AGENT_DASHBOARD_20260930')
from research_value import order
lines=Path(sys.argv[1]).read_text().splitlines(keepends=True)
for line in order(lines)[0]:
    sys.stdout.write(line)
