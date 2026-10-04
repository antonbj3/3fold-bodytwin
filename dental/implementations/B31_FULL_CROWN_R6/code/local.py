from dental_release.paths import expand as _release_expand
import os, sys
from pathlib import Path
LANE_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(LANE_ROOT.parent / 'PROOF_LANE_FULL_CROWN_R5/code'))
from construct_b import *
R6 = LANE_ROOT
D6 = Path(_release_expand('@DENTAL_WORK_ROOT@/PROOF_LANE_FULL_CROWN_R6'))
D6.mkdir(exist_ok=True)
threadpool_limits(limits=1)

def state6(phase, gate, nxt):
    save(R6 / 'CURRENT_WORK_STATE.json', dict(lane='PROOF_LANE-full-crown-r6', phase=phase, latest_gate=gate, next_operation=nxt, updated_utc=now(), review_state='PENDING_INDEPENDENT_REVIEW'))

def first():
    return next(((a, b) for (a, b) in inputs() if a['key'] == '079905ebf9504544_molar'))
