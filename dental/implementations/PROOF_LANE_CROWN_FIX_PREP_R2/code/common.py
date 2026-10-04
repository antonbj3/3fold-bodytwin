from dental_release.paths import expand as _release_expand
import os, sys
from pathlib import Path
SOURCE_ROOT = Path(__file__).resolve().parents[1]
ROOT = Path(os.environ.get('DENT_R2_RUN_ROOT', SOURCE_ROOT))
PREV = Path(os.environ.get('DENT_R2_PARENT_ROOT', _release_expand('@DENTAL_IMPLEMENTATIONS@/PROOF_LANE_CROWN_FIX_PREP')))
os.environ['DENT_FIX_RUN_ROOT'] = str(ROOT)
os.environ['DENT_FIX_DATA_ROOT'] = os.environ.get('DENT_R2_DATA_ROOT', _release_expand('@DENTAL_WORK_ROOT@/PROOF_LANE_CROWN_FIX_PREP_R2'))
sys.path.insert(0, str(PREV / 'code'))
from core import *
OLD_D = PARENT_D / 'PROOF_LANE_CROWN_FIX_PREP'
for n in ['exact_clearance', 'exact_surface']:
    if not (D / n).exists():
        (D / n).symlink_to(OLD_D / n)

def state(phase, latest, nxt):
    dump(R / 'CURRENT_WORK_STATE.json', dict(lane='PROOF_LANE-crown-fix-prep-r2', updated_utc=now(), phase=phase, latest_gate=latest, next_operation=nxt, review_state='PENDING_INDEPENDENT_REVIEW', data_root=D))

def inputs():
    return read(PREV / 'INPUTS.json')['records']

def native(rec):
    info = read(PREV / 'raw' / (rec['key'] + '_R2_SUPPORTS.json'))['native']
    return (*meshread(info['path']), info)
