"""Immutable review binding, so later demo timing updates do not invalidate it."""
from common import *
import shutil

def run():
    from integrity import verify
    verification = verify()
    result = read(ROOT / 'results.json')
    if not result['headline_results']['all_faults_rejected']:
        raise ValueError('Do not seal failed implementation controls')
    if result['headline_results']['blind_images'] != 32:
        raise ValueError('Unexpected rater package coverage')
    dest = ROOT / 'review_artifacts'
    if dest.exists():
        raise ValueError('Immutable review artifacts already exist')
    dest.mkdir()
    files = ['results.json', 'LEADERBOARD.md', 'raw/R1.json', 'raw/ITEMS.json', 'raw/R1B_BOUNDARY.json', 'raw/R3.json', 'raw/FUNCTIONAL_ROWS.json', 'raw/PHYSICAL_ROWS.json.gz', 'raw/R5_CONTRAST.json', 'raw/R6_EXTERNAL_SUPPORT.json', 'raw/CONTROLS.json', 'raw/CONTACT_NORMS.json', 'raw/X60_TABLE.json', 'raw/LAB_PACKAGE.json', 'raw/RUN_COST.json', 'raw/EXTERNAL_REPLAY.json', 'raw/FINAL_INPUT_VERIFICATION.json']
    artifacts = []
    for name in files:
        target = dest / Path(name).name
        shutil.copy2(ROOT / name, target)
        artifacts.append(dict(path=str(target.relative_to(ROOT)), sha256=sha(target), bytes=target.stat().st_size))
    code = {str(p.relative_to(ROOT)): sha(p) for p in sorted((ROOT / 'gencad_bench_v5').glob('*.py'))}
    for name in ['run_all.sh', 'replay_external.sh']:
        code[name] = sha(ROOT / name)
    freeze(ROOT / 'REVIEW_SNAPSHOT.json', dict(claim_type='capability', review_state='PENDING_INDEPENDENT_REVIEW', purpose='Stable retrospective result binding; does not assert graph experiment readiness', result=result, artifacts=artifacts, code_sha256=code, input_verification=verification, external_model_lock_sha256=sha(ROOT / 'raw/EXTERNAL_MODEL_LOCK.json'), frozen_image_package_sha256=sha(ROOT / 'FROZEN_IMAGE_PACKAGE.json'), input_lock_sha256=sha(ROOT / 'INPUT_LOCK.json')))
    return dict(status='SEALED_PENDING_REVIEW', path=str(ROOT / 'REVIEW_SNAPSHOT.json'), sha256=sha(ROOT / 'REVIEW_SNAPSHOT.json'), artifacts=len(artifacts))
if __name__ == '__main__':
    print(run())
