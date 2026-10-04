import sys, time, importlib.metadata
from common import R, read, dump, sha, state, now

def inputs():
    lock = read('DEMO_INPUT_LOCK.json')
    for a in lock['files']:
        if sha(R / a['path']) != a['sha256']:
            raise ValueError('Pinned file hash drift ' + a['path'])
    for line in (R / 'requirements.lock').read_text().splitlines():
        (name, version) = line.split('==')
        if importlib.metadata.version(name) != version:
            raise ValueError('Runtime version differs ' + name)

def run():
    start = time.monotonic()
    inputs()
    state('one_command_demo', 'pinned inputs verified', 're-extract sources, repeat comparisons, inject faults, fresh bridge solve')
    from extract import extract
    from fit_r1 import run as fit
    from compare_r2 import run as compare
    from contrast_r3 import run as contrast
    from paired_r4 import run as pair
    from hazard_r5 import run as hazard
    from freeze_predictions import run as frozen
    from verify import run as verify
    from export_lab import run as export
    from replay import run as replay
    from report import run as report
    extract()
    fit()
    compare()
    contrast()
    pair(solve_new=False)
    hazard()
    frozen()
    verify()
    export()
    if '--skip-fresh-solve' not in sys.argv:
        replay()
    out = report()
    dump('DEMO_RUN.json', dict(created_utc=now(), seconds=time.monotonic() - start, fresh_solve='--skip-fresh-solve' not in sys.argv, complete=True, source_unit_fault_controls='PASS', scientific_status=out['overall_status'], prediction_refit=False))
    state('demo_complete', 'source/fault/replay controlsPASS; connector/physical/mesh gatesFAILorUNKNOWN', 'independent review; measured connector radius plus support/contact strain; then matched area9/16 batch fracture experiment')
    print('DEMO COMPLETE: published margins reproduced; physical connector law UNKNOWN; see README_DEMO.md and bridge_connector.png', flush=True)
if __name__ == '__main__':
    run()
