"""Clean-copy reproduction from this lane's small code/inputs only."""
from pathlib import Path
import tempfile, shutil, subprocess, json, time
from common import ROOT, save

def main():
    t0 = time.perf_counter()
    with tempfile.TemporaryDirectory(prefix='X51-clean-', dir=ROOT / 'raw') as td:
        p = Path(td)
        for folder in ['code', 'inputs']:
            dest = p / folder
            dest.mkdir()
            for f in (ROOT / folder).iterdir():
                if f.is_file():
                    shutil.copy2(f, dest / f.name)
        for pattern in ['run_all.sh', 'INPUT_LOCK.json', 'PREREG*.json', 'PREREG*.sha256', 'FROZEN_PREDICTIONS*.json', 'FROZEN_PREDICTIONS*.sha256']:
            for f in ROOT.glob(pattern):
                shutil.copy2(f, p / f.name)
        for folder in ['raw', 'history', 'figures']:
            (p / folder).mkdir()
        run = subprocess.run([str(p / 'run_all.sh')], cwd=p, capture_output=True, text=True, timeout=60)
        (ROOT / 'raw/CLEAN_REPLAY.stdout').write_text(run.stdout)
        (ROOT / 'raw/CLEAN_REPLAY.stderr').write_text(run.stderr)
        out = {'exit_code': run.returncode, 'wall_s': time.perf_counter() - t0, 'source_only': 'own lane code, frozen small inputs and contracts', 'network_fetch': False}
        if run.returncode == 0:
            j = json.loads((p / 'results.json').read_text())
            out['controls'] = j['controls']
            out['source_load_replay'] = next((c['detail'] for c in j['R1']['checks'] if c['name'] == 'X25_CENSORED_LOAD_REPLAY'))
            out['figure_generated'] = (p / 'figures/LAB_ACCEPTANCE_DRIFT.png').exists()
            out['physical_unknown_preserved'] = j['outcomes']['empirical_alarm_day'] == 'UNKNOWN'
        save('raw/CLEAN_REPLAY.json', out)
        print(json.dumps(out, indent=2))
        if run.returncode:
            raise SystemExit(run.returncode)
if __name__ == '__main__':
    main()
