#!/usr/bin/env python3
"""Launch one public-data cloud lane from the already trusted research checkout."""
import datetime as dt
import fcntl
import json
import os
from pathlib import Path
import pty
import re
import select
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
REPO = Path('external_research_path')
COORDINATOR = 'local_config_path/bin/coordinator'
PROMPT = (
    'Read BRIEF.md and PREREG.md. Carry out this bounded task using only public data '
    'or the bundled code. Do not use private repositories, subagents, paid extra usage, '
    'or Git push. Write RESULTS.md, results.json, and small reproducible code. '
    'Keep the computation bounded to 20 minutes; then stop remaining sweeps and report '
    'completed cases plus PARTIAL/UNKNOWN for unfinished cases. '
    'Commit the result locally on your session branch. In your FINAL message print each '
    'result file verbatim, in UTF-8, as =====FILE <path>===== followed by its contents '
    'and =====END FILE=====. Use relative paths; include RESULTS.md and results.json; '
    'total output at most 200 KB. State source URLs/DOIs, uncertainty, counterexamples, '
    'and whether each preregistered criterion passed. Do not claim validation from model output alone.'
)


def git(*args):
    return subprocess.run(['git', *args], cwd=REPO, check=True, capture_output=True, text=True).stdout.strip()


def main():
    if len(sys.argv) != 3:
        raise SystemExit('usage: launch_cloud_lane.py CLOUD-<ID> <bundle_dir>')
    lane_id, bundle_arg = sys.argv[1:]
    bundle = Path(bundle_arg).resolve()
    if not re.fullmatch(r'CLOUD-[A-Z0-9-]+', lane_id) or bundle.name != lane_id:
        raise SystemExit('lane id must match CLOUD-<ID> bundle directory')
    if not bundle.is_dir() or not all((bundle / n).is_file() for n in ('BRIEF.md', 'PREREG.md')):
        raise SystemExit('bundle must contain BRIEF.md and PREREG.md')
    files = [p for p in bundle.rglob('*') if '.git' not in p.relative_to(bundle).parts]
    if any(p.is_symlink() for p in files):
        raise SystemExit('bundle may not contain symlinks')
    size = sum(p.stat().st_size for p in files if p.is_file())
    if size > 20_000_000:
        raise SystemExit('bundle exceeds 20 MB')
    with (ROOT / 'launcher.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        intent = ROOT / (lane_id + '.intent.json')
        if intent.exists():
            raise SystemExit('existing launch intent; inspect before retrying')
        if git('status', '--porcelain'):
            raise SystemExit('trusted checkout has local changes')
        git('checkout', 'main')
        baseline = git('rev-parse', 'HEAD')
        git('checkout', '-B', 'lane/' + lane_id, baseline)
        for item in REPO.iterdir():
            if item.name == '.git':
                continue
            if item.is_dir() and not item.is_symlink():
                shutil.rmtree(item)
            else:
                item.unlink()
        for item in bundle.iterdir():
            if item.name == '.git':
                continue
            target = REPO / item.name
            if item.is_dir():
                shutil.copytree(item, target)
            else:
                shutil.copy2(item, target)
        git('add', '-A')
        subprocess.run(['git', '-c', 'user.name=BodyTwin Cloud',
                        '-c', 'user.email=bodytwin-cloud@localhost', 'commit', '-qm',
                        'Preregister ' + lane_id], cwd=REPO, check=True)
        commit = git('rev-parse', 'HEAD')
        record = {'lane': lane_id, 'time': dt.datetime.now(dt.timezone.utc).isoformat(),
                  'requested_model': 'coordinator', 'status': 'launch_intent',
                  'bundle_bytes': size, 'branch': 'lane/' + lane_id,
                  'baseline_commit': baseline, 'bundle_commit': commit}
        intent.write_text(json.dumps(record) + '\n')
        env = os.environ.copy()
        for key in ('ANTHROPIC_API_KEY', 'ANTHROPIC_AUTH_TOKEN'):
            env.pop(key, None)
        env.update(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
        cmd = [COORDINATOR, '--model', 'coordinator', '--effort', 'high', '--cloud', PROMPT]
        child, fd = pty.fork()
        if child == 0:
            os.chdir(REPO)
            os.execve(COORDINATOR, cmd, env)
        transcript = ''
        with (ROOT / (lane_id + '.launch.log')).open('w') as log:
            while True:
                readable, _, _ = select.select([fd], [], [], 60)
                if not readable:
                    continue
                try:
                    chunk = os.read(fd, 16384)
                except OSError:
                    break
                if not chunk:
                    break
                s = chunk.decode(errors='replace')
                log.write(s)
                log.flush()
                transcript = (transcript + s)[-65536:]
                match = re.findall(r'https://coordinator\.ai/code/(session_[A-Za-z0-9]+)', transcript)
                if match and 'session_id' not in record:
                    record['session_id'] = match[-1]
                    record['session_url'] = 'https://coordinator.ai/code/' + match[-1]
                    intent.write_text(json.dumps(record) + '\n')
                print(s, end='', flush=True)
        _, status = os.waitpid(child, 0)
        os.close(fd)
        record['launcher_exit_code'] = os.waitstatus_to_exitcode(status)
        record['status'] = 'submitted' if record.get('session_id') and record['launcher_exit_code'] == 0 else 'uncertain'
        intent.write_text(json.dumps(record) + '\n')
        if record.get('session_id'):
            with (ROOT / 'receipts.jsonl').open('a') as out:
                fcntl.flock(out, fcntl.LOCK_EX)
                out.write(json.dumps(record) + '\n')
                out.flush()
                os.fsync(out.fileno())
        git('checkout', 'main')
        raise SystemExit(record['launcher_exit_code'] or (0 if record.get('session_id') else 1))


if __name__ == '__main__':
    main()
