"""OpenCode on OVH with an explicit account (A/B/C). Copy of SPACE_SWARM_MULTIKEY_20260923/run_profile.py,
adapted: all profiles are isolated under /opt/agents/profiles/<X> (auth.json only, chmod 600)."""
import argparse
import json
import os
from pathlib import Path

ROOT = Path('/opt/agents')
MODELS = {'reserve_worker': 'opencode-go/reserve_worker-v4.1-flash', 'swarm': 'opencode/space-swarm-free'}


def environment(label, model, web=False):
    profile = ROOT / 'profiles' / label
    account_paths = {folder: profile / folder for folder in ('config', 'data', 'state', 'cache')}
    auth = account_paths['data'] / 'opencode/auth.json'
    provider, model_id = model.split('/', 1)
    if provider not in json.loads(auth.read_text()):
        raise SystemExit('Requested provider is not configured in this profile.')
    env = os.environ.copy()
    for key in list(env):
        if key.endswith('_API_KEY') or key.startswith('OPENCODE_'):
            env.pop(key)
    for variable, folder in [('XDG_CONFIG_HOME', 'config'), ('XDG_DATA_HOME', 'data'), ('XDG_STATE_HOME', 'state'), ('XDG_CACHE_HOME', 'cache')]:
        path = account_paths[folder]
        path.mkdir(parents=True, exist_ok=True, mode=0o700)
        env[variable] = str(path)
    config = {'model': model, 'small_model': model, 'default_agent': 'build',
              'enabled_providers': [provider], 'provider': {provider: {'whitelist': [model_id]}},
              'autoupdate': False, 'share': 'disabled',
              'permission': {'task': 'deny', 'external_directory': 'deny', 'webfetch': 'allow' if web else 'deny', 'websearch': 'allow' if web else 'deny'}}
    env.update(OPENCODE_DISABLE_PROJECT_CONFIG='true', OPENCODE_CONFIG_CONTENT=json.dumps(config),
               OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', NUMEXPR_NUM_THREADS='1',
               PYTHONPATH='/opt/bt/vendor')
    return env


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('profile', choices=('A', 'B', 'C'))
    parser.add_argument('model', choices=tuple(MODELS))
    parser.add_argument('--dir', required=True)
    parser.add_argument('--title', default='isolated_profile_job')
    parser.add_argument('prompt')
    args = parser.parse_args()
    model = MODELS[args.model]
    # Web only for packets with no internal data (marker file ALLOW_WEB in the job directory)
    env = environment(args.profile, model, web=os.path.exists(os.path.join(args.dir, 'ALLOW_WEB')))
    os.chdir(args.dir)
    executable = '/opt/agents/.opencode/bin/opencode'
    os.execve(executable, [executable, 'run', '--pure', '--auto', '--agent', 'build',
                           '--model', model, '--dir', str(Path.cwd()), '--format', 'json',
                           '--title', args.title, args.prompt], env)
