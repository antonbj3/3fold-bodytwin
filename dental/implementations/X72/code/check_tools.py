from common import *
import subprocess, time

def main():
    paths = ['vendor/1216032/tests/test_contact_port_v1.py', 'vendor/6e2defd/tests/test_contact_port_adapter.py', 'vendor/e590aa3/tests/test_farkas_witness.py', 'vendor/7abf765/tests/test_normal_cone.py', 'vendor/95ec803/tests/test_family_sweep.py']
    cmd = [sys.executable, '-m', 'pytest', '-q', '-rs', *paths]
    env = dict(os.environ, PYTEST_DISABLE_PLUGIN_AUTOLOAD='1', PYTHONPATH=str(ROOT / 'vendor'))
    t = time.perf_counter()
    p = subprocess.run(cmd, cwd=ROOT, env=env, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    (ROOT / 'raw/UPSTREAM_TESTS.txt').write_text(p.stdout)
    write(ROOT / 'raw/UPSTREAM_TESTS.json', dict(command=cmd, env_overrides=dict(PYTEST_DISABLE_PLUGIN_AUTOLOAD='1', PYTHONPATH='vendor'), exit_code=p.returncode, output=p.stdout, wall_s=time.perf_counter() - t))
    assert p.returncode == 0, p.stdout
    print(p.stdout.splitlines()[-1])
if __name__ == '__main__':
    main()
