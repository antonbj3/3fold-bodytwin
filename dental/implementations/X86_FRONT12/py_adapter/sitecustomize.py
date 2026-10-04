"""Adapt archived bind destinations at launch; installed binaries/modules stay untouched."""
from dental_release.paths import expand as _release_expand
import os, subprocess, sys, shutil
from pathlib import Path
_original = subprocess.Popen.__init__
_batch4 = Path(_release_expand('@DENTAL_WORK_ROOT@/DEMO48_PACKAGE_BATCH4')).resolve()
_template = _batch4 / 'GENCAD_V3'
_work = _batch4 / 'work'
_original_copy2 = shutil.copy2

def copy_with_budget(src, dst, *a, **kw):
    source = Path(src).resolve()
    destination = Path(dst).resolve()
    if source.is_relative_to(_template) and destination.is_relative_to(_work):
        required = int(os.environ.get('DENTAL_X86_REGULAR_COPY_BYTES', '0'))
        available = int(os.environ.get('DENTAL_X86_REMAINING_INTERMEDIATE_BYTES', '3000000000'))
        if required >= available:
            raise RuntimeError(f'RESOURCE_LIMIT: original regular-file copy requires {required} bytes; remaining own3GB budget {available} bytes; no numerical comparison reached')
    return _original_copy2(src, dst, *a, **kw)
shutil.copy2 = copy_with_budget

def adapt(args):
    if not isinstance(args, (list, tuple)) or not args or Path(str(args[0])).name != 'bwrap':
        return args
    argv = list(args)
    hidden = []
    changes = []
    i = 1
    inherited_root = any((argv[j] in ['--bind', '--ro-bind'] and argv[j + 1:j + 3] == ['/', '/'] for j in range(len(argv) - 2)))
    while i < len(argv):
        if argv[i] == '--tmpfs' and i + 1 < len(argv):
            hidden.append(argv[i + 1].rstrip('/'))
            i += 2
        elif argv[i] in ['--bind', '--ro-bind', '--bind-try', '--ro-bind-try', '--dev-bind', '--dev-bind-try'] and i + 2 < len(argv):
            dst = argv[i + 2]
            if inherited_root and dst.startswith('/') and (not any((dst == h or dst.startswith(h + '/') for h in hidden))):
                resolved = str(Path(dst).resolve())
                if resolved != dst:
                    argv[i + 2] = resolved
                    changes.append((dst, resolved))
            i += 3
        else:
            i += 1
    if changes:
        print('Archived bind destinations resolved:', changes, file=sys.stderr, flush=True)
    return argv

def init(self, args, *a, **kw):
    return _original(self, adapt(args), *a, **kw)
subprocess.Popen.__init__ = init
