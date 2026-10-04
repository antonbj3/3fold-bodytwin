"""Load a source kernel with isolated legacy module names."""
from contextlib import contextmanager
import importlib.util
from pathlib import Path
import sys
import ast

ROOT = Path(__file__).resolve().parents[1]

def functions(demo, relative, names, namespace):
    """Execute specified unmodified function definitions without source acquisition."""
    path=ROOT/'implementations'/demo/relative
    tree=ast.parse(path.read_text())
    definitions=[node for node in tree.body if isinstance(node,(ast.FunctionDef,ast.AsyncFunctionDef)) and node.name in names]
    if {node.name for node in definitions}!=set(names):raise ValueError('Missing specified source function')
    scope=dict(namespace)
    exec(compile(ast.Module(body=definitions,type_ignores=[]),str(path),'exec'),scope)
    return scope

@contextmanager
def kernel(demo, relative, dependencies=()):
    path = ROOT/'implementations'/demo/relative
    original = list(sys.path)
    saved = dict(sys.modules)
    sys.path.insert(0, str(path.parent))
    # Source operators originally used short sibling names; do not cross-contaminate demos.
    for sibling in path.parent.glob('*.py'):
        sys.modules.pop(sibling.stem, None)
    try:
        spec = importlib.util.spec_from_file_location('release_'+demo+'_'+path.stem, path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        yield module
    finally:
        sys.path[:] = original
        for key in list(sys.modules):
            if key not in saved:
                loaded = getattr(sys.modules[key], '__file__', '') or ''
                if str(ROOT/'implementations') in loaded:
                    del sys.modules[key]
        for key, value in saved.items():
            if key not in sys.modules or getattr(sys.modules[key], '__file__', '') != getattr(value, '__file__', ''):
                sys.modules[key] = value
