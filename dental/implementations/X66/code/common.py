from pathlib import Path
import json, hashlib, datetime, ast, types
ROOT = Path(__file__).resolve().parents[1]
DENTAL = ROOT.parents[1]

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()

def put(p, d):
    p = ROOT / p
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(d, ensure_ascii=False, indent=2, allow_nan=False, default=lambda x: x.item() if hasattr(x, 'item') else str(x)) + '\n')

def read(p):
    return json.loads((ROOT / p).read_text())

def old_functions(path, names, ns):
    s = Path(path).read_text()
    tree = ast.parse(s)
    body = [n for n in tree.body if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef)) and n.name in names]
    assert {n.name for n in body} == set(names)
    mod = types.SimpleNamespace(**ns)
    exec(compile(ast.Module(body=body, type_ignores=[]), str(path), 'exec'), mod.__dict__)
    return mod
