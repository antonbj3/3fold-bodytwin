"""Find checks that cannot fail: comparisons where both sides are the same expression.

Usage: python3 find_vacuous.py [max_hits] [root_directory]

The graph lane (anton-4d, 2/10) measured that this is the most common and hardest form to see, because the line
LOOKS LIKE a comparison: abs(x - x) < tol is always true. Four instances were in BodyTwin.

Searches results/*/**.py for:
  1. SELF_DIFF   - a - a, a == a, a / a within the same expression (identical operands)
  2. SELF_ABS    - abs(<expression> - <same expression>) < tol
Reports file, line and the expression. Makes no change; this is a search tool.
Uses the AST, not regular expressions, so comments and strings give no false hits.
"""
import ast, sys, json
from pathlib import Path

# Root: argv[2] if given, otherwise BodyTwin's results. The graph lane (2/10) needed it movable in order to run
# the tool against the graph engine and copied the file instead; now no copy is needed.
R = Path(sys.argv[2]) if len(sys.argv) > 2 else Path('results')


def same(a, b):
    try:
        return ast.dump(a) == ast.dump(b)
    except Exception:
        return False


# Words that make a self-comparison DANGEROUS: the field promises reproduction, invariance or acceptance.
# Without them, x == x is usually the NaN idiom and unit/unit a dimension check - both legitimate.
PROMISE = ('reproduc', 'matche', 'matched', 'invarian', 'agree', 'consistent', 'verif', 'check',
           'pass', 'holds', 'identical', 'unchanged', 'preserv', 'equal', 'confirm', 'validat')


def _promise_near(tree, lineno):
    """True if an assignment or dict key on the same line names a promise."""
    for n in ast.walk(tree):
        if getattr(n, 'lineno', None) != lineno and getattr(n, 'lineno', None) != lineno - 1:
            continue
        if isinstance(n, ast.Constant) and isinstance(n.value, str) and any(w in n.value.lower() for w in PROMISE):
            return n.value[:60]
        if isinstance(n, ast.Name) and any(w in n.id.lower() for w in PROMISE):
            return n.id
        if isinstance(n, ast.Attribute) and any(w in n.attr.lower() for w in PROMISE):
            return n.attr
    return None


def scan(path):
    try:
        tree = ast.parse(path.read_text(errors='ignore'))
    except (SyntaxError, ValueError, OSError):
        return []
    out = []
    for node in ast.walk(tree):
        hit = None
        if isinstance(node, ast.BinOp) and isinstance(node.op, (ast.Sub, ast.Div)) and same(node.left, node.right):
            hit = 'SELF_' + type(node.op).__name__.upper()
        elif isinstance(node, ast.Compare) and len(node.comparators) == 1 and same(node.left, node.comparators[0]):
            hit = 'SELF_COMPARE'
        if hit:
            promise = _promise_near(tree, node.lineno)
            if promise:
                out.append((node.lineno, hit, f'[{promise}] ' + ast.unparse(node)[:140]))
    return out


def main():
    limit = int(sys.argv[1]) if len(sys.argv) > 1 else 0
    rows = []
    for p in sorted(R.glob('*/**/*.py')):
        for line, kind, src in scan(p):
            rows.append({'file': str(p.relative_to(R)), 'line': line, 'kind': kind, 'src': src})
            if limit and len(rows) >= limit:
                break
        if limit and len(rows) >= limit:
            break
    print(json.dumps({'found': len(rows), 'rows': rows}, ensure_ascii=False, indent=1))


if __name__ == '__main__':
    main()
