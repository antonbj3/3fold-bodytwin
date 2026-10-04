"""Experimental exact conflict API. Mathematical certificates, caller-bound physics.

LP rows: Ax <= b, free variables. Cone rows: Ax=b, x in one SOC3 or orthant.
A claim is a group of rows; minimum is claim cardinality, conditional on domain.
No prose/sign-to-equation inference. No floating point may certify a decision.
"""
from __future__ import annotations
from fractions import Fraction as Q
from decimal import Decimal
from itertools import combinations
import hashlib
import json
import time
from ._farkas_check import check_dual, check_primal, rational, dot
__all__ = ['farkas_witness', 'verify_certificate', 'adapt_claim_federation', 'adapt_bodytwin', 'adapt_g3_pair', 'load_exact_json']

class ContractError(ValueError):
    pass

class BudgetError(ValueError):
    pass

def exact(v):
    if isinstance(v, Decimal):
        v = str(v)
    q = rational(v)
    if max(q.numerator.bit_length(), q.denominator.bit_length()) > 256:
        raise ContractError('EXACT_INPUT_BIT_BUDGET')
    return q

def load_exact_json(text):
    """Preserve decimal JSON tokens, rather than blessing rounded Python floats."""
    return json.loads(text, parse_float=Decimal, parse_constant=lambda x: (_ for _ in ()).throw(ContractError('NONFINITE_JSON')))

def _json(v):
    if isinstance(v, (Q, Decimal)):
        return str(v)
    if isinstance(v, dict):
        return {str(k): _json(x) for (k, x) in v.items()}
    if isinstance(v, (list, tuple)):
        return [_json(x) for x in v]
    if isinstance(v, float):
        raise ContractError('FLOAT_INPUT_NOT_EXACT')
    if v is None or isinstance(v, (str, int, bool)):
        return v
    raise ContractError('UNSUPPORTED_JSON_TYPE')

def _hash(v):
    return hashlib.sha256(json.dumps(_json(v), sort_keys=True, separators=(',', ':'), ensure_ascii=False).encode()).hexdigest()

def _unit(u):
    if not isinstance(u, str) or not u:
        raise ContractError('MISSING_UNIT')
    return u
_UNITS = {'1': ('1', Q(1)), '%': ('1', Q(1, 100)), 'M': ('M', Q(1)), 'mM': ('M', Q(1, 1000)), 'uM': ('M', Q(1, 10 ** 6)), 'µM': ('M', Q(1, 10 ** 6)), 'nM': ('M', Q(1, 10 ** 9)), 'V': ('V', Q(1)), 'mV': ('V', Q(1, 1000)), 'Pa': ('Pa', Q(1)), 'kPa': ('Pa', Q(1000)), 'N': ('N', Q(1))}

def _scale(src, dst):
    if src == dst:
        return Q(1)
    (a, b) = (_UNITS.get(src), _UNITS.get(dst))
    if a is None or b is None or a[0] != b[0]:
        raise ContractError('INCOMPATIBLE_OR_UNSUPPORTED_UNIT_CONVERSION')
    return a[1] / b[1]

def _source(s):
    if not isinstance(s, dict) or not (s.get('id') or (s.get('path') and s.get('sha256'))):
        raise ContractError('MISSING_SOURCE')

def _parse(model):
    if not isinstance(model, dict) or model.get('schema') != 'farkas-claims/v1':
        raise ContractError('MISSING_EXPLICIT_CONSTRAINT_MODEL')
    _json(model)
    (vs, cs) = (model.get('variables'), model.get('claims'))
    if not isinstance(vs, list) or not vs or (not isinstance(cs, list)) or (not cs):
        raise ContractError('EMPTY_OR_MISSING_VARIABLES_OR_CLAIMS')
    ids = [v.get('id') for v in vs]
    if any((not isinstance(i, str) or not i for i in ids)) or len(set(ids)) != len(ids):
        raise ContractError('INVALID_VARIABLE_IDS')
    if len(vs) > 16:
        raise ContractError('VARIABLE_BUDGET')
    for v in vs:
        if not v.get('quantity'):
            raise ContractError('MISSING_QUANTITY')
        _unit(v.get('unit'))
    scope = model.get('scope')
    if not isinstance(scope, dict) or not scope:
        raise ContractError('MISSING_COMMON_SCOPE')
    kind = model.get('domain', 'inequality')
    if kind not in ('inequality', 'orthant', 'soc3') or (kind == 'soc3' and len(vs) != 3):
        raise ContractError('UNSUPPORTED_CONE_DOMAIN')
    if kind == 'soc3' and len({v['unit'] for v in vs}) != 1:
        raise ContractError('SOC_COORDINATES_MUST_SHARE_UNIT_AFTER_EXPLICIT_SCALING')
    (names, groups) = (set(), [])
    for c in cs:
        cid = c.get('id')
        if not isinstance(cid, str) or not cid or cid in names:
            raise ContractError('INVALID_OR_DUPLICATE_CLAIM_ID')
        names.add(cid)
        _source(c.get('source'))
        if c.get('scope') != scope:
            raise ContractError('INCOMPATIBLE_SCOPE')
        if c.get('survival_status') == 'REJECTED_CONVENTION':
            raise ContractError('REJECTED_CLAIM_IS_NOT_A_HARD_CONSTRAINT')
        rows = []
        if 'interval' in c:
            if kind != 'inequality' or c.get('variable') not in ids:
                raise ContractError('INTERVAL_REQUIRES_FREE_LINEAR_VARIABLE')
            j = ids.index(c['variable'])
            s = _scale(_unit(c.get('unit')), vs[j]['unit'])
            (lo, hi) = c['interval']
            (lo, hi) = (exact(lo) * s, exact(hi) * s)
            if lo > hi:
                raise ContractError('INVALID_INTERVAL')
            a = [Q(0)] * len(vs)
            a[j] = Q(1)
            rows = [(tuple(a), hi, 'upper'), (tuple((-q for q in a)), -lo, 'lower')]
        else:
            if not c.get('constraints'):
                raise ContractError('MISSING_EXPLICIT_CONSTRAINTS')
            for (k, r) in enumerate(c['constraints']):
                coeff = r.get('coefficients', {})
                unit = _unit(r.get('unit'))
                if not isinstance(coeff, dict) or set(coeff) - set(ids):
                    raise ContractError('UNKNOWN_VARIABLE_IN_ROW')
                a = tuple((exact(coeff.get(i, 0)) for i in ids))
                cu = r.get('coefficient_units', {})
                for (j, i) in enumerate(ids):
                    if a[j]:
                        expected = '1' if unit == vs[j]['unit'] else unit + '/' + vs[j]['unit']
                        if cu.get(i, '1') != expected:
                            raise ContractError('COEFFICIENT_UNIT_CONTRACT_MISMATCH')
                b = exact(r.get('rhs'))
                rel = r.get('relation')
                if kind != 'inequality' and rel != '=':
                    raise ContractError('CONE_SUPPORTS_EQUALITY_ROWS_ONLY')
                if rel not in ('<=', '>=', '='):
                    raise ContractError('UNSUPPORTED_RELATION')
                if kind != 'inequality' or rel in ('<=', '='):
                    rows.append((a, b, str(k) + ':upper'))
                if kind == 'inequality' and rel in ('>=', '='):
                    rows.append((tuple((-q for q in a)), -b, str(k) + ':lower'))
        groups.append({'id': cid, 'source': c['source'], 'rows': rows})
    return {'kind': kind, 'variables': vs, 'groups': groups, 'scope': scope, 'input_sha256': _hash(model)}

def _rows(p, subset):
    (A, b, labels) = ([], [], [])
    for i in subset:
        for (a, v, label) in p['groups'][i]['rows']:
            A.append(a)
            b.append(v)
            labels.append({'claim_id': p['groups'][i]['id'], 'row': label})
    return (A, b, labels)

def _simplify(rows, limit):
    out = {}
    for (a, b, y) in rows:
        if not any(a):
            if b < 0:
                return [(a, b, y)]
            continue
        s = abs(next((q for q in a if q)))
        key = tuple((q / s for q in a))
        if key not in out or b / s < out[key][1]:
            out[key] = (key, b / s, tuple((q / s for q in y)))
    if len(out) > limit:
        raise BudgetError('FOURIER_MOTZKIN_ROW_BUDGET')
    if any((max(q.numerator.bit_length(), q.denominator.bit_length()) > 4096 for (a, b, y) in out.values() for q in (*a, b, *y))):
        raise BudgetError('EXACT_ARITHMETIC_BIT_BUDGET')
    return list(out.values())

def _eliminate(rows, limit):
    pos = [r for r in rows if r[0][-1] > 0]
    neg = [r for r in rows if r[0][-1] < 0]
    if len(pos) * len(neg) + len(rows) > limit * 4:
        raise BudgetError('FOURIER_MOTZKIN_PAIR_BUDGET')
    out = [(a[:-1], b, y) for (a, b, y) in rows if a[-1] == 0]
    for (a, b, y) in pos:
        for (c, d, z) in neg:
            (s, t) = (a[-1], -c[-1])
            out.append((tuple((x / s + w / t for (x, w) in zip(a[:-1], c[:-1]))), b / s + d / t, tuple((x / s + w / t for (x, w) in zip(y, z)))))
    return _simplify(out, limit)

def _fm(A, b, n, limit=4096):
    m = len(A)
    rows = _simplify([(tuple(a), v, tuple((Q(i == j) for j in range(m)))) for (i, (a, v)) in enumerate(zip(A, b))], limit)
    history = []
    for dim in range(n, 0, -1):
        bad = next((r for r in rows if not any(r[0]) and r[1] < 0), None)
        if bad:
            return ('INCONSISTENT', bad[2])
        history.append(rows)
        rows = _eliminate(rows, limit)
    bad = next((r for r in rows if r[1] < 0), None)
    if bad:
        return ('INCONSISTENT', bad[2])
    x = []
    for rows in reversed(history):
        (lo, hi) = ([], [])
        for (a, b, _) in rows:
            if a[-1]:
                v = (b - dot(a[:-1], x)) / a[-1]
                (lo if a[-1] < 0 else hi).append(v)
        (l, h) = (max(lo) if lo else None, min(hi) if hi else None)
        if l is not None and h is not None and (l > h):
            raise ContractError('INTERNAL_BACKSUBSTITUTION_FAILURE')
        x.append(l if l is not None else min(Q(0), h) if h is not None else Q(0))
    return ('CONSISTENT', tuple(x))

def _cone(A, b, n, kind):
    """Untrusted numerical discovery; only frozen exact checker decides."""
    try:
        import clarabel
        import numpy as np
        import scipy.sparse as sp
        settings = clarabel.DefaultSettings()
        settings.verbose = False
        settings.max_threads = 4
        settings.max_iter = 100
        AA = sp.csc_matrix(np.vstack((np.array(A, dtype=float), -np.eye(n))))
        bb = np.r_[np.array(b, dtype=float), np.zeros(n)]
        cone = clarabel.SecondOrderConeT(n) if kind == 'soc3' else clarabel.NonnegativeConeT(n)
        sol = clarabel.DefaultSolver(sp.csc_matrix((n, n)), np.zeros(n), AA, bb, [clarabel.ZeroConeT(len(A)), cone], settings).solve()
        for den in (10 ** 4, 10 ** 6, 10 ** 8):
            y = [Q(float(q)).limit_denominator(den) for q in sol.z[:len(A)]]
            if check_dual(kind, A, b, y).status == 'NEJ':
                return ('INCONSISTENT', y)
            x = [Q(float(q)).limit_denominator(den) for q in sol.x]
            if check_primal(kind, A, b, x).status == 'JA':
                return ('CONSISTENT', x)
        return ('UNKNOWN', 'NO_EXACT_PRIMAL_OR_STRICT_DUAL; NUMERICAL_STATUS=' + str(sol.status))
    except ImportError:
        return ('UNKNOWN', 'OPTIONAL_CLARABEL_UNAVAILABLE')
    except (ArithmeticError, ValueError, RuntimeError) as e:
        return ('UNKNOWN', 'CONE_DISCOVERY_FAILED:' + str(e))

def _cert(p, subset, status, witness):
    (A, b, labels) = _rows(p, subset)
    return {'kind': p['kind'], 'A': _json(A), 'b': _json(b), 'row_labels': labels, 'claim_ids': [p['groups'][i]['id'] for i in subset], 'status': status, 'witness': _json(witness)}

def _check(p, subset, c):
    (A, b, _) = _rows(p, subset)
    if c.get('status') not in ('CONSISTENT', 'INCONSISTENT'):
        return False
    if c != _cert(p, subset, c['status'], c['witness']):
        return False
    if not A:
        return c['status'] == 'CONSISTENT' and len(c['witness']) == len(p['variables']) and all((exact(q) == 0 for q in c['witness']))
    checker = check_dual if c['status'] == 'INCONSISTENT' else check_primal
    return checker(p['kind'], A, b, c['witness']).status == ('NEJ' if c['status'] == 'INCONSISTENT' else 'JA')

def _unknown(reason):
    return {'status': 'UNKNOWN', 'reason': reason, 'minimal_subset': None, 'dual_weights': [], 'suggested_discriminating_measurement': {'status': 'UNKNOWN', 'reason': reason}, 'certificate': None, 'minimality': 'UNKNOWN'}

def farkas_witness(claims, *, max_claims=12, max_oracle_calls=4096, row_budget=4096):
    """Return JSON-able status/minimum/weights/OED/certificate. No side effects."""
    start = time.perf_counter()
    calls = 0
    try:
        if isinstance(claims, dict) and claims.get('adapter_unknown'):
            return _unknown(claims['adapter_unknown'])
        p = _parse(claims)
        N = len(p['groups'])
        cache = {}
        if N > max_claims:
            return _unknown('CLAIM_COUNT_BUDGET')
        if min(max_claims, max_oracle_calls, row_budget) < 1:
            return _unknown('INVALID_BUDGET')

        def oracle(sub):
            nonlocal calls
            sub = tuple(sub)
            if sub in cache:
                return cache[sub]
            if calls >= max_oracle_calls:
                raise BudgetError('MINIMUM_SEARCH_ORACLE_BUDGET')
            calls += 1
            (A, b, _) = _rows(p, sub)
            if not A:
                (status, w) = ('CONSISTENT', [Q(0)] * len(p['variables']))
            elif p['kind'] == 'inequality':
                (status, w) = _fm(A, b, len(p['variables']), row_budget)
            else:
                (status, w) = _cone(A, b, len(p['variables']), p['kind'])
            if status == 'UNKNOWN':
                cache[sub] = (status, w)
            else:
                c = _cert(p, sub, status, w)
                if not _check(p, sub, c):
                    raise ContractError('INTERNAL_CERTIFICATE_FAILED')
                cache[sub] = (status, c)
            return cache[sub]
        full = tuple(range(N))
        (status, c) = oracle(full)
        if status == 'UNKNOWN':
            return _unknown(c)
        if status == 'CONSISTENT':
            cert = {'schema': 'farkas-certificate/v1', 'input_sha256': p['input_sha256'], 'decision': c, 'minimum_proof': None}
            out = _unknown('EXACT_PRIMAL_WITNESS')
            out.update(status=status, certificate=cert, minimal_subset=[], minimality='NOT_APPLICABLE')
        else:
            (core, dc) = (full, c)
            minimum = None
            blockers = []
            try:
                for k in range(1, N + 1):
                    unresolved = []
                    found = None
                    for sub in combinations(range(N), k):
                        (st, cc) = oracle(sub)
                        if st == 'UNKNOWN':
                            unresolved.append(list(sub))
                        elif st == 'INCONSISTENT':
                            found = (sub, cc)
                            break
                    if found is not None:
                        (core, dc) = found
                        if not blockers:
                            smaller = [oracle(sub)[1] for sub in combinations(range(N), k - 1)]
                            if all((s['status'] == 'CONSISTENT' for s in smaller)):
                                minimum = {'cardinality': k, 'feasible_k_minus_1': smaller}
                        break
                    blockers.extend(unresolved)
            except BudgetError as e:
                blockers.append(str(e))
            cert = {'schema': 'farkas-certificate/v1', 'input_sha256': p['input_sha256'], 'decision': dc, 'minimum_proof': minimum}
            weights = []
            for (label, weight) in zip(dc['row_labels'], dc['witness']):
                if exact(weight):
                    weights.append({**label, 'weight': str(exact(weight)), 'source': p['groups'][next((i for i in core if p['groups'][i]['id'] == label['claim_id']))]['source']})
            out = {'status': 'INCONSISTENT', 'reason': 'EXACT_STRICT_DUAL_CERTIFICATE', 'minimal_subset': dc['claim_ids'] if minimum else None, 'conflict_subset': dc['claim_ids'], 'minimality': 'CARDINALITY_MINIMUM' if minimum else 'UNKNOWN', 'minimum_blockers': blockers, 'dual_weights': weights, 'certificate': cert, 'suggested_discriminating_measurement': _oed(p, core, claims.get('measurements', []), row_budget) if minimum else {'status': 'UNKNOWN', 'reason': 'MINIMUM_NOT_CERTIFIED'}}
        out['cost'] = {'oracle_calls': calls, 'total_query_wall_s': time.perf_counter() - start, 'fit': 'NOT_APPLICABLE', 'physical_validation': 'UNMEASURED'}
        return out
    except (ValueError, TypeError, KeyError, IndexError, AttributeError, ZeroDivisionError, OverflowError) as e:
        return _unknown(str(e))

def verify_certificate(claims, certificate):
    """Rebind full input and independently check primal/dual and cardinality proof."""
    try:
        p = _parse(claims)
        c = certificate
        if c['schema'] != 'farkas-certificate/v1' or c['input_sha256'] != p['input_sha256']:
            return False
        ids = [g['id'] for g in p['groups']]
        d = c['decision']
        sub = tuple((ids.index(i) for i in d['claim_ids']))
        if len(set(sub)) != len(sub) or not _check(p, sub, d):
            return False
        proof = c.get('minimum_proof')
        if proof is None:
            return d['status'] == 'INCONSISTENT' or sub == tuple(range(len(ids)))
        k = proof['cardinality']
        if d['status'] != 'INCONSISTENT' or k != len(sub) or (not 1 <= k <= len(ids)):
            return False
        expected = set(combinations(range(len(ids)), k - 1))
        seen = set()
        for pc in proof['feasible_k_minus_1']:
            ss = tuple((ids.index(i) for i in pc['claim_ids']))
            if ss not in expected or ss in seen or pc['status'] != 'CONSISTENT' or (not _check(p, ss, pc)):
                return False
            seen.add(ss)
        return seen == expected
    except (ValueError, TypeError, KeyError, IndexError, AttributeError, ZeroDivisionError, OverflowError):
        return False

def _project(A, b, m, limit):
    n = len(m)
    lifted = [(Q(0), *a) for a in A]
    lifted += [(Q(1), *(-q for q in m)), (Q(-1), *m)]
    rhs = [*b, Q(0), Q(0)]
    size = len(lifted)
    rows = _simplify([(a, v, tuple((Q(i == j) for j in range(size)))) for (i, (a, v)) in enumerate(zip(lifted, rhs))], limit)
    for _ in range(n):
        rows = _eliminate(rows, limit)
    (lower, upper) = (None, None)
    for (a, v, _) in rows:
        if a[0] > 0:
            val = v / a[0]
            upper = val if upper is None else min(upper, val)
        elif a[0] < 0:
            val = v / a[0]
            lower = val if lower is None else max(lower, val)
        elif v < 0:
            raise ContractError('REPAIR_MODEL_IS_INCONSISTENT')
    return (lower, upper)

def _oed(p, core, measurements, limit):
    if not isinstance(measurements, list):
        return {'status': 'UNKNOWN', 'reason': 'MEASUREMENTS_MUST_BE_A_LIST'}
    unknown = {'status': 'UNKNOWN', 'reason': 'NO_EXPLICIT_MEASUREMENT_CONTRACT'}
    if p['kind'] != 'inequality':
        return {'status': 'UNKNOWN', 'reason': 'CONE_OED_PROJECTION_NOT_SUPPORTED'}
    if len(core) < 2:
        return {'status': 'UNKNOWN', 'reason': 'FEWER_THAN_TWO_REPAIR_HYPOTHESES'}
    ids = [v['id'] for v in p['variables']]
    (candidates, reasons) = ([], [])
    for m in measurements:
        if not isinstance(m, dict):
            reasons.append({'measurement_id': None, 'reason': 'MEASUREMENT_MUST_BE_AN_OBJECT'})
            continue
        try:
            if not isinstance(m.get('id'), str) or not m['id'] or (not m.get('quantity')) or (not m.get('where')) or (m.get('scope') != p['scope']):
                raise ContractError('INCOMPLETE_MEASUREMENT_LOCATION_OR_SCOPE')
            unit = _unit(m.get('unit'))
            coeff = m.get('coefficients', {})
            if not isinstance(coeff, dict):
                raise ContractError('MEASUREMENT_COEFFICIENTS_MUST_BE_AN_OBJECT')
            if set(coeff) - set(ids):
                raise ContractError('MEASUREMENT_UNKNOWN_VARIABLE')
            a = tuple((exact(coeff.get(i, 0)) for i in ids))
            for (j, v) in enumerate(p['variables']):
                if a[j] and v['unit'] != unit:
                    raise ContractError('MEASUREMENT_UNIT_MISMATCH')
            (err, res, cost) = (exact(m['error_bound']), exact(m['resolution']), exact(m['cost']))
            if err < 0 or res <= 0 or cost <= 0 or (not any(a)):
                raise ContractError('INVALID_MEASUREMENT_PRECISION_OR_COST')
            bands = []
            for removed in core:
                sub = tuple((i for i in core if i != removed))
                (A, b, _) = _rows(p, sub)
                (lo, hi) = _project(A, b, a, limit)
                if lo is None or hi is None:
                    raise ContractError('UNBOUNDED_REPAIR_PREDICTION')
                bands.append({'removed_claim': p['groups'][removed]['id'], 'prediction': [str(lo), str(hi)], 'observable_band': [str(lo - err), str(hi + err)]})
            gaps = []
            for (x, y) in combinations(bands, 2):
                (l, h) = map(Q, x['observable_band'])
                (u, v) = map(Q, y['observable_band'])
                gaps.append(max(u - h, l - v))
            gap = min(gaps)
            if gap <= 0:
                raise ContractError('REPAIR_PREDICTIONS_OVERLAP_AT_GIVEN_NOISE')
            candidates.append({'status': 'CERTIFIED_DISCRIMINATING_UNDER_MODEL', 'measurement_id': m['id'], 'quantity': m['quantity'], 'unit': unit, 'where': _json(m['where']), 'scope': _json(p['scope']), 'error_bound': str(err), 'resolution': str(res), 'cost': str(cost), 'hypotheses': bands, 'min_observable_gap': str(gap), 'score': str(gap / res / cost), 'optimality_scope': 'MAX_SCORE_AMONG_SUPPLIED_CERTIFIABLY_DISJOINT_LINEAR_MEASUREMENTS', 'assumptions': 'Exactly one claim in this core is removed; other claims outside the core are excluded. Caller binds model and bounded sensor error to reality.'})
        except (ValueError, KeyError, TypeError, IndexError, AttributeError, ArithmeticError) as e:
            reasons.append({'measurement_id': m.get('id'), 'reason': str(e)})
    if candidates:
        chosen = max(candidates, key=lambda c: (Q(c['score']), c['measurement_id']))
        chosen['unusable_candidates'] = reasons
        return chosen
    if reasons:
        unknown.update(reason='NO_CERTIFIABLY_DISJOINT_MEASUREMENT', candidates=reasons)
    return unknown

def _adapter_guard(func):
    from functools import wraps

    @wraps(func)
    def wrapped(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except (ValueError, TypeError, KeyError, AttributeError, IndexError) as e:
            return {'adapter_unknown': 'MALFORMED_ADAPTER_INPUT:' + str(e)}
    return wrapped

@_adapter_guard
def adapt_claim_federation(graph):
    """Use RAW graph claims before Federation.add_graph discards extension fields."""
    cs = graph.get('claims', [])
    if not cs or any(('farkas_claim' not in c for c in cs)) or (not graph.get('farkas_context')):
        return {'adapter_unknown': 'FEDERATION_SIGN_AND_VALIDITY_ARE_NOT_HARD_CONSTRAINTS'}
    ctx = graph['farkas_context']
    converted = []
    gid = graph.get('graph_id')
    if not isinstance(gid, str) or not gid:
        return {'adapter_unknown': 'MISSING_FEDERATION_GRAPH_ID'}
    for c in cs:
        item = dict(c['farkas_claim'])
        item['id'] = gid + ':' + c['id']
        converted.append(item)
    return {**ctx, 'schema': 'farkas-claims/v1', 'claims': converted, 'federation_binding': {'graph_id': gid, 'raw_graph_sha256': _hash(graph)}}

@_adapter_guard
def adapt_bodytwin(node):
    """Native prose does not supply a mathematical model. Preserve node lineage."""
    if not isinstance(node, dict) or not node.get('farkas_model'):
        return {'adapter_unknown': 'BODYTWIN_NODE_HAS_NO_EXPLICIT_CONSTRAINT_MODEL', 'node_id': node.get('id') if isinstance(node, dict) else None}
    return {**node['farkas_model'], 'native_binding': {'node_id': node['id'], 'node_sha256': _hash(node), 'native_node_fully_resolved': False}}

@_adapter_guard
def adapt_g3_pair(edge):
    """Strictly adapt existing G3 atomic sides, preserving method/regime/quantity.

Point values are nominal model assertions, not manufactured confidence bands.
A proof would be conditional on treating these exact nominal values as hard.
"""
    sides = edge.get('evidence', {}).get('sides', [])
    if len(sides) != 2:
        return {'adapter_unknown': 'G3_EDGE_HAS_NO_TWO_TYPED_SIDES'}
    if any((s.get('survival_status') == 'REJECTED_CONVENTION' for s in sides)):
        return {'adapter_unknown': 'G3_REJECTED_CONVENTION_NOT_A_VALID_CLAIM'}
    fields = ('quantity', 'population', 'regime', 'measurement_method')
    if any((f not in s for f in fields for s in sides)):
        return {'adapter_unknown': 'G3_MISSING_SCOPE_FIELDS'}
    if any((sides[0][f] != sides[1][f] for f in fields)):
        return {'adapter_unknown': 'G3_DIFFERENT_QUANTITY_POPULATION_REGIME_OR_METHOD'}
    scope = {f: sides[0][f] for f in fields if f != 'quantity'}
    claims = []
    for (i, s) in enumerate(sides):
        value = s['value']
        interval = value if isinstance(value, list) else [value, value]
        claims.append({'id': edge['contradiction_node_id'] + ':' + str(i), 'variable': 'q', 'unit': s['unit'], 'interval': interval, 'scope': scope, 'source': s['source']})
    return {'schema': 'farkas-claims/v1', 'variables': [{'id': 'q', 'quantity': sides[0]['quantity'], 'unit': sides[0]['unit']}], 'scope': scope, 'claims': claims, 'nominal_value_interpretation': 'EXACT_MODEL_ASSERTION_NOT_EMPIRICAL_UNCERTAINTY'}
