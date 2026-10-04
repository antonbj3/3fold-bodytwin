"""Dependency-free PORT v1 types; exact transport is not physical certification.

All lengths have an explicit unit. Signed gaps are surface clearance (+ apart).
UNKNOWN and MULTIPLE are distinct. Importing a witness does not replay it.
"""
from dataclasses import dataclass, fields, is_dataclass
from enum import Enum
from fractions import Fraction as Q
from hashlib import sha256
import json, math

class Unit(str, Enum):
    M = 'm'
    MM = 'mm'

class Status(str, Enum):
    OPEN = 'OPEN'
    CLOSED = 'CLOSED'
    UNIQUE = 'UNIQUE'
    MULTIPLE = 'MULTIPLE'
    EMPTY = 'EMPTY'
    UNKNOWN = 'UNKNOWN'

class Assurance(str, Enum):
    MODEL_BOUND = 'MODEL_BOUND'
    CONDITIONAL = 'CONDITIONAL'
    NUMERICAL = 'NUMERICAL'
    POINT = 'POINT'

def status(value):
    aliases = {'UNCERTAIN': 'UNKNOWN', 'UNCERTAIN': 'UNKNOWN', 'UNCERTAIN': 'UNKNOWN', 'AMBIGUOUS': 'MULTIPLE', 'UNIQUE': 'UNIQUE', 'BRANCH_SET': 'MULTIPLE'}
    return Status(aliases.get(value, value))

def exact(x):
    if isinstance(x, bool):
        raise ValueError('bool is not a measurement')
    if isinstance(x, Q):
        return x
    if isinstance(x, (int, str)):
        return Q(x)
    if isinstance(x, float) and math.isfinite(x):
        return Q(x)
    raise ValueError('finite integer, rational string, Fraction or binary float required')

def encode(x):
    if isinstance(x, Q):
        return str(x)
    if isinstance(x, Enum):
        return x.value
    if is_dataclass(x):
        return {f.name: encode(getattr(x, f.name)) for f in fields(x)}
    if isinstance(x, dict):
        return {str(k): encode(v) for (k, v) in x.items()}
    if isinstance(x, (tuple, list)):
        return [encode(v) for v in x]
    if isinstance(x, float):
        if not math.isfinite(x):
            raise ValueError('nonfinite JSON')
        return str(exact(x))
    if x is None or type(x) in (str, int, bool):
        return x
    raise ValueError('unsupported JSON value')

def canonical(x):
    return json.dumps(encode(x), sort_keys=True, separators=(',', ':'), ensure_ascii=False, allow_nan=False)

def digest(x):
    return sha256(canonical(x).encode()).hexdigest()

def required(x):
    if not isinstance(x, str) or not x.strip():
        raise ValueError('nonempty identity/scope required')

def hash_required(x):
    if not isinstance(x, str) or len(x) != 64 or any((c not in '0123456789abcdef' for c in x)):
        raise ValueError('full lowercase SHA-256 required')

@dataclass(frozen=True)
class Source:
    id: str
    sha256: str
    scope: str

    def __post_init__(self):
        required(self.id)
        hash_required(self.sha256)
        required(self.scope)

@dataclass(frozen=True)
class Gap:
    lower: Q
    upper: Q
    unit: Unit
    source: Source
    sigma: Q | None = None
    assurance: Assurance = Assurance.NUMERICAL
    condition: str = ''
    blocked: bool = False

    def __post_init__(self):
        object.__setattr__(self, 'unit', Unit(self.unit))
        object.__setattr__(self, 'assurance', Assurance(self.assurance))
        for k in ('lower', 'upper'):
            object.__setattr__(self, k, exact(getattr(self, k)))
        if self.lower > self.upper:
            raise ValueError('unordered gap interval')
        if self.sigma is not None:
            object.__setattr__(self, 'sigma', exact(self.sigma))
            if self.sigma < 0:
                raise ValueError('negative sigma')
        if not isinstance(self.source, Source) or type(self.blocked) is not bool:
            raise ValueError('typed source and boolean completion guard required')
        if self.assurance == Assurance.CONDITIONAL:
            required(self.condition)

    def to(self, unit):
        unit = Unit(unit)
        ratio = Q(1) if unit == self.unit else Q(1, 1000) if unit == Unit.M else Q(1000)
        return Gap(self.lower * ratio, self.upper * ratio, unit, self.source, None if self.sigma is None else self.sigma * ratio, self.assurance, self.condition, self.blocked)

    def classify(self, threshold=0):
        t = exact(threshold)
        return Status.UNKNOWN if self.blocked else Status.CLOSED if self.upper < t else Status.OPEN if self.lower > t else Status.UNKNOWN

@dataclass(frozen=True)
class ContactCandidate:
    id: str
    participants: tuple[str, str]
    frame: str
    gap: Gap
    normal: tuple[Q, Q, Q] | None = None

    def __post_init__(self):
        required(self.id)
        required(self.frame)
        object.__setattr__(self, 'participants', tuple(self.participants))
        if len(self.participants) != 2 or self.participants[0] == self.participants[1]:
            raise ValueError('two distinct explicit participants required')
        for p in self.participants:
            required(p)
        if not isinstance(self.gap, Gap):
            raise ValueError('Gap required')
        if self.normal is not None:
            n = tuple((exact(x) for x in self.normal))
            if len(n) != 3 or sum((x * x for x in n)) != 1:
                raise ValueError('exact unit normal required')
            object.__setattr__(self, 'normal', n)

@dataclass(frozen=True)
class Witness:
    kind: str
    payload_json: str
    input_sha256: str
    source: Source

    def __post_init__(self):
        required(self.kind)
        hash_required(self.input_sha256)
        if not isinstance(self.source, Source):
            raise ValueError('witness provenance required')
        if canonical(json.loads(self.payload_json)) != self.payload_json:
            raise ValueError('canonical witness JSON required')

@dataclass(frozen=True)
class Branch:
    id: str
    state_json: str
    witnesses: tuple[Witness, ...]

    def __post_init__(self):
        required(self.id)
        if canonical(json.loads(self.state_json)) != self.state_json:
            raise ValueError('canonical state JSON required')
        object.__setattr__(self, 'witnesses', tuple(self.witnesses))
        if not self.witnesses or any((not isinstance(w, Witness) for w in self.witnesses)):
            raise ValueError('every branch needs a witness')

@dataclass(frozen=True)
class BranchSet:
    branches: tuple[Branch, ...]
    complete: bool
    scope: str
    coverage: Witness | None = None
    reason: str = ''

    def __post_init__(self):
        object.__setattr__(self, 'branches', tuple(self.branches))
        required(self.scope)
        if type(self.complete) is not bool or any((not isinstance(b, Branch) for b in self.branches)):
            raise ValueError('typed branches and boolean completeness required')
        if len({b.id for b in self.branches}) != len(self.branches):
            raise ValueError('duplicate branch IDs')
        if self.complete and (not isinstance(self.coverage, Witness)):
            raise ValueError('complete claim requires coverage witness')

    @property
    def status(self):
        if not self.complete:
            return Status.UNKNOWN
        return Status.EMPTY if not self.branches else Status.UNIQUE if len(self.branches) == 1 else Status.MULTIPLE

@dataclass(frozen=True)
class AffineGap:
    """g=constant+sum(coeff*latent); domains in Port, all quantities in metres."""
    candidate_id: str
    constant_m: Q
    coefficients: tuple[tuple[str, Q], ...]

    def __post_init__(self):
        required(self.candidate_id)
        object.__setattr__(self, 'constant_m', exact(self.constant_m))
        c = tuple(((key, exact(v)) for (key, v) in self.coefficients))
        for (key, v) in c:
            required(key)
        if len({k for (k, v) in c}) != len(c):
            raise ValueError('duplicate latent variable')
        object.__setattr__(self, 'coefficients', c)

    def bounds(self, box):
        lo = hi = self.constant_m
        for (key, c) in self.coefficients:
            (a, b) = box[key]
            v = sorted((c * a, c * b))
            lo += v[0]
            hi += v[1]
        return (lo, hi)

@dataclass(frozen=True)
class Port:
    producer: str
    consumer: str
    contacts: tuple[ContactCandidate, ...]
    branch_set: BranchSet
    latent_box: tuple[tuple[str, Q, Q], ...] = ()
    affine_gaps: tuple[AffineGap, ...] = ()
    evidence: tuple[Witness, ...] = ()
    physical_status: Status = Status.UNKNOWN
    schema: str = 'threefold-contact-port/v1'

    def __post_init__(self):
        required(self.producer)
        required(self.consumer)
        object.__setattr__(self, 'evidence', tuple(self.evidence))
        if any((not isinstance(w, Witness) for w in self.evidence)):
            raise ValueError('typed evidence required')
        if self.schema != 'threefold-contact-port/v1':
            raise ValueError('PORT version unsupported')
        if not isinstance(self.branch_set, BranchSet):
            raise ValueError('BranchSet required')
        object.__setattr__(self, 'contacts', tuple(self.contacts))
        if any((not isinstance(c, ContactCandidate) for c in self.contacts)):
            raise ValueError('typed contacts required')
        if len({c.id for c in self.contacts}) != len(self.contacts):
            raise ValueError('duplicate candidate IDs')
        object.__setattr__(self, 'physical_status', Status(self.physical_status))
        box = tuple(((k, exact(a), exact(b)) for (k, a, b) in self.latent_box))
        if len({k for (k, a, b) in box}) != len(box) or any((a > b for (k, a, b) in box)):
            raise ValueError('ordered unique latent bounds required')
        for (k, a, b) in box:
            required(k)
        object.__setattr__(self, 'latent_box', box)
        object.__setattr__(self, 'affine_gaps', tuple(self.affine_gaps))
        candidates = {c.id: c for c in self.contacts}
        domains = {k: (a, b) for (k, a, b) in box}
        if len({a.candidate_id for a in self.affine_gaps}) != len(self.affine_gaps):
            raise ValueError('duplicate affine candidate')
        for a in self.affine_gaps:
            if not isinstance(a, AffineGap) or a.candidate_id not in candidates:
                raise ValueError('unbound affine gap')
            try:
                (lo, hi) = a.bounds(domains)
            except KeyError as exc:
                raise ValueError('missing latent domain') from exc
            g = candidates[a.candidate_id].gap.to(Unit.M)
            if not g.lower <= lo <= hi <= g.upper:
                raise ValueError('affine gap escapes marginal')
        if self.physical_status != Status.UNKNOWN:
            raise ValueError('physical certification not supported by PORT v1')

    def to_json(self):
        return canonical(self)

    @classmethod
    def from_json(cls, text):

        def strict(d, keys):
            if set(d) != set(keys):
                raise ValueError('missing/unknown PORT fields')
        d = json.loads(text)
        strict(d, [f.name for f in fields(cls)])

        def source(s):
            strict(s, ['id', 'sha256', 'scope'])
            return Source(**s)

        def gap(g):
            strict(g, [f.name for f in fields(Gap)])
            return Gap(**{**g, 'source': source(g['source'])})

        def witness(w):
            strict(w, [f.name for f in fields(Witness)])
            return Witness(**{**w, 'source': source(w['source'])})

        def branch(b):
            strict(b, [f.name for f in fields(Branch)])
            return Branch(**{**b, 'witnesses': tuple((witness(w) for w in b['witnesses']))})
        cs = []
        for c in d['contacts']:
            strict(c, [f.name for f in fields(ContactCandidate)])
            cs.append(ContactCandidate(**{**c, 'gap': gap(c['gap'])}))
        bs = d['branch_set']
        strict(bs, [f.name for f in fields(BranchSet)])
        bs = BranchSet(**{**bs, 'branches': tuple((branch(b) for b in bs['branches'])), 'coverage': None if bs['coverage'] is None else witness(bs['coverage'])})
        aff = []
        for a in d['affine_gaps']:
            strict(a, [f.name for f in fields(AffineGap)])
            aff.append(AffineGap(**a))
        result = cls(**{**d, 'contacts': tuple(cs), 'branch_set': bs, 'affine_gaps': tuple(aff), 'evidence': tuple((witness(w) for w in d['evidence']))})
        if json.dumps(d, sort_keys=True, separators=(',', ':'), ensure_ascii=False, allow_nan=False) != result.to_json():
            raise ValueError('wire numbers must be canonical rational strings')
        return result

def replay_witnesses(branch_set, replayers):
    """Replay input binding and every witness. No replay means no admission.

    replayers[kind](payload) must verify coverage as well as local existence.
    Identity of the replayer/code is an obligation in the caller's manifest.
    """
    ws = [w for b in branch_set.branches for w in b.witnesses]
    if branch_set.coverage is not None:
        ws.append(branch_set.coverage)
    if not ws:
        return False
    for w in dict.fromkeys(ws):
        data = json.loads(w.payload_json)
        if digest(data) != w.input_sha256 or w.kind not in replayers:
            return False
        if not replayers[w.kind](data):
            return False
    return True
