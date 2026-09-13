"""Explicit three-state coupling evidence; abstention never means acceptance."""
from collections.abc import Mapping
from dataclasses import dataclass


@dataclass(frozen=True)
class EvidenceVerdict:
    state: str
    missing: tuple[str, ...]

    @property
    def accepted(self):
        return self.state == 'PASS'


@dataclass(frozen=True)
class CouplingEvidence:
    names: tuple[str, ...]

    def __post_init__(self):
        if (type(self.names) is not tuple or not self.names
                or any(type(k) is not str or not k for k in self.names)
                or len(set(self.names)) != len(self.names)):
            raise ValueError('Declare unique nonempty comparison names')

    def assess(self, comparisons):
        if not isinstance(comparisons, Mapping) or set(comparisons) - set(self.names):
            raise ValueError('Undeclared comparison inventory')
        if any(v is not None and type(v) is not bool for v in comparisons.values()):
            raise ValueError('Comparisons must be literal booleans or null')
        missing = tuple(k for k in self.names if comparisons.get(k) is None)
        if missing:
            return EvidenceVerdict('ABSTAIN', missing)
        return EvidenceVerdict('PASS' if all(comparisons.values()) else 'FAIL', ())


def self_test():
    c = CouplingEvidence(('a', 'b'))
    sample = {'a': True, 'b': None}
    before = sample.copy()
    cases = ({}, sample, {'a': None, 'b': None}, {'a': False, 'b': None})
    refused = 0
    for value in ([], {'x': True}, {'a': 1}, {'a': 'true'}, {'b': 0.0}):
        try:
            c.assess(value)
        except ValueError:
            refused += 1
    declarations = 0
    for names in ((), ['a'], ('a', 'a'), ('',), (1,)):
        try:
            CouplingEvidence(names)
        except ValueError:
            declarations += 1
    return dict(positive_negative=c.assess({'a': True, 'b': True}).accepted and c.assess({'a': True, 'b': False}).state == 'FAIL',
                abstentions=all(c.assess(v).state == 'ABSTAIN' and not c.assess(v).accepted for v in cases),
                refusals=refused == 5 and declarations == 5,
                nonmutation=sample == before,
                repeat=[c.assess(v) for v in cases] == [c.assess(v) for v in cases])


if __name__ == '__main__':
    import json
    gates = self_test()
    print(json.dumps({'gates': gates, 'overall_pass': all(gates.values())}, sort_keys=True))
    raise SystemExit(0 if all(gates.values()) else 2)
