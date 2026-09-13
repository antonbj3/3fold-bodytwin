"""Explicit numerical gate-record contracts without inference or artifact writes."""
from dataclasses import dataclass
from collections.abc import Mapping


@dataclass(frozen=True)
class GateVerdict:
    accepted: bool
    required_passed: int
    required_total: int
    optional_failed: tuple[str, ...]


@dataclass(frozen=True)
class GateContract:
    names: tuple[str, ...]
    required: tuple[str, ...]
    overall_key: str

    def __post_init__(self):
        for values in (self.names,self.required):
            if type(values) is not tuple or not values or any(type(v) is not str or not v for v in values) or len(set(values))!=len(values):
                raise ValueError('Nonempty unique explicit gate tuples required')
        if not set(self.required)<=set(self.names):
            raise ValueError('Required gates must belong to the full inventory')
        if type(self.overall_key) is not str or not self.overall_key or self.overall_key=='gates':
            raise ValueError('Explicit separate overall key required')

    def assess(self,record):
        if not isinstance(record,Mapping):raise ValueError('Gate record must be a mapping')
        gates=record.get('gates')
        if not isinstance(gates,Mapping) or set(gates)!=set(self.names):
            raise ValueError('Exact declared gate inventory required')
        if any(type(v) is not bool for v in gates.values()) or type(record.get(self.overall_key)) is not bool:
            raise ValueError('Literal boolean gates and overall value required')
        passed=sum(gates[k] for k in self.required)
        accepted=passed==len(self.required)
        if record[self.overall_key]!=accepted:
            raise ValueError('Overall value contradicts required gates')
        optional=tuple(sorted(k for k in self.names if k not in self.required and not gates[k]))
        return GateVerdict(accepted,passed,len(self.required),optional)


def selftest():
    import copy
    contract=GateContract(('a','b','bonus'),('a','b'),'accepted')
    positive={'gates':{'a':True,'b':True,'bonus':True},'accepted':True}
    negative={'gates':{'a':True,'b':False,'bonus':True},'accepted':False}
    optional={'gates':{'a':True,'b':True,'bonus':False},'accepted':True}
    before=copy.deepcopy((positive,negative,optional))
    first=[contract.assess(r) for r in (positive,negative,optional)]
    second=[contract.assess(r) for r in (positive,negative,optional)]
    refused=0
    bad=[None,{},dict(positive,gates={'a':True,'b':True}),dict(positive,gates={'a':True,'b':True,'bonus':True,'extra':True}),dict(positive,gates={'a':1,'b':True,'bonus':True}),dict(positive,accepted=1),dict(negative,accepted=True)]
    for record in bad:
        try:contract.assess(record)
        except ValueError:refused+=1
    declarations=0
    for args in [((),('a',),'accepted'),(('a',),(),'accepted'),(('a','a'),('a',),'accepted'),(('a',),('b',),'accepted'),(('a',),('a',),'gates'),(['a'],('a',),'accepted')]:
        try:GateContract(*args)
        except ValueError:declarations+=1
    return dict(known_outcomes=[r.accepted for r in first]==[True,False,True],
                optional_failure=first[2].optional_failed==('bonus',),
                full_repeat=first==second,nonmutation=before==(positive,negative,optional),
                record_refusals=refused==7,declaration_refusals=declarations==6)


if __name__=='__main__':
    import json
    gates=selftest();print(json.dumps({'gates':gates}))
    raise SystemExit(0 if all(gates.values()) else 2)
