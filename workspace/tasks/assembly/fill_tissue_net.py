#!/usr/bin/env python3
"""Fill the declared-but-empty tissue constraint net, in the canonical schema, additively.

Why this script and not a new file. I first wrote `data/CONSTRAINT_NET_TISSUE.json`, which would have
been a fifth coordinate system. The graph lane measured the actual state on 2026-10-03: the slot already
exists and is empty — `CONSTRAINT_NETS.json` here is 21 bytes holding `{"bodytwin": {}}`, and dental's
229 KB file carries vehicle, projector, factory and the stress map while its own `bodytwin` key is also
empty. An empty declared slot is a different thing from "not built"; it is a gap waiting for exactly
this material. So this writes into that key, in the schema the existing nets use
(`id`, `between`, `constraint`, `status`, `evidence`), and nothing is moved.

Why a constraint net at all. This is the layer the operator means by the graph: quantities in quantified
relation, each edge tightened by a measurement rather than asserted, with status TIGHT, OPEN or UNKNOWN.
The 3950-node MECHANISM graph is not that layer — 3674 of its nodes have an empty `depends_on`, so it is
a well-sourced flat list, and dependency-based priority cannot mean anything over it.

Every edge below carries a number this session read from the producing run's own artefact file. An edge
with no verified number is not written; the relation is recorded as UNKNOWN and says what is missing.
"""
from __future__ import annotations

import json
from pathlib import Path

W = Path('.')
TARGET = W / 'CONSTRAINT_NETS.json'
SOURCE = W / 'data/CONSTRAINT_NET_TISSUE.json'


def main() -> int:
    draft = json.loads(SOURCE.read_text())
    existing = json.loads(TARGET.read_text()) if TARGET.exists() else {}
    if existing.get('bodytwin'):
        print('bodytwin key is not empty; refusing to overwrite. Merge by hand.')
        return 1

    edges = []
    for i, e in enumerate(draft['edges'], 1):
        ev = e['evidence']
        edge = {
            'id': f"T-E{i}-{e['from']}-vs-{e['to']}",
            'between': sorted({e['from'], e['to']}),
            'constraint': e['relation'],
            'status': e['status'],
            'evidence': f"{ev['file']} :: {ev['key']} = {ev['value']}",
        }
        # Optional fields only when present: a generator cannot be unpacked with **, which is what my
        # first version tried.
        if 'owner' in e:
            edge['owner'] = e['owner']
        if 'note' in e:
            edge['note'] = e['note']
        edges.append(edge)

    existing['bodytwin'] = {
        'tissue_constraint_net': {
            'seed': ('The operator asked for the graph layer: variables in quantified relation with '
                     'edge status, not bookkeeping. This net was assembled 2026-10-03 from the night\'s '
                     'own measurements after the graph lane established that this slot was declared and '
                     'empty.'),
            'objective': ('Say, per pair of quantities, whether the relation is constrained by a '
                          'measurement (TIGHT), measured to be insufficient or undetermined with the gap '
                          'as the number (OPEN), or asserted in a model with nothing bearing on it '
                          '(UNKNOWN).'),
            'schema_note': ('Same shape as vehicle_constraints and projector_constraints: variables with '
                            'unit and level, edges with between, constraint, status and evidence. '
                            'Additive — nothing was moved and no other net was touched.'),
            'source_files': [str(SOURCE.relative_to(W))],
            'variables': draft['variables'],
            'edges': edges,
            'stress_points': draft['stress_points'],
            'counts': {
                'variables': len(draft['variables']),
                'edges': len(edges),
                'TIGHT': sum(1 for e in edges if e['status'] == 'TIGHT'),
                'OPEN': sum(1 for e in edges if e['status'] == 'OPEN'),
                'UNKNOWN': sum(1 for e in edges if e['status'] == 'UNKNOWN'),
            },
            'review_state': 'PENDING_INDEPENDENT_REVIEW',
        }
    }
    TARGET.write_text(json.dumps(existing, ensure_ascii=False, indent=1))
    c = existing['bodytwin']['tissue_constraint_net']['counts']
    print(f"  variables {c['variables']}  edges {c['edges']}  "
          f"TIGHT {c['TIGHT']}  OPEN {c['OPEN']}  UNKNOWN {c['UNKNOWN']}")
    print(f"  written into the declared slot: {TARGET.relative_to(W)}")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
