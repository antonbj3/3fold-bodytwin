"""Write each brief from the node's REGION in the net, not from a one-hop edge list.

The operator's correction: a neighbourhood is a volume, not a list. One hop shows what a quantity is
directly tied to; two hops show the path through which it could unlock something it does not touch.
A job that only sees its own edge can only look for a second reading of its own number, which is the
narrow thing I kept asking for.

It also exposed a defect in the net rather than in the brief. `stromal_index` is NAMED in the
constraint on T-E3 -- "thickness, index, hydration AND surface allocation together leave 0.3056 D
undetermined" -- but it was not an endpoint of that edge, so a neighbourhood walk returned zero edges
for it. Being named in the text of a constraint is not the same as being connected by it, and only
generating a brief from the graph made that visible. The endpoint is now added.
"""
from __future__ import annotations

import json
from pathlib import Path

W = Path('')
NET = W / 'CONSTRAINT_NETS.json'
HOPS = 2

NODE_FOR = {
    'CORNEAL-THICKNESS': 'stromal_thickness',
    'LENS-POSITION': 'lens_position',
    'POSTERIOR-CORNEAL-CYLINDER': 'posterior_corneal_astigmatism',
    'STROMAL-INDEX': 'stromal_index',
    'IOL-POWER-LABEL': 'iol_label_power',
    'AXIAL-LENGTH': 'axial_length',
}


def region(edges: list[dict], start: str, hops: int = HOPS):
    """The induced subgraph within `hops` of start: nodes with their distance, and the edges inside it."""
    dist = {start: 0}
    frontier = {start}
    for h in range(1, hops + 1):
        nxt = set()
        for e in edges:
            if any(v in frontier for v in e['between']):
                for v in e['between']:
                    if v not in dist:
                        dist[v] = h
                        nxt.add(v)
        frontier = nxt
        if not frontier:
            break
    inside = [e for e in edges if all(v in dist for v in e['between'])]
    return dist, inside


def main() -> int:
    net = json.loads(NET.read_text())['bodytwin']['tissue_constraint_net']
    edges = net['edges']
    for tag, node in NODE_FOR.items():
        dist, inside = region(edges, node)
        direct = [e for e in inside if node in e['between']]
        further = [e for e in inside if node not in e['between']]
        openish = [e for e in inside if e['status'] in ('OPEN', 'UNKNOWN')]

        def fmt(e):
            other = ' + '.join(v for v in e['between'] if v != node) or node
            return f"- **{other}** [{e['status']}] {e['constraint']}"

        body = [f'# {node}', '',
                f'The node is in our constraint network. Here is the region around it: '
                f'{len(dist)} nodes within two steps, {len(inside)} edges between them, '
                f'{len(openish)} of them undetermined. All numbers are our own measurements.', '',
                '## Directly on the node']
        body += [fmt(e) for e in direct] or ['- (none)']
        if further:
            body += ['', '## One step away — this is where the path continues']
            body += [f"- **{' + '.join(e['between'])}** [{e['status']}] {e['constraint'][:190]}"
                     for e in further]
        body += ['',
                 'TIGHT = the relation established with a number. OPEN = undetermined, the gap is given as a number. '
                 'UNKNOWN = we do not know.', '',
                 'What would move any of the edges in this region? It need not be the node '
                 'itself — a measurement two steps away can unlock it. Tell us what you find, how it '
                 'differs from ours, and which edge it concerns.', '',
                 'Give DOI or PMID. Status PENDING_INDEPENDENT_REVIEW.', '']
        out = W / 'results' / f'BT-2ND-{tag}' / 'BRIEF.md'
        out.write_text('\n'.join(body))
        print(f'  {tag:28s} {len(dist):>2d} nodes, {len(inside):>2d} edges '
              f'({len(direct)} direct, {len(openish)} undetermined)')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
