from common import *
from collections import defaultdict

def topology(faces):
    edges = defaultdict(list)
    adj = defaultdict(set)
    for (i, face) in enumerate(faces):
        for (a, b) in zip(face, face[1:] + face[:1]):
            edges[tuple(sorted((a, b)))].append(i)
    for fs in edges.values():
        for i in fs:
            adj[i].update(fs)
    unseen = set(range(len(faces)))
    components = 0
    while unseen:
        components += 1
        stack = [unseen.pop()]
        while stack:
            q = adj[stack.pop()] & unseen
            unseen.difference_update(q)
            stack.extend(q)
    boundary = [e for (e, fs) in edges.items() if len(fs) == 1]
    ba = defaultdict(set)
    for (a, b) in boundary:
        ba[a].add(b)
        ba[b].add(a)
    nodes = set(ba)
    cycles = 0
    while nodes:
        cycles += 1
        stack = [nodes.pop()]
        while stack:
            q = ba[stack.pop()] & nodes
            nodes.difference_update(q)
            stack.extend(q)
    manifold = bool(boundary and all((len(x) == 2 for x in ba.values())) and all((len(x) <= 2 for x in edges.values())))
    return dict(components=components, section_components=cycles, boundary_edges=len(boundary), manifold_section=manifold, all_pass=components == 1 and cycles == 1 and manifold)

def clipped(v, f, d, h):
    s = [sum((F(float(a)) * b for (a, b) in zip(p, d))) - h for p in v]
    assert all((x != 0 for x in s)), 'Exact plane/source-vertex degeneracy'
    crossed = {}
    faces = []
    for face in f:
        p = []
        for (a, b) in zip(face, np.roll(face, -1)):
            (a, b) = (int(a), int(b))
            if s[a] < 0:
                p.append(a)
            if (s[a] < 0) != (s[b] < 0):
                edge = tuple(sorted((a, b)))
                if edge not in crossed:
                    crossed[edge] = len(v) + len(crossed)
                p.append(crossed[edge])
        if len(p) >= 3:
            faces.append(p)
    return faces

def run():
    rows = []
    faults = []
    for part in read(R / 'EXPORTS.json')['rows']:
        rec = next((x for x in inputs() if x['key'] == part['key']))
        (v, f, info) = native(rec)
        d = [F(x) for x in part['full_preparation_formula']['d_exact']]
        h = F(part['full_preparation_formula']['h_exact'])
        polys = clipped(v, f, d, h)
        t = topology(polys)
        assert t['all_pass']
        if not faults:
            i = next((i for (i, p) in enumerate(polys) if all((x < len(v) for x in p))))
            remove = topology(polys[:i] + polys[i + 1:])
            extra = topology(polys + [[10 ** 8, 10 ** 8 + 1, 10 ** 8 + 2]])
            faults = [dict(name='removed_retained_facet_rejected', passed=not remove['all_pass'], observed=remove), dict(name='isolated_lower_surface_rejected', passed=not extra['all_pass'], observed=extra)]
            assert all((x['passed'] for x in faults))
        rows.append(dict(key=part['key'], material=part['material'], source_path=info['path'], source_sha256=sha(info['path']), d_exact=part['full_preparation_formula']['d_exact'], h_exact=str(h), clipped_facets=len(polys), topology=t, resolution='PER_SURFACE_REGION', scope='Exact clipping incidence from nonzero rational plane signs, on the inherited closed embedded native T. One section cycle caps one connected lower solid. No rounded intersection coordinate used.'))
    out = dict(claim_type='capability', rows=rows, controls=faults, all_pass=all((x['topology']['all_pass'] for x in rows)))
    dump(R / 'RESULTS_LOWER_TOPOLOGY.json', out)
    print('Exact lower-native topology', len(rows), 'passed;2fault controls passed', flush=True)
    return out
if __name__ == '__main__':
    run()
