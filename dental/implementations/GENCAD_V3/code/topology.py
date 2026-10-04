"""Face-fan topology without altering source triangle geometry."""
import collections, numpy as np

def split_fans(v, f):
    v = np.asarray(v)
    f = np.asarray(f, int)
    edges = collections.defaultdict(list)
    incident = collections.defaultdict(list)
    for (i, face) in enumerate(f):
        for (j, a) in enumerate(face):
            incident[int(a)].append((i, j))
        for (a, b) in zip(face, np.roll(face, -1)):
            edges[tuple(sorted((int(a), int(b))))].append(i)
    if any((len(x) > 2 for x in edges.values())):
        raise ValueError('nonmanifold source edge cannot be split by vertex fans')
    out = np.empty_like(f)
    verts = []
    for (a, corners) in incident.items():
        remaining = {i for (i, j) in corners}
        while remaining:
            todo = [min(remaining)]
            component = set()
            while todo:
                i = todo.pop()
                if i not in remaining:
                    continue
                remaining.remove(i)
                component.add(i)
                for b in f[i]:
                    if int(b) != a:
                        todo.extend((j for j in edges[tuple(sorted((a, int(b))))] if j in remaining))
            index = len(verts)
            verts.append(v[a])
            for (i, j) in corners:
                if i in component:
                    out[i, j] = index
    vv = np.asarray(verts)
    if not np.array_equal(vv[out], v[f]):
        raise ValueError('source triangle moved during topology split')
    directed = collections.defaultdict(list)
    for (i, face) in enumerate(out):
        for (a, b) in zip(face, np.roll(face, -1)):
            directed[tuple(sorted((int(a), int(b))))].append((i, 1 if a < b else -1))
    adjacency = collections.defaultdict(list)
    for pairs in directed.values():
        if len(pairs) == 2:
            ((i, di), (j, dj)) = pairs
            adjacency[i].append((j, -di * dj))
            adjacency[j].append((i, -di * dj))
    signs = {}
    for start in range(len(out)):
        if start in signs:
            continue
        signs[start] = 1
        todo = [start]
        while todo:
            i = todo.pop()
            for (j, factor) in adjacency[i]:
                want = signs[i] * factor
                if j in signs:
                    if signs[j] != want:
                        raise ValueError('nonorientable source patch')
                else:
                    signs[j] = want
                    todo.append(j)
    flipped = np.array([i for (i, s) in signs.items() if s < 0], int)
    if len(flipped):
        out[flipped] = out[flipped][:, ::-1]
    before = np.sort(v[f], axis=1)
    after = np.sort(vv[out], axis=1)
    if not np.array_equal(before, after):
        raise ValueError('native face changed')
    return (vv, out, dict(input_vertices=len(v), output_vertices=len(vv), flipped_faces=len(flipped), source_coordinate_change_mm=0.0))
