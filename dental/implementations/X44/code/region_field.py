"""K09: native-tetra region partition and lazy signed distances; no voxel approximation.
Distances are piecewise-linear source-geometry distances in mm. Source identity is exact;
float64 numeric queries are tested, not formally interval-certified.
"""
from dental_release.paths import expand as _release_expand
import sys, time, zipfile, json, itertools
import numpy as np
import meshio, trimesh
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from rtree import index
from common import *
sys.path.insert(0, str(ROOT / 'cells/physics'))
from fe_reference_io import read_msh2, tet_volumes, tooth_components
SOURCE = Path(_release_expand('@DENTAL_DATA_ROOT@/geometry/Open-Full-Jaw/dataset/@DENTAL_CASE_ID@.zip'))

def source_files():
    D.mkdir(parents=True, exist_ok=True)
    assert footprint() < 3000000000
    z = zipfile.ZipFile(SOURCE)
    out = {}
    for ext in ['msh', 'vtk']:
        n = _release_expand(f'@DENTAL_CASE_ID@/output/mandible/volumetric_meshes/mandible_volumetric_mesh.{ext}')
        p = D / Path(n).name
        if not p.exists():
            assert footprint() + z.getinfo(n).file_size < 3000000000
            with z.open(n) as f, p.open('wb') as g:
                import shutil
                shutil.copyfileobj(f, g)
        out[ext] = p
    return out

def interfaces(V, T, tag):
    F = np.concatenate([T[:, [1, 2, 3]], T[:, [0, 3, 2]], T[:, [0, 1, 3]], T[:, [0, 2, 1]]])
    own = np.tile(np.arange(len(T)), 4)
    key = np.sort(F, axis=1)
    order = np.lexsort(key.T[::-1])
    key = key[order]
    F = F[order]
    own = own[order]
    same = np.all(key[1:] == key[:-1], axis=1)
    paired = np.flatnonzero(same)
    assert not np.any(same[:-1] & same[1:]), 'non-manifold tet faces'
    unpaired = np.ones(len(F), bool)
    unpaired[paired] = False
    unpaired[paired + 1] = False
    rows = []
    for r in (1, 2, 3):
        a = paired[tag[own[paired]] != tag[own[paired + 1]]]
        choose = np.concatenate([np.flatnonzero(unpaired & (tag[own] == r)), a[tag[own[a]] == r], (a + 1)[tag[own[a + 1]] == r]])
        f = F[choose].copy()
        o = own[choose]
        tri = V[f]
        n = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
        inward = np.einsum('ij,ij->i', n, V[T[o]].mean(1) - tri.mean(1)) > 0
        f[inward] = f[inward][:, [0, 2, 1]]
        rows.append((f, o))
    both = paired[tag[own[paired]] != tag[own[paired + 1]]]
    (a, b) = (own[both], own[both + 1])
    f = F[both]
    area = np.linalg.norm(np.cross(V[f[:, 1]] - V[f[:, 0]], V[f[:, 2]] - V[f[:, 0]]), axis=1) / 2
    return (rows, dict(a=a, b=b, f=f, area=area))

class RegionField:

    def __init__(self, V, T, tag, ec, pc, faces, pdl_tooth):
        (self.V, self.T, self.tag, self.ec, self.pc) = (V, T, tag, ec, pc)
        self.body = np.zeros(len(T), np.int32)
        self.body[tag == 1] = ec[tag == 1] + 1
        for c in np.unique(pc[tag == 2]):
            self.body[pc == c] = 1000 + pdl_tooth[int(c)] + 1
        self.body[tag == 3] = 2000
        self.mesh = {r: trimesh.Trimesh(V, f[0], process=False) for (r, f) in zip((1, 2, 3), faces)}
        for (r, m) in self.mesh.items():
            assert m.is_watertight, ('open native material boundary', r)
        xyz = V[T]
        lo = xyz.min(1)
        hi = xyz.max(1)
        pr = index.Property()
        pr.dimension = 3
        self.spatial = index.Index(((i, tuple(np.r_[lo[i], hi[i]]), None) for i in range(len(T))), properties=pr)

    def owners(self, P):
        reg = np.zeros(len(P), np.int32)
        body = np.zeros(len(P), np.int32)
        witness = np.full(len(P), -1, np.int64)
        for (j, p) in enumerate(P):
            ids = np.array(list(self.spatial.intersection(tuple(np.r_[p - 1e-10, p + 1e-10]))), dtype=int)
            if not len(ids):
                continue
            X = self.V[self.T[ids]]
            a = np.transpose(X[:, 1:] - X[:, 0, None], (0, 2, 1))
            b = np.linalg.solve(a, (p - X[:, 0])[..., None])[..., 0]
            good = np.all(b >= -1e-09, axis=1) & (b.sum(1) <= 1 + 1e-09)
            ix = ids[good]
            if not len(ix):
                continue
            if len(np.unique(self.tag[ix])) > 1:
                raise ValueError('set-valued boundary owner requires sided query')
            q = ix[0]
            reg[j] = self.tag[q]
            body[j] = self.body[q]
            witness[j] = q
        return (reg, body, witness)

    def query(self, P):
        (reg, body, witness) = self.owners(P)
        sdf = []
        for r in (1, 2, 3):
            (_, d, _) = trimesh.proximity.closest_point(self.mesh[r], P)
            sdf.append(np.where(reg == r, -d, d))
        return dict(region=reg, body=body, witness=witness, sdf_mm=np.array(sdf).T)

def run():
    prereg = verify('PREREG_K09.json')
    st = time.perf_counter()
    paths = source_files()
    (V, T, tag) = read_msh2(paths['msh'], cache=False)
    v0 = time.perf_counter()
    other = meshio.read(paths['vtk'])
    Tv = other.cells_dict['tetra']
    tags = other.cell_data_dict
    print('VTK label keys', list(tags), flush=True)
    if not tags:
        raise ValueError('published VTK has no independent region labels')
    candidates = [(k, v['tetra']) for (k, v) in tags.items() if 'tetra' in v and set(np.unique(v['tetra'])) == {1, 2, 3}]
    assert candidates, [(k, np.unique(v['tetra']).tolist()) for (k, v) in tags.items() if 'tetra' in v]
    (key, tv) = candidates[0]
    vol = np.abs(tet_volumes(V, T))
    volv = np.abs(tet_volumes(other.points, Tv))
    (ec, _) = tooth_components(T, tag, len(V), 1)
    (pc, _) = tooth_components(T, tag, len(V), 2)
    (faces, inter) = interfaces(V, T, tag)
    pdl_tooth = {}
    assignment = []
    for c in np.unique(pc[tag == 2]):
        (ai, bi) = (inter['a'], inter['b'])
        sel = (tag[ai] == 1) & (pc[bi] == c) | (tag[bi] == 1) & (pc[ai] == c)
        toothid = np.where(tag[ai[sel]] == 1, ec[ai[sel]], ec[bi[sel]])
        weights = inter['area'][sel]
        u = np.unique(toothid)
        areas = np.array([weights[toothid == i].sum() for i in u])
        assert len(u) == 1, ('PDL body touches multiple teeth', int(c), u.tolist())
        pdl_tooth[int(c)] = int(u[0])
        assignment.append(dict(pdl_component=int(c), tooth_component=int(u[0]), interface_area_mm2=float(areas[0]), resolution='PER_SURFACE_REGION'))
    field = RegionField(V, T, tag, ec, pc, faces, pdl_tooth)
    rng = np.random.default_rng(44)
    ix = np.concatenate([rng.choice(np.flatnonzero(tag == r), 300, replace=False) for r in (1, 2, 3)])
    P = V[T[ix]].mean(1)
    lo = V.min(0)
    hi = V.max(0)
    out = rng.uniform(lo - 2, hi + 2, (60, 3))
    out[:, 2] = hi[2] + 2
    P = np.vstack([P, out])
    expected = np.r_[tag[ix], np.zeros(60, int)]
    q = field.query(P)
    sign = []
    dist = []
    for r in (1, 2, 3):
        inside = field.mesh[r].contains(P)
        sign.append(int(np.sum(inside != (expected == r))))
        pi = P[(r - 1) * 300:(r - 1) * 300 + 3]
        brute = []
        for p in pi:
            (_, d, _) = trimesh.proximity.closest_point_naive(field.mesh[r], p[None])
            brute.append(d[0])
        dist.append(float(np.max(np.abs(np.abs(q['sdf_mm'][(r - 1) * 300:(r - 1) * 300 + 3, r - 1]) - brute))))
    regions = []
    for r in (1, 2, 3):
        (f, _) = faces[r - 1]
        a = np.linalg.norm(np.cross(V[f[:, 1]] - V[f[:, 0]], V[f[:, 2]] - V[f[:, 0]]), axis=1).sum() / 2
        t = Tv[tv == r]
        ff = np.concatenate([t[:, [1, 2, 3]], t[:, [0, 3, 2]], t[:, [0, 1, 3]], t[:, [0, 2, 1]]])
        keys = np.sort(ff, axis=1)
        (uu, inds, cnt) = np.unique(keys, axis=0, return_index=True, return_counts=True)
        sf = ff[inds[cnt == 1]]
        av = np.linalg.norm(np.cross(other.points[sf[:, 1]] - other.points[sf[:, 0]], other.points[sf[:, 2]] - other.points[sf[:, 0]]), axis=1).sum() / 2
        vv = float(vol[tag == r].sum())
        vr = float(volv[tv == r].sum())
        regions.append(dict(region=r, tets=int((tag == r).sum()), volume_mm3=vv, reference_volume_mm3=vr, volume_relative_error=abs(vv - vr) / vr, boundary_area_mm2=float(a), reference_boundary_area_mm2=float(av), area_relative_error=abs(a - av) / av, resolution='PER_SURFACE_REGION'))
    mismatch = int((q['region'] != expected).sum())
    mut = q['region'].copy()
    mut[expected == 2] = 1
    mut_rejected = bool(np.any(mut != expected))
    gates = dict(volume=all((x['volume_relative_error'] <= 1e-08 for x in regions)), area=all((x['area_relative_error'] <= 1e-08 for x in regions)), owners=mismatch == 0, sign=sum(sign) == 0, distance=max(dist) <= 1e-08, mutation_rejected=mut_rejected, pdl_mapping=len(pdl_tooth) == int(ec.max() + 1))
    np.savez_compressed(D / 'region_field.npz', V=V, T=T, tag=tag, body=field.body, ec=ec, pc=pc, **{f'faces_{r}': faces[r - 1][0] for r in (1, 2, 3)})
    np.savez_compressed(L / 'raw/K09_PROBES.npz', points_mm=P, expected_region=expected, **q)
    report = dict(chain='K09', claim_type='capability', outcome='PASS_SOURCE_GEOMETRY' if all(gates.values()) else 'FAIL', resolution='PER_POINT', timescale='SIMULTANEOUS', external_referent=prereg['external_referent'], regions=regions, gates=gates, ownership_mismatch_count=mismatch, sign_mismatch_counts=sign, distance_max_errors_mm=dist, body_counts=dict(teeth=int(ec.max() + 1), pdl=int(pc.max() + 1), bone=1), pdl_tooth_interfaces=assignment, rejection=dict(candidate_tets=len(T), rejected=0, reason='All native valid tets retained; boundary owner queried with sided semantics'), source_inputs=[dict(path=str(p), sha256=sha(p)) for p in paths.values()], control='Independent published VTK parse + surface ray parity + direct all-triangle distance; same information parity, no method superiority', cost=cost(st), uncertainty='Float64 source piecewise-linear geometry only; patient/histology accuracy UNKNOWN. Full jaw field is lazy, not a voxel grid; no edit/remesh/FE stability claim.')
    write(L / 'raw/K09.json', report)
    state('K09_COMPLETE', gates, 'K13 compile source material and conditional PDL')
    print(json.dumps(dict(chain='K09', gates=gates, cost=report['cost'])), flush=True)
if __name__ == '__main__':
    run()
