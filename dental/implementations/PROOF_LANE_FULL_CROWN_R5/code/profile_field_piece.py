from construct_b import *
t = time.perf_counter()
(r, rr) = inputs()[1]
m = npz(r['mesh_path'])
tri = m['vertices'][m['faces'][m['roles'] == 1]]
(v, f) = compact(m['vertices'], m['faces'][m['roles'] == 1][:, ::-1])
mesh0 = trimesh.Trimesh(v, f, process=False)
rim = v[loops(mesh0)[0]]
xyz = np.random.default_rng(29).uniform(v.min(0), v.max(0), (65536, 3))
t0 = time.perf_counter()
(d, i, q, n) = igl.signed_distance(xyz, v, f, igl.SIGNED_DISTANCE_TYPE_PSEUDONORMAL)
t1 = time.perf_counter()
gap = boundary_gap(q, rim)
t2 = time.perf_counter()
print(json.dumps(dict(signed_distance65536_s=t1 - t0, boundary_gap65536_s=t2 - t1, triangles=len(f), rim_vertices=len(rim)), indent=2))
