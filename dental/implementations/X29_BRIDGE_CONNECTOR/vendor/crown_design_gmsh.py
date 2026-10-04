"""crown_design_gmsh.py — tetmeshing of a closed STL area with gmsh (run with /usr/bin/python 3 , where gmsh 4.15 is located).

Usage: /usr/bin/python 3 crown_design_gmsh.py <in.stl> <out.npz> <size_max_mm> [ size_min_mm ]
Surface triangles are retained (discreet surface); internal tetraeders grow towards size_max. Output : V ( N , 3 ) [ mm ], T (M, 4 ) 0-based.
"""
import sys
import numpy as np
import gmsh
if __name__ == '__main__':
    (stl, out, smax) = (sys.argv[1], sys.argv[2], float(sys.argv[3]))
    smin = float(sys.argv[4]) if len(sys.argv) > 4 else 0.0
    gmsh.initialize()
    gmsh.option.setNumber('General.Terminal', 0)
    gmsh.option.setNumber('General.NumThreads', 4)
    gmsh.merge(stl)
    surfs = gmsh.model.getEntities(2)
    loop = gmsh.model.geo.addSurfaceLoop([s[1] for s in surfs])
    gmsh.model.geo.addVolume([loop])
    gmsh.model.geo.synchronize()
    gmsh.option.setNumber('Mesh.MeshSizeMax', smax)
    if smin > 0:
        gmsh.option.setNumber('Mesh.MeshSizeMin', smin)
    gmsh.option.setNumber('Mesh.MeshSizeExtendFromBoundary', 1)
    gmsh.option.setNumber('Mesh.Algorithm3D', 1)
    gmsh.option.setNumber('Mesh.Optimize', 1)
    gmsh.model.mesh.generate(3)
    (nt, nc, _) = gmsh.model.mesh.getNodes()
    V = nc.reshape(-1, 3)
    idx = {int(t): i for (i, t) in enumerate(nt)}
    (et, etags, enodes) = gmsh.model.mesh.getElements(3)
    T = None
    for (k, t) in enumerate(et):
        if t == 4:
            T = np.array([idx[int(n)] for n in enodes[k]], np.int64).reshape(-1, 4)
    gmsh.finalize()
    used = np.unique(T)
    remap = -np.ones(len(V), np.int64)
    remap[used] = np.arange(len(used))
    np.savez(out, V=V[used], T=remap[T])
    print('tets', len(T), 'nodes', len(used))
