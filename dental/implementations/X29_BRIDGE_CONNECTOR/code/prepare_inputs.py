"""Pin only the required SDF and the actual existing FE operator, not a source tree."""
from dental_release.paths import expand as _release_expand
from pathlib import Path
import shutil, json, subprocess
from common import R, X1B, sha, dump, read, now

def run():
    items = [(X1B / 'inputs/geometry/D1_model.npz', 'inputs/D1_model.npz'), (X1B / 'inputs/geometry/D1_grid.json', 'inputs/D1_grid.json')]
    W = Path(_release_expand('@DENTAL_INPUT_ROOT@/workspace'))
    for (name, src) in [('crown_design_fe.py', W / 'cells/design/crown_design_fe.py'), ('crown_design_gmsh.py', W / 'cells/design/crown_design_gmsh.py'), ('crown_fit_geometry.py', W / 'cells/manufacturing/crown_fit_geometry.py'), ('ccx_highres.py', W / 'cells/solvers/ccx_highres.py')]:
        items.append((src, 'vendor/' + name))
    manifest = []
    for (p, local) in items:
        q = R / local
        q.parent.mkdir(exist_ok=True, parents=True)
        if q.exists() and sha(q) != sha(p):
            raise ValueError('Pinned input changed ' + local)
        shutil.copyfile(p, q)
        manifest.append(dict(source=str(p), local=local, sha256=sha(q), bytes=q.stat().st_size))
    solver = X1B / 'vendor/solver/ccx'
    manifest.append(dict(source=str(solver), local=None, sha256=sha(solver), bytes=solver.stat().st_size, role='shared_local_solver_readonly'))
    dump('INPUT_MANIFEST.json', dict(created_utc=now(), items=manifest, solver=str(solver), solver_library_path=str(solver.parent / 'lib'), geometry_license='STS-Tooth3D CC BY4.0 DOI10.5281/zenodo.10597292; tooth-derived model fixture, repeated same tooth, not measured bridge', source_grid_sha256=sha(X1B / 'inputs/geometry/D1_grid.json')))
    print('Pinned input MB', sum((p['bytes'] for p in manifest if p['local'])) / 1000000.0)
if __name__ == '__main__':
    run()
