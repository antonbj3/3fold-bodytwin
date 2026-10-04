from common import *
import importlib.metadata
import construction, elliptic, loop_crown, height_crown, offset_crown, score_assembly, export_prep, path_certificate, continuous_film

def run():
    versions = {x: importlib.metadata.version(x) for x in ['numpy', 'scipy', 'trimesh', 'matplotlib', 'threadpoolctl']}
    modules = sorted({str(Path(m.__file__).resolve()) for m in list(sys.modules.values()) if getattr(m, '__file__', None) and ('/dental/' in m.__file__ or '/dental_sol_night/' in m.__file__) and Path(m.__file__).is_file()})
    binaries = [D / 'exact_mesh_clearance', D / 'exact_projected_loop', D / 'exact_projection_inside', D / 'exact_gap_upper', D / 'exact_surface', D / 'exact_clearance', PARENT_D / 'PROOF_LANE_FULL_CROWN_R6/intersections']
    dump(R / 'DEPENDENCIES.json', dict(python=sys.executable, python_version=sys.version, versions=versions, imported_local_modules=[dict(path=p, sha256=sha(p)) for p in modules], binaries=[dict(path=str(p.resolve()), sha256=sha(p)) for p in binaries], scope='Pinned local replay environment. Parent code/data referenced read-only; no copied source trees. No guarantee of portability to another software stack.'))
    files = set()
    for pat in ['PREREG*', 'DECOMPOSITION*', 'FROZEN_PREDICTIONS*', 'RESULTS*.json', 'EXPORTS.json', 'DEPENDENCIES.json', 'run_all.sh', 'build*.sh']:
        files.update((p for p in R.glob(pat) if p.is_file()))
    files.update((p for p in (R / 'code').rglob('*') if p.is_file() and p.suffix in ['.py', '.cpp']))
    files.update((p for p in (R / 'exports').glob('*.json') if p.is_file()))
    files.update((p for p in (R / 'raw').glob('*.json') if p.name.endswith(('_CERT.json', '_PATH.json', '_FILM.json', '_FILM_UPPER.json'))))
    files.update([R / 'raw/FRESH_RECOMPUTE.json', R / 'raw/EXACT_CARRIER_AUDIT.json', R / 'PREPARATION_F_CERTIFICATES.json'])
    files.update((p for p in D.rglob('*') if p.is_file() and (not p.is_symlink()) and ('recompute_' not in str(p.relative_to(D))) and ('_REPLAY' not in p.name)))
    files.update((Path(p) for p in modules))
    files.update(binaries)
    for rec in inputs():
        files.update([Path(rec['private_path']), Path(rec['public_path']), Path(rec['outer']['mesh_path'])])
        (_, _, info) = native(rec)
        files.add(Path(info['path']))
    for p in (R / 'raw').glob('*_CERT.json'):
        c = read(p)
        for k in ['support_path', 'query_path']:
            if k in c:
                files.add(Path(c[k]))
    rows = [dict(path=str(p), bytes=p.stat().st_size, sha256=sha(p)) for p in sorted(files)]
    dump(R / 'MANIFEST.json', dict(created_utc=now(), files=rows, scope='Immutable evidence and local code/data dependencies. Live state, final narrative, resource receipts, graph feedback and replay output are excluded to avoid circular hashes.'))
    print('Manifest files', len(rows), 'bytes', sum((x['bytes'] for x in rows)))
if __name__ == '__main__':
    run()
