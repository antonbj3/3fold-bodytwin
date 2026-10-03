"Geometry exports to JSON, OpenSim and signed-distance fields."
from __future__ import annotations

import json
import re
import subprocess
from pathlib import Path

import numpy as np

from . import identity as I

HERE = Path(__file__).resolve().parent


def point_ids(res):
    """Map of stable ids -> (position mm, 3x3 cov) for every point entity of the instance."""
    reg = res.instance.reg
    part = reg.part
    vals = res.instance.values()
    out = {}
    for e in reg.entities.values():
        if e.kind in ('landmark', 'attachment', 'derived'):
            key = e.name if e.kind == 'landmark' else ('att:' + e.name if e.kind == 'attachment' else None)
            if e.kind == 'derived':
                key = 'HC6' if e.name == 'HJC' else None
            cov = None
            if res.points is not None and key in res.points:
                cov = res.points[key][1]
            out[e.id] = (np.asarray(vals[e.base], float), cov)
    return out


def sanitize(i):
    return 'id_' + re.sub(r'[^A-Za-z0-9_]', '_', i)


# ------------------------------------------------------------------ JSON
def to_json(res, path):
    pts = point_ids(res)
    vals = res.instance.values()
    reg = res.instance.reg
    d = dict(schema='n7a.geomgr.instance/0.1', subject=res.subject, frame=res.frame, units=res.instance.units,
             side=res.instance.side, registry_manifest_sha256=reg.manifest_hash(), registry=reg.manifest(),
             accepted=res.accepted, certificate=res.certificate, kappa=res.kappa, op_log=res.instance.log,
             points={k: dict(xyz_mm=v[0].tolist(), cov_mm2=None if v[1] is None else np.asarray(v[1]).tolist())
                     for k, v in pts.items()},
             regions={e.id: dict(n_vertices=len(e.indices), centroid_mm=vals[e.base]['centroid'].tolist(),
                                 vertex_ids_sha256=I._h(list(map(int, e.indices))))
                      for e in reg.entities.values() if e.kind == 'region'},
             frames={e.id: dict(origin_mm=vals[e.base]['origin'].tolist(), axes_rows=vals[e.base]['axes'].tolist())
                     for e in reg.entities.values() if e.kind == 'frame'},
             constructs={e.id: dict(xyz_mm=np.asarray(vals[e.base]).tolist(), rule=e.rule)
                         for e in reg.entities.values() if e.kind == 'construct'},
             wraps={e.id: dict(centre_mm=vals[e.base]['c'].tolist(), axis=vals[e.base]['a'].tolist(),
                               radius_mm=vals[e.base]['R'], s0_mm=vals[e.base]['s0'].tolist(), length_mm=vals[e.base]['L'],
                               rule=e.rule) for e in reg.entities.values() if e.kind == 'wrap'},
             fields={k: dict(kind=v['kind'], units=v['units'], source=v['source'],
                             n_finite=int(np.isfinite(v['values']).sum()),
                             values_sha256=I._h([None if not np.isfinite(x) else round(float(x), 9) for x in v['values']]))
                     for k, v in getattr(res.instance, 'fields', {}).items()},
             sigma_mm_per_vertex=dict(median=float(np.median(res.sigma_mm)), p95=float(np.percentile(res.sigma_mm, 95)))
             if res.sigma_mm is not None else None)
    Path(path).write_text(json.dumps(d, indent=0, default=float))
    return d


def ids_in_json(path):
    d = json.loads(Path(path).read_text())
    return set(d['points']) | set(d['regions']) | set(d['frames']) | set(d.get('constructs', {})) | set(d.get('wraps', {}))


# ------------------------------------------------------------------ OpenSim
def _isb(res):
    vals = res.instance.values()
    f = vals[f'{res.instance.reg.part}/frame/ISB']
    return f['origin'], f['axes']


def to_opensim(res, outdir, name=None, density_kg_m3=1900.0):
    """Minimal OpenSim 4 model: ground + body <part> welded at the ISB frame; mesh + one marker per point id.
    Mass/inertia from the closed mesh at a declared homogeneous bone density (not segment mass)."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    name = name or f'{res.subject}_{res.instance.reg.part}'
    o, Ax = _isb(res)
    V = (res.instance.V - o) @ Ax.T / 1000.0            # m in ISB frame
    F = res.instance.F
    obj = outdir / f'{name}.obj'
    with open(obj, 'w') as f:
        f.write('# geomgr N7a export, metres, ISB femur frame\n')
        np.savetxt(f, V, fmt='v %.9f %.9f %.9f')
        np.savetxt(f, F + 1, fmt='f %d %d %d')
    # mass properties (divergence theorem)
    t = V[F]
    vol6 = np.einsum('ij,ij->i', t[:, 0], np.cross(t[:, 1], t[:, 2]))
    vol = vol6.sum() / 6
    com = (vol6[:, None] * t.sum(1)).sum(0) / (24 * vol)
    mass = density_kg_m3 * vol
    part = res.instance.reg.part
    markers = []
    idmap = {}
    for i, (p, _) in sorted(point_ids(res).items()):
        q = (p - o) @ Ax.T / 1000.0
        nm = sanitize(i)
        idmap[nm] = i
        markers.append(f'''        <Marker name="{nm}">
          <socket_parent_frame>/bodyset/{part}</socket_parent_frame>
          <location>{q[0]:.12g} {q[1]:.12g} {q[2]:.12g}</location>
          <fixed>true</fixed>
        </Marker>''')
    # inertia about COM: coarse (point-mass on vertices scaled to mass) -> declared approximate
    Vc = V - com
    I3 = (mass / len(V)) * ((Vc ** 2).sum() * np.eye(3) - Vc.T @ Vc)
    osim = f'''<?xml version="1.0" encoding="UTF-8" ?>
<OpenSimDocument Version="40000">
  <Model name="{name}">
    <credits>geomgr N7a (internal). Body frame = ISB femur frame of the instance. Mass = {density_kg_m3} kg/m3 x mesh volume (bone only, not segment mass); inertia approximate.</credits>
    <length_units>meters</length_units>
    <force_units>N</force_units>
    <gravity>0 -9.8066500000000005 0</gravity>
    <BodySet name="bodyset">
      <objects>
        <Body name="{part}">
          <attached_geometry>
            <Mesh name="{part}_geom">
              <mesh_file>{obj.name}</mesh_file>
            </Mesh>
          </attached_geometry>
          <mass>{mass:.12g}</mass>
          <mass_center>{com[0]:.12g} {com[1]:.12g} {com[2]:.12g}</mass_center>
          <inertia>{I3[0,0]:.9g} {I3[1,1]:.9g} {I3[2,2]:.9g} {I3[0,1]:.9g} {I3[0,2]:.9g} {I3[1,2]:.9g}</inertia>
        </Body>
      </objects>
    </BodySet>
    <JointSet name="jointset">
      <objects>
        <WeldJoint name="ground_{part}">
          <socket_parent_frame>ground_offset</socket_parent_frame>
          <socket_child_frame>{part}_offset</socket_child_frame>
          <frames>
            <PhysicalOffsetFrame name="ground_offset">
              <socket_parent>/ground</socket_parent>
              <translation>0 0 0</translation>
              <orientation>0 0 0</orientation>
            </PhysicalOffsetFrame>
            <PhysicalOffsetFrame name="{part}_offset">
              <socket_parent>/bodyset/{part}</socket_parent>
              <translation>0 0 0</translation>
              <orientation>0 0 0</orientation>
            </PhysicalOffsetFrame>
          </frames>
        </WeldJoint>
      </objects>
    </JointSet>
    <MarkerSet name="markerset">
      <objects>
{chr(10).join(markers)}
      </objects>
    </MarkerSet>
  </Model>
</OpenSimDocument>
'''
    path = outdir / f'{name}.osim'
    path.write_text(osim)
    (outdir / f'{name}.idmap.json').write_text(json.dumps(idmap, indent=0))
    expected = {sanitize(i): ((p - o) @ Ax.T / 1000.0).tolist() for i, (p, _) in point_ids(res).items()}
    return dict(osim=str(path), obj=str(obj), n_markers=len(markers), mass_kg=float(mass), expected=expected)


OSIM_CHECK = r'''
import json, sys, opensim as osim
out = {}
for p in sys.argv[1:]:
    try:
        m = osim.Model(p)
        s = m.initSystem()
        ms = m.getMarkerSet()
        d = {}
        for i in range(ms.getSize()):
            mk = ms.get(i)
            loc = mk.get_location()
            d[mk.getName()] = [loc.get(0), loc.get(1), loc.get(2)]
        b = m.getBodySet().get(0)
        g = b.get_attached_geometry(0)
        out[p] = dict(ok=True, n_markers=ms.getSize(), markers=d, body=b.getName(), mass=b.getMass(),
                      mesh=osim.Mesh.safeDownCast(g).get_mesh_file(), n_bodies=m.getBodySet().getSize())
    except Exception as e:
        out[p] = dict(ok=False, error=str(e))
print(json.dumps(out))
'''


def check_opensim(exports):
    """Load each .osim in OpenSim 4.6 (one subprocess for all) and compare marker locations."""
    from ._local import MSK_PY
    paths = [e['osim'] for e in exports]
    r = subprocess.run(['nice', '-n', '19', str(MSK_PY), '-c', OSIM_CHECK, *paths], capture_output=True, text=True,
                       env={'OMP_NUM_THREADS': '1', 'PATH': '/usr/bin:/bin'}, timeout=1800)
    if r.returncode != 0:
        return dict(ok=False, stderr=r.stderr[-2000:])
    got = json.loads(r.stdout.strip().splitlines()[-1])
    rows = []
    for e in exports:
        g = got[e['osim']]
        if not g.get('ok'):
            rows.append(dict(osim=e['osim'], ok=False, error=g.get('error')))
            continue
        exp = e['expected']
        lost = sorted(set(exp) - set(g['markers']))
        dev = max(float(np.max(np.abs(np.array(g['markers'][k]) - np.array(v)))) for k, v in exp.items() if k in g['markers'])
        rows.append(dict(osim=e['osim'], ok=not lost and dev <= 1e-9, n_markers=g['n_markers'], lost=lost,
                         max_marker_dev_m=dev, mesh=g['mesh'], mass=g['mass']))
    return dict(ok=all(x['ok'] for x in rows), rows=rows)












# ------------------------------------------------------------------ SDF via the field engine
def to_sdf(items, outdir, pitch=2.0):
    """items: list of (name, V, F). One subprocess (field venv, CPU warp) for all meshes."""
    from ._local import FIELD_PY, STAGING
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    inp = outdir / '_sdf_input.npz'
    np.savez(inp, names=np.array([n for n, _, _ in items]), **{f'V_{k}': V for k, (_, V, _) in enumerate(items)},
             **{f'F_{k}': F for k, (_, _, F) in enumerate(items)})
    r = subprocess.run(['nice', '-n', '19', str(FIELD_PY), str(HERE / '_sdf_worker.py'), str(inp), str(outdir),
                        str(pitch), str(STAGING)], capture_output=True, text=True, timeout=7200,
                       env={'OMP_NUM_THREADS': '2', 'CUDA_VISIBLE_DEVICES': '', 'PATH': '/usr/bin:/bin',
                            'HOME': str(Path.home())})
    inp.unlink(missing_ok=True)
    if r.returncode != 0:
        return dict(ok=False, stderr=r.stderr[-3000:])
    return json.loads((outdir / 'sdf_report.json').read_text())








