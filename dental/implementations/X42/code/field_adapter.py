"""A concrete mesh port for the field engine; does not initialize any GPU.

All inputs in the benchmark-local site frame, mm. V3's roof constraints certify
the graph patches, not the enclosing vertical rim/full clinical restoration.
"""
import numpy as np
from generator import conditional, branch
from legacy.geometry import shell, export_stl

def generate_mesh(task, shape_model, contact_model=None, beta=0.0):
    d = conditional(task, shape_model) if contact_model is None else branch(task, shape_model, contact_model, beta)
    if d['status'] != 'DESIGN':
        return d
    outer = np.asarray(d['outer_vertices'])
    inner = np.asarray(d['inner_vertices'])
    (v, f) = shell(task['xy'], outer[:, 2], inner[:, 2], task['faces'])
    return {'status': 'DESIGN', 'vertices': v, 'faces': f, 'units': 'mm', 'frame': task['frame'], 'source_resolution': 'PER_POINT', 'timescale_to_manufacturing': 'HANDOVER', 'roof_contract': 'Same V3 wall/film/continuous-obstacle constraints; vertical rim and full manufacturing UNKNOWN', 'outer_vertices': outer, 'inner_vertices': inner}

def export_research_mesh(port, path):
    if port['status'] != 'DESIGN':
        raise ValueError('refused design has no export')
    export_stl(path, port['vertices'], port['faces'])

def prepare_field(port, pitch, lo, *, winding_library, ptx, edt_library):
    """Native engine port after explicitly supplying its libraries/hardware lease.

Caller adds dental/references/field_engine/src/field_engine to sys.path.
Return handle belongs to this thread; caller must close it. This GPU-backed
bridge is documented, not exercised by the CPU-only X42 experiment.
 """
    import trimesh
    mesh = trimesh.Trimesh(port['vertices'], port['faces'], process=False)
    if not mesh.is_watertight or not mesh.is_winding_consistent:
        raise ValueError('field port refuses open or inconsistent dental shell')
    from mesh_field_rt_v1 import PreparedMeshField
    return PreparedMeshField(port['vertices'], port['faces'], pitch, lo, winding_library=winding_library, ptx=ptx, edt_library=edt_library)

def cpu_field(port, pitch, lo, method='raypar_vindning'):
    import trimesh
    mesh = trimesh.Trimesh(port['vertices'], port['faces'], process=False)
    if not mesh.is_watertight or not mesh.is_winding_consistent:
        raise ValueError('field port refuses open or inconsistent dental shell')
    import faltkarna_v1_mesh_to_sdf as backend
    diag = {}
    result = backend.surface_raster_and_flood(port['vertices'], port['faces'], pitch, lo, metod=method, diag_ut=diag, ytkorrektion=True)
    return (result, diag)
