"""Preserved unsuccessful positional sidecar construction."""
import datetime, json
from pathlib import Path
import numpy as np
import trimesh
R = Path(__file__).resolve().parents[1]

def positional_transfer(old_labels, new_faces):
    return old_labels.copy()

def main():
    m = trimesh.load_mesh(R / 'inputs/exports/D1/crown.stl', process=False)
    tri = np.asarray(m.vertices)[m.faces]
    label = np.where(tri.mean(1)[:, 2] >= np.median(tri.mean(1)[:, 2]), 'transport_test_A', 'transport_test_B')
    order = np.arange(len(tri))[::-1]
    imported = trimesh.Trimesh(m.vertices, m.faces[order], process=False)
    same = bool(np.allclose(np.sort(m.triangles_center, axis=0), np.sort(imported.triangles_center, axis=0), atol=0, rtol=0))
    wrong = positional_transfer(label, imported.faces)
    oracle = label[order]
    x = dict(round='R1', claim_type='capability', outcome='FAIL', faces=len(tri), geometry_unchanged=same, volume_relative_difference=abs(m.volume - imported.volume) / abs(m.volume), wrong_region_facets=int((wrong != oracle).sum()), region_disagreement_fraction=float(np.mean(wrong != oracle)), resolution='PER_SURFACE_REGION', region_names_are_synthetic=True, cause='Face indices are not stable under ordinary facet-order conversion', next_operation='Replace index identity by bijective oriented physical-triangle correspondence')
    (R / 'rounds/R1.json').write_text(json.dumps(x, indent=2) + '\n')
    (R / 'HANDOFF_R1.md').write_text(f"R1 executed: geometry unchanged={same}, region mismatch={x['region_disagreement_fraction']:.6f}; FAIL. No threshold changed. Next construction uses oriented physical triangles to recover original regions and refuses ambiguity. Region probe is synthetic; independent geometry import is trimesh.\n")
    (R / 'CURRENT_WORK_STATE.json').write_text(json.dumps(dict(lane='X38-export-gate', phase='R1_DECIDED_R2_BUILDING', updated_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(), latest_gate=x, next_operation='Run geometry-bound R2 successor and statistical consumers'), indent=2) + '\n')
    print(json.dumps(x))
if __name__ == '__main__':
    main()
