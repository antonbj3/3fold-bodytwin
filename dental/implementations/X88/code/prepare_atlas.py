"""Export anatomy context only from the corrected X12 atlas; no disease join."""
import csv
import hashlib
import json
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT.parent / 'LANE_X12_PULPY3D/raw/R8_teeth.csv'

def main():
    rows = list(csv.DictReader(SOURCE.open()))
    valid = [r for r in rows if r['domain_truncated'] == 'False']
    eligible = [r for r in valid if int(r['fdi']) in (36, 37, 46, 47)]
    exported = []
    for r in eligible:
        exported.append(dict(atlas_case=r['case'], canonical_fdi=int(r['fdi']), source_pulpy_fdi=int(r['source_pulpy_fdi']), semantic_map=r['semantic_map'], voxel_mm=float(r['voxel_mm']), horn_boundary_min_mm=float(r['horn_boundary_min_mm']), geometry_resolution='PER_TOOTH', CBCT_periapical_disease_category=None, conversion_probability=None, join_kind='ANATOMY_CONTEXT_ONLY', same_subject_as_RCT=False))
    result = dict(source=str(SOURCE), source_sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(), source_rows=len(rows), domain_valid_rows=len(valid), exported_molar_rows=len(exported), unique_molar_cases=len(set((r['atlas_case'] for r in exported))), excluded=dict(truncated=len(rows) - len(valid), different_FDI=len(valid) - len(eligible)), clinical_join_rejected=len(exported), clinical_join_rejection_fraction=1.0, rejection_reason='Atlas anatomy has no linked preop periapical category, haemostasis or treatment outcome; no shared subject IDs.', earlier_brief_355_cases='Outdated. Corrected atlas accepts 283 cases, with 282 nontruncated cases; use R8 and preserve semantic_map.', license=dict(ToothFairy2='CC BY-SA 4.0', Pulpy3D_pulp_extension='UNKNOWN explicit license in reviewed sources', provenance='LANE_X12_PULPY3D/sources/LICENSE_EVIDENCE.json; no imaging copied/distributed'), rows=exported)
    (ROOT / 'raw/ATLAS_CONTEXT.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for (k, v) in result.items() if k != 'rows'}, indent=2))
if __name__ == '__main__':
    main()
