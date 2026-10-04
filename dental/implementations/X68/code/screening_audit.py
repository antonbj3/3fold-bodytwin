"""Bounded screening accounting; deferred papers are not declared negative."""
from dental_release.paths import expand as _release_expand
import json, pathlib, hashlib, xml.etree.ElementTree as ET
ROOT = pathlib.Path(__file__).resolve().parents[1]
CORPUS = pathlib.Path(_release_expand('@DENTAL_CORPUS_ROOT@/europepmc/fulltext'))
REASONS = {'PMC10314359': 'Torque outcome tables are images, no unique final-bore numerical profile recovered from local XML. Four protocols include tap/cortical/conical operations. Deferred extraction, not evidence of absent effect.', 'PMC10531925': 'Measured nominal drill and manufacturer contact areas are tabulated, but arm torque outcomes are plotted; 80 Ncm motor ceiling censors insertion. Not included in frozen complete-case regression.', 'PMC10366104': 'Measured torque is DRILLING torque during osteotomy, not implant insertion torque. Wrong compared quantity.', 'PMC10245014': 'Tool mechanics and drill torque, not paired measured implant insertion torque. Wrong compared quantity.', 'PMC10058042': 'Prescribed 30/45/55 Ncm insertion conditions for surface-damage experiment; not independent response observations with measured density/bore profile.', 'PMC10253638': 'Porcine tibia, ordinal density III-IV, no measured absolute density or numerical final-bore profile. Individual IT quantized in 5Ncm steps; different population than prereg PU inference.', 'PMC10329581': 'Human RCT; ordinal bone quality and tabulated drilling prescriptions; torque is <40/>=40 categories with source unit N. No matched continuous IT/cortex/geometry rows. Do not convert N into Ncm silently.', 'PMC10214678': 'Post-extraction socket/apical-depth experiment: evolving regional contact is the manipulated variable. Nominal final bore does not describe the socket/contact state. Not silently pooled with full-length bore arms.', 'PMC10605303': 'Porcine ribs; density inclusion criterion <750 reported HU, not paired calibrated density; numerical pitch1.35mm and regional thread heights0.4/0.5mm ARE reported. Actual final osteotomy profile/cortex thickness absent from quantitative tables. Cannot fit thread covariate to R1 PU complete cases.'}

def main():
    manifest = json.loads((ROOT / 'raw/source_manifest.json').read_text())
    screened = json.loads((ROOT / 'raw/screening_titles.json').read_text())
    records = []
    for r in screened:
        pmc = r['pmc']
        v = dict(pmc=pmc, title=r['title'], locator=r['doi'])
        if pmc in manifest:
            v.update(status='EXTRACTED_SOURCE', reason='Typed arm or separately retained endpoint evidence; see raw/arms.json and raw/other_evidence.json')
        elif pmc in REASONS:
            v.update(status='EXCLUDED_OR_DEFERRED_COMPLETE_CASE', reason=REASONS[pmc])
        elif (ROOT / 'raw/table_previews' / f'{pmc}.txt').exists():
            v.update(status='DEFERRED_TABLE_PREVIEW_ONLY', reason='Preview produced; paired full input contract not fully audited. Not a negative measurement.')
        else:
            v.update(status='DEFERRED_METADATA_ONLY', reason='Keyword/title flags only; not a systematic fulltext exclusion.')
        records.append(v)
    additional = []
    for (pmc, reason) in REASONS.items():
        p = CORPUS / (pmc + '.xml')
        x = ET.parse(p)
        additional.append(dict(pmc=pmc, locator=str(p), sha256=hashlib.sha256(p.read_bytes()).hexdigest(), doi=next((a.text for a in x.findall('.//article-id') if a.get('pub-id-type') == 'doi')), reason=reason))
    counts = {k: sum((r['status'] == k for r in records)) for k in sorted({r['status'] for r in records})}
    arms = json.loads((ROOT / 'raw/extraction_summary.json').read_text())
    out = dict(keyword_candidates=200, candidate_source_status_counts=counts, source_status_fractions={k: v / 200 for (k, v) in counts.items()}, records=records, additional_source_checks=additional, legacy_small_XML_source_verifications=5, arm_dropout_fraction=arms['rejected_arms'] / arms['candidate_arms'], complete_case_arm_counts=arms, claim='Bounded targeted extraction, not complete retrieval or absence proof for 200 candidate articles', measured_thread_geometry_in_separately_retained_ex_vivo_source='PMC10605303 §2.2; one study family, not three fully paired families')
    (ROOT / 'SCREENING_AUDIT.json').write_text(json.dumps(out, indent=2, ensure_ascii=False) + '\n')
    print(json.dumps(dict(source_status_counts=counts, arm_dropout_fraction=out['arm_dropout_fraction']), indent=2))
if __name__ == '__main__':
    main()
