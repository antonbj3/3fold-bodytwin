"""Lossless report-field extraction. Lexical numeric mentions stay UNVERIFIED.

Only compared_quantity is parsed, never our_value or a derived The swarm result.
"""
from pathlib import Path
import json, re, collections
P = Path(__file__).resolve().parent
rows = [json.loads(l) for l in (P / 'ALL_EXTERNAL_REFERENTS.jsonl').read_text().splitlines()]
verified = {x['candidate_id']: x for x in map(json.loads, (P / 'VERIFICATIONS.jsonl').read_text().splitlines())}
pattern = re.compile('(?<![\\w.])([-+−]?\\d+(?:[.,]\\d+)?(?:\\s*(?:±|\\+/-|\\+\\s*/\\s*-)\\s*\\d+(?:[.,]\\d+)?)?)\\s*(GPa|MPa|Ncm|N\\s*cm|mm\\^?2|mm²|mm2|mm|µm|μm|um|ue|microstrain|°C|degC|K|rpm|Hz|%|percent|W/mK|N)(?![\\w])', re.I)
out = []
for r in rows:
    x = r['external_referent']
    q = str(x.get('compared_quantity', ''))
    loc = str(x.get('locator', ''))
    mentions = [{'printed_claim': m.group(0), 'numeric_lexeme': m.group(1), 'unit_lexeme': m.group(2), 'span_in_compared_quantity': [m.start(), m.end()], 'nearby_quantity_context': q[max(0, m.start() - 80):m.end() + 80], 'status': 'UNVERIFIED_REPORT_PROSE; may be hypothesis, model or measurement'} for m in pattern.finditer(q)]
    identifiers = {'doi_mentions': re.findall('10\\.\\d{4,9}/[^\\s,;\\)\\]}]+', loc), 'PMID_mentions': re.findall('PMID\\s*[:=]?\\s*(\\d+)', loc, re.I), 'PMC_mentions': re.findall('PMC\\d+', loc), 'urls': re.findall('https?://[^\\s\\)\\]}]+', loc)}
    row = {'candidate_id': r['candidate_id'], 'path': r['path'], 'source_sha256': r['source_sha256'], 'rank': r['rank'], 'rank_score': r['rank_score'], 'target_id': r['target_id'], 'source': x.get('source', x.get('citation', 'UNKNOWN_AT_REPORT_LEVEL')), 'source_identifiers': identifiers, 'locator': x.get('locator'), 'compared_quantity': x.get('compared_quantity'), 'external_value': x.get('value', x.get('referent_value')), 'unit': x.get('unit', x.get('units')), 'population': x.get('population'), 'protocol': x.get('protocol'), 'kind': x.get('kind'), 'refutes_us': x.get('refutes_us'), 'lexical_measurement_candidates': mentions, 'original_external_referent': x, 'verification': verified.get(r['candidate_id']), 'status': 'REPORT_EXTRACTION_ONLY; numeric admission requires independent primary observation', 'missing_metadata_rule': 'UNKNOWN when no explicit source field; full locator/compared prose retained. Lexical mentions do not inherit measurement status.'}
    out.append(row)
(P / 'EXTRACTED_FACIT_CANDIDATES.jsonl').write_text(''.join((json.dumps(x, ensure_ascii=False) + '\n' for x in out)))
(P / 'EXTRACTION_DIAGNOSTICS.json').write_text(json.dumps({'referents': len(out), 'lexical_numeric_unit_mentions': sum((len(x['lexical_measurement_candidates']) for x in out)), 'explicit_field_counts': {k: sum((x[k] is not None for x in out)) for k in ['external_value', 'unit', 'population', 'protocol']}, 'ranking': 'Heuristic: primary identifier20, table/figure/abstract locator10, unit mention8, refutes flag5, populated our_value3. Inspection priority, not optimality proof.', 'resolution_level': 'PHENOMENOLOGICAL', 'missing_population_protocol': 'Never inferred from our_value'}, indent=2) + '\n')
print('lossless extracted referents', len(out))
