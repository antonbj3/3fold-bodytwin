"""Read scoped review disclosures; the signed-off release anchor remains trusted."""
import re

def code_eligible(record):
    return (record.get('classification') in ("ACCEPTED",'HOLDS_WITH_CORRECTION')
            and record.get('action') in ('ADD','ALREADY_PRESENT')
            and all(isinstance(record.get(k),str) and re.fullmatch('[0-9a-f]{64}',record[k])
                    for k in ('result_sha256','review_sha256')))

def count_sufficiency(records):
    """Identical entry count can conceal one inadmissible code scope."""
    import copy
    a=copy.deepcopy(records);b=copy.deepcopy(records)
    b[-1]['classification']='UNKNOWN'
    summaries=[len(a),len(b)]
    eligible=[sum(map(code_eligible,a)),sum(map(code_eligible,b))]
    return dict(summary_entry_counts=summaries,identity_error=summaries[0]-summaries[1],
                downstream_eligible_counts=eligible,downstream_difference=eligible[0]-eligible[1],
                minimum_extension='Retain per-entry scope, acceptance class and exact result/review binding; counts alone do not certify coverage',
                resolution='Software ledger; no physical evidence',timescale='SIMULTANEOUS')
