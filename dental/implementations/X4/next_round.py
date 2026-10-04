import json
from freeze import ROOT, freeze, write, state, now

def prepare_r2():
    r = json.loads((ROOT / 'RESULTS_R1.json').read_text())
    write(ROOT / 'HANDOFF_R1.md.json', dict(round='R1', outcome=r['outcome'], gates=r['gates'], methods=r['methods'], source_semantics='Post virtual, Original defective', next_construction='R2 expert inferior rail plus finite-thickness cut constraint inside search', updated_utc=now()))
    (ROOT / 'HANDOFF_R1.md').write_text('# R1 decided\n\n' + r['outcome'] + '; ' + str(r['n_available']) + '/118 available.\n\nDP vs independent exhaustive chord solver ties at numerical precision. R1 absolute 3mm and complete-cohort gates remain frozen. See RESULTS_R1.json and its per-case/raw hashes. Post is virtual; Original is defective.\n\nNext changes the information contract and geometric operation: supply the expert inferior rail, then constrain finite cylinder miter extents inside search. It is assisted CAD planning, not blind completion.\n')
    freeze('R2', dict(prediction=['Pre', 'expert Post inferior rail and common cranial alignment'], evaluation=['remaining Post surface', 'Original defective anatomy'], blind_reconstruction=False), 'expert inferior rail -> thickness-constrained segment partition -> mitered graft; same target curve given to actual upstream OsteoOpt geometric simplifier', dict(median_primary_p95_mm_max=3.0, all_cases_available=True, paired_gain_vs_R1_mm_min=1.0, solver_max_absolute_disagreement_mm=1e-08, all_miters_non_crossing=True, every_gate_fault_injection=True))
    p = ROOT / 'PREREG_R2.json'
    reg = json.loads(p.read_text())
    if 'additional_contract' not in reg:
        reg['additional_contract'] = dict(candidate='all partitions with 1-3 segments, minimum length and finite cylinder miter extent constraints', strongest_control='independent cross-product chord cost matrix + independently enumerated same feasible partitions', open_competitor='unmodified pinned OsteoOpt PolylineSimplifier.java, same curve, min length10mm and segment budget; full BO/ArtiSynth NOT_RUN', forbidden_claim='gain in blind Pre-only anatomical reconstruction', paired_R1_case_ids=[c['case'] for c in r['cases'] if c['numerically_available']], validation_quantity='full expert added surface, not only fitted inferior rail', information_acquisition_cost='complete Post component extraction + registration + inferior curve extraction charged; not a free landmark')
        write(p, reg)
        from freeze import sha
        (ROOT / 'PREREG_R2.sha256').write_text(sha(p) + '  ' + p.name + '\n')
    state('R2_preregistered', 'execute expert-assisted planning with equal-information external geometric competitor', r['outcome'])
if __name__ == '__main__':
    prepare_r2()
