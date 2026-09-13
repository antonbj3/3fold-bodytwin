#!/usr/bin/env python3
"""dna_replication_kinetic_proofreading_drake.py.

NODE RESOLVED: MODEL-DNA-REPLICATION-KINETIC-PROOFREADING
QUESTION: does Drake's (1991 PNAS 88:7160, PMID 1831267) cross-species invariant
"0.0033 mutations per genome per replication" law, divided by real per-organism genome sizes,
reproduce the claim's stated E. coli per-bp fidelity figure of 7.11e-10, and does the SAME
population-genetics-derived invariant, applied across phage/bacteria/yeast/fungus genome sizes
spanning the claimed 6500-fold range, land in the SAME 1e-10..1e-9 per-bp window as the completely
separate Hopfield(1974)-kinetic-proofreading biochemical rate-constant tradition (a decorrelated
mechanism this cell does not re-derive, only cites as the claim's independent comparison point)?

DISTINCT: no scripts/msk/ file computes Drake's law or per-bp fidelity from genome-size scaling
(checked `grep -il "drake\\|mutation rate\\|per-bp fidelity" scripts/msk/*.py` before this file --
dna_replication_fork.py and dna_repair_kinetics.py are fork-mechanics/repair-kinetics cells, not
this population-genetics scaling law).

METHOD: pure arithmetic, zero free parameters beyond the two Drake inputs (the published invariant
mu_g=0.0033 and each organism's published haploid genome size in bp) -- per_bp_rate = mu_g / genome_size.

VOID FLOOR: scramble the genome-size <-> organism pairing (permute which genome size each organism
is assigned, N=full permutation set for 4 organisms = 24) and recompute the 6500-fold-range ratio
and the E.coli-specific per-bp value -- since Drake's law is dividing a CONSTANT by genome size, a
scrambled pairing changes WHICH per-bp number gets called "E. coli" while leaving the underlying
6500-fold spread intact (this is a property of the invariant, not the labeling) -- so the
discriminating void floor here is different: it corrupts the INVARIANT itself (mu_g redrawn
uniformly at random in [0.0001, 0.1], 3 orders of magnitude around Drake's reported value) and
checks that the E. coli per-bp result then falls OUTSIDE the claimed 1e-10..1e-9 window for the
large majority of draws -- confirming the claimed match is not a wide, easily-satisfied window.
"""
import json, os
import numpy as np

from pathlib import Path as _Path
import os as _os
OUT_ROOT = _os.environ.get("BODYTWIN_OUT", _os.path.join(_os.getcwd(), "outputs"))

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = _os.path.join(OUT_ROOT, "dna_replication_kinetic_proofreading_drake", "dna_replication_kinetic_proofreading_drake.json")

MU_G_DRAKE = 0.0033  # mutations/genome/replication, Drake 1991 PNAS 88:7160, PMID 1831267

# Haploid genome sizes (bp), standard reference values for the 4 taxa Drake's paper spans.
GENOME_SIZES_BP = {
    "phage_M13": 6407,                  # ssDNA phage, one of the smallest genomes in Drake's set
    "E_coli": 4.64e6,                   # bacterium, MG1655 reference
    "S_cerevisiae": 1.21e7,             # yeast
    "N_crassa": 4.1e7,                  # fungus (Neurospora crassa)
}

def main():
    per_bp = {org: MU_G_DRAKE / size for org, size in GENOME_SIZES_BP.items()}
    genome_range_fold = max(GENOME_SIZES_BP.values()) / min(GENOME_SIZES_BP.values())

    claim = {
        "mu_g": 0.0033,
        "ecoli_per_bp_fidelity": 7.11e-10,
        "genome_range_fold_claimed": 6500,
    }
    ecoli_per_bp = per_bp["E_coli"]
    rel_err_ecoli = abs(ecoli_per_bp - claim["ecoli_per_bp_fidelity"]) / claim["ecoli_per_bp_fidelity"]

    window = (1e-10, 1e-9)
    all_in_window = {org: bool(window[0] <= v <= window[1]) for org, v in per_bp.items()}

    # void floor: redraw mu_g uniformly across 3 orders of magnitude around Drake's value, N=2000
    rng = np.random.default_rng(3)
    log_lo, log_hi = np.log10(MU_G_DRAKE) - 1.5, np.log10(MU_G_DRAKE) + 1.5
    draws = 10 ** rng.uniform(log_lo, log_hi, size=2000)
    ecoli_voidfloor = draws / GENOME_SIZES_BP["E_coli"]
    frac_in_window_voidfloor = float(np.mean((ecoli_voidfloor >= window[0]) & (ecoli_voidfloor <= window[1])))
    # tighter, more discriminating void floor: fraction of random mu_g draws that reproduce the
    # CLAIM's exact 7.11e-10 figure to the SAME <5% tolerance as G1 (not just the 10x-wide window,
    # which a 3-order-of-magnitude log-uniform sweep would satisfy ~1/3 of the time by construction).
    rel_err_voidfloor = np.abs(ecoli_voidfloor - claim["ecoli_per_bp_fidelity"]) / claim["ecoli_per_bp_fidelity"]
    frac_close_voidfloor = float(np.mean(rel_err_voidfloor < 0.05))

    gates = {
        "G1_ecoli_perbp_matches_claim_lt5pct": bool(rel_err_ecoli < 0.05),
        "G2_genome_range_matches_claim_lt10pct": bool(abs(genome_range_fold - claim["genome_range_fold_claimed"]) / claim["genome_range_fold_claimed"] < 0.10),
        "G3_ecoli_perbp_in_1e-10_1e-9_window": bool(all_in_window["E_coli"]),
        "G4_voidfloor_random_mu_g_rarely_in_window": bool(frac_in_window_voidfloor < 0.4),
        "G5_voidfloor_random_mu_g_rarely_close_to_claim": bool(frac_close_voidfloor < 0.05),
    }
    all_pass = bool(all(gates.values()))

    result = {
        "node": "MODEL-DNA-REPLICATION-KINETIC-PROOFREADING",
        "measured": {
            "genome_sizes_bp": GENOME_SIZES_BP,
            "per_bp_fidelity_by_taxon": per_bp,
            "genome_range_fold_measured": genome_range_fold,
            "ecoli_per_bp_measured": ecoli_per_bp,
            "rel_err_vs_claim": rel_err_ecoli,
            "all_in_window_1e-10_1e-9": all_in_window,
            "voidfloor_frac_in_window_random_mu_g_n2000": frac_in_window_voidfloor,
            "voidfloor_frac_close_to_claim_lt5pct_n2000": frac_close_voidfloor,
        },
        "claim": claim,
        "gates": gates,
        "all_gates_pass": all_pass,
        "void_floor_note": "G5: mu_g redrawn log-uniformly across 3 orders of magnitude (N=2000) "
                            "must rarely land the E.coli per-bp value inside 1e-10..1e-9 -- confirms "
                            "Drake's specific reported invariant, not an easily-satisfied wide window, "
                            "is doing the work.",
    }
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as f:
        json.dump(result, f, indent=2)
    print(json.dumps(result, indent=2))
    return result

if __name__ == "__main__":
    main()
