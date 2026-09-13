#!/usr/bin/env python3
"""genetic_code_error_minimization.py.

NODE RESOLVED: MODEL-GENETIC-CODE-ERROR-MINIMIZATION
QUESTION: does the standard genetic code's structure (deterministic facts about the fixed
codon table -- no experimental data, no Monte Carlo needed) reproduce the three numbers the node's
claim asserts: (1) wobble-redundancy 68.9% synonymous at codon position 3, (2) eta^2 clustering by
2nd-position nucleotide in the 0.578-0.756 range (using a physicochemical amino-acid distance, here
Grantham distance -- a standard, independently published table, not fit to this task), and (3) a
2-opt local-search re-optimization of "how much better could the code be" landing near the claimed
QAP-certified global optimum 3.489 (claim: 2-opt converges to 3.4886).

DISTINCT from any MSK cell: this is the only cell touching the genetic code table itself (no sibling
stem in scripts/msk/ as of  -- checked via `ls scripts/msk | grep -i codon`).

METHOD: (1)+(2) are exact deterministic computations off the fixed standard codon table + a
published Grantham(1974) distance matrix -- no randomness. (3) is a real 2-opt local-search
optimization over code-table permutations, seeded from the standard code, minimizing
sum_{codon pairs differing by 1 point mutation} (Grantham distance between their two amino acids),
normalized the same way the claim's error-cost metric is defined (mean effect per 1-mutation
neighbor, in units of the STANDARD code's value =1.0, so a fully-random code averages ~higher and
the "better" codes score LOWER cost -- claim's ratio 3.489 is (random-code mean cost)/(best-found
cost); we reproduce that ratio definition explicitly below).

VOID FLOOR: scrambled control 1 = randomly permuted codon-to-amino-acid table (100 draws) for the
position-3 wobble and eta^2 checks -- MUST destroy the wobble asymmetry (position 3 no longer
distinguished from positions 1/2) and collapse eta^2 toward 0 (no structured clustering by any
arbitrarily-labeled "position"). Void floor 2 for the optimization: 2-opt run from a RANDOM start
permutation (not the standard code) must NOT reliably reach a cost within the claimed few-percent
of the standard-code-seeded run -- if it does, the "standard code sits near-optimal" claim is not
distinguishing structure from a generic local-search artifact.
"""
import json, os, random
import numpy as np

from pathlib import Path as _Path
import os as _os
OUT_ROOT = _os.environ.get("BODYTWIN_OUT", _os.path.join(_os.getcwd(), "outputs"))

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = _os.path.join(OUT_ROOT, "genetic_code_error_minimization", "genetic_code_error_minimization.json")

BASES = "TCAG"
CODONS = [a + b + c for a in BASES for b in BASES for c in BASES]

# Standard genetic code (DNA-style codons), '*' = stop.
STD_TABLE = {}
_aas = ("F F L L L L L L I I I M V V V V S S S S P P P P T T T T A A A A "
        "Y Y * * H H Q Q N N K K D D E E C C * W R R R R S S R R G G G G").split()
for cod, aa in zip(CODONS, _aas):
    STD_TABLE[cod] = aa

# Grantham (1974) physicochemical distance matrix (published, amino-acid pairwise; symmetric).
# Subset sufficient to cover all 20 aa + stop treated as max-distance sentinel.
GRANTHAM = {
 ('SER','ARG'):110,('SER','LEU'):145,('SER','PRO'):74,('SER','THR'):58,('SER','ALA'):99,
 ('SER','VAL'):124,('SER','GLY'):56,('SER','ILE'):142,('SER','PHE'):155,('SER','TYR'):144,
 ('SER','CYS'):112,('SER','HIS'):89,('SER','GLN'):68,('SER','ASN'):46,('SER','LYS'):121,
 ('SER','ASP'):65,('SER','GLU'):80,('SER','MET'):135,('SER','TRP'):177,
 ('ARG','LEU'):102,('ARG','PRO'):103,('ARG','THR'):71,('ARG','ALA'):112,('ARG','VAL'):96,
 ('ARG','GLY'):125,('ARG','ILE'):97,('ARG','PHE'):97,('ARG','TYR'):77,('ARG','CYS'):180,
 ('ARG','HIS'):29,('ARG','GLN'):43,('ARG','ASN'):86,('ARG','LYS'):26,('ARG','ASP'):96,
 ('ARG','GLU'):54,('ARG','MET'):91,('ARG','TRP'):101,
 ('LEU','PRO'):98,('LEU','THR'):92,('LEU','ALA'):96,('LEU','VAL'):32,('LEU','GLY'):138,
 ('LEU','ILE'):5,('LEU','PHE'):22,('LEU','TYR'):36,('LEU','CYS'):198,('LEU','HIS'):99,
 ('LEU','GLN'):113,('LEU','ASN'):145,('LEU','LYS'):107,('LEU','ASP'):172,('LEU','GLU'):138,
 ('LEU','MET'):15,('LEU','TRP'):61,
 ('PRO','THR'):38,('PRO','ALA'):27,('PRO','VAL'):68,('PRO','GLY'):42,('PRO','ILE'):95,
 ('PRO','PHE'):114,('PRO','TYR'):110,('PRO','CYS'):169,('PRO','HIS'):77,('PRO','GLN'):76,
 ('PRO','ASN'):91,('PRO','LYS'):103,('PRO','ASP'):108,('PRO','GLU'):93,('PRO','MET'):87,
 ('PRO','TRP'):147,
 ('THR','ALA'):58,('THR','VAL'):69,('THR','GLY'):59,('THR','ILE'):89,('THR','PHE'):103,
 ('THR','TYR'):92,('THR','CYS'):149,('THR','HIS'):47,('THR','GLN'):42,('THR','ASN'):65,
 ('THR','LYS'):78,('THR','ASP'):85,('THR','GLU'):65,('THR','MET'):81,('THR','TRP'):128,
 ('ALA','VAL'):64,('ALA','GLY'):60,('ALA','ILE'):94,('ALA','PHE'):113,('ALA','TYR'):112,
 ('ALA','CYS'):195,('ALA','HIS'):86,('ALA','GLN'):91,('ALA','ASN'):111,('ALA','LYS'):106,
 ('ALA','ASP'):126,('ALA','GLU'):107,('ALA','MET'):84,('ALA','TRP'):148,
 ('VAL','GLY'):109,('VAL','ILE'):29,('VAL','PHE'):50,('VAL','TYR'):55,('VAL','CYS'):192,
 ('VAL','HIS'):84,('VAL','GLN'):96,('VAL','ASN'):133,('VAL','LYS'):97,('VAL','ASP'):152,
 ('VAL','GLU'):121,('VAL','MET'):21,('VAL','TRP'):88,
 ('GLY','ILE'):135,('GLY','PHE'):153,('GLY','TYR'):147,('GLY','CYS'):159,('GLY','HIS'):98,
 ('GLY','GLN'):87,('GLY','ASN'):80,('GLY','LYS'):127,('GLY','ASP'):94,('GLY','GLU'):98,
 ('GLY','MET'):127,('GLY','TRP'):184,
 ('ILE','PHE'):21,('ILE','TYR'):33,('ILE','CYS'):198,('ILE','HIS'):94,('ILE','GLN'):109,
 ('ILE','ASN'):149,('ILE','LYS'):102,('ILE','ASP'):168,('ILE','GLU'):134,('ILE','MET'):10,
 ('ILE','TRP'):61,
 ('PHE','TYR'):22,('PHE','CYS'):205,('PHE','HIS'):100,('PHE','GLN'):116,('PHE','ASN'):158,
 ('PHE','LYS'):102,('PHE','ASP'):177,('PHE','GLU'):140,('PHE','MET'):28,('PHE','TRP'):40,
 ('TYR','CYS'):194,('TYR','HIS'):83,('TYR','GLN'):99,('TYR','ASN'):143,('TYR','LYS'):85,
 ('TYR','ASP'):160,('TYR','GLU'):122,('TYR','MET'):36,('TYR','TRP'):37,
 ('CYS','HIS'):174,('CYS','GLN'):154,('CYS','ASN'):139,('CYS','LYS'):202,('CYS','ASP'):154,
 ('CYS','GLU'):170,('CYS','MET'):196,('CYS','TRP'):215,
 ('HIS','GLN'):24,('HIS','ASN'):68,('HIS','LYS'):32,('HIS','ASP'):81,('HIS','GLU'):40,
 ('HIS','MET'):87,('HIS','TRP'):115,
 ('GLN','ASN'):46,('GLN','LYS'):53,('GLN','ASP'):61,('GLN','GLU'):29,('GLN','MET'):101,
 ('GLN','TRP'):130,
 ('ASN','LYS'):94,('ASN','ASP'):23,('ASN','GLU'):42,('ASN','MET'):142,('ASN','TRP'):174,
 ('LYS','ASP'):101,('LYS','GLU'):56,('LYS','MET'):95,('LYS','TRP'):110,
 ('ASP','GLU'):45,('ASP','MET'):160,('ASP','TRP'):181,
 ('GLU','MET'):126,('GLU','TRP'):152,
 ('MET','TRP'):67,
}
AA3 = {'F':'PHE','L':'LEU','I':'ILE','M':'MET','V':'VAL','S':'SER','P':'PRO','T':'THR','A':'ALA',
       'Y':'TYR','H':'HIS','Q':'GLN','N':'ASN','K':'LYS','D':'ASP','E':'GLU','C':'CYS','W':'TRP',
       'R':'ARG','G':'GLY'}

def grantham(aa1, aa2):
    if aa1 == aa2:
        return 0.0
    if aa1 == '*' or aa2 == '*':
        return 215.0  # sentinel: worst-case distance (stop codon = catastrophic, not a substitution)
    a, b = AA3[aa1], AA3[aa2]
    key = (a, b) if (a, b) in GRANTHAM else (b, a)
    return float(GRANTHAM.get(key, 100.0))  # 100 = table mean fallback for the few untabulated pairs

def wobble_synonymous_fraction(table):
    """Fraction of position-3 point mutations (same codon, base changed at position 3 only) that
    are SYNONYMOUS (same amino acid) -- the claim's 68.9% number."""
    same = 0
    total = 0
    for cod in CODONS:
        aa0 = table[cod]
        for b in BASES:
            if b == cod[2]:
                continue
            mut = cod[0] + cod[1] + b
            total += 1
            if table[mut] == aa0:
                same += 1
    return same / total

def eta2_by_position(table, pos):
    """eta^2 (variance-explained) of pairwise Grantham distance-for-single-point-mutations,
    grouped by whether the two codons SHARE the nucleotide at `pos` (0-indexed) -- the claim's
    'clustering by 2nd-position nucleotide' check (pos=1)."""
    groups = {}
    all_d = []
    for cod in CODONS:
        aa0 = table[cod]
        for p in range(3):
            for b in BASES:
                if b == cod[p]:
                    continue
                mut = cod[:p] + b + cod[p+1:]
                d = grantham(aa0, table[mut])
                key = cod[pos] if p != pos else None  # group by nucleotide AT pos, for mutations NOT at pos
                if p == pos:
                    continue
                groups.setdefault(cod[pos], []).append(d)
                all_d.append(d)
    grand_mean = np.mean(all_d)
    ss_tot = sum((d - grand_mean) ** 2 for d in all_d)
    ss_between = 0.0
    for g, vals in groups.items():
        m = np.mean(vals)
        ss_between += len(vals) * (m - grand_mean) ** 2
    return ss_between / ss_tot if ss_tot > 0 else 0.0

def mean_1mut_cost(table):
    """Mean Grantham distance over all single-point-mutation neighbor pairs -- the code's
    'error cost'. Lower = more error-minimizing (mutations tend to land on similar amino acids)."""
    ds = []
    for cod in CODONS:
        aa0 = table[cod]
        for p in range(3):
            for b in BASES:
                if b == cod[p]:
                    continue
                mut = cod[:p] + b + cod[p+1:]
                ds.append(grantham(aa0, table[mut]))
    return float(np.mean(ds))

def two_opt_optimize(start_table, iters=4000, seed=0):
    """2-opt local search over the codon-to-amino-acid ASSIGNMENT (permute which 1 of the 21
    tokens -- 20 aa + stop -- each of the 64 codons' amino-acid LABEL is, holding codon sequences
    fixed; swap two amino-acid labels across the whole table and accept if mean cost improves).
    This mirrors 'how much better could the code be' by re-optimizing the codon<->aa MAPPING."""
    rng = random.Random(seed)
    table = dict(start_table)
    aa_list = sorted(set(table.values()))
    cur_cost = mean_1mut_cost(table)
    for _ in range(iters):
        a, b = rng.sample(aa_list, 2)
        codons_a = [c for c, aa in table.items() if aa == a]
        codons_b = [c for c, aa in table.items() if aa == b]
        # swap labels a<->b everywhere
        for c in codons_a:
            table[c] = b
        for c in codons_b:
            table[c] = a
        new_cost = mean_1mut_cost(table)
        if new_cost < cur_cost:
            cur_cost = new_cost
        else:
            for c in codons_a:
                table[c] = a
            for c in codons_b:
                table[c] = b
    return table, cur_cost

def scramble_table(seed):
    rng = random.Random(seed)
    aas = list(STD_TABLE.values())
    rng.shuffle(aas)
    return dict(zip(CODONS, aas))

def main():
    std_wobble = wobble_synonymous_fraction(STD_TABLE)
    std_eta2 = eta2_by_position(STD_TABLE, pos=1)
    std_cost = mean_1mut_cost(STD_TABLE)

    # void floor 1: random-relabeled codon tables, N=100
    rand_wobbles, rand_eta2s, rand_costs = [], [], []
    for s in range(100):
        t = scramble_table(s)
        rand_wobbles.append(wobble_synonymous_fraction(t))
        rand_eta2s.append(eta2_by_position(t, pos=1))
        rand_costs.append(mean_1mut_cost(t))
    rand_wobble_mean = float(np.mean(rand_wobbles))
    rand_eta2_mean = float(np.mean(rand_eta2s))
    rand_cost_mean = float(np.mean(rand_costs))

    # 2-opt seeded FROM the standard code
    _, std_seeded_cost = two_opt_optimize(STD_TABLE, iters=3000, seed=1)
    ratio_std_seeded = rand_cost_mean / std_seeded_cost

    # void floor 2: 2-opt seeded from a RANDOM start (avg of 5 seeds)
    rand_start_costs = []
    for s in range(5):
        _, c = two_opt_optimize(scramble_table(1000 + s), iters=3000, seed=2000 + s)
        rand_start_costs.append(c)
    rand_start_cost_mean = float(np.mean(rand_start_costs))
    ratio_rand_seeded = rand_cost_mean / rand_start_cost_mean

    claim = {
        "wobble_synonymous_fraction": 0.689,
        "eta2_range": [0.578, 0.756],
        "qap_optimum_ratio": 3.489,
        "two_opt_ratio": 3.4886,
    }
    gates = {
        "G1_wobble_matches_claim": bool(abs(std_wobble - claim["wobble_synonymous_fraction"]) < 0.03),
        "G2_wobble_voidfloor_destroyed": bool(abs(rand_wobble_mean - std_wobble) > 0.05),
        "G3_eta2_in_claimed_range": bool(claim["eta2_range"][0] - 0.05 <= std_eta2 <= claim["eta2_range"][1] + 0.05),
        "G4_eta2_voidfloor_collapses": bool(rand_eta2_mean < 0.15),
        "G5_std_seeded_ratio_near_claim": bool(abs(ratio_std_seeded - claim["two_opt_ratio"]) / claim["two_opt_ratio"] < 0.30),
        "G6_voidfloor_random_start_worse_or_far": bool(ratio_rand_seeded < ratio_std_seeded * 0.7 or abs(ratio_rand_seeded - ratio_std_seeded) / ratio_std_seeded > 0.15),
    }
    all_pass = bool(all(gates.values()))
    result = {
        "node": "MODEL-GENETIC-CODE-ERROR-MINIMIZATION",
        "measured": {
            "std_wobble_synonymous_fraction": std_wobble,
            "std_eta2_pos2": std_eta2,
            "std_mean_1mut_cost": std_cost,
            "random_relabel_wobble_mean_n100": rand_wobble_mean,
            "random_relabel_eta2_mean_n100": rand_eta2_mean,
            "random_relabel_cost_mean_n100": rand_cost_mean,
            "two_opt_from_standard_code_cost": std_seeded_cost,
            "ratio_random_mean_over_std_seeded_2opt": ratio_std_seeded,
            "two_opt_from_random_start_cost_mean_n5": rand_start_cost_mean,
            "ratio_random_mean_over_randomstart_2opt": ratio_rand_seeded,
        },
        "claim": claim,
        "gates": gates,
        "all_gates_pass": all_pass,
        "void_floor_note": "G2/G4/G6 are the void floors (scrambled codon-table relabeling, N=100, "
                            "and 2-opt from a random start, N=5 seeds) -- they must FAIL to reproduce "
                            "the standard-code structure, else the gate is not discriminating.",
    }
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as f:
        json.dump(result, f, indent=2)
    print(json.dumps(result, indent=2))
    return result

if __name__ == "__main__":
    main()
