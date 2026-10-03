# LANE_XSEED_REVIEW

Independent review, by a stronger model, of the swarm’s (The swarm/swarm_worker) results on Anton’s X bookmark seeds. Findings from weaker models do not count as "holds" or "bug" until a stronger model has verified them. Result directory `results/LANE_XSEED_REVIEW/`.

## Att granska (results/BT-XSEED-*)

1. **BT-XSEED-3B1-F1 (highest priority)** — heart transplant, step 2: φ (recipient-cell fraction) from 59 EPIC rs probes. Claims: (a) precision σ(φ)=0,0077 against requirement 0,073 (9,4× margin); (b) "allele-specific methylation (ASM)" makes φ_hat = φ·f with bias 0,10–0,17, 83 % of the original confounder, so the design still cannot distinguish rejuvenation from cell replacement. **Coordinator’s doubt**: the rs probes on 450K/EPIC are SNP genotyping probes (designed for sample identity, not CpG) — ASM should not affect them in the way modelled. Go to the Illumina manifest and the literature on rs-probe signal (and on SNP-under-probe effects), rerun the parent’s code (`results/BT-XSEED-3B1-F1/src`) and determine: is (a) correct? is (b) a real effect, a modelling error, or does it apply to probes other than rs probes? Deliver a verdict with numbers and what it means for the seed’s design.
2. **BT-XSEED-3B1** — the identity g_hat = g + φ (the coordinator has checked the algebraic core; review assumptions).
3. Sample the negatives: **3B4 myelin** (storage 2–4 J versus 1,15 MJ/day), **3B5 photons/electrons** (flow change 2,9e6× below the detection floor), **4B2 federated calibration** (sufficient statistics = exact pooling; leakage in per-client Fisher), **2B3 infarct material** (placement optimum outside the border zone). For each: does the decisive number hold on recomputation? Is there an error that reverses the conclusion?

## Leverans

REVIEW_TABLE.json: job → claim → verdict (CONFIRMED / WRONG / PARTIAL / UNDETERMINED) → recomputed number → file:line. RESULTS.md for readers. If 3B1-F1(b) is wrong: write what the corrected design would give and whether the seed’s step 2 then holds.

Inga interna data. Allt PENDING_INDEPENDENT_REVIEW.
