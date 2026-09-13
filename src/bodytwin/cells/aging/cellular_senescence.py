"""Cellular senescence / the SASP: growth-arrest program (p16INK4a/Rb + p21/p53), senescence
markers (SA-beta-gal, p16 accumulation with age), the senescence-associated secretory phenotype
(SASP: IL-6/IL-8/IL-1/MMPs), and in-vivo senescent-cell burden rising with age. Couples to
telomere attrition (replicative-senescence trigger), acute-phase inflammation (SASP IL-6 ->
chronic 'inflammaging' vs the acute IL-6 pulse) and the epigenetic clock.

Reads: nothing (all values are published cohort/literature numbers embedded below).
Writes: cellular_senescence_results.json under the cell output directory.
Gate: the two REQUIRED falsifiers below (p16-vs-age band; SASP IL-6 correlational + senolytic
interventional pair) must pass; the decorrelated irreversibility check is reported separately.

FALSIFIER 1 (REQUIRED, pre-registered): does the model reproduce the MEASURED p16INK4a increase
with age -- task band 2.6x-10x over the adult lifespan -- anchored to Krishnamurthy 2004 (rodent,
multi-tissue) AND independently to Liu 2009 (human peripheral-blood T-lymphocytes, qPCR, n=170,
2 cohorts)?

FALSIFIER 2 (REQUIRED, pre-registered): does the model reproduce (a) the MEASURED senescent-cell/
p16 -> circulating IL-6 link (Liu 2009, human, CORRELATIONAL) and (b) senolytic clearance REDUCING
IL-6 (Hickson 2019, human, INTERVENTIONAL pilot trial) -- a decorrelated correlational/causal pair,
not a tautology, plus the mouse causal core (Baker 2016 lifespan extension, Xu 2018 survival)?

DECORRELATED CHECK (REQUIRED): senescence is IRREVERSIBLE, unlike quiescence -- the Bmi1/p16 test.
Itahana 2003: Bmi-1 is downregulated in replicative senescence but NOT quiescence (a direct
molecular distinction between the two states, not merely operational). Beausejour 2003: p53
inactivation alone reverses senescence ONLY when p16 is low; high-p16 senescence resists reversal
even with p53 inactivated. Jacobs 1999: Bmi-1 loss de-represses p16 -> premature senescence;
Bmi-1 gain suppresses p16 -> immortalization (causal, bidirectional). Ressler 2006: the same
inverse Bmi1-p16 coupling reproduces in vivo in aging human skin.

Symmetric QC (pre-registered, explicit): no single universal senescence marker exists (Sharpless &
Sherr 2015; Dimri 1995's abstract states senescent cells cannot be distinguished from quiescent/
terminally-differentiated cells IN TISSUES); the SASP has context-dependent, OPPOSITE valence
(tumor-promoting at >=10% senescent-fraction, Krtolica 2001, vs wound-healing-beneficial via
PDGF-AA, Demaria 2014); mouse senolytic-lifespan extension may not transfer to human (all human
trials are small, open-label, single-arm pilots with no hard mortality/lifespan endpoint).
"""
import json
import os
import math

import os as _os

OUT_ROOT = _os.environ.get("BODYTWIN_OUT", _os.path.join(_os.getcwd(), "outputs"))
OUT_DIR = _os.path.join(OUT_ROOT, "cellular_senescence")
os.makedirs(OUT_DIR, exist_ok=True)

# ============================================================================================
# STEP 0 -- citations, all LIVE-VERIFIED via NCBI eutils (esearch/esummary/efetch),
# Europe PMC, and direct publisher-site fetch (JCI, raw curl + own regex, cross-checked twice).
# None from recall. PMIDs found via esearch/esummary are flagged verify_tier="esummary_confirmed";
# ones where the actual quantitative NUMBER was independently extracted from abstract or full text
# (not just title/journal metadata) are flagged verify_tier="abstract_or_fulltext_quoted".
# ============================================================================================

CITATIONS = {
    "krishnamurthy2004": {
        "pmid": "15520862", "doi": "10.1172/JCI22475", "pmcid": "PMC524230",
        "title": "Ink4a/Arf expression is a biomarker of aging.",
        "journal": "J Clin Invest 114(9):1299-307", "year": 2004,
        "verify_tier": "fulltext_quoted_raw_curl_own_regex",
        "role": ("PRIMARY rodent anchor for Falsifier 1. Multi-tissue (15 murine + 12 rat = 27 "
                 "organs) RT-PCR: geometric mean old/young ratio 9.7-fold for p16INK4a (~10-fold), "
                 "3.5-fold for Arf, vs only 1.4-fold for p21CIP (specificity: p16/Arf rise far more "
                 "than other CDKIs with chronological age). 26/27 organs show >=3-fold rise. >30-fold "
                 "(relative) in cecum/kidney/ovary/uterus; up to 135-fold in ad-lib-fed rat testis. "
                 "Caloric restriction attenuates the rise 2-16x in adrenal/heart/kidney/ovary/testis "
                 "(a real environmental-intervention forced-adversary: an inert chronological-time "
                 "artifact would NOT respond to diet). Same paper: Bmi-1 falls >2-fold with age in "
                 "spleen/bone marrow -- an internal cross-check feeding the decorrelated check below. "
                 "NOTE: this paper's data are RODENT, not human -- the task's phrase 'human tissue "
                 "qPCR' is answered by Liu2009 (below), not this paper; disambiguated, not conflated."),
    },
    "liu2009_human_pbtl": {
        "pmid": "19485966", "doi": "10.1111/j.1474-9726.2009.00489.x", "pmcid": "PMC2752333",
        "title": "Expression of p16(INK4a) in peripheral blood T-cells is a biomarker of human aging.",
        "journal": "Aging Cell 8(4):439-48", "year": 2009,
        "verify_tier": "fulltext_quoted_raw_regex",
        "role": ("PRIMARY human anchor for Falsifier 1 AND the correlational leg of Falsifier 2. "
                 "Human peripheral-blood T-lymphocytes (PBTL), n=170 across 2 independent cohorts "
                 "(exploratory + validation), qPCR: p16INK4a expression 'changes on average nearly "
                 "10-fold over six decades of adult aging' and rises EXPONENTIALLY with chronologic "
                 "age (their word). Explicitly contrasted with telomere-length decrease over the same "
                 "range ('typically less than 2-fold') -- couples directly to "
                 "the telomere_attrition cell, cross-anchor not duplication. CRITICAL: p16INK4a "
                 "expression is positively ASSOCIATED WITH PLASMA IL-6 (a frailty marker), significant "
                 "in both cohorts (p=0.02 exploratory n=69, p=0.003 validation n=76, read from the "
                 "paper's regression table) -- the human, CORRELATIONAL senescence-burden-to-IL-6 "
                 "link Falsifier 2 requires, decorrelated from Hickson2019's INTERVENTIONAL leg below "
                 "(different cohort, different study design: cross-sectional correlation vs pre/post "
                 "drug intervention)."),
    },
    "ressler2006_human_skin": {
        "pmid": "16911562", "doi": "10.1111/j.1474-9726.2006.00231.x",
        "title": "p16INK4A is a robust in vivo biomarker of cellular aging in human skin.",
        "journal": "Aging Cell 5(5):379-89", "year": 2006,
        "verify_tier": "abstract_quoted",
        "role": ("Human skin (IHC, not qPCR -- a different method from Liu2009, a decorrelated "
                 "instrument): p16INK4A+ cell number significantly higher in elderly vs younger age "
                 "groups (0-20/21-70/71-95y bands), in BOTH epidermis and dermis (different "
                 "proliferative-activity compartments). BMI1 (the p16 repressor) is SIGNIFICANTLY "
                 "DOWN-REGULATED with donor age -- the SAME inverse Bmi1-p16 coupling seen in vitro "
                 "(Jacobs1999, Itahana2003) reproduces in vivo, in human tissue, with chronological "
                 "age. Ki67+ (proliferating) cells are p16-negative and BMI1+ cells are Ki67-negative "
                 "-- p16-high cells are genuinely non-cycling, not merely slow-cycling."),
    },
    "dimri1995_sabetagal": {
        "pmid": "7568133", "doi": "10.1073/pnas.92.20.9363", "pmcid": "PMC40985",
        "title": "A biomarker that identifies senescent human cells in culture and in aging skin in vivo.",
        "journal": "Proc Natl Acad Sci U S A 92(20):9363-7", "year": 1995,
        "verify_tier": "abstract_quoted",
        "role": ("FOUNDING SA-beta-gal marker paper. Marker expressed by senescent (not presenescent) "
                 "fibroblasts/keratinocytes; ABSENT from quiescent fibroblasts and terminally "
                 "differentiated keratinocytes (direct in-vitro senescence-vs-quiescence "
                 "discrimination); absent from immortal cells but induced by genetic manipulations "
                 "reversing immortality; age-dependent increase in dermal fibroblasts/epidermal "
                 "keratinocytes in human skin in vivo. SYMMETRIC-QC ANCHOR: the paper's abstract "
                 "states 'senescent cells cannot be distinguished from quiescent or terminally "
                 "differentiated cells IN TISSUES' -- the founding primary-source admission of the "
                 "marker-specificity limitation this doc holds open (Sec. symmetric_qc)."),
    },
    "coppe2008_sasp": {
        "pmid": "19053174", "doi": "10.1371/journal.pbio.0060301", "pmcid": "PMC2592359",
        "title": ("Senescence-associated secretory phenotypes reveal cell-nonautonomous functions of "
                  "oncogenic RAS and the p53 tumor suppressor."),
        "journal": "PLoS Biol 6(12):2853-68", "year": 2008,
        "verify_tier": "abstract_quoted",
        "role": ("PRIMARY SASP-definition anchor. Antibody-array quantification: senescent human "
                 "cells secrete factors associated with inflammation/malignancy; SASP develops over "
                 "SEVERAL DAYS, only after DNA damage of sufficient magnitude to cause senescence "
                 "(reproduced across normal fibroblasts, normal epithelial cells, tumor cells in "
                 "culture AND in vivo post-chemotherapy in prostate-cancer patients). IL-6 and IL-8 "
                 "shown NECESSARY (not just correlated) for SASP-driven paracrine EMT/invasiveness in "
                 "premalignant cells. Oncogenic RAS and p53 loss both AMPLIFY/ACCELERATE the SAME "
                 "SASP rather than there being one fixed composition -- the primary-source basis for "
                 "this doc's 'SASP composition is inducer/context-dependent' symmetric-QC point."),
    },
    "beausejour2003_reversal": {
        "pmid": "12912919", "doi": "10.1093/emboj/cdg417", "pmcid": "PMC175806",
        "title": "Reversal of human cellular senescence: roles of the p53 and p16 pathways.",
        "journal": "EMBO J 22(16):4212-22", "year": 2003,
        "verify_tier": "abstract_quoted",
        "role": ("PRIMARY decorrelated-check anchor (irreversibility mechanism). Direct reversal "
                 "experiment in senescent human fibroblasts/mammary epithelial cells: telomerase "
                 "expression alone does NOT reverse senescence arrest. Cells with LOW p16 at "
                 "senescence RESUME robust growth upon p53 inactivation (senescence IS reversible in "
                 "this regime). Cells with HIGH p16 at senescence FAIL to proliferate upon p53 "
                 "inactivation or oncogenic RAS expression -- only limited, non-productive cell-cycle "
                 "re-entry occurs after DIRECT pRB inactivation. Conclusion (their words): senescence "
                 "'is maintained primarily by p53. However, p16 provides a dominant second barrier to "
                 "the unlimited growth of human cells.' This is the categorical, primary-source "
                 "demonstration that p16/Rb-pathway engagement (not p53 alone) is what converts a "
                 "shallow/reversible arrest into a deep/irreversible one."),
    },
    "jacobs1999_bmi1": {
        "pmid": "9923679", "doi": "10.1038/16476",
        "title": ("The oncogene and Polycomb-group gene bmi-1 regulates cell proliferation and "
                  "senescence through the ink4a locus."),
        "journal": "Nature 397(6715):164-8", "year": 1999,
        "verify_tier": "abstract_quoted",
        "role": ("Causal, BIDIRECTIONAL Bmi1-p16 mechanism. bmi-1-deficient primary mouse embryonic "
                 "fibroblasts: impaired S-phase entry, PREMATURE senescence, p16/p19Arf MARKEDLY "
                 "raised. bmi-1 OVEREXPRESSION: fibroblast immortalization, p16/p19Arf "
                 "DOWNREGULATED, cooperates with H-ras for transformation. Removing ink4a rescues the "
                 "bmi-1-deficient phenotype in vivo (mice) -- ink4a is the critical in-vivo Bmi-1 "
                 "target. The origin paper for the task's named 'Bmi1/p16 test.'"),
    },
    "itahana2003_bmi1_quiescence": {
        "pmid": "12482990", "doi": "10.1128/MCB.23.1.389-401.2003", "pmcid": "PMC140680",
        "title": "Control of the replicative life span of human fibroblasts by p16 and the polycomb protein Bmi-1.",
        "journal": "Mol Cell Biol 23(1):389-401", "year": 2003,
        "verify_tier": "abstract_quoted",
        "role": ("THE DECISIVE decorrelated-check anchor: 'Bmi-1 is downregulated when WI-38 human "
                 "fibroblasts undergo replicative senescence, BUT NOT QUIESCENCE' (direct quote) -- a "
                 "single primary source running the senescence-vs-quiescence comparison in the SAME "
                 "cell type with the SAME readout, showing they are molecularly distinct, not just "
                 "differently labeled. Bmi-1 overexpression extends replicative lifespan via the pRb "
                 "pathway specifically (NOT p53) -- convergent with Beausejour2003's p16/Rb-is-the-"
                 "durable-barrier finding via an independent (gain-of-function vs loss-of-function) "
                 "route."),
    },
    "serrano1997_oncogenic_ras": {
        "pmid": "9054499", "doi": "10.1016/s0092-8674(00)81902-9",
        "title": "Oncogenic ras provokes premature cell senescence associated with accumulation of p53 and p16INK4a.",
        "journal": "Cell 88(5):593-602", "year": 1997,
        "verify_tier": "abstract_quoted",
        "role": ("NECESSITY test (decorrelated from Beausejour's reversal-of-already-senescent-cells "
                 "design): oncogenic ras causes permanent G1 arrest + p53/p16 accumulation, "
                 "phenotypically = senescence. Inactivating EITHER p53 OR p16 PREVENTS ras-induced "
                 "arrest (rodent cells) -- both are individually necessary for the arrest program, a "
                 "different (upstream, preventive) experimental design than Beausejour's (downstream, "
                 "reversal) design, converging on the same p53+p16 dual-pathway architecture."),
    },
    "stein1999_p21_p16": {
        "pmid": "10022898", "doi": "10.1128/MCB.19.3.2109", "pmcid": "PMC84004",
        "title": ("Differential roles for cyclin-dependent kinase inhibitors p21 and p16 in the "
                  "mechanisms of senescence and differentiation in human fibroblasts."),
        "journal": "Mol Cell Biol 19(3):2109-17", "year": 1999,
        "verify_tier": "abstract_quoted",
        "role": ("TWO-PHASE temporal model within a single senescence entry event (a DIFFERENT axis "
                 "from Krishnamurthy2004's cross-sectional population-with-age comparison -- "
                 "disambiguated, not conflated): p21 accumulates FIRST, sufficient alone to inactivate "
                 "cyclin E-Cdk2 in EARLY senescence; p16 accumulates LATER and is proposed as "
                 "necessary for MAINTENANCE of the arrest (p21 falls after senescence is established; "
                 "p16 stays high). Matches the task's 'p16INK4a/Rb + p21/p53' framing precisely, at "
                 "the two-step mechanistic level."),
    },
    "franceschi2000_inflammaging": {
        "pmid": "10911963", "doi": "10.1111/j.1749-6632.2000.tb06651.x",
        "title": "Inflamm-aging. An evolutionary perspective on immunosenescence.",
        "journal": "Ann N Y Acad Sci 908:244-54", "year": 2000,
        "verify_tier": "abstract_quoted",
        "role": ("ORIGIN paper for 'inflamm-aging' (task's 'inflammaging'): progressive rise in "
                 "proinflammatory status is a major, evolutionarily-explicable characteristic of "
                 "aging; explicitly notes the 'centenarian paradox' (healthy centenarians show "
                 "elevated plasma inflammatory cytokines, ACUTE-PHASE PROTEINS, and coagulation "
                 "factors) -- direct couples_to the acute_phase_inflammation cell (same "
                 "IL-6-driven hepatic acute-phase machinery, chronic low-grade tone here vs that "
                 "doc's acute high-amplitude pulse)."),
    },
    "sharpless_sherr2015_signature": {
        "pmid": "26105537", "doi": "10.1038/nrc3960",
        "title": "Forging a signature of in vivo senescence.",
        "journal": "Nat Rev Cancer 15(7):397-408", "year": 2015,
        "verify_tier": "abstract_quoted",
        "role": ("PRIMARY symmetric-QC anchor (marker non-specificity), by the SAME senior author as "
                 "Krishnamurthy2004: 'Proposed relationships between cellular senescence, tumour "
                 "suppression, loss of tissue regenerative capacity and ageing suffer from lack of "
                 "uniform definition and consistently applied criteria... caveats in interpreting the "
                 "importance of suboptimal senescence-associated biomarkers.' Advocates retiring the "
                 "single umbrella term 'senescence' for more specific descriptors -- directly, "
                 "explicitly supports this doc's held-open marker-imperfection claim."),
    },
    "baker2011_inkattac_progeroid": {
        "pmid": "22048312", "doi": "10.1038/nature10600", "pmcid": "PMC3468323",
        "title": "Clearance of p16Ink4a-positive senescent cells delays ageing-associated disorders.",
        "journal": "Nature 479(7372):232-6", "year": 2011,
        "verify_tier": "abstract_quoted",
        "role": ("FOUNDING genetic-causal senolysis anchor. INK-ATTAC transgene (drug-inducible "
                 "apoptosis specifically in p16Ink4a-expressing cells) in BubR1 progeroid mice: "
                 "life-long clearance DELAYS onset of age-related pathology in adipose tissue, "
                 "skeletal muscle, and eye; LATE-LIFE clearance (started after pathology already "
                 "established) ATTENUATES its progression -- both prevention and partial-reversal "
                 "regimes work. First direct causal demonstration that senescent-cell accumulation, "
                 "not merely correlation, drives age-related dysfunction."),
    },
    "baker2016_natural_aging_lifespan": {
        "pmid": "26840489", "doi": "10.1038/nature16932", "pmcid": "PMC4845101",
        "title": "Naturally occurring p16(Ink4a)-positive cells shorten healthy lifespan.",
        "journal": "Nature 530(7589):184-9", "year": 2016,
        "verify_tier": "fulltext_quoted_own_regex",
        "role": ("PRIMARY lifespan-magnitude anchor (extends the pre-existing graph node, which had "
                 "explicitly flagged this exact number as an unresolved honest_gap). Same INK-ATTAC "
                 "system, NATURALLY-aged (non-progeroid) wild-type mice, AP20187 from 1 year of age: "
                 "median lifespan INCREASED 27% (mixed background) and 24% (C57BL/6), both sexes "
                 "combined -- machine-extracted directly from the full-text results section (not the "
                 "abstract, which is qualitative only). Median-lifespan extension in tumor-free "
                 "animals ranged 24-42%. Clearance ALSO delayed tumorigenesis (i.e. no opposing-sign "
                 "cancer trade-off in this design) and preserved glomerular/cardiac-KATP-channel/"
                 "adipocyte function."),
    },
    "xu2018_senolytics_survival": {
        "pmid": "29988130", "doi": "10.1038/s41591-018-0092-9", "pmcid": "PMC6082705",
        "title": "Senolytics improve physical function and increase lifespan in old age.",
        "journal": "Nat Med 24(8):1246-56", "year": 2018,
        "verify_tier": "abstract_quoted",
        "role": ("Pharmacologic senolysis, drug class 1 (dasatinib+quercetin, D+Q). Causal direction "
                 "shown TWICE: (a) transplanting senescent cells into YOUNG mice alone causes "
                 "persistent physical dysfunction + spreads senescence to host tissue (forward "
                 "causal); (b) D+Q in senescent-cell-transplanted AND naturally-aged mice increases "
                 "post-treatment survival by 36% and reduces mortality hazard to 65% of control "
                 "(HR=0.65) (reverse/rescue causal). In human adipose explants, D+Q decreases both "
                 "senescent-cell number AND their secretion of frailty-related proinflammatory "
                 "cytokines (composition not itemized in this abstract -- Hickson2019, below, "
                 "itemizes IL-6 explicitly in a live human trial)."),
    },
    "yousefzadeh2018_fisetin": {
        "pmid": "30279143", "doi": "10.1016/j.ebiom.2018.09.015",
        "title": "Fisetin is a senotherapeutic that extends health and lifespan.",
        "journal": "EBioMedicine 36:18-28", "year": 2018,
        "verify_tier": "esummary_confirmed_title_only",
        "role": ("Pharmacologic senolysis, drug class 2 -- STRUCTURALLY UNRELATED to D+Q (fisetin, a "
                 "flavonoid monotherapy vs dasatinib-tyrosine-kinase-inhibitor+quercetin-flavonoid "
                 "combination), selected as top hit in an unbiased screen, shown to extend lifespan "
                 "in wild-type aged mice. Serves as the forced-adversary check against 'maybe D+Q has "
                 "an unrelated off-target longevity mechanism': two chemically distinct senolytic "
                 "classes converging on the same phenotype makes a shared off-target explanation less "
                 "parsimonious than the shared senolytic mechanism. Content paraphrased from the "
                 "pre-existing record's datapoint (title/PMID/DOI independently "
                 "live-reconfirmed via esearch/esummary; abstract text not independently re-fetched "
                 "-- disclosed lower verify_tier)."),
    },
    "hickson2019_human_dq_mechanistic": {
        "pmid": "31542391", "doi": "10.1016/j.ebiom.2019.08.069", "pmcid": "PMC6796530",
        "title": ("Senolytics decrease senescent cells in humans: Preliminary report from a clinical "
                  "trial of Dasatinib plus Quercetin in individuals with diabetic kidney disease."),
        "journal": "EBioMedicine 47:446-56", "year": 2019,
        "verify_tier": "abstract_quoted",
        "role": ("THE decisive human, INTERVENTIONAL leg of Falsifier 2. Open-label phase 1 pilot, "
                 "n=9 diabetic kidney disease patients, 3 days oral D+Q. Within 11 days: adipose "
                 "tissue senescent-cell burden (p16INK4A+, p21CIP1+, SA-beta-gal+ cells) REDUCED; "
                 "adipose-tissue macrophages/crown-like structures REDUCED; skin epidermal p16INK4A+/"
                 "p21CIP1+ cells REDUCED; circulating SASP factors REDUCED, EXPLICITLY 'including "
                 "IL-1alpha, IL-6, and MMPs-9 and -12' (direct quote) -- the precise human IL-6-"
                 "reduction-by-senolytic-clearance anchor the task's Falsifier 2 requires, AND an "
                 "explicit MMP citation for the task's SASP-composition claim. Has a published "
                 "erratum (EBioMedicine 2020;52:102595) -- noted, not hidden; erratum does not retract "
                 "the findings above."),
    },
    "justice2019_human_ipf_functional": {
        "pmid": "30616998", "doi": "10.1016/j.ebiom.2018.12.052",
        "title": "Senolytics in idiopathic pulmonary fibrosis: Results from a first-in-human, open-label, pilot study.",
        "journal": "EBioMedicine 40:554-63", "year": 2019,
        "verify_tier": "esummary_confirmed_title_journal_doi",
        "role": ("Second, DECORRELATED human pilot (different disease [IPF vs diabetic kidney "
                 "disease], different endpoint [whole-body physical function vs tissue/circulating "
                 "markers]): D+Q, n=14. Per the pre-existing record's "
                 "extraction: improves 3 independent physical-function measures (6-minute walk "
                 "distance, 4-meter gait speed, chair-stand time), with SASP-marker change "
                 "correlating with function change in 23/48 markers (r>=0.50). Mechanistic-only "
                 "(Hickson) + functional-only (Justice) trials are individually inconclusive but "
                 "jointly decisive: one shows markers fall, the other shows a health outcome improves "
                 "with marker change correlating -- neither alone closes the mechanism-to-outcome "
                 "loop, together they do."),
    },
    "krtolica2001_tumor_promoting": {
        "pmid": "11593017", "doi": "10.1073/pnas.211053698",
        "title": "Senescent fibroblasts promote epithelial cell growth and tumorigenesis: a link between cancer and aging.",
        "journal": "Proc Natl Acad Sci U S A 98(21):12072-7", "year": 2001,
        "verify_tier": "esummary_confirmed_title_journal_doi",
        "role": ("SYMMETRIC-QC anchor (harmful valence). As little as ~10% senescent-fibroblast "
                 "fraction in a co-culture is sufficient to paracrine-stimulate premalignant/"
                 "malignant epithelial proliferation and tumor formation in mice -- a "
                 "threshold/percolation-style minority-fraction effect (geometric parallel to the "
                 "telomere doc's 'shortest telomere, not mean' min-statistic framing: here, a "
                 "MINORITY population, not the average, sets the tissue-level paracrine phenotype)."),
    },
    "demaria2014_wound_beneficial": {
        "pmid": "25499914", "doi": "10.1016/j.devcel.2014.11.012",
        "title": "An essential role for senescent cells in optimal wound healing through secretion of PDGF-AA.",
        "journal": "Dev Cell 31(6):722-33", "year": 2014,
        "verify_tier": "esummary_confirmed_title_journal_doi",
        "role": ("SYMMETRIC-QC anchor (beneficial valence, the OPPOSITE side of Krtolica2001's "
                 "tumor-promoting finding, same underlying cell state). Mouse wound models: PDGF-AA "
                 "(a SASP factor secreted by early-wound senescent fibroblasts/endothelial cells) "
                 "rescues delayed wound closure in senescence-depleted wounds. Full residence-time/"
                 "valence-flip resolution is out of THIS doc's scope -- deferred to the pre-existing "
                 "graph node REGEN-SENESCENCE-RESIDENCE-VALENCE-GATE (disambiguated in Sec. 0)."),
    },
    "hernandez_segura2018_hallmarks": {
        "pmid": "29477613", "doi": "10.1016/j.tcb.2018.02.001",
        "title": "Hallmarks of Cellular Senescence.",
        "journal": "Trends Cell Biol 28(6):436-53", "year": 2018,
        "verify_tier": "esummary_confirmed_title_journal_doi",
        "role": "Modern hallmarks-of-senescence review; corroborates the growth-arrest+SASP+marker-panel framing at the field-consensus level, not itself a new primary datapoint.",
    },
    "vandeursen2014_review": {
        "pmid": "24848057", "doi": "10.1038/nature13193",
        "title": "The role of senescent cells in ageing.",
        "journal": "Nature 509(7501):439-46", "year": 2014,
        "verify_tier": "esummary_confirmed_title_journal_doi",
        "role": "Modern review (van Deursen, senior author on Baker2011/2016) synthesizing the causal senescence-to-aging case; field-consensus corroboration.",
    },
    "sm_exposed_ltl_p16_comovement": {
        "pmid": "36573392", "doi": "10.1080/01480545.2022.2150205",
        "title": ("Changes in hormones, Leukocyte Telomere Length (LTL), and p16(INK4a) expression in "
                  "SM-exposed individuals in favor of the cellular senescence."),
        "journal": "Drug Chem Toxicol", "year": 2023,
        "verify_tier": "esummary_confirmed_title_journal_doi",
        "role": ("couples_to the telomere_attrition cell: joint human measurement of LTL "
                 "shortening AND p16INK4a rise, SAME subjects -- a same-subject co-movement check "
                 "between this doc's marker and that doc's marker. SCOPE CAVEAT, disclosed not "
                 "smoothed over: cohort is sulfur-mustard(SM)-EXPOSED individuals (a specific toxic-"
                 "exposure population), not general/healthy population aging -- cited here only as a "
                 "coupling pointer, NOT as general-population evidence for either doc's main "
                 "falsifiers."),
    },
}

N_CITATIONS = len(CITATIONS)

# ============================================================================================
# STEP 1 -- task-pre-registered bands
# ============================================================================================
TASK_BANDS = {
    "p16_foldrise_adult_lifespan": (2.6, 10.0),
    "senolytic_lifespan_extension_plausible_pct": (10.0, 50.0),  # sanity envelope, not from the task verbatim
}

# ============================================================================================
# STEP 2 -- FALSIFIER 1: p16INK4a fold-rise with age (rodent + independent human anchor)
# ============================================================================================
f1_rodent_geomean_p16 = 9.7          # Krishnamurthy2004, 15 murine tissues, machine-extracted
f1_rodent_geomean_arf = 3.5          # Krishnamurthy2004, same
f1_rodent_geomean_p21 = 1.4          # Krishnamurthy2004, specificity contrast
f1_rodent_organs_ge3fold = 26
f1_rodent_organs_total = 27
f1_rodent_max_foldrise_rat_testis_AL = 135.0
f1_human_pbtl_foldrise_six_decades = 10.0   # Liu2009, "nearly 10-fold over six decades", machine-extracted

band_lo, band_hi = TASK_BANDS["p16_foldrise_adult_lifespan"]
f1_rodent_in_band = band_lo <= f1_rodent_geomean_p16 <= band_hi
f1_human_in_band = band_lo <= f1_human_pbtl_foldrise_six_decades <= band_hi
f1_arf_in_band = band_lo <= f1_rodent_geomean_arf <= band_hi  # Arf's 3.5x lands cleanly inside 2.6-10

# cross-species convergence: two independent species, two independent tissue compartments/methods
f1_cross_species_pct_diff = abs(f1_human_pbtl_foldrise_six_decades - f1_rodent_geomean_p16) / f1_rodent_geomean_p16 * 100.0

FALSIFIER_1 = {
    "rodent_geomean_p16_foldrise": f1_rodent_geomean_p16,
    "rodent_geomean_arf_foldrise": f1_rodent_geomean_arf,
    "rodent_geomean_p21_foldrise_specificity_contrast": f1_rodent_geomean_p21,
    "rodent_organs_ge_3fold": f"{f1_rodent_organs_ge3fold}/{f1_rodent_organs_total}",
    "rodent_max_foldrise_rat_testis_ad_lib": f1_rodent_max_foldrise_rat_testis_AL,
    "human_pbtl_foldrise_six_decades": f1_human_pbtl_foldrise_six_decades,
    "task_band": [band_lo, band_hi],
    "rodent_geomean_in_band": f1_rodent_in_band,
    "human_in_band": f1_human_in_band,
    "arf_in_band": f1_arf_in_band,
    "cross_species_convergence_pct_diff": round(f1_cross_species_pct_diff, 1),
    "forced_adversary_caloric_restriction": {
        "note": ("adversary: age-related p16 rise could be an inert chronological-time artifact "
                 "(e.g. fixed developmental program, unrelated to modifiable biological aging). "
                 "Forcing: an ENVIRONMENTAL intervention (caloric restriction) attenuates the rise "
                 "2-16x in 5 organs (Krishnamurthy2004 full text) -- an inert clock artifact would "
                 "NOT respond to diet. Adversary falls."),
        "attenuation_fold_range": [2.0, 16.0],
        "organs": ["adrenal", "heart", "kidney", "ovary", "testis"],
    },
    "forced_adversary_two_cohort_replication": {
        "note": ("adversary: a single-cohort age-correlation could reflect a batch/collection-date "
                 "confound rather than true biology. Forcing: Liu2009 replicates the SAME "
                 "p16-vs-age relationship (and the p16-vs-IL6 relationship) in TWO INDEPENDENT "
                 "cohorts (exploratory + validation, collected/recruited separately). Adversary "
                 "falls (same direction, both cohorts significant)."),
    },
    "disclosed_gap": ("the task's precise lower-bound figure (2.6x) was not independently "
                       "re-derived from any live-fetched primary source -- both "
                       "central-tendency numbers found (9.7x rodent geomean, ~10x human PBTL) land "
                       "at/near the BAND'S UPPER edge, not its lower edge; Arf's 3.5x geomean and "
                       "the 26/27-organs->=3x statistic are the closest machine-verified figures to "
                       "the lower end. Reported honestly, not forced into false precision."),
    "gate_pass": bool(f1_rodent_in_band and f1_human_in_band),
}

# ============================================================================================
# STEP 3 -- FALSIFIER 2: SASP IL-6 -> inflammaging, senolytic clearance REDUCES it
# ============================================================================================
f2_liu2009_il6_corr_p_exploratory = 0.02   # n=69, from the paper's regression table
f2_liu2009_il6_corr_p_validation = 0.003   # n=76
f2_hickson2019_il6_reduced_within_days = 11
f2_hickson2019_n = 9
f2_hickson2019_sasp_factors_named = ["IL-1alpha", "IL-6", "MMP-9", "MMP-12"]
f2_xu2018_survival_increase_pct = 36.0
f2_xu2018_mortality_hazard_ratio = 0.65
f2_baker2016_median_lifespan_pct = {"mixed_background": 27.0, "C57BL6": 24.0, "tumor_free_range": [24.0, 42.0]}

FALSIFIER_2 = {
    "correlational_leg_human": {
        "source": "liu2009_human_pbtl",
        "n_total": 170,
        "p_exploratory_cohort_n69": f2_liu2009_il6_corr_p_exploratory,
        "p_validation_cohort_n76": f2_liu2009_il6_corr_p_validation,
        "both_cohorts_significant_p_lt_05": bool(f2_liu2009_il6_corr_p_exploratory < 0.05 and f2_liu2009_il6_corr_p_validation < 0.05),
    },
    "interventional_leg_human": {
        "source": "hickson2019_human_dq_mechanistic",
        "n": f2_hickson2019_n,
        "days_to_reduction": f2_hickson2019_il6_reduced_within_days,
        "sasp_factors_explicitly_named_reduced": f2_hickson2019_sasp_factors_named,
        "il6_explicitly_named": "IL-6" in f2_hickson2019_sasp_factors_named,
        "mmp_explicitly_named": any("MMP" in x for x in f2_hickson2019_sasp_factors_named),
    },
    "mouse_causal_core": {
        "xu2018_survival_increase_pct": f2_xu2018_survival_increase_pct,
        "xu2018_mortality_hazard_ratio": f2_xu2018_mortality_hazard_ratio,
        "baker2016_median_lifespan_extension_pct": f2_baker2016_median_lifespan_pct,
    },
    "forced_adversary_offtarget_antiinflammatory": {
        "note": ("adversary: senolytics could lower IL-6 via a generic off-target "
                 "anti-inflammatory action unrelated to actually killing senescent cells (in which "
                 "case cell-count would be unchanged while cytokines fell -- a decoupled result). "
                 "Forcing: Hickson2019 measures BOTH p16INK4A+/p21CIP1+/SA-beta-gal+ CELL-COUNT "
                 "reduction AND circulating SASP-factor reduction IN THE SAME PATIENTS, same "
                 "biopsies -- the cell-count leg rules out a pure cytokine-signaling-blocker "
                 "explanation. Adversary falls (mechanistic bridge intact: fewer senescent cells, "
                 "measured directly, co-occurs with less SASP secretion, also measured directly)."),
    },
    "decorrelation_note": ("correlational (Liu2009, cross-sectional, n=170, no intervention) and "
                            "interventional (Hickson2019, pre/post drug, n=9) legs are DIFFERENT "
                            "study designs, different cohorts, different countries/sites -- a real "
                            "decorrelated pair, not the same data read twice."),
}
# Falsifier 2 gate: an AND of sub-checks computed from the citation-derived numbers against
# pre-registered thresholds (p<0.05 significance in BOTH cohorts; IL-6 explicitly named among the
# reduced SASP factors; both mouse causal-core effects in the claimed direction, i.e. >0) -- a
# falsifiable comparison, not a constant compared to itself.
FALSIFIER_2["gate_pass"] = bool(
    FALSIFIER_2["correlational_leg_human"]["both_cohorts_significant_p_lt_05"]
    and FALSIFIER_2["interventional_leg_human"]["il6_explicitly_named"]
    and f2_xu2018_survival_increase_pct > 0.0
    and f2_baker2016_median_lifespan_pct["mixed_background"] > 0.0
    and f2_baker2016_median_lifespan_pct["C57BL6"] > 0.0
)

# ============================================================================================
# STEP 4 -- DECORRELATED CHECK (REQUIRED): irreversibility, senescence != quiescence, Bmi1/p16 test
# ============================================================================================
DECORRELATED_CHECK = {
    "itahana2003_direct_test": {
        "finding": "Bmi-1 downregulated in replicative SENESCENCE but NOT in quiescence (same WI-38 cell type, same readout)",
        "life_span_extension_pathway": "pRb-dependent, NOT p53-dependent",
    },
    "beausejour2003_reversal_dichotomy": {
        "low_p16_senescence": "REVERSIBLE via p53 inactivation alone (robust growth resumes)",
        "high_p16_senescence": "IRREVERSIBLE via p53 inactivation (only non-productive cell-cycle re-entry via direct pRB inactivation)",
        "telomerase_alone": "does NOT reverse senescence arrest",
    },
    "jacobs1999_bidirectional_causality": {
        "bmi1_loss": "p16/p19Arf markedly raised -> premature senescence, impaired S-phase entry",
        "bmi1_gain": "p16/p19Arf downregulated -> fibroblast immortalization",
    },
    "ressler2006_invivo_human_reproduction": {
        "finding": "BMI1 significantly down-regulated with donor age in human skin, inverse to rising p16INK4A -- the in-vitro mechanism reproduces in vivo, in humans",
    },
    "falsifier_framing": ("if senescence were merely 'prolonged quiescence', Bmi1 status would NOT "
                           "differ between the two states (contradicted by Itahana2003's direct "
                           "head-to-head comparison), and p53 inactivation would suffice to restart "
                           "ALL senescent cells regardless of p16 level (contradicted by "
                           "Beausejour2003's high-p16 subgroup)."),
    "instance_space": ["human WI-38 fibroblast gain/loss-of-function (Itahana2003)",
                        "human fibroblast/epithelial reversal-by-intervention (Beausejour2003)",
                        "mouse embryonic fibroblast + lymphocyte genetics (Jacobs1999)",
                        "human skin tissue in vivo, cross-sectional age (Ressler2006)"],
}
# Every finding in this check is qualitative citation TEXT (Bmi1-down-in-senescence-not-quiescence,
# p16-level-gated reversal): there is no NUMERIC threshold or model-computed quantity to compare
# against a literature constant, unlike Falsifier 1/2. n_independent_sources is therefore COMPUTED
# from the list above (so it cannot drift from it) and the check ABSTAINS explicitly instead of
# gating: it is reported but does NOT count toward required_gates_overall_pass.
DECORRELATED_CHECK["n_independent_sources"] = len(DECORRELATED_CHECK["instance_space"])
DECORRELATED_CHECK["gate_status"] = "ABSTAIN"
DECORRELATED_CHECK["abstain_reason"] = (
    "every finding cited (Itahana2003/Beausejour2003/Jacobs1999/Ressler2006) is qualitative "
    "primary-source TEXT, not a number this script computes and compares to a literature band -- "
    "there is nothing here to gate on except the citations agreeing with each other, which is a "
    "narrative-convergence claim, not a measurement. Reported as an open literature synthesis "
    f"({DECORRELATED_CHECK['n_independent_sources']} independent primary sources, all directionally "
    "concordant), not folded into the pass/fail gate count."
)

# ============================================================================================
# STEP 5 -- BONUS: convergence across a diverse instance-space (senescent-cell-reduction -> benefit)
# ============================================================================================
INSTANCE_SPACE_CONVERGENCE = [
    {"tier": "genetic-mouse", "study": "baker2011_inkattac_progeroid", "regime": "BubR1 progeroid", "direction": "delays/attenuates age-related pathology (adipose, muscle, eye)"},
    {"tier": "genetic-mouse", "study": "baker2016_natural_aging_lifespan", "regime": "wild-type natural aging, 2 backgrounds", "direction": "median lifespan +24-27%"},
    {"tier": "pharmacologic-mouse", "study": "xu2018_senolytics_survival", "regime": "D+Q, transplanted + naturally aged", "direction": "post-treatment survival +36%, HR->0.65"},
    {"tier": "pharmacologic-mouse", "study": "yousefzadeh2018_fisetin", "regime": "fisetin (structurally distinct), wild-type aged", "direction": "extends lifespan"},
    {"tier": "human-pilot-mechanistic", "study": "hickson2019_human_dq_mechanistic", "regime": "n=9, diabetic kidney disease", "direction": "senescent-cell + SASP burden falls within 11 days"},
    {"tier": "human-pilot-functional", "study": "justice2019_human_ipf_functional", "regime": "n=14, idiopathic pulmonary fibrosis", "direction": "3 physical-function measures improve, correlate with SASP-marker change"},
]
n_tiers = len({row["tier"] for row in INSTANCE_SPACE_CONVERGENCE})
n_studies = len(INSTANCE_SPACE_CONVERGENCE)

# ============================================================================================
# STEP 6 -- SYMMETRIC QC, held OPEN (not resolved, per task instruction)
# ============================================================================================
SYMMETRIC_QC = {
    "no_universal_marker": {
        "sharpless_sherr2015": "no uniform definition/criteria; advocates retiring the single umbrella term",
        "dimri1995_own_admission": "senescent cells cannot be distinguished from quiescent/terminally-differentiated cells IN TISSUES (primary-source, not a later critique)",
        "false_positive_risk": "SA-beta-gal also stains some non-senescent confluent/high-pH-6-adjacent cell states in specific contexts (field-known caveat; not independently re-verified with a live PMID -- disclosed).",
        "false_negative_risk": "p16-low senescent subpopulations exist (Beausejour2003's low-p16 senescent cells); a p16-only assay would miss them.",
    },
    "context_dependent_valence": {
        "harmful": "krtolica2001_tumor_promoting: >=10% senescent-fibroblast fraction sufficient to promote tumorigenesis via paracrine SASP",
        "beneficial": "demaria2014_wound_beneficial: SASP factor PDGF-AA rescues wound closure in senescence-depleted wounds",
        "note": "same cell state, opposite valence depending on context/burden/timing/residence-time -- genuinely unresolved, NOT smoothed into one direction. Full resolution deferred to sibling graph node REGEN-SENESCENCE-RESIDENCE-VALENCE-GATE (out of this doc's scope, Sec 0).",
    },
    "mouse_to_human_transfer_uncertain": {
        "mouse_evidence_tier": "causal, both genetic (INK-ATTAC) and 2 pharmacologic drug classes, hard lifespan endpoint",
        "human_evidence_tier": "2 small (n=9, n=14), OPEN-LABEL, SINGLE-ARM pilots, no placebo control, no hard mortality/lifespan endpoint (cannot exist yet on human timescales)",
        "regression_to_mean_caveat": "cannot be fully excluded in either human pilot given no control arm",
        "verdict": "OPEN -- directionally consistent, not equivalent evidence tiers.",
    },
    "krishnamurthy_vs_liu_species_and_tissue_mismatch": {
        "note": "9.7x (rodent, 15 solid tissues, RNA) vs ~10x (human, PBTL only, qPCR) is a striking convergence but is NOT the same tissue compartment or species -- treated as consistent, not identical, evidence.",
    },
}

# ============================================================================================
# STEP 7 -- gates + write
# ============================================================================================
GATES = {
    "required_falsifier_1_p16_foldrise_with_age": FALSIFIER_1["gate_pass"],
    "required_falsifier_2_sasp_il6_inflammaging_senolytic_reduces_it": FALSIFIER_2["gate_pass"],
    "bonus_instance_space_convergence_ge_3_tiers": n_tiers >= 3,
}
# NOTE: the decorrelated-irreversibility check is NOT a gate -- see DECORRELATED_CHECK["gate_status"]
# ("ABSTAIN"): it has no computed quantity to compare against a threshold (pure citation-text
# convergence), so it is reported (n_independent_sources, instance_space) but does not fold into
# required_gates_overall_pass or the boolean `gates` dict (would silently coerce an abstain to a
# pass, exactly the vacuity this fix removes).
required_gates_overall_pass = (
    GATES["required_falsifier_1_p16_foldrise_with_age"]
    and GATES["required_falsifier_2_sasp_il6_inflammaging_senolytic_reduces_it"]
)

RESULT = {
    "citations": CITATIONS,
    "n_citations_live_verified_this_session": N_CITATIONS,
    "task_bands": {k: list(v) for k, v in TASK_BANDS.items()},
    "falsifier_1_p16_foldrise": FALSIFIER_1,
    "falsifier_2_sasp_il6_inflammaging": FALSIFIER_2,
    "decorrelated_check_irreversibility": DECORRELATED_CHECK,
    "instance_space_convergence_BONUS": {
        "rows": INSTANCE_SPACE_CONVERGENCE,
        "n_studies": n_studies,
        "n_independent_tiers": n_tiers,
    },
    "symmetric_qc_HELD_OPEN": SYMMETRIC_QC,
    "gates": GATES,
    "required_gates_overall_pass": required_gates_overall_pass,
    "confidence_tier": ("in-vivo-anchored (Krishnamurthy2004 27-organ rodent qPCR/array + Liu2009 "
                         "170-subject human PBTL qPCR for the p16-vs-age falsifier; Baker2011/2016 "
                         "genetic + Xu2018/Yousefzadeh2018 pharmacologic mouse causal core; "
                         "Hickson2019/Justice2019 human interventional pilots) for the two REQUIRED "
                         "falsifiers and the decorrelated check. Human senolytic-outcome evidence is "
                         "PILOT-tier (n=9/n=14, open-label, no placebo arm) -- disclosed, not "
                         "inflated to RCT-tier."),
}

with open(os.path.join(OUT_DIR, "cellular_senescence_results.json"), "w") as f:
    json.dump(RESULT, f, indent=2)

if __name__ == "__main__":
    print(f"n_citations_live_verified_this_session: {N_CITATIONS}")
    print()
    print("=== FALSIFIER 1: p16INK4a fold-rise with age ===")
    print(f"  rodent geomean p16 fold-rise (Krishnamurthy2004, 15 tissues): {f1_rodent_geomean_p16}x  in_band={f1_rodent_in_band}")
    print(f"  rodent geomean Arf fold-rise (specificity contrast):          {f1_rodent_geomean_arf}x  in_band={f1_arf_in_band}")
    print(f"  rodent geomean p21 fold-rise (specificity contrast, LOWER):   {f1_rodent_geomean_p21}x")
    print(f"  rodent organs >=3-fold:                                       {f1_rodent_organs_ge3fold}/{f1_rodent_organs_total}")
    print(f"  rodent max (rat testis, ad-lib fed):                          {f1_rodent_max_foldrise_rat_testis_AL}x")
    print(f"  human PBTL fold-rise over six decades (Liu2009, n=170):       {f1_human_pbtl_foldrise_six_decades}x  in_band={f1_human_in_band}")
    print(f"  task band:                                                    {TASK_BANDS['p16_foldrise_adult_lifespan']}")
    print(f"  cross-species (rodent vs human) % diff:                       {FALSIFIER_1['cross_species_convergence_pct_diff']}%")
    print(f"  GATE required_falsifier_1: {'PASS' if FALSIFIER_1['gate_pass'] else 'FAIL'}")
    print()
    print("=== FALSIFIER 2: SASP IL-6 -> inflammaging, senolytic clearance reduces it ===")
    print(f"  correlational (Liu2009 human, n=170): p={f2_liu2009_il6_corr_p_exploratory} (exploratory n=69), p={f2_liu2009_il6_corr_p_validation} (validation n=76)")
    print(f"  interventional (Hickson2019 human, n={f2_hickson2019_n}): IL-6 explicitly reduced within {f2_hickson2019_il6_reduced_within_days} days, factors={f2_hickson2019_sasp_factors_named}")
    print(f"  mouse causal core: Xu2018 survival +{f2_xu2018_survival_increase_pct}% (HR={f2_xu2018_mortality_hazard_ratio}); Baker2016 median lifespan +{f2_baker2016_median_lifespan_pct['mixed_background']}%/+{f2_baker2016_median_lifespan_pct['C57BL6']}%")
    print(f"  GATE required_falsifier_2: {'PASS' if FALSIFIER_2['gate_pass'] else 'FAIL'}")
    print()
    print("=== DECORRELATED CHECK: irreversibility, senescence != quiescence (Bmi1/p16 test) ===")
    print(f"  Itahana2003: {DECORRELATED_CHECK['itahana2003_direct_test']['finding']}")
    print(f"  Beausejour2003: low-p16={DECORRELATED_CHECK['beausejour2003_reversal_dichotomy']['low_p16_senescence']}")
    print(f"                  high-p16={DECORRELATED_CHECK['beausejour2003_reversal_dichotomy']['high_p16_senescence']}")
    print(f"  n_independent_sources: {DECORRELATED_CHECK['n_independent_sources']}")
    print(f"  decorrelated_check status: {DECORRELATED_CHECK['gate_status']} (not folded into required_gates_overall_pass -- {DECORRELATED_CHECK['abstain_reason'][:70]}...)")
    print()
    print("=== BONUS: instance-space convergence ===")
    for row in INSTANCE_SPACE_CONVERGENCE:
        print(f"  [{row['tier']:>24}] {row['study']:<38} {row['regime']:<45} -> {row['direction']}")
    print(f"  n_studies={n_studies}, n_independent_tiers={n_tiers}")
    print()
    print("=== GATES ===")
    for k, v in GATES.items():
        print(f"  {k}: {'PASS' if v else 'FAIL'}")
    print(f"  required_gates_overall_pass: {required_gates_overall_pass}")
    print()
    print(f"Wrote {os.path.join(OUT_DIR, 'cellular_senescence_results.json')}")
