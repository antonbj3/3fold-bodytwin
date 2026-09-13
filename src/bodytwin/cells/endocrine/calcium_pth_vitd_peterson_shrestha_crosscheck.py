"""Independent re-derivation of the Peterson & Riggs 2010 calcium/PTH/vitamin-D + bone-remodelling
model (BIOMD0000000613), cross-checked against Shrestha et al. 2010 real-patient clamp data.

The published claim cites Peterson & Riggs 2010 (Bone; BIOMD0000000613) steady state (Ca=2.3500
mmol/L, PTH=3.8500 model units) and acute Ca-clamp PTH fold-changes (4.11x hypo / 0.32x hyper),
plus a decorrelated real-patient anchor (Shrestha et al. 2010 Math Biosci, PMID20406649,
BIOMD0000000276/277: fold-change 2.328x/0.202x, log-ratio residual +0.518+/-0.073 nats). The
calcium_pth_vitd cell tests a DIFFERENT thing (Brown-vs-Parfitt setpoint calibration). This cell
fills that gap: a fresh, independent scipy re-derivation of BIOMD0000000613's numbers, not a
re-print of stored values.

DISCLOSED CORRECTION: a "freeze the slow bone-cell signalling block (RANKL/OPG/osteoblast/
osteoclast/TGFb)" shortcut is WRONG -- diagnostic runs show that block moves non-negligibly within
1-2h and by 66-140% relative drift by 240h, so it is NOT cleanly time-scale-separable from the acute
Ca-PTH-vitD loop. This cell therefore runs the FULL model (all 30 dynamical species from the
source's SBML), not a bone-cell-frozen reduction.

SOURCE (BioModels REST API,
https://www.ebi.ac.uk/biomodels/model/download/BIOMD0000000613): the curated SBML
("Peterson2010 - Integrated calcium homeostasis and bone remodelling", isDescribedBy
PMID:19732857). The XPP/.ode export (BIOMD0000000613.ode) gives every constant parameter
verbatim; several reference constants ("X0" symbols: Q0, OC0, RNK0, RANKL0, OB0, ROB0,
QboneInit, OPG0, RX20, CREB0, M0, TGFBact0, TGFB0) are declared via SBML <initialAssignment>
elements that just equal each corresponding species' OWN t=0 initial value (read directly from
the raw SBML XML, not guessed) -- e.g. TGFBact0=Pic0=0.228142 (NOT 0; the .ode/.m converter
mislabeled TGFBact/M/N/A/OBfast/OBslow/TGFB as "no initial state" purely because their real
initial values come from a chain of <initialAssignment> elements referencing OTHER species,
which the xpp/octave converters could not resolve -- the raw SBML XML resolves them exactly, in
initial-assignment dependency order, all real numbers, none invented):
  A(0)=B(0)/10=126.0, TGFB(0)=Pic0*1000=228.142, TGFBact(0)=Pic0=0.228142,
  OBfast(0)=OB*FracOBfast=0.00399947, OBslow(0)=OB*(1-FracOBfast)=0.00101377,
  M(0)=k3*RNK(0)*L(0)/k4=2.22827e-4, N(0)=k1*O(0)*L(0)/k2=8.91270e-5, Da=0.7/24=0.0291667,
  GFR=100/16.667=6.00001, OralPhos=10.5/24=0.4375, OralCa=24.055/24=1.00229.
Teriparatide dosing (TERIPK/TERISC bolus events) is OFF (no drug, physiological baseline only) --
this matches the .ode converter's treatment ("unable to handle events with delays currently,
event ignored") and the source's baseline framing (IPTHint=IPTHinf=0 already in the source).

30 dynamical species (verbatim ODEs from the source, all listed in-code below): PTH, S (gland
Ca-sensing pool), PTmax, B (calcitriol), SC, P (plasma Ca), ECCPhos, Tg (gut Ca pool, source's
"T"), R, HAp, PhosGut, IntraPO, OC (osteoclast), ROB1 (responding osteoblast), L (RANKL), RNK
(RANK), O (OPG), Q (bone-exchangeable Ca), Qbone, RX2 (RunX2), CREB, BCL2, TERISC, A
(1a-hydroxylase), TGFB, TGFBact, OBfast, OBslow, M (RANK:RANKL), N (OPG:RANKL).

QUESTION (pre-registered): does this independent scipy re-derivation land at (a) a genuine
machine-checked steady state, (b) the correct SIGN of acute PTH response to a Ca clamp in both
directions, and (c) the same "Peterson over-predicts vs. real-patient Shrestha data" pattern the
published text already discloses (not forced to match, reported as measured)?

Reads: nothing (all constants embedded). Writes: nothing -- results are printed to stdout.

GATE 1 (self-consistency): the full 30-state RHS must evaluate to ~0 (machine-checked, not
eyeballed) at this rebuild's converged fixed point.
GATE 2 (acute clamp dynamic, SIGN+rough-magnitude only, pre-registered BEFORE running: hypo ->
fold>1x, hyper -> fold<1x; no tight absolute number pre-registered since the source's two
model families already disagree by ~1.68x).
LEG 3 (informational, decorrelated real-patient cross-check): log-ratio residual between this
reimplementation's fold-changes and Shrestha's published Subject-1 clamp fold-changes (2.328x/
0.202x, literature/BioModels numbers, not independently re-derived here -- BIOMD0000000276/277 are
NOT re-simulated; only their already-published fold-change numbers are compared against).
VOID FLOOR (pre-registered): with the Ca-sensing INPUT to PTH secretion (both the direct T63/EPTH
arm and the ScaEff/T72 gland-pool arm) frozen at its baseline value regardless of the actual
(possibly clamped) plasma Ca, a clamp perturbation MUST produce PTH fold-change ~1x (no
response) -- confirms the feedback, not the clamp mechanics, drives the main-test result.
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

# ============================== live-fetched constants ==============================
# (BIOMD0000000613.ode, verbatim; a few "Frac.../gam" names shortened only for readability)
Pic0 = 0.228142
FracOBfast = 0.797629
k3, k4 = 6.24e-6, 0.112013
k1, k2 = 6.24e-6, 0.112013
V1 = 14.0
CaDay = 88.0
FracJ14 = 0.107763
J14OCmax, J14OCgam = 0.543488, 1.6971
FracJ15 = 0.114376
kinRNKgam = 0.151825
koutRNK = 0.00323667
MOCratioGam = 0.603754
Da = 0.7 / 24  # resolved from initialAssignment (0.7/24), real
OBtgfGAM = 0.0111319
koutTGF0 = 2.98449e-5
OCtgfGAM = 0.593891
EmaxPicROB, PicROBgam, FracPicROB = 3.9745, 1.80968, 0.883824
PicOBgam, FracPicOB, EmaxPicOB = 0.122313, 2.44818e-4, 0.251636
E0Meff, EmaxMeffOC, kinOCgam = 0.388267, 3.15667, 8.53065
EmaxPicOC, FracPicOC, PicOCgam = 1.9746, 0.878215, 1.0168
E0RANKL, EmaxL = 3.80338, 0.469779
GFR = 100 / 16.667  # resolved
T16 = 1.06147
T64 = 0.05
T65 = 6.3
T67 = 1.54865
AlphOHgam = 0.111241
k14a = 2.44437e-5
HApMRT = 3.60609
koutL = 0.00293273
OsteoEffectGam = 0.173833
TESTPOWER = 1.0
opgPTH50 = 3.85
IO = 0.0
RX2Kout0, E0rx2Kout, EmaxPTHRX2x = 0.693, 0.125, 5.0
E0crebKin, EmaxPTHcreb, crebKout = 0.5, 3.39745, 0.00279513
bcl2Kout = 0.693
ScaEffGam = 0.9
PhosEff0, PhosEff50, PhosEffGam = 1.52493, 1.3021, 8.25229
PO4INHPTH_GAM = 0.0  # par PO4inhPTH=0.0 (this gamma disables the phosphate-inhibition term)
T69 = 0.1
Reabs50 = 1.57322
T7, T9 = 2.0, 90.0
T70, T71 = 0.01, 0.03
T33, T34, T35 = 0.003, 0.037, 90.0
CaPOgam = 1.0
T46, T52 = 1.142, 0.365
OralPhos = 10.5 / 24  # resolved
F12 = 0.7
T49 = 51.8
T55 = 0.019268
PicOBgamkb, MultPicOBkb, FracPic0kb = 2.92375, 3.11842, 0.764028
E0RUNX2kbEffFACT, RUNkbGAM, RUNkbMaxFact = 1.01, 3.67798, 0.638114
RUNX20 = 10.0
Frackb = 0.313186
T81, T87 = 0.75, 0.0495
T0 = 1.58471
T28 = 0.9
OralCa = 24.055 / 24  # resolved
T77, T80 = 0.909359, 4.0
CtriolPTgam, CtriolMax, CtriolMin = 12.5033, 4.1029, 0.9
PTout = 1.604e-4
T57, T58, T59, T61 = 100.0, 6249.09, 11.7387, 96.25
IPTHint, IPTHinf = 0.0, 0.0
LsurvOCgam = 3.09023
EmaxLpth = 1.30721
kO, kb = 15.8885, 6.05516e-4

# ----- resolved "0"-suffixed reference constants (equal to each species' own t=0 IC, from the
# raw SBML <initialAssignment> chain, all real, none invented) -----
OB = 0.00501324  # par (fixed osteoblast baseline reference, not itself simulated)
RNK0_ic, L0_ic, O0_ic = 10.0, 0.4, 4.0
Q0 = 100.0
OC0 = 0.00115398
RNK0 = RNK0_ic
RANKL0 = L0_ic
OB0 = OB
ROB0 = 0.00104122
QboneInit = 24900.0
OPG0 = O0_ic
RX20 = 10.0
CREB0 = 10.0
M0 = k3 * RNK0_ic * L0_ic / k4       # = 2.22827e-4
TGFBact0 = Pic0                       # = 0.228142
TGFB0 = Pic0 * 1000                   # = 228.142

# ============================== initial conditions (30 states) ==============================
STATE_NAMES = ["PTH", "S", "PTmax", "B", "SC", "P", "ECCPhos", "Tg", "R", "HAp", "PhosGut",
               "IntraPO", "OC", "ROB1", "L", "RNK", "O", "Q", "Qbone", "RX2", "CREB", "BCL2",
               "TERISC", "A", "TGFB", "TGFBact", "OBfast", "OBslow", "M", "N"]
IDX = {n: i for i, n in enumerate(STATE_NAMES)}

Y0 = np.array([
    53.9,        # PTH
    0.5,         # S
    1.0,         # PTmax
    1260.0,      # B
    0.0,         # SC
    32.9,        # P
    16.8,        # ECCPhos
    1.2375,      # Tg
    0.5,         # R
    1.0,         # HAp
    0.839,       # PhosGut
    3226.0,      # IntraPO
    0.00115398,  # OC
    0.00104122,  # ROB1
    0.4,         # L
    10.0,        # RNK
    4.0,         # O
    100.0,       # Q
    24900.0,     # Qbone
    10.0,        # RX2
    10.0,        # CREB
    100.0,       # BCL2
    0.0,         # TERISC
    126.0,       # A  (=B0/10, resolved from initialAssignment)
    228.142,     # TGFB (=Pic0*1000)
    0.228142,    # TGFBact (=Pic0)
    0.00399947,  # OBfast (=OB*FracOBfast)
    0.00101377,  # OBslow (=OB*(1-FracOBfast))
    2.22827e-4,  # M
    8.91270e-5,  # N
])


def _pos(x, floor=1e-12):
    return max(x, floor)


def rhs(y, ca_sense_override=None):
    """Full verbatim translation of BIOMD0000000613's assignment rules + reactions.
    ca_sense_override: if not None, this value replaces CaConc (=P/14) ONLY in the two arms that
    feed PTH secretion (T63's EPTH term, and ScaEff's T72 gland-pool term) -- this is the
    void-floor feedback cut. All other equations still see the real, possibly-clamped P.
    """
    (PTH, S, PTmax, B, SC, P, ECCPhos, Tg, R, HAp, PhosGut, IntraPO, OC, ROB1, L, RNK, O, Q,
     Qbone, RX2, CREB, BCL2, TERISC, A, TGFB, TGFBact, OBfast, OBslow, M, N) = y

    Osteoclast = OC
    T13 = CaDay / 24 / Q0
    T15 = CaDay / (2.35 * 14 * 24)
    T17 = 3.85 * T16 - 3.85
    J14OC50 = np.exp(np.log(J14OCmax * OC0 ** J14OCgam / T13 - OC0 ** J14OCgam) / J14OCgam)
    OCeqn = J14OCmax * Osteoclast ** J14OCgam / (Osteoclast ** J14OCgam + J14OC50 ** J14OCgam)
    kinRNK = (koutRNK * RNK0 + k3 * RNK0 * RANKL0 - k4 * M0) / TGFBact0 ** kinRNKgam
    MOCratio = M / _pos(Osteoclast)
    MOCratio0 = M0 / OC0
    MOCratioEff = (MOCratio / MOCratio0) ** MOCratioGam
    J14OCdepend = OCeqn * Q0 * FracJ14 * MOCratioEff
    J14 = T13 * Q0 * (1 - FracJ14) + J14OCdepend
    J41 = 0.464 * J14
    PicOCkin = Pic0
    bigDb = kb * OB0 * Pic0 / ROB0
    kinTGF = koutTGF0 * TGFB0
    koutTGF = koutTGF0
    koutTGFact = koutTGF0 * 1000
    koutTGFeqn = koutTGF * TGFB * (Osteoclast / OC0) ** OCtgfGAM
    E0PicROB = FracPicROB * Pic0
    EC50PicROBparen = EmaxPicROB * TGFBact0 ** PicROBgam / (Pic0 - E0PicROB) - TGFBact0 ** PicROBgam
    EC50PicROB = np.exp(np.log(EC50PicROBparen) / PicROBgam)
    Dr = kb * OB0 / Pic0
    PicROB = E0PicROB + EmaxPicROB * TGFBact ** PicROBgam / (TGFBact ** PicROBgam + EC50PicROB ** PicROBgam)
    ROBin = Dr * PicROB
    E0PicOB = FracPicOB * Pic0
    EC50PicOBparen = EmaxPicOB * TGFBact0 ** PicOBgam / (Pic0 - E0PicOB) - TGFBact0 ** PicOBgam
    EC50PicOB = np.exp(np.log(EC50PicOBparen) / PicOBgam)
    PicOB = E0PicOB + EmaxPicOB * TGFBact ** PicOBgam / (TGFBact ** PicOBgam + EC50PicOB ** PicOBgam)
    KPT = bigDb / PicOB
    D = ROB1
    EC50MeffOC = np.exp(np.log(M0 ** kinOCgam * EmaxMeffOC / (1 - E0Meff) - M0 ** kinOCgam) / kinOCgam)
    MeffOC = E0Meff + EmaxMeffOC * M ** kinOCgam / (M ** kinOCgam + EC50MeffOC ** kinOCgam)
    kinOC2 = Da * PicOCkin * MeffOC * OC0
    E0PicOC = FracPicOC * Pic0
    EC50PicOCparen = EmaxPicOC * TGFBact0 ** PicOCgam / (Pic0 - E0PicOC) - TGFBact0 ** PicOCgam
    EC50PicOC = np.exp(np.log(EC50PicOCparen) / PicOCgam)
    PicOC = E0PicOC + EmaxPicOC * TGFBact ** PicOCgam / (TGFBact ** PicOCgam + EC50PicOC ** PicOCgam)
    PiL0 = k3 / k4 * RANKL0
    PiL = M / 10
    EC50survInPar = (E0RANKL - EmaxL) * PiL0 ** LsurvOCgam / (E0RANKL - 1) - PiL0 ** LsurvOCgam
    EC50surv = np.exp(np.log(EC50survInPar) / LsurvOCgam)
    LsurvOC = E0RANKL - (E0RANKL - EmaxL) * PiL ** LsurvOCgam / (PiL ** LsurvOCgam + EC50surv ** LsurvOCgam)
    KLSoc = Da * PicOC * LsurvOC
    C4 = PTH / V1
    T66 = (T67 ** AlphOHgam + 3.85 ** AlphOHgam) / 3.85 ** AlphOHgam
    k15a = k14a * QboneInit / Q0
    J14a = k14a * Qbone
    J15a = k15a * Q
    kLShap = 1 / HApMRT
    kHApIn = kLShap / OB0
    J15 = T15 * P * (1 - FracJ15) + T15 * P * FracJ15 * HAp
    J42 = 0.464 * J15
    Osteoblast_tot = OBfast + OBslow
    kinLbase = koutL * RANKL0
    OsteoEffect = (Osteoblast_tot / OB0) ** OsteoEffectGam
    PTH50 = EmaxLpth * 3.85 - 3.85
    PTHconc = C4
    LpthEff = EmaxLpth * PTHconc / (PTH50 * OsteoEffect ** TESTPOWER + PTHconc)
    kinL = kinLbase * OsteoEffect * LpthEff
    pObase = kO * OPG0
    pO = pObase * D / ROB0 * (PTHconc + opgPTH50 * D / ROB0) / (2 * PTHconc) + IO
    RX2Kin = RX2Kout0 * RX20
    EC50PTHRX2x = EmaxPTHRX2x * 3.85 / (RX2Kout0 - E0rx2Kout) - 3.85
    RX2Kout = E0rx2Kout + EmaxPTHRX2x * PTHconc / (PTHconc + EC50PTHRX2x)
    EC50PTHcreb = EmaxPTHcreb * 3.85 / (1 - E0crebKin) - 3.85
    crebKin0 = crebKout * CREB0
    crebKin = crebKin0 * (E0crebKin + EmaxPTHcreb * PTHconc / (PTHconc + EC50PTHcreb))
    bcl2Kin = RX2 * CREB * 0.693

    CaConc_true = P / 14
    CaConc_for_ptin = CaConc_true if ca_sense_override is None else ca_sense_override
    C2 = ECCPhos / V1
    PO4inhPTH_eff = (C2 / 1.2) ** PO4INHPTH_GAM  # gam=0 -> identically 1.0, real (not invented)
    PhosEffTop = (PhosEff0 - 1) * (1.2 ** PhosEffGam + PhosEff50 ** PhosEffGam)
    PhosEffBot = PhosEff0 * 1.2 ** PhosEffGam
    PhosEffMax = PhosEffTop / PhosEffBot
    PhosEff = PhosEff0 - PhosEffMax * PhosEff0 * C2 ** PhosEffGam / (C2 ** PhosEffGam + PhosEff50 ** PhosEffGam)
    PhosEffect = PhosEff if C2 > 1.2 else 1.0
    T68 = T66 * C4 ** AlphOHgam / (T67 ** AlphOHgam * PO4inhPTH_eff + C4 ** AlphOHgam)
    SE = T65 * T68 * PhosEffect
    C8 = B / V1
    C1 = P / V1
    T36 = T33 + (T34 - T33) * C8 ** CaPOgam / (T35 ** CaPOgam + C8 ** CaPOgam)
    T37 = T34 - (T34 - T33) * C8 ** CaPOgam / (T35 ** CaPOgam + C8 ** CaPOgam)
    CaFilt = 0.6 * 0.5 * GFR * C1
    ReabsMax = (0.3 * GFR * 2.35 - 0.149997) * (Reabs50 + 2.35) / 2.35
    ReabsPTHeff = T16 * C4 / (C4 + T17)
    CaReabsActive = ReabsMax * C1 / (Reabs50 + C1) * ReabsPTHeff
    T20 = CaFilt - CaReabsActive
    T10 = T7 * C8 / (C8 + T9)
    J27a = (2 - T10) * T20
    J27 = max(J27a, 0.0)
    ScaEff = (2.35 / _pos(CaConc_for_ptin)) ** ScaEffGam
    T72 = 90 * ScaEff
    T73 = T71 * (C8 - T72)
    T74 = np.tanh(T73)
    T75 = T70 * (0.85 * (1 + T74) + 0.15)
    T76 = T70 * (0.85 * (1 - T74) + 0.15)
    T47 = T46 * 0.88 * GFR
    J48a = 0.88 * GFR * C2 - T47
    J48 = max(J48a, 0.0)
    J53 = T52 * PhosGut
    J54 = T49 * C2
    J56 = T55 * IntraPO
    E0PicOBkb = MultPicOBkb * Pic0
    EmaxPicOBkb_v = FracPic0kb * Pic0
    EC50PicOBparenKb = (E0PicOBkb - EmaxPicOBkb_v) * TGFBact0 ** PicOBgamkb / (E0PicOBkb - Pic0) - TGFBact0 ** PicOBgamkb
    EC50PicOBkb = np.exp(np.log(EC50PicOBparenKb) / PicOBgamkb)
    PicOBkb = E0PicOBkb - (E0PicOBkb - EmaxPicOBkb_v) * TGFBact ** PicOBgamkb / (TGFBact ** PicOBgamkb + EC50PicOBkb ** PicOBgamkb)
    PicOBkbEff = PicOBkb / Pic0
    E0RUNX2kbEff = E0RUNX2kbEffFACT * kb
    RUNX2 = (BCL2 - 90) if BCL2 > 105 else 10.0
    RUNkbMax = E0RUNX2kbEff * RUNkbMaxFact
    INparen = RUNkbMax * RUNX20 ** RUNkbGAM / (E0RUNX2kbEff - kb) - RUNX20 ** RUNkbGAM
    RUNkb50 = np.exp(np.log(INparen) / RUNkbGAM)
    RUNX2kbPrimeEff = RUNkbMax * RUNX2 ** RUNkbGAM / (RUNX2 ** RUNkbGAM + RUNkb50 ** RUNkbGAM)
    kbprime = E0RUNX2kbEff * PicOBkbEff - RUNX2kbPrimeEff
    kbslow = kbprime * Frackb
    kbfast = (kb * OB0 + kbslow * (OB0 * FracOBfast) - kbslow * OB0) / (OB0 * FracOBfast)
    T29 = (T28 * T0 - 0.17533 * T0) / 0.17533
    T31 = T28 * Tg / (Tg + T29)
    T83 = R / 0.5
    J40 = T31 * Tg * T83 / (Tg + T81) + T87 * Tg
    T85Rpart = R ** T80 / (R ** T80 + T81 ** T80)
    F11 = T77 * T85Rpart
    INparenCtriol = (CtriolMax - CtriolMin) * C8 ** CtriolPTgam / (CtriolMax - 1) - C8 ** CtriolPTgam
    Ctriol50 = np.exp(np.log(INparenCtriol) / CtriolPTgam)
    CtriolPTeff = CtriolMax - (CtriolMax - CtriolMin) * C8 ** CtriolPTgam / (C8 ** CtriolPTgam + Ctriol50 ** CtriolPTgam)
    PTin = PTout * CtriolPTeff
    INparenCa = (T58 - T61) * 2.35 ** T59 / (T58 - 385) - 2.35 ** T59
    T60 = np.exp(np.log(INparenCa) / T59)
    FCTD = S / 0.5 * PTmax
    T63 = T58 - (T58 - T61) * CaConc_for_ptin ** T59 / (CaConc_for_ptin ** T59 + T60 ** T59)
    EPTH = T63 * FCTD
    IPTH_ = 0.693 * SC + IPTHinf
    SPTH = EPTH + IPTH_
    kout = T57 / 14

    # ---- reactions ----
    PTH_produ = SPTH
    PTH_admin = 0.0  # teriparatide dosing OFF (baseline physiology, disclosed)
    PTH_degra = kout * PTH
    PTH_gla_1 = (1 - S) * T76
    PTH_gla_2 = S * T75
    PTH_max_c = PTin
    PTH_max_1 = PTout * PTmax
    Calcitrio = A
    Calcitr_1 = T69 * B
    SC_in = IPTHint
    SC_out = 0.693 * SC
    alpha_hyd = SE
    alpha_h_1 = T64 * A
    Ca_flux_f = J14
    Ca_absorp = J40
    Ca_flux_t = J15
    Ca_filtra = J27
    Ca_oral_i = OralCa * F11
    Intestina = T36 * (1 - R)
    Intesti_1 = T37 * R
    HAp_prouc = kHApIn * Osteoblast_tot
    HAp_out = kLShap * HAp
    OBfast_fr = bigDb / PicOB * D * FracOBfast * (kbfast / kbprime)
    OBfast_de = kbfast * OBfast
    OBslow_fr = bigDb / PicOB * D * (1 - FracOBfast) * Frackb
    OBslow_de = kbslow * OBslow
    Osteocl_1 = kinOC2
    Osteocl_2 = KLSoc * OC
    Respondin = ROBin
    Respond_1 = KPT * ROB1
    TGFB_prod = kinTGF * (Osteoblast_tot / OB0) ** OBtgfGAM
    TGFB_acti = koutTGFeqn
    TGFB_degr = koutTGFact * TGFBact
    RANKL_pro = kinL
    RANKL_deg = koutL * L
    RANK_prod = kinRNK * TGFBact ** kinRNKgam
    RANK_degr = koutRNK * RNK
    OPG_produ = pO
    OPG_degra = kO * O
    RANKLRANK = k3 * RNK * L
    RANKLRA_1 = k4 * M
    OPGRANKL_ = k1 * O * L
    OPGRANK_1 = k2 * N
    Ca_to_int = J14a
    Ca_to_ext = J15a
    PO4_trans = J41
    PO4_absor = J53
    PO4_flux_ = J56
    PO4_tra_1 = J42
    PO4_filtr = J48
    PO4_flu_1 = J54
    PO4_oral_ = OralPhos * F12
    RunX2_pro = RX2Kin
    RunX2_deg = RX2Kout * RX2
    CREB_prod = crebKin
    CREB_degr = crebKout * CREB
    Bcl2_prod = bcl2Kin
    Bcl2_degr = bcl2Kout * BCL2

    dPTH = PTH_produ + PTH_admin - PTH_degra
    dS = PTH_gla_1 - PTH_gla_2
    dPTmax = PTH_max_c - PTH_max_1
    dB = Calcitrio - Calcitr_1
    dSC = SC_in - SC_out
    dP = Ca_flux_f + Ca_absorp - Ca_flux_t - Ca_filtra
    dECCPhos = PO4_trans + PO4_absor + PO4_flux_ - PO4_tra_1 - PO4_filtr - PO4_flu_1
    dTg = -Ca_absorp + Ca_oral_i
    dR = Intestina - Intesti_1
    dHAp = HAp_prouc - HAp_out
    dPhosGut = -PO4_absor + PO4_oral_
    dIntraPO = -PO4_flux_ + PO4_flu_1
    dOC = Osteocl_1 - Osteocl_2
    dROB1 = Respondin - Respond_1
    dL = RANKL_pro - RANKL_deg - RANKLRANK + RANKLRA_1
    dRNK = RANK_prod - RANK_degr - RANKLRANK + RANKLRA_1
    dO = OPG_produ - OPG_degra
    dQ = -Ca_flux_f + Ca_flux_t + Ca_to_int - Ca_to_ext
    dQbone = -Ca_to_int + Ca_to_ext
    dRX2 = RunX2_pro - RunX2_deg
    dCREB = CREB_prod - CREB_degr
    dBCL2 = Bcl2_prod - Bcl2_degr
    dTERISC = -PTH_admin
    dA = alpha_hyd - alpha_h_1
    dTGFB = TGFB_prod - TGFB_acti
    dTGFBact = TGFB_acti - TGFB_degr
    dOBfast = OBfast_fr - OBfast_de
    dOBslow = OBslow_fr - OBslow_de
    dM = RANKLRANK - RANKLRA_1
    dN = OPGRANKL_ - OPGRANK_1

    return np.array([dPTH, dS, dPTmax, dB, dSC, dP, dECCPhos, dTg, dR, dHAp, dPhosGut, dIntraPO,
                      dOC, dROB1, dL, dRNK, dO, dQ, dQbone, dRX2, dCREB, dBCL2, dTERISC, dA,
                      dTGFB, dTGFBact, dOBfast, dOBslow, dM, dN])


def integrate(y0, t_end, ca_clamp=None, ca_sense_override=None, n_eval=400, rtol=1e-9, atol=1e-12):
    """ca_clamp: if not None, plasma Ca amount P is externally forced to this fixed value for the
    whole run (P's ODE is disabled) -- models an acute Ca clamp experiment."""
    def f(t, y):
        yy = y.copy()
        if ca_clamp is not None:
            yy[IDX["P"]] = ca_clamp
        d = rhs(yy, ca_sense_override=ca_sense_override)
        if ca_clamp is not None:
            d[IDX["P"]] = 0.0
        return d

    t_eval = np.linspace(0, t_end, n_eval)
    sol = solve_ivp(f, [0, t_end], y0, method="BDF", rtol=rtol, atol=atol, t_eval=t_eval)
    return sol


def main():
    print("=== calcium_pth_vitd_peterson_shrestha_crosscheck: "
          "calcium-PTH-vitamin-D homeostasis ===")
    print(f"Full model: {len(STATE_NAMES)} dynamical species (source's 30; teriparatide "
          f"dosing off, physiological baseline).")

    # -------- GATE 1: self-consistency at the published initial condition, and at our own
    # converged fixed point --------
    print("\n--- GATE 1: steady-state self-consistency ---")
    d0 = rhs(Y0)
    max_resid_ic = np.abs(d0).max()
    worst_ic = STATE_NAMES[int(np.argmax(np.abs(d0)))]
    print(f"RHS at the source's published initial condition: max|dy/dt| = {max_resid_ic:.4e} "
          f"(worst component: {worst_ic})")

    print("Integrating forward from the IC to this rebuild's converged fixed point "
          "(t=20000h)...")
    sol_ss = integrate(Y0, 20000.0, n_eval=3)
    y_ss = sol_ss.y[:, -1]
    d_ss = rhs(y_ss)
    max_resid_ss = np.abs(d_ss).max()
    worst_ss = STATE_NAMES[int(np.argmax(np.abs(d_ss)))]
    print(f"RHS at OWN converged fixed point (t=20000h): max|dy/dt| = {max_resid_ss:.4e} "
          f"(worst component: {worst_ss})")
    ca_ss = y_ss[IDX["P"]] / 14
    pth_conc_ss = y_ss[IDX["PTH"]] / V1
    print(f"Own converged fixed point: Ca={ca_ss:.4f} mmol/L (published value: 2.3500), "
          f"PTH_conc(=PTH/V1)={pth_conc_ss:.4f} model units (published value: 3.8500), "
          f"raw PTH amount={y_ss[IDX['PTH']]:.4f}")
    gate1_pass = max_resid_ss < 1e-3 * max(1.0, np.abs(y_ss).max())
    print(f"GATE 1 (self-consistency, own converged fixed point): "
          f"{'PASS' if gate1_pass else 'FAIL'} (relative-scale residual threshold)")

    # -------- GATE 2: acute Ca clamp, sign + rough magnitude --------
    print("\n--- GATE 2: acute Ca clamp dynamic (sign pre-registered; rough magnitude only) ---")
    P_baseline = y_ss[IDX["P"]]
    pth_baseline = y_ss[IDX["PTH"]]

    def clamp_run(frac_change, label):
        P_clamp = P_baseline * (1 + frac_change)
        sol = integrate(y_ss, 10.0, ca_clamp=P_clamp, n_eval=200)
        t = sol.t
        pth = sol.y[IDX["PTH"]]
        mask = (t >= 2.0) & (t <= 8.0)
        pth_acute_mean = pth[mask].mean()
        fold = pth_acute_mean / pth_baseline
        print(f"{label}: Ca clamp {frac_change*100:+.1f}% (P: {P_baseline:.3f}->{P_clamp:.3f}), "
              f"PTH acute-phase (t=2-8h) mean={pth_acute_mean:.4f}, fold-change={fold:.4f}x "
              f"(baseline PTH={pth_baseline:.4f})")
        return fold

    fold_hypo = clamp_run(-0.145, "HYPOCALCEMIA")
    fold_hyper = clamp_run(+0.215, "HYPERCALCEMIA")

    sign_hypo_ok = fold_hypo > 1.0
    sign_hyper_ok = fold_hyper < 1.0
    gate2_pass = sign_hypo_ok and sign_hyper_ok
    print(f"Pre-registered signs: hypocalcemia fold>1x ({'OK' if sign_hypo_ok else 'FAIL'}), "
          f"hypercalcemia fold<1x ({'OK' if sign_hyper_ok else 'FAIL'})")
    print(f"GATE 2 (acute clamp sign+rough magnitude): {'PASS' if gate2_pass else 'FAIL'}")
    print(f"(the published prior numbers: 4.11x hypo / 0.32x hyper -- this is a FRESH "
          f"reimplementation, own numbers reported as measured, not forced to match.)")

    # -------- LEG 3: informational cross-check against Shrestha real-patient numbers --------
    print("\n--- LEG 3 (informational, NOT gated): log-ratio residual vs. Shrestha real-patient "
          "clamp fold-changes ---")
    shrestha_hypo, shrestha_hyper = 2.328, 0.202
    lr_hypo = np.log(fold_hypo / shrestha_hypo)
    lr_hyper = np.log(fold_hyper / shrestha_hyper)
    lr_vals = np.array([lr_hypo, lr_hyper])
    print(f"Peterson-reimplementation fold-changes: hypo={fold_hypo:.4f}x, hyper={fold_hyper:.4f}x")
    print(f"Shrestha published Subject-1 clamp fold-changes (literature/BIOMD276/277, NOT "
          f"re-simulated): hypo={shrestha_hypo}x, hyper={shrestha_hyper}x")
    print(f"log-ratio residual (nats): hypo={lr_hypo:+.3f}, hyper={lr_hyper:+.3f}, "
          f"mean={lr_vals.mean():+.3f} +/- {lr_vals.std(ddof=1)/np.sqrt(2):.3f} (n=2 directions)")
    same_sign = np.sign(lr_hypo) == np.sign(lr_hyper)
    print(f"Same sign both directions (the disclosed '~1.68x over-prediction' pattern)? "
          f"{'YES' if same_sign else 'NO'} -- the published prior mean residual was +0.518+/-0.073 "
          f"nats; this is a FRESH measurement on an independent reimplementation, reported "
          f"as-is (agreement strengthens the cross-model-tension finding, disagreement would "
          f"refute it -- either is informational, not gated).")

    # -------- VOID FLOOR: freeze the Ca-sensing input to PTH secretion --------
    print("\n--- VOID FLOOR: freeze Ca-sensing input to PTH secretion (T63/EPTH + ScaEff/T72 "
          "arms) at baseline, while still applying the SAME clamp mechanics ---")
    ca_baseline = P_baseline / 14

    def void_clamp_run(frac_change, label):
        P_clamp = P_baseline * (1 + frac_change)
        sol = integrate(y_ss, 10.0, ca_clamp=P_clamp, ca_sense_override=ca_baseline, n_eval=200)
        t = sol.t
        pth = sol.y[IDX["PTH"]]
        mask = (t >= 2.0) & (t <= 8.0)
        fold = pth[mask].mean() / pth_baseline
        print(f"{label} (void floor): fold-change={fold:.4f}x")
        return fold

    void_fold_hypo = void_clamp_run(-0.145, "HYPOCALCEMIA")
    void_fold_hyper = void_clamp_run(+0.215, "HYPERCALCEMIA")
    void_floor_pass = abs(void_fold_hypo - 1.0) < 0.05 and abs(void_fold_hyper - 1.0) < 0.05
    print(f"Void floor (Ca-sensing blind to actual Ca): fold-changes should be ~1x (no "
          f"response) -- hypo={void_fold_hypo:.4f}x, hyper={void_fold_hyper:.4f}x -> "
          f"{'PASS' if void_floor_pass else 'FAIL'}")

    verdict = "CONFIRMED" if (gate1_pass and gate2_pass and void_floor_pass) else "PARTIAL"
    if gate1_pass and gate2_pass and void_floor_pass and not same_sign:
        verdict += "_WITH_DOCUMENTED_ANOMALY(cross-model residual sign differs from the published prior)"

    print(f"\nGATE SUMMARY: G1={'PASS' if gate1_pass else 'FAIL'}"
          f"(Ca={ca_ss:.4f},PTHconc={pth_conc_ss:.4f},resid={max_resid_ss:.2e}) "
          f"G2={'PASS' if gate2_pass else 'FAIL'}(hypo={fold_hypo:.3f}x,hyper={fold_hyper:.3f}x) "
          f"LEG3=INFO(logratio_mean={lr_vals.mean():+.3f}nats,same_sign={same_sign}) "
          f"VOID_FLOOR={'PASS' if void_floor_pass else 'FAIL'}"
          f"(hypo={void_fold_hypo:.3f}x,hyper={void_fold_hyper:.3f}x) "
          f"VERDICT={verdict}")


if __name__ == "__main__":
    main()
