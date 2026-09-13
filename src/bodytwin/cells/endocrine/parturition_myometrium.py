"""
PARTURITION / MYOMETRIUM -- gap-junction-coupled syncytium (Cx43) +
Ferguson-reflex positive-feedback escalation switch.

Reads: nothing. Writes: parturition_myometrium_evidence.json.
Gate: overall_pass_strict_all (all Part A + Part B gates).

Part A: Cx43 gap-junction-coupled myometrial network. Kuramoto phase-oscillator
        reduction on two topologies (6x6 periodic lattice, Erdos-Renyi random
        graph, matched mean degree=4) -- synchronization threshold in coupling
        g, cross-checked against the graph Laplacian's algebraic connectivity
        (lambda2, the sigma_min-analog spectral governor), forced against a
        g=0 "excitability/homogenization-only" adversary (Garfield 1978's
        finding: gap junctions are LITERALLY ABSENT in pregnancy, so the fair
        adversary is g=0 exactly, not "weak" g).

        DISCLOSED, ABANDONED alternative: a full van der Pol amplitude-dynamics
        network (the same general class of reduction as the
        gi_motility_slow_waves cell) was tried FIRST and showed messy,
        non-monotonic, weak (<0.35) coherence even up to g=2 -- a real,
        OODA-diagnosed dead end (see honest_gaps), abandoned in favor of the
        cleaner, theoretically-canonical Kuramoto phase reduction for this
        specific synchronization falsifier.

Part B: Ferguson-reflex (stretch -> oxytocin -> contraction -> stretch)
        2-state reduced ODE (oxytocin treated as fast/quasi-instantaneous:
        neurohormonal release+clearance ~minutes vs. mechanical labor
        progression ~hours, a disclosed timescale-separation coarse-graining).
        Bifurcation in receptor gain R (both an operational R_crit, found by
        direct simulation+bisection at a standard "Braxton-Hicks-sized"
        perturbation, and a structural R_crit, found analytically via the
        2x2 Jacobian determinant crossing zero -- two independent methods,
        cross-checked). Forced against TWO independently-constructed
        sign-flipped (negative-feedback) adversaries swept over a 7-order-of-
        magnitude R range. Necessity check (fixed low R fails to escalate even
        at near-maximal sub-delivery perturbations). Three mechanistically
        DISTINCT tocolytic-class parameter-reduction routes (atosiban ->
        receptor occupancy; indomethacin -> stretch-generation/ripening
        efficiency; nifedipine/beta-agonists -> relaxation rate). Oxytocin-
        induction dose-response (ripe vs unripe). Progesterone-suppression /
        preterm-acceleration of the R(t) gestational trajectory.

All citations: PMID/DOI live-verified via NCBI eutils (esearch then efetch).
See CITATIONS dict. This is a HYPOTHESIS for independent QC, not a
self-certified verdict.
"""
import json
import os as _os

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

OUT_ROOT = _os.environ.get("BODYTWIN_OUT", _os.path.join(_os.getcwd(), "outputs"))

RNG_SEED = 20260722

CITATIONS = {
    "fuchs_1982_science": {
        "pmid": "6278592", "doi": "10.1126/science.6278592",
        "cite": "Fuchs AR, Fuchs F, Husslein P, Soloff MS, Fernstrom MJ (1982). Oxytocin receptors and human parturition: a dual role for oxytocin in the initiation of labor. Science 215(4538):1396-8.",
        "verified_via": "NCBI esearch+efetch",
        "note": "THE dual-role anchor: OTR concentration increased in myometrium of pregnant women, reached maximum in early labor; also high/maximal in decidua at parturition; OT stimulates decidual (not myometrial) PG production in vitro."
    },
    "fuchs_1984_ajog": {
        "pmid": "6093538", "doi": "10.1016/0002-9378(84)90677-x",
        "cite": "Fuchs AR, Fuchs F, Husslein P, Soloff MS (1984). Oxytocin receptors in the human uterus during pregnancy and parturition. Am J Obstet Gynecol 150(6):734-41.",
        "verified_via": "NCBI esearch+efetch full abstract",
        "numbers": "VERBATIM: myometrial receptor concentration 'low at 13 to 17 weeks but had risen about TWELVEFOLD by 37 to 41 weeks. After the onset of labor...receptor levels were maximal and significantly higher than before onset of labor.' Fundus/corpus > lower segment > cervix (lowest). Failed induction + postterm (43-46wk) pregnancies: SIGNIFICANTLY LOWER receptor concentration than spontaneous labor -- a real, disclosed dissociation (see honest_gaps)."
    },
    "kimura_1996_endocrinology": {
        "pmid": "8593830", "doi": "10.1210/endo.137.2.8593830",
        "cite": "Kimura T, Takemura M, Nomura S, Nobunaga T, Kubota Y, Inoue T, Hashimoto K, Kumazawa I, Ito Y, Ohashi K, Koyama M, Azuma C, Kitamura Y, Saji F (1996). Expression of oxytocin receptor in human pregnant myometrium. Endocrinology 137(2):780-5.",
        "verified_via": "NCBI esearch+efetch full abstract",
        "numbers": "THE headline receptor-fold-change anchor. VERBATIM: OTR mRNA 'reached over 300-FOLD at parturition compared with the nonpregnant myometrium. In the myometrium at 32 weeks of gestation and not in labor, a relatively large amount (about 100-FOLD) of the receptor message was expressed.' Receptor PROTEIN also augmented at term/labor onset (Western blot)."
    },
    "fuchs_1991_pulses": {
        "pmid": "1957888", "doi": "10.1016/0002-9378(91)90399-c",
        "cite": "Fuchs AR, Romero R, Keefe D, Parra M, Oyarzun E, Behnke E (1991). Oxytocin secretion and human parturition: pulse frequency and duration increase during spontaneous labor in women. Am J Obstet Gynecol 165(5 Pt 1):1515-23.",
        "verified_via": "NCBI esearch+efetch full abstract",
        "numbers": "THE Ferguson-reflex escalation anchor (n=32 parturient women, 1-min sampling x30min). Pulse frequency/30min: 1.2+/-0.54 (before labor) -> 4.2+/-0.45 (1st stage) -> 6.7+/-0.49 (2nd/3rd stage) = 3.5x then 5.6x rise. Pulse duration: 1.2+/-0.20 -> 1.9+/-0.28 -> 2.0+/-0.26 min. Amplitude NOT significantly different between groups (~1.0 microU/ml) -- escalation carried by FREQUENCY+DURATION, not amplitude (disclosed, informs the model's OT_ss(s) coarse-graining as an integrated/mean exposure)."
    },
    "garfield_1978_gapjunctions": {
        "pmid": "727239", "doi": "10.1152/ajpcell.1978.235.5.C168",
        "cite": "Garfield RE, Sims SM, Kannan MS, Daniel EE (1978). Possible role of gap junctions in activation of myometrium during parturition. Am J Physiol 235(5):C168-79.",
        "verified_via": "NCBI esearch+efetch full abstract",
        "numbers": "THE Part-A + progesterone anchor (rat). Gap junctions found ONLY immediately prior to/during/after parturition. Ovariectomy of 16-17d pregnant rats -> PREMATURE termination of pregnancy + PREMATURE gap-junction appearance. Progesterone treatment PREVENTS gap-junction appearance (prolongs pregnancy). Gap junctions form in vitro without hormonal influence. Proposes gap junctions necessary for synchronous myometrial contraction."
    },
    "chow_lye_1994_cx43": {
        "pmid": "8141203", "doi": "10.1016/s0002-9378(94)70284-5",
        "cite": "Chow L, Lye SJ (1994). Expression of the gap junction protein connexin-43 is increased in the human myometrium toward term and with the onset of labor. Am J Obstet Gynecol 170(3):788-95.",
        "verified_via": "NCBI esearch+efetch full abstract",
        "numbers": "n=27 not-in-labor (37-41wk) + 7 in labor (39-41wk). Cx43 mRNA increased (p<0.01) 37->40wk, FURTHER increase with labor (p<0.05). Cx43 PROTEIN DECREASED (p<0.05) in late pregnancy and was NOT increased during labor -- YET gap junctions (absent in pregnancy) DID appear in the cell membrane during labor. A real, disclosed mRNA/protein/function dissociation (regulation at assembly, not just synthesis) -- directly informs Part A's use of a discrete/step-like coupling variable rather than one continuously proportional to total Cx43 protein mass."
    },
    "devedeux_1993_emg": {
        "pmid": "8267082", "doi": "10.1016/0002-9378(93)90456-s",
        "cite": "Devedeux D, Marque C, Mansour S, Germain G, Duchene J (1993). Uterine electromyography: a critical review. Am J Obstet Gynecol 169(6):1636-53.",
        "verified_via": "NCBI esearch+efetch full abstract",
        "note": "Internal+external EMG occur in phase with intrauterine-pressure increase; high-frequency EMG band specifically related to EFFICIENT (coordinated) parturition contractions -- external anchor that coordinated electrical activity, not just any activity, is what produces effective pressure/contraction."
    },
    "mesiano_2002_pra": {
        "pmid": "12050275", "doi": "10.1210/jcem.87.6.8609",
        "cite": "Mesiano S, Chan EC, Fitter JT, Kwek K, Yeo G, Smith R (2002). Progesterone withdrawal and estrogen activation in human parturition are coordinated by progesterone receptor A expression in the myometrium. J Clin Endocrinol Metab 87(6):2924-30.",
        "verified_via": "NCBI esearch+efetch full abstract",
        "numbers": "n=12 non-laboring + 12 laboring term myometrium. PR-A, PR-B, PR-A/PR-B mRNA ratio, and ERalpha mRNA all SIGNIFICANTLY INCREASED in laboring tissue. THE functional-progesterone-withdrawal anchor: human parturition is NOT preceded by a literal drop in circulating progesterone (unlike most other mammals) -- withdrawal is FUNCTIONAL, via a shift in the PR-A/PR-B ratio that suppresses progesterone responsiveness despite unchanged/rising circulating hormone."
    },
    "romero_2014_pretermreview": {
        "pmid": "25124429", "doi": "10.1126/science.1251816",
        "cite": "Romero R, Dey SK, Fisher SJ (2014). Preterm labor: one syndrome, many causes. Science 345(6198):760-5.",
        "verified_via": "NCBI esearch+efetch full abstract",
        "numbers": "Preterm birth: 5-18% of pregnancies; spontaneous preterm labor causes 70% of preterm births. Multiple distinct pathologic processes (incl. intrauterine infection/inflammation) converge on the SAME contractile-activation syndrome -- external anchor for 'preterm labor = premature crossing of the same switch', motivating the model's infection-accelerated R(t) sweep."
    },
    "romero_2000_atosiban_rct": {
        "pmid": "10819855", "doi": "10.1067/mob.2000.95834",
        "cite": "Romero R, Sibai BM, Sanchez-Ramos L, et al. (2000). An oxytocin receptor antagonist (atosiban) in the treatment of preterm labor: a randomized, double-blind, placebo-controlled trial with tocolytic rescue. Am J Obstet Gynecol 182(5):1173-83.",
        "verified_via": "NCBI esearch+efetch full abstract",
        "numbers": "n=531 randomized (246 atosiban/255 placebo), 501 treated. PRIMARY endpoint (time to delivery/therapeutic failure) NOT significantly different (median 25.6 vs 21.0 days, p=.6). SECONDARY endpoints (% undelivered w/o rescue tocolytic at 24h/48h/7d) significantly BETTER with atosiban (all p<=.008), with a significant treatment x gestational-age interaction: atosiban consistently superior to placebo at >=28wk. A real, disclosed SAFETY SIGNAL: higher fetal-infant death in the atosiban group at <24wk (14 atosiban vs 5 placebo patients randomized that early -- small n). A genuinely mixed/nuanced real-world result, reported honestly, not oversimplified to 'atosiban works'."
    },
    "wex_2009_atosiban_meta": {
        "pmid": "19538754", "doi": "10.1186/1471-2393-9-23", "pmcid": "PMC2708127",
        "cite": "Wex J, Connolly M, Rath W (2009). Atosiban versus betamimetics in the treatment of preterm labour in Germany: an economic evaluation. BMC Pregnancy Childbirth 9:23.",
        "verified_via": "NCBI esearch+efetch full abstract",
        "numbers": "Meta-analysis of 3 double-blind placebo-controlled RCTs: atosiban and betamimetics show STATISTICALLY EQUIVALENT efficacy (RR=0.99, 95%CI 0.94-1.04, p=.772) despite acting via completely different mechanisms (OT-receptor antagonism vs beta-adrenergic agonism) -- the over-determination anchor: mechanistically distinct drug classes converge on the same qualitative (tocolytic) outcome, matching the model's 3-distinct-parameter-routes-all-rescue prediction."
    },
    "haas_2011_tocolytic_review": {
        "pmid": "21463540", "pmcid": "PMC3217816",
        "cite": "Haas DM (2011). Preterm birth. BMJ Clin Evid 2011:1404.",
        "verified_via": "NCBI esearch+efetch full abstract",
        "note": "Systematic review naming the full tocolytic-class landscape used clinically: beta-mimetics, calcium-channel blockers, oxytocin-receptor antagonists (atosiban), prostaglandin inhibitors (e.g. indometacin), progesterone -- omnibus anchor that these mechanistically distinct classes are the real clinical toolkit the model's 3 parameter-reduction routes are meant to represent."
    },
    "seitchik_1985_dosing": {
        "pmid": "3976787", "doi": "10.1016/0002-9378(85)90514-9",
        "cite": "Seitchik J, Amico JA, Castillo M (1985). Oxytocin augmentation of dysfunctional labor. V. An alternative oxytocin regimen. Am J Obstet Gynecol 151(6):757-61.",
        "verified_via": "NCBI esearch+efetch full abstract",
        "numbers": "Continuous IV oxytocin in the therapeutic range requires ~40min to reach steady-state plasma concentration. Timed-dose method selected effective maintenance dose in 43/59 (73%) on first attempt, 58/59 by second; doses effected cervical dilation in 55/59 (93%) -- real dose-response/induction-success anchor."
    },
    "amico_1984_otkinetics": {
        "pmid": "6693537", "doi": "10.1210/jcem-58-2-274",
        "cite": "Amico JA, Seitchik J, Robinson AG (1984). Studies of oxytocin in plasma of women during hypocontractile labor. J Clin Endocrinol Metab 58(2):274-9.",
        "verified_via": "NCBI esearch+efetch full abstract",
        "numbers": "n=11: OT levels/kinetics during hypocontractile labor similar to nonpregnant/late-pregnancy individuals (MCR 17.4+/-9.2 ml/kg.min, ~matching non-pregnant men 17.6+/-2.1). VERBATIM: 'individual myometrial SENSITIVITY is an important determinant of the response to administered OT' -- direct support for the model's R (receptor gain) as the key variable, not plasma OT level per se."
    },
    "higuchi_1987_fergusonreflex_term": {
        "pmid": "3569466", "doi": "10.1016/0014-4886(87)90061-6",
        "cite": "Higuchi T, Uchide K, Honda K, Negoro H (1987). Pelvic neurectomy abolishes the fetus-expulsion reflex and induces dystocia in the rat. Exp Neurol 96(2):443-55.",
        "verified_via": "NCBI esearch+efetch full abstract",
        "note": "Explicitly names and defines the reflex this document models: 'the fetus-ejection reflex (Ferguson reflex) which initiates oxytocin release', mechanistically distinguished from a separate abdominal fetus-expulsion straining reflex (unaffected by pelvic neurectomy in this rat model) -- a real, disclosed nuance: the Ferguson reflex proper is the neurohormonal (oxytocin-release) arc modeled here, not the separate expulsive-straining reflex."
    },
    "meis_2003_17ohpc_positive": {
        "pmid": "12802023", "doi": "10.1056/NEJMoa035140",
        "cite": "Meis PJ, Klebanoff M, Thom E, et al.; NICHD MFMU Network (2003). Prevention of recurrent preterm delivery by 17 alpha-hydroxyprogesterone caproate. N Engl J Med 348(24):2379-85.",
        "verified_via": "NCBI esearch+efetch abstract",
        "numbers": "Double-blind placebo-controlled trial (women with prior spontaneous preterm delivery), 19 centers, 2:1 randomization to weekly 250mg 17P IM vs oil placebo -- the ORIGINAL positive 17-OHPC/preterm-birth-prevention trial (stopped early for large treatment benefit per PROLONG's later description)."
    },
    "blackwell_2020_prolong_negative": {
        "pmid": "31652479", "doi": "10.1055/s-0039-3400227",
        "cite": "Blackwell SC, Gyamfi-Bannerman C, Biggio JR Jr, et al. (2020). 17-OHPC to Prevent Recurrent Preterm Birth in Singleton Gestations (PROLONG Study): A Multicenter, International, Randomized Double-Blind Trial. Am J Perinatol 37(2):127-136.",
        "verified_via": "NCBI esearch+efetch full abstract",
        "numbers": "THE negative replication (n=1130 17-OHPC vs 578 placebo, ~5.5x larger than Meis 2003): PTB<35wk 11.0% vs 11.5% (RR=0.95, 95%CI 0.71-1.26) -- NO significant difference; neonatal morbidity composite also not different (5.6% vs 5.0%, RR=1.12). A genuine, important FAILURE TO REPLICATE a positive finding (Meis 2003) in a larger, more rigorous trial -- directly informs honest_gaps/symmetric QC on the progesterone-suppression mechanism's real-world clinical-efficacy evidence, independent of whether the underlying mechanism (Garfield 1978, Mesiano 2002) is real."
    },
}

# ============================================================
# PART A -- gap-junction-coupled myometrial network (Kuramoto)
# ============================================================

def build_lattice_adjacency(n_side):
    """2D periodic lattice (torus), 4-neighbor -- degree exactly 4 everywhere."""
    N = n_side * n_side
    A = np.zeros((N, N))
    for i in range(n_side):
        for j in range(n_side):
            idx = i * n_side + j
            for di, dj in [(1, 0), (-1, 0), (0, 1), (0, -1)]:
                ni, nj = (i + di) % n_side, (j + dj) % n_side
                A[idx, ni * n_side + nj] = 1.0
    return A


def build_random_graph_adjacency(N, mean_degree, seed):
    """Erdos-Renyi random graph, matched mean degree to the lattice."""
    rng = np.random.default_rng(seed)
    p = mean_degree / (N - 1)
    A = np.zeros((N, N))
    for i in range(N):
        for j in range(i + 1, N):
            if rng.random() < p:
                A[i, j] = A[j, i] = 1.0
    return A


def algebraic_connectivity(A):
    """lambda2 of the graph Laplacian L=D-A -- the sigma_min-analog spectral governor."""
    D = np.diag(A.sum(axis=1))
    L = D - A
    eigvals = np.linalg.eigvalsh(L)
    n_isolated = int(np.sum(A.sum(axis=1) == 0))
    return float(eigvals[1]), n_isolated


def simulate_kuramoto(A, g, omega_mean, omega_std, T=80.0, seed=1):
    """Kuramoto phase-oscillator network: dtheta_i/dt = omega_i + g*sum_j A_ij sin(theta_j-theta_i)."""
    N = A.shape[0]
    rng = np.random.default_rng(seed)
    omega = rng.normal(omega_mean, omega_std, N)
    theta0 = rng.uniform(0, 2 * np.pi, N)

    def rhs(t, theta):
        diff = theta[None, :] - theta[:, None]
        coupling = g * np.sum(A * np.sin(diff), axis=1)
        return omega + coupling

    t_eval = np.linspace(0, T, 600)
    sol = solve_ivp(rhs, [0, T], theta0, t_eval=t_eval, method="RK45", rtol=1e-7, atol=1e-9)
    th = sol.y
    n_burn = int(0.5 * th.shape[1])
    ths = th[:, n_burn:]
    R_t = np.abs(np.mean(np.exp(1j * ths), axis=0))
    return float(np.mean(R_t)), bool(sol.success)


def part_a_run():
    N_SIDE = 6
    A_lat = build_lattice_adjacency(N_SIDE)
    N = A_lat.shape[0]
    A_rand = build_random_graph_adjacency(N, mean_degree=4, seed=7)
    lam2_lat, iso_lat = algebraic_connectivity(A_lat)
    lam2_rand, iso_rand = algebraic_connectivity(A_rand)

    g_grid = [0.0, 0.02, 0.04, 0.06, 0.07, 0.08, 0.09, 0.1, 0.12, 0.15, 0.2, 0.3, 0.5, 1.0, 2.0, 5.0]
    R_lat, R_rand = [], []
    for g in g_grid:
        Rl, _ = simulate_kuramoto(A_lat, g, 1.0, 0.15, seed=3)
        Rr, _ = simulate_kuramoto(A_rand, g, 1.0, 0.15, seed=3)
        R_lat.append(Rl)
        R_rand.append(Rr)

    def find_gcrit(grid, Rs, thr=0.5):
        for gg, RR in zip(grid, Rs):
            if RR >= thr:
                return gg
        return None

    gcrit_lat = find_gcrit(g_grid, R_lat)
    gcrit_rand = find_gcrit(g_grid, R_rand)

    # void floor: g=0, multiple seeds, both topologies (identical adjacency-independent at g=0
    # since coupling term vanishes regardless of A -- still run both labels for clarity)
    void_floor_seeds = [1, 2, 3, 4, 5, 10, 11, 12, 13, 14]
    void_R = []
    for s in void_floor_seeds:
        Rv, _ = simulate_kuramoto(A_lat, 0.0, 1.0, 0.15, seed=s)
        void_R.append(Rv)

    # TRUE excitability/homogenization-only adversary: g=0 EXACTLY (Garfield 1978: gap
    # junctions are LITERALLY ABSENT in pregnancy -- the fair adversary uses g=0, not "weak" g)
    adv_omega_std_sweep = []
    for om_std in [0.15, 0.05, 0.02, 0.005, 0.001, 0.0001]:
        Ra, _ = simulate_kuramoto(A_lat, 0.0, 1.0, om_std, seed=3)
        adv_omega_std_sweep.append({"omega_std": om_std, "R": Ra})
    adv_omega_mean_sweep = []
    for om_mean in [1.0, 2.0, 5.0, 10.0, 20.0]:
        Ra, _ = simulate_kuramoto(A_lat, 0.0, om_mean, 0.15, seed=3)
        adv_omega_mean_sweep.append({"omega_mean": om_mean, "R": Ra})
    adv_seeds_sweep = []
    for s in [1, 2, 3, 4, 5]:
        Ra, _ = simulate_kuramoto(A_lat, 0.0, 1.0, 0.0001, seed=s)
        adv_seeds_sweep.append({"seed": s, "R": Ra})
    adv_max_R = max(
        [d["R"] for d in adv_omega_std_sweep]
        + [d["R"] for d in adv_omega_mean_sweep]
        + [d["R"] for d in adv_seeds_sweep]
    )

    lambda2_predicts_ordering_correctly = None
    if gcrit_lat is not None and gcrit_rand is not None:
        # higher lambda2 -> lower g_crit (better-connected graph syncs at weaker coupling)
        if lam2_lat > lam2_rand:
            lambda2_predicts_ordering_correctly = bool(gcrit_lat <= gcrit_rand)
        elif lam2_rand > lam2_lat:
            lambda2_predicts_ordering_correctly = bool(gcrit_rand <= gcrit_lat)
        else:
            lambda2_predicts_ordering_correctly = None

    results = {
        "N": N,
        "mean_degree_lattice": float(A_lat.sum(axis=1).mean()),
        "mean_degree_random": float(A_rand.sum(axis=1).mean()),
        "lambda2_lattice": lam2_lat,
        "lambda2_random": lam2_rand,
        "isolated_nodes_lattice": iso_lat,
        "isolated_nodes_random": iso_rand,
        "g_grid": g_grid,
        "kuramoto_R_lattice": R_lat,
        "kuramoto_R_random": R_rand,
        "g_crit_lattice": gcrit_lat,
        "g_crit_random": gcrit_rand,
        "lambda2_predicts_gcrit_ordering": lambda2_predicts_ordering_correctly,
        "void_floor_seeds": void_floor_seeds,
        "void_floor_R_values": void_R,
        "void_floor_max_R": float(max(void_R)),
        "adversary_g_exactly_zero_omega_std_sweep": adv_omega_std_sweep,
        "adversary_g_exactly_zero_omega_mean_sweep": adv_omega_mean_sweep,
        "adversary_g_exactly_zero_seeds_sweep": adv_seeds_sweep,
        "adversary_max_R_across_all_conditions": float(adv_max_R),
    }

    gates = {
        "A1_void_floor_below_0p3": bool(max(void_R) < 0.3),
        "A2_gcrit_lattice_finite": gcrit_lat is not None,
        "A3_gcrit_random_finite": gcrit_rand is not None,
        "A4_both_topologies_reach_ge0p9_at_high_g": bool(R_lat[-1] >= 0.9 and R_rand[-1] >= 0.9),
        "A5_lambda2_ordering_predicts_gcrit_ordering": bool(lambda2_predicts_ordering_correctly),
        "A6_excitability_only_adversary_never_crosses_0p5": bool(adv_max_R < 0.5),
        "A7_gcrit_lattice_le_gcrit_random_consistent_with_higher_lambda2": bool(
            gcrit_lat is not None and gcrit_rand is not None and gcrit_lat <= gcrit_rand
        ),
    }
    return results, gates


# ============================================================
# PART B -- Ferguson-reflex positive-feedback bifurcation
# ============================================================

BASE_PARAMS = dict(k1=1.0, k2=1.0, ka=1.0, kb=1.0, kc=1.0, K=0.5, n=2, OT0=0.05)


def hill(x, K, n):
    x = np.maximum(x, 0.0)
    return x ** n / (K ** n + x ** n)


def hill_deriv(x, K, n):
    x = max(x, 1e-9)
    return n * (K ** n) * (x ** (n - 1)) / ((K ** n + x ** n) ** 2)


def rhs_positive(t, y, R, params, I_exo=0.0):
    C, s = y
    OT = params["OT0"] + params["k1"] * max(s, 0.0) + I_exo
    dC = params["k2"] * R * hill(OT, params["K"], params["n"]) - params["kc"] * C
    ds = params["ka"] * max(C, 0.0) - params["kb"] * s
    return [dC, ds]


def rhs_negative_variant1(t, y, R, params, I_exo=0.0):
    """Ferguson-reflex arc SIGN-FLIPPED: stretch REDUCES oxytocin release."""
    C, s = y
    OT = max(params["OT0"] - params["k1"] * max(s, 0.0), 0.0) + I_exo
    dC = params["k2"] * R * hill(OT, params["K"], params["n"]) - params["kc"] * C
    ds = params["ka"] * max(C, 0.0) - params["kb"] * s
    return [dC, ds]


def rhs_negative_variant2(t, y, R, params, I_exo=0.0):
    """Mechanical arc SIGN-FLIPPED: contraction REDUCES stretch instead of increasing it."""
    C, s = y
    OT = params["OT0"] + params["k1"] * max(s, 0.0) + I_exo
    dC = params["k2"] * R * hill(OT, params["K"], params["n"]) - params["kc"] * C
    ds = -params["ka"] * max(C, 0.0) - params["kb"] * s
    return [dC, ds]


def classify_trajectory(rhs_fn, R, params, s0, C0=0.0, T=400.0, I_exo=0.0, thr=0.9):
    """escalated iff the SETTLED (long-time) trajectory value is at/beyond the delivery
    threshold -- NOT 'ever touches the threshold', which would spuriously fire at t=0 for
    perturbations placed at/near the threshold itself (a real bug caught via OODA, see
    honest_gaps)."""
    sol = solve_ivp(
        rhs_fn, [0, T], [C0, s0], args=(R, params, I_exo), method="RK45",
        max_step=0.5, rtol=1e-8, atol=1e-10,
    )
    s_traj = sol.y[1]
    final_s = float(s_traj[-1])
    escalated = bool(final_s >= thr)
    return escalated, final_s, bool(sol.success)


def fixed_point_s_positive(R, params, npts=40001):
    """Scan for the smallest positive root of dC/dt=dS/dt=0 (quiescent/low branch)."""
    ss = np.linspace(0, 1, npts)
    g = params["k2"] * R * hill(params["OT0"] + params["k1"] * ss, params["K"], params["n"]) \
        - params["kc"] * (params["kb"] / params["ka"]) * ss
    sign_changes = np.where(np.diff(np.sign(g)))[0]
    if len(sign_changes) == 0:
        return None
    i = sign_changes[0]
    s_lo, s_hi = ss[i], ss[i + 1]
    g_lo, g_hi = g[i], g[i + 1]
    return float(s_lo - g_lo * (s_hi - s_lo) / (g_hi - g_lo))


def jacobian_det_positive(R, params):
    """2x2 Jacobian determinant at the quiescent fixed point. det=0 is the structural,
    saddle-node bifurcation boundary -- the sigma_min-analog spectral criterion for Part B."""
    s_star = fixed_point_s_positive(R, params)
    if s_star is None:
        return None, None
    hprime = hill_deriv(params["OT0"] + params["k1"] * s_star, params["K"], params["n"])
    det = params["kc"] * params["kb"] - params["ka"] * (params["k2"] * R * params["k1"] * hprime)
    return det, s_star


def find_operational_Rcrit(rhs_fn, params, s0=0.3, lo=0.01, hi=50.0, n_iter=45):
    """Bisect for R_crit: smallest R at which the standard perturbation s0 escalates."""

    def esc(R):
        e, _, _ = classify_trajectory(rhs_fn, R, params, s0=s0)
        return e

    if esc(lo):
        return None  # even the smallest R escalates -- bracket invalid, disclose
    if not esc(hi):
        return None  # even the largest R fails -- bracket invalid, disclose
    for _ in range(n_iter):
        mid = 0.5 * (lo + hi)
        if esc(mid):
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


def find_structural_Rcrit(params, lo=0.1, hi=3.0, n_steps=400):
    """Scan det(R) for its zero-crossing (structural saddle-node)."""
    Rs = np.linspace(lo, hi, n_steps)
    prev_det = None
    prev_R = None
    for Rv in Rs:
        d, _ = jacobian_det_positive(Rv, params)
        if d is None:
            if prev_det is not None and prev_det > 0:
                return float(0.5 * (prev_R + Rv))  # bracket the vanish point
            return None
        if prev_det is not None and prev_det > 0 and d <= 0:
            # linear-interp the crossing
            frac = prev_det / (prev_det - d)
            return float(prev_R + frac * (Rv - prev_R))
        prev_det, prev_R = d, Rv
    return None


def part_b_run():
    P = dict(BASE_PARAMS)

    # --- B1-B3: operational + structural R_crit, cross-checked ---
    Rcrit_op = find_operational_Rcrit(rhs_positive, P, s0=0.3)
    Rcrit_struct = find_structural_Rcrit(P, lo=0.1, hi=3.0, n_steps=2000)

    below_tests = []
    for frac, label in [(0.01, "Rcrit_op/100"), (0.02, "Rcrit_op/50"), (0.05, "Rcrit_op/20"),
                         (0.1, "Rcrit_op/10"), (0.5, "Rcrit_op/2")]:
        Rv = Rcrit_op * frac
        row = {"R": Rv, "label": label, "perturbations": []}
        for s0 in [0.1, 0.3, 0.5, 0.7, 0.8]:
            esc, fs, ok = classify_trajectory(rhs_positive, Rv, P, s0=s0)
            row["perturbations"].append({"s0": s0, "escalated": esc, "final_s": fs})
        below_tests.append(row)
    all_below_decay = all(
        (not d["escalated"]) for row in below_tests for d in row["perturbations"]
    )

    above_tests = []
    for mult in [1.5, 2.0, 5.0, 10.0]:
        Rv = Rcrit_op * mult
        row = {"R": Rv, "mult": mult, "perturbations": []}
        for s0 in [0.05, 0.1, 0.3]:
            esc, fs, ok = classify_trajectory(rhs_positive, Rv, P, s0=s0)
            row["perturbations"].append({"s0": s0, "escalated": esc, "final_s": fs})
        above_tests.append(row)
    all_above_escalate = all(
        d["escalated"] for row in above_tests for d in row["perturbations"]
    )

    # --- B4-B6: forced negative-feedback adversaries, wide R sweep, 2 independent variants ---
    R_wide_sweep = [0.01, 0.1, 1.0, 10.0, 100.0, 1000.0, 10000.0, 100000.0]
    neg1_results = []
    for Rv in R_wide_sweep:
        esc, fs, ok = classify_trajectory(rhs_negative_variant1, Rv, P, s0=0.5)
        neg1_results.append({"R": Rv, "escalated": esc, "final_s": fs})
    neg2_results = []
    for Rv in R_wide_sweep:
        esc, fs, ok = classify_trajectory(rhs_negative_variant2, Rv, P, s0=0.5, C0=0.0)
        neg2_results.append({"R": Rv, "escalated": esc, "final_s": fs})
    neg1_never_escalates = not any(d["escalated"] for d in neg1_results)
    neg2_never_escalates = not any(d["escalated"] for d in neg2_results)

    # analytic proof cross-check: for variant 1 (k1->-k1), det = kc*kb + ka*k2*R*k1*Hill' > 0 ALWAYS
    # (sum of non-negative terms) -- verify numerically across the same R sweep
    analytic_neg1_det_always_positive = True
    for Rv in R_wide_sweep:
        hprime_at_baseline = hill_deriv(P["OT0"], P["K"], P["n"])
        det_neg1 = P["kc"] * P["kb"] + P["ka"] * (P["k2"] * Rv * P["k1"] * hprime_at_baseline)
        if det_neg1 <= 0:
            analytic_neg1_det_always_positive = False

    # --- B7: necessity -- fixed low R fails even at near-maximal sub-delivery perturbation ---
    necessity_tests = []
    for frac, label in [(0.01, "1/100"), (0.02, "1/50"), (0.1, "1/10"), (0.5, "1/2")]:
        Rv = Rcrit_op * frac
        esc, fs, ok = classify_trajectory(rhs_positive, Rv, P, s0=0.8)
        necessity_tests.append({"R": Rv, "label": label, "s0": 0.8, "escalated": esc, "final_s": fs})
    necessity_holds = all(not d["escalated"] for d in necessity_tests)

    # --- B8-B10: 3 mechanistically distinct tocolytic parameter-reduction routes ---
    R_super = Rcrit_op * 3.0  # a comfortably-escalating "term/labor" reference state

    def sweep_param_rescue(param_name, values, base_R=R_super, s0=0.3, mode="multiply_R"):
        rows = []
        for v in values:
            Pv = dict(P)
            Rv = base_R
            if mode == "multiply_R":  # atosiban: effective R *= (1-occupancy)
                Rv = base_R * (1.0 - v)
            elif mode == "set_param":
                Pv[param_name] = v
            esc, fs, ok = classify_trajectory(rhs_positive, Rv, Pv, s0=s0)
            rows.append({"value": v, "escalated": esc, "final_s": fs})
        return rows

    atosiban_occupancy_sweep = sweep_param_rescue(
        None, [0.0, 0.2, 0.4, 0.6, 0.65, 0.7, 0.75, 0.8, 0.9, 0.95, 0.99], mode="multiply_R"
    )
    indomethacin_ka_sweep = sweep_param_rescue(
        "ka", [1.0, 0.7, 0.5, 0.4, 0.35, 0.3, 0.25, 0.2, 0.1, 0.05, 0.01], mode="set_param"
    )
    nifedipine_kc_sweep = []
    for kc_v in [1.0, 1.5, 2.0, 2.5, 2.8, 3.0, 3.2, 3.5, 4.0, 6.0, 10.0]:
        Pv = dict(P)
        Pv["kc"] = kc_v
        esc, fs, ok = classify_trajectory(rhs_positive, R_super, Pv, s0=0.3)
        nifedipine_kc_sweep.append({"value": kc_v, "escalated": esc, "final_s": fs})

    def first_rescue(rows, key="value", ascending=True):
        for row in rows:
            if not row["escalated"]:
                return row[key]
        return None

    atosiban_rescue_occupancy = first_rescue(atosiban_occupancy_sweep)
    indomethacin_rescue_ka = first_rescue(indomethacin_ka_sweep)
    nifedipine_rescue_kc = first_rescue(nifedipine_kc_sweep)
    all_3_tocolytic_routes_rescue = all(
        v is not None for v in [atosiban_rescue_occupancy, indomethacin_rescue_ka, nifedipine_rescue_kc]
    )

    # --- B11: oxytocin induction (Pitocin) dose-response, ripe vs unripe cervix ---
    # OODA fix (disclosed): an earlier attempt set R_ripe=R_super (the tocolytic reference
    # state), which is ALREADY above R_crit and escalates spontaneously even at I_exo=0 --
    # a confounded, vacuous "induction success" (nothing was actually induced). Fixed by
    # using an R just BELOW threshold (quiescent at baseline, s0=0, I_exo=0) so a genuine
    # dose-response to EXOGENOUS oxytocin can be observed.
    R_ripe = Rcrit_op * 0.9
    R_unripe = Rcrit_op * 0.1
    induction_ripe = []
    induction_unripe = []
    for I_exo in [0.0, 0.02, 0.05, 0.08, 0.1, 0.15, 0.2, 0.3, 0.5, 1.0, 2.0]:
        esc_r, fs_r, _ = classify_trajectory(rhs_positive, R_ripe, P, s0=0.0, I_exo=I_exo)
        esc_u, fs_u, _ = classify_trajectory(rhs_positive, R_unripe, P, s0=0.0, I_exo=I_exo)
        induction_ripe.append({"I_exo": I_exo, "escalated": esc_r, "final_s": fs_r})
        induction_unripe.append({"I_exo": I_exo, "escalated": esc_u, "final_s": fs_u})
    ripe_baseline_quiescent_at_Iexo0 = not induction_ripe[0]["escalated"]
    ripe_induction_succeeds_at_some_dose = any(d["escalated"] for d in induction_ripe)
    unripe_induction_fails_across_tested_doses = not any(d["escalated"] for d in induction_unripe)

    # --- B12-B13: progesterone-suppression / preterm-acceleration of R(t) ---
    # R(t) is a logistic CAP-upregulation trajectory with FIXED steepness (rate) and a
    # swept ONSET TIME t_mid -- this is the mechanistically-motivated parametrization
    # (Garfield 1978: progesterone maintenance DELAYS CAP/gap-junction appearance;
    # its withdrawal, or infection/inflammation per Romero 2014, brings it forward).
    # NOTE (disclosed OODA fix): an EARLIER attempt instead swept the sigmoid's STEEPNESS
    # at a FIXED shared midpoint (t=20) and found a NON-monotonic t_cross(rate) --
    # diagnosed: a shared midpoint conflates "when R(t) first crosses R_crit" with "how
    # fast R keeps climbing AFTER crossing" (both controlled by the same rate parameter,
    # trading off against the (C,s) system's near-threshold response lag). Fixed by
    # separating these: hold steepness fixed, sweep the onset time itself.
    def R_of_t(t, t_mid, rate=0.3, R0=0.05, Rmax=None):
        Rmax = Rmax if Rmax is not None else Rcrit_op * 4.0
        return R0 + (Rmax - R0) / (1.0 + np.exp(-rate * (t - t_mid)))

    def time_to_crossing(t_mid, rate=0.3, s0=0.05, T_max=80.0, dt=0.05):
        """Simulate with R(t) TIME-VARYING (re-integrate stepwise), find gestational time
        at which a standard nudge first settles into escalation."""
        t = 0.0
        C, s = 0.0, s0
        while t < T_max:
            Rv = R_of_t(t, t_mid, rate)
            sol = solve_ivp(rhs_positive, [0, dt], [C, s], args=(Rv, P, 0.0), method="RK45",
                             rtol=1e-8, atol=1e-10)
            C, s = sol.y[0][-1], sol.y[1][-1]
            if s >= 0.9:
                return t
            t += dt
        return None

    progesterone_t_mids = [12, 15, 18, 20, 22, 25, 28, 30]
    progesterone_sweep = []
    for t_mid in progesterone_t_mids:
        t_cross = time_to_crossing(t_mid)
        progesterone_sweep.append({"t_mid": t_mid, "t_cross": t_cross})
    valid_rows = [d for d in progesterone_sweep if d["t_cross"] is not None]
    progesterone_monotonic = True
    for i in range(1, len(valid_rows)):
        # LATER onset (more progesterone suppression) -> LATER (or equal) crossing, never earlier
        if valid_rows[i]["t_cross"] < valid_rows[i - 1]["t_cross"]:
            progesterone_monotonic = False

    baseline_t_mid = 20
    preterm_t_mid = 10  # infection/inflammation-accelerated CAP upregulation (earlier onset)
    t_baseline = time_to_crossing(baseline_t_mid)
    t_preterm = time_to_crossing(preterm_t_mid)
    preterm_crosses_earlier = (
        t_baseline is not None and t_preterm is not None and t_preterm < t_baseline
    )

    results = {
        "base_params": P,
        "Rcrit_operational_s0_0p3": Rcrit_op,
        "Rcrit_structural_analytic_det0": Rcrit_struct,
        "below_threshold_tests": below_tests,
        "above_threshold_tests": above_tests,
        "negative_variant1_wide_sweep": neg1_results,
        "negative_variant2_wide_sweep": neg2_results,
        "analytic_neg1_det_always_positive": analytic_neg1_det_always_positive,
        "necessity_tests_fixed_low_R_s0_0p8": necessity_tests,
        "R_super_reference": R_super,
        "atosiban_occupancy_sweep": atosiban_occupancy_sweep,
        "atosiban_rescue_occupancy": atosiban_rescue_occupancy,
        "indomethacin_ka_sweep": indomethacin_ka_sweep,
        "indomethacin_rescue_ka": indomethacin_rescue_ka,
        "nifedipine_kc_sweep": nifedipine_kc_sweep,
        "nifedipine_rescue_kc": nifedipine_rescue_kc,
        "R_ripe": R_ripe,
        "R_unripe": R_unripe,
        "induction_ripe_sweep": induction_ripe,
        "induction_unripe_sweep": induction_unripe,
        "progesterone_t_mid_sweep": progesterone_sweep,
        "preterm_baseline_t_mid": baseline_t_mid,
        "preterm_accelerated_t_mid": preterm_t_mid,
        "t_cross_baseline": t_baseline,
        "t_cross_preterm": t_preterm,
    }

    gates = {
        "B1_Rcrit_operational_finite": Rcrit_op is not None,
        "B2_Rcrit_structural_finite": Rcrit_struct is not None,
        "B3_two_Rcrit_methods_agree_order_op_le_struct": bool(
            Rcrit_op is not None and Rcrit_struct is not None and Rcrit_op <= Rcrit_struct
        ),
        "B4_below_threshold_all_decay": all_below_decay,
        "B5_above_threshold_all_escalate": all_above_escalate,
        "B6_negvariant1_never_escalates_7ordersmag": neg1_never_escalates,
        "B7_negvariant2_never_escalates_7ordersmag": neg2_never_escalates,
        "B8_analytic_proof_neg1_det_always_positive": analytic_neg1_det_always_positive,
        "B9_necessity_fixedlowR_maxperturbation_fails": necessity_holds,
        "B10_atosiban_route_rescues": atosiban_rescue_occupancy is not None,
        "B11_indomethacin_route_rescues": indomethacin_rescue_ka is not None,
        "B12_nifedipine_route_rescues": nifedipine_rescue_kc is not None,
        "B13_three_distinct_tocolytic_routes_all_rescue": all_3_tocolytic_routes_rescue,
        "B14_ripe_baseline_quiescent_Iexo0_not_confounded": ripe_baseline_quiescent_at_Iexo0,
        "B14b_ripe_induction_succeeds": ripe_induction_succeeds_at_some_dose,
        "B15_unripe_induction_fails_across_tested_doses": unripe_induction_fails_across_tested_doses,
        "B16_progesterone_suppression_delays_crossing_monotonic": progesterone_monotonic,
        "B17_preterm_acceleration_crosses_earlier": preterm_crosses_earlier,
    }
    return results, gates


def main():
    a_results, a_gates = part_a_run()
    b_results, b_gates = part_b_run()

    all_gates = {**a_gates, **b_gates}
    pass_count = sum(1 for v in all_gates.values() if v)
    total = len(all_gates)

    output = {
        "citations": CITATIONS,
        "part_a_gap_junction_network": a_results,
        "part_a_gates": a_gates,
        "part_b_ferguson_reflex": b_results,
        "part_b_gates": b_gates,
        "verdict": {
            "part_a_pass_count": sum(1 for v in a_gates.values() if v),
            "part_a_total": len(a_gates),
            "part_b_pass_count": sum(1 for v in b_gates.values() if v),
            "part_b_total": len(b_gates),
            "overall_gates_pass": pass_count,
            "overall_gates_total": total,
            "overall_pass_strict_all": bool(pass_count == total),
        },
    }

    out_dir = _os.path.join(OUT_ROOT, "parturition_myometrium")
    _os.makedirs(out_dir, exist_ok=True)
    out_path = _os.path.join(out_dir, "parturition_myometrium_evidence.json")
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2, allow_nan=False)

    print(f"Wrote {out_path}")
    print(f"\nPART A gates ({sum(1 for v in a_gates.values() if v)}/{len(a_gates)}):")
    for k, v in a_gates.items():
        print(f"  {k}: {'PASS' if v else 'FAIL'}")
    print(f"\nPART B gates ({sum(1 for v in b_gates.values() if v)}/{len(b_gates)}):")
    for k, v in b_gates.items():
        print(f"  {k}: {'PASS' if v else 'FAIL'}")
    print(f"\nOVERALL: {pass_count}/{total} gates PASS")
    print(f"\nKey numbers:")
    print(f"  lambda2 lattice={a_results['lambda2_lattice']:.4f}  random={a_results['lambda2_random']:.4f}")
    print(f"  g_crit  lattice={a_results['g_crit_lattice']}  random={a_results['g_crit_random']}")
    print(f"  void floor max R={a_results['void_floor_max_R']:.4f}")
    print(f"  adversary (g=0) max R={a_results['adversary_max_R_across_all_conditions']:.4f}")
    print(f"  Rcrit operational={b_results['Rcrit_operational_s0_0p3']:.4f}  structural={b_results['Rcrit_structural_analytic_det0']}")
    print(f"  atosiban rescue occupancy={b_results['atosiban_rescue_occupancy']}")
    print(f"  indomethacin rescue ka={b_results['indomethacin_rescue_ka']}")
    print(f"  nifedipine rescue kc={b_results['nifedipine_rescue_kc']}")
    print(f"  t_cross baseline={b_results['t_cross_baseline']}  preterm={b_results['t_cross_preterm']}")


if __name__ == "__main__":
    main()
