# LANE_CAPACITY_FROM_TMAO

Bind the transporter capacity with the substance whose extraction ratio allows it, and check against a thermodynamic ceiling. Results folder `results/LANE_CAPACITY_FROM_TMAO/`. Form A.

## Why this lane, and why right now

A swarm job today delivered a certificate that redirects the entire capacity chain: **No renal clearance measurement can provide an upper limit on transporter turnover, except for substances whose extraction ratio is sufficiently far below one.**

- **Metformin, E = 0,698: separation fails.** Capacity range is unbounded upwards with probability 0,113.
- **TMAO, E = 0,368: separation succeeds.** The interval is finite.
- **Kreatinin, E = 0,200: lyckas.**

The entire 243× gap was built on metformin, thus the only one of three substances whose clearance is demonstrably unable to bind capacity. That explains why the chain gave unreasonable demands: it asked for a quantity that is saturated with flow.

At the same time, there is a held-out TMAO-ankare that no capacity chain has used: urine TMAO-clearance **219 ± 78 mL/min** against creatinine **119 ± 21** and urea **55 ± 14** in SAME individuals (Hai et al., PLoS ONE 2015;10:e0143731). The ratio 219/119 = 1,840 is net secretion that filtration cannot carry, and it is measured paired — which is unusual and valuable, because pooled-versus-paired two of our conclusions today.

## Do this

1. **Calculate the capacity limit from TMAO**, not from metformin. Using the paired measurement, propagate the dispersion and specify its nature per source. Supply a FINITE range for the carrier's turnover or for the capacity per area, and explicitly say what quantity you are committing to.
2. **The regime must be declared before the figures are compared.** This refuted our biggest claim today: a demand was calculated as flow = k_cat × N, i.e. saturated regime, and was compared against a measured value in linear regime where clearance = V_max/K_m. The difference is K_m/C and was a factor over one hundred. For each term, indicate whether C is far below, near, or above K_m, using the form required by the regime.
3. **OIndependent control that does not go via flow at all:** a concentration ratio at steady state is set by thermodynamics and not by delivery. Derive the intracellular accumulation allowed by membrane potential and pH gradient for a cation, and compare with published accumulation in cells expressing the transporter alone and together with the efflux transporter, respectively. That limit holds regardless of turnover, which is the whole point of having it.
4. **Test if the two paths are compatible.** If the flow-separated TMAO limit and the thermodynamic ceiling give incompatible intervals, at least one assumption is wrong — name which one and what would determine. Two independent paths that happen to match is a stronger result than one path that does.

## Strongest control and falsifier

- **Control:** the metformin chain as it stands, i.e. the one that gave 243× and as the certificate shows cannot bind the capacity. The profit should be a finite interval where it gave an unlimited.
- **Forger:** if the TMAO limit also becomes unlimited upwards when the dispersion is propagated, the certificate separation does not hold in practice for our anchor, and then it must be said — it would be an important negative about a result we just posted.
- **Prohibited:** to use metformin clearance as capacity limit; to mix saturated and linear regime in the same ratio; pooling amount from one cohort with activity from another (that artifact was 2,17–7,76× in another lane today); to treat an estimating equation value as a measured filtration.

## Deliverable

`PORT.json` with the capacity range out of TMAO, the thermodynamic ceiling, and a verdict if the two are compatible. Plus a row about which of our existing renal conclusions rest on metformin and thus need to be recalculated. Narrow follow-ups in FOLLOWUPS.json with external_referent complete — no template repeated per node.

No internal data, no patient data. Everything PENDING_INDEPENDENT_REVIEW.
