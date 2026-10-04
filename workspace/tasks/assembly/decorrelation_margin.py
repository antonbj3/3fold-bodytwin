#!/usr/bin/env python3
"""Replace the failed probability output with a measured decorrelation bracket and a margin.

Why this and not a better probability. The forecast-verification cell scored both surgical decisions as
probabilities and the result was negative: skill +0.038 on implant power and -0.085 on the toric
cylinder, i.e. worse than quoting the cohort rate. The Murphy decomposition said where it failed --
resolution 0.0034 -- which is the signature of collapsing independent error components into one
estimated scale, so the forecast cannot tell patients apart. The operator pointed at how the graph
engine is meant to handle uncertainty, and the graph lane named the module: decorrelation is treated
there as a certifiable resource that must be MEASURED rather than assumed, because a flag nobody
measures fails open.

So the question changes from "what is the probability" to "do the error components actually decorrelate,
and what does that assumption buy". That is answerable with the data on disk and it is falsifiable.

The decomposition, and why each row sums to the real error. A stage decomposition is only honest if the
stages add up to the error that was actually observed, otherwise the covariance is of something else.

  implant power, in dioptres of refraction at the spectacle plane:
    grid   = (power on the 0.5 D manufacturing grid - bias-corrected recommendation) * dR/dP
    chain  = (bias-corrected recommendation - hindsight-correct power)               * dR/dP
    sum    = the twin's actual miss in dioptres, exactly

  toric cylinder, in dioptres of cylinder:
    posterior_information = prediction with measured posterior - prediction with population estimate
    chain_residual        = prediction with population estimate - measured outcome
    sum                   = the residual cylinder actually observed, exactly

The third error source is instrument repeatability, and it is NOT in these matrices. It is measured over
300 eyes from a different cohort and instrument, with no overlap with the 20 and 69 surgical eyes, so
there is no unit definition under which it could be stacked into the same M x K array. Stacking it
silently would be the fail-open the module warns about. It is therefore reported separately as a
declared component with its own scale and the population it was estimated on, and the bracket is given
both with and without it so the cost of that choice is visible rather than hidden.

A caution carried from the graph lane: min_units is 20, so the implant-power verdict sits exactly on the
floor by construction and is borderline whatever it says. The toric case at 69 is not.

WITHDRAWN, and why it is the result rather than a setback. The toric decomposition is certified
DECORRELATED with rho_hat -0.7992 and a margin of 1.1847 that lies outside the module's own bracket of
0.4719 to 0.6613, which is what sent me to the mechanism: the two toric stages are (measured - population) and (population - outcome), so the
population prediction appears in both with opposite signs. It carries 1.285 times the variance of stage
one on its own, so the anticorrelation is INDUCED BY THE CONSTRUCTION and nothing about the errors was
measured. The certification is withdrawn.

The general statement, which is the useful part: an additive stage decomposition of a total that was
observed only once is not identifiable. Any two stages summing to a fixed observed error share a term by
construction, and their correlation is an artefact of where the split was placed. The implant-power case
escapes this precisely because the grid-quantisation stage is computed from the 0.5 D manufacturing grid
rather than from the outcome, so there is no shared term -- and its correlation is 0.026, which is what
an honest near-zero looks like next to a constructed -0.80. That is also exactly the module's own stated
precondition, real per-unit per-stage error SAMPLES, and it means the toric case needs a second
independent measurement per eye before any margin can be certified for it.

CORRECTION to my own reasoning, from the graph lane, 2026-10-03. I read the margin lying outside its
bracket as proof that the input was incompatible with the model. That inference is wrong and the number
was right: net-negative covariance shrinks 1'Sigma1 by cancellation, so a margin above the RSS endpoint
is legitimate arithmetic. What the bracket excursion actually signals is that the function left its own
documented regime, since RSS and SUM are documented as the rho=0 and rho=1 endpoints. The withdrawal
stands, but on the identifiability ground alone -- the stages were not independently measured -- and not
because the margin was impossible. The verifier now returns ANTICORRELATED-CHECK-CONSTRUCTION for this
case, keeping RSS because it stays conservative while setting decorrelation_verified false because
nothing about independence was shown.

A bound that came out of their reproduction and that limits where this failure can occur at all: with K
stages the minimum attainable equicorrelation is -1/(K-1), so at K=8 nothing below -0.143 exists and the
anticorrelated branch is unreachable. My decomposition used K=2, which is the only reason -0.80 was
available to me.
"""
from __future__ import annotations

import json
import statistics
import sys
from pathlib import Path

import numpy as np

W = Path('.')
ENGINE = Path('the public staging tree/3fold-graph-engine/src')
OUT = W / 'results/ASSEMBLY_DECORRELATION_MARGIN'
SPEC_D = 0.5

sys.path.insert(0, str(ENGINE))
from graph_engine.stage_decorrelation_verifier import (        # noqa: E402
    covariance_aware_margin, from_samples_margin, heavy_tail_index, stage_decorrelation_verifier)


def stages_power() -> tuple[np.ndarray, list[str]]:
    rows = json.loads((W / 'results/ASSEMBLY_IOL_DECISION/DECISION_V1.json').read_text())['rows']
    out = []
    for r in rows:
        s = r['sensitivity_dR_dP']
        grid = (r['recommended_power_on_manufacturing_grid_D']
                - r['recommended_power_bias_corrected_D']) * s
        chain = (r['recommended_power_bias_corrected_D'] - r['hindsight_correct_power_D']) * s
        out.append([grid, chain])
    return np.asarray(out), ['grid_quantisation', 'chain_model_error']


def stages_toric() -> tuple[np.ndarray, list[str]]:
    rows = json.loads((W / 'results/ASSEMBLY_TORIC_DECISION/TORIC_V1.json').read_text())['rows']
    out = []
    for r in rows:
        meas, pop = r['measured_posterior']['predicted_cyl_D'], r['population_estimate']['predicted_cyl_D']
        out.append([meas - pop, pop - r['measured_cyl_D']])
    return np.asarray(out), ['posterior_information', 'chain_residual']


def instrument_component() -> dict:
    """The declared third component: measured, but on a different cohort and instrument."""
    src = json.loads((W / 'results/LANE_CORNEA_SHAPE/r3/REPEATABILITY_V1.json').read_text())
    sds = [r['KrSEQ']['sample_sd_D'] for r in src['records']]
    return {
        'component': 'instrument_repeatability',
        'scale_sd_D': round(statistics.mean(sds), 5),
        'scale_sd_D_p95': round(sorted(sds)[int(0.95 * len(sds))], 5),
        'eyes_it_was_measured_on': len(sds),
        'overlap_with_surgical_eyes': 0,
        'why_not_in_the_matrix': ('a different cohort and instrument with no overlap, so no unit '
                                  'definition stacks it into the same M x K array; stacking it silently '
                                  'would be the fail-open this module exists to prevent'),
    }


def repair_response(X: np.ndarray, stage_names: list[str]) -> dict:
    """Does IMPROVING a stage improve the system? Under negative correlation it need not.

    A swarm audit (BT-FW48-AUTO-4fa8c5a121ab31) refuted its own brief's discriminator -- a linear
    fit of residual against bin width, which read R^2 = 0.021 on a 6-point ladder and 0.630 on a
    9-point ladder for the same arm -- and replaced it with the response of the spread to a
    per-operator repair. Measured there, reported spreads understated the truth by 1.90x to 2.74x
    because opposite-sign errors cancelled. The same test belongs here, because this module's
    favourable margins come from exactly that: the toric chain is certified on rho_hat = -0.796.

    So for each stage, set its error to zero -- the best possible repair -- and recompute the
    combined standard deviation. A margin that rests on cancellation gets WORSE when a stage is
    repaired, and a consumer who reads the margin as a bound that survives improvement is misled.
    The interior optimum for a two-stage chain is sigma_j = -rho * sigma_other, which is also
    reported, because below it further improvement costs rather than pays.
    """
    sd_now = float(np.std(X.sum(axis=1), ddof=1))
    out = {'combined_sd_D': round(sd_now, 5), 'per_stage_repair': {}}
    for j, nm in enumerate(stage_names):
        Y = X.copy()
        Y[:, j] = 0.0
        sd_rep = float(np.std(Y.sum(axis=1), ddof=1))
        out['per_stage_repair'][nm] = {
            'combined_sd_after_full_repair_D': round(sd_rep, 5),
            'factor_vs_now': round(sd_rep / sd_now, 4),
            'repair_helps': bool(sd_rep < sd_now)}
    if X.shape[1] == 2:
        s0, s1 = (float(x) for x in X.std(axis=0, ddof=1))
        rho = float(np.corrcoef(X[:, 0], X[:, 1])[0, 1])
        out['two_stage_interior_optimum'] = {
            'rho': round(rho, 4),
            'optimal_sd_stage0_D': round(max(-rho * s1, 0.0), 5),
            'optimal_sd_stage1_D': round(max(-rho * s0, 0.0), 5),
            'stage0_above_optimum_factor': round(s0 / (-rho * s1), 4) if rho < 0 else None,
            'stage1_above_optimum_factor': round(s1 / (-rho * s0), 4) if rho < 0 else None}
    out['independent_combined_sd_D'] = round(float(np.sqrt((X.std(axis=0, ddof=1) ** 2).sum())), 5)
    out['cancellation_removes_fraction'] = round(1.0 - sd_now / out['independent_combined_sd_D'], 4)
    return out


def analyse(name: str, X: np.ndarray, stage_names: list[str], extra_sd: float | None = None) -> dict:
    if extra_sd is not None:
        # The declared component enters as an independent column with the stated scale, drawn once with a
        # fixed seed. This is a DECLARED prior component, not a measurement of these eyes.
        rng = np.random.default_rng(0)
        X = np.hstack([X, rng.normal(0.0, extra_sd, size=(X.shape[0], 1))])
        stage_names = stage_names + ['instrument_repeatability_declared']

    row_sum = X.sum(axis=1)
    v = stage_decorrelation_verifier(X, rho_tol=0.15, min_units=20, n_boot=400, seed=0)
    cam = covariance_aware_margin(X, SPEC_D, ci_pctl=99.0, n_boot=300, seed=0, min_units=20)
    fsm = from_samples_margin(X, SPEC_D, alpha=0.05, kurtosis_flag=True, min_units=20)
    hti = heavy_tail_index(X, frac=0.15, hill_thresh=0.35)
    return {
        'decision': name,
        'units': int(X.shape[0]),
        'stages': stage_names,
        'per_stage_sd_D': [round(float(s), 5) for s in X.std(axis=0, ddof=1)],
        'row_sum_equals_observed_error': {
            'max_abs_row_sum_D': round(float(np.max(np.abs(row_sum))), 5),
            'mean_abs_row_sum_D': round(float(np.mean(np.abs(row_sum))), 5)},
        'spec_D': SPEC_D,
        'decorrelation_verifier': v,
        # The verifier's rule is one-sided: it passes when the correlation's CI UPPER bound is below
        # the tolerance, so a strongly ANTI-correlated chain passes a gate named "decorrelation".
        # That is how rho_hat = -0.796 was certified. Stated here rather than left to be inferred.
        'decorrelation_gate_is_one_sided': True,
        'repair_response': repair_response(X, stage_names),
        'covariance_aware_margin': cam,
        'from_samples_margin': fsm,
        'heavy_tail': hti,
    }


def main() -> int:
    inst = instrument_component()
    sd = inst['scale_sd_D']
    results = []
    for name, fn in (('implant_power', stages_power), ('toric_cylinder', stages_toric)):
        X, names = fn()
        results.append(analyse(name, X, names))
        results.append(analyse(name + '_with_declared_instrument', X, names, extra_sd=sd))

    summary = {
        'question': 'do the error components decorrelate, and what does assuming it buy',
        'replaces': ('a Brier-scored probability whose resolution was 0.0034, i.e. one estimated scale '
                     'shared by all patients'),
        'instrument_component_declared_not_stacked': inst,
        'claim_type': 'capability',
        'control': ('the two explicit endpoints the module computes: spec/sum(sd) at rho=1 and '
                    'spec/sqrt(sum sd^2) at rho=0; the bracket between them is the measured value of '
                    'the independence assumption rather than an assumed one'),
        'caution': ('min_units is 20, so the implant-power verdict sits on the floor by construction and '
                    'is borderline whatever it reports; the toric case at 69 units is not'),
        'results': results,
        'review_state': 'PENDING_INDEPENDENT_REVIEW',
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / 'DECORRELATION_V1.json').write_text(json.dumps(summary, indent=1, ensure_ascii=False, default=str))

    for r in results:
        v, c = r['decorrelation_verifier'], r['covariance_aware_margin']
        print(f"  {r['decision']:42s} M={r['units']:3d} K={len(r['stages'])}")
        print(f"    rho_hat {v.get('rho_hat')}  {v.get('verdict')} -> {v.get('composition')}  "
              f"decorrelation_verified={v.get('decorrelation_verified')}")
        print(f"    margin {c.get('margin')}  bracket rho=0 {c.get('rss_margin')} .. "
              f"rho=1 {c.get('sum_margin')}  ({c.get('verdict')})")
        print(f"    rowsum|max {r['row_sum_equals_observed_error']['max_abs_row_sum_D']} D  "
              f"heavy tail {r['heavy_tail'].get('verdict', r['heavy_tail'])}")
    print(f"  written: {OUT / 'DECORRELATION_V1.json'}")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
