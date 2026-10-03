# PREREG — BT-HX-Q154

## Scope and hypothesis

The substrate is D-glucose and the cell type is an endothelial cell in a brain capillary. The model represents a network-equivalent membrane patch for 1,0 g brain tissue, with separate pools for luminal glucose, free endothelial glucose, receptor-bound glucose, vesicular glucose, abluminal glucose, and lactabouten.

Hypothesis H1 is that the unidirectional glucose transfer at basal plasma concentration is mainly determined by a saturated luminal transporter (GLUT1/GLUT3-liknande), while passive diffusion is small. H2 is that active efflux and receptor-directed transcytosis are separate, direction-dependent pathways. H3 is that endothelial metabolism is a possible net loss from the intact glucose pool and yields a separate lactate pool. The null model H0 is the same geometry and concentrations but without saturated, vesicular and metabolic pathways.

## Previous work

There was no model or result file from a previously aborted session in this directory. `agent.log` only showed directory checking and reading of `BRIEF.md`; no node IDs or reusable BodyTwin files were available to reuse. The external paths mentioned in `NIGHT_PREAMBLE.md` do not exist in this container environment. Therefore, the work is based only on the question and published literature.

## Builds on

- `BRIEF.md`, `inputs/NIGHT_PREAMBLE.md` and `inputs/QUESTION.md`.
- `agent.log`: no previously deployed node or executable.
- Hertz et al. (1981), DOI `10.1172/JCI110073`, verified via PMC-page.

## Not redone

- No internal BodyTwin data, restricted model data runs, external solver runtime or private results.
- No literature review; the source is used only as a numerical comparison.
- No 3D geometry or substrate-specific efflux identification; these are the next resolution steps.

## Predicted quantity and frozen reference

Primary prediction is the time average over the last 600 s of the total luminal glucose flux to the endothelium, including passive luminal diffusion, luminal transporter influx, and receptor endocytosis:

`J_u = J_pass,L->E + J_carrier,L->E + J_receptor,endo`, i µmol g⁻¹ min⁻¹.

The reference is the following published measurement:

- Hertz, M. M.; Paulson, O. B.; Barry, D. I.; Christiansen, J. S.; Svendsen, P. A. (1981), “Insulin increases glucose transfer across the blood-brain barrier in man”, *Journal of Clinical Investigation* 67(3):597–604.
- DOI: `10.1172/JCI110073`.
- Aabstract value before insulin: **0,46 µmol g⁻¹ min⁻¹** unidirectional `[14C]D-glukos` flow; the extraction is given as 14 %. Under insulin is given 0,66 µmol g⁻¹ min⁻¹, which is not the baseline.
- Web verified page: https://pmc.ncbi.nlm.nih.gov/articles/PMC370607/ (article title, author, journal, DOI and abstract checked 2026-09-25).

## Frozen acceptance criterion

The model is accepted for primary prediction of

`abs(J_u_model - 0,46) / 0,46 <= 0,30`

with `J_u_model` mean value of `[t_end-600, t_end]`, and if all pools are non-negative, the numerical mass balance is at most `1e-10 mol/s` and `test_model.py` passes the analytical limiting case. This is a numerical implementation and order of magnitude criterion, not a proof that each mechanism is unique.

## Controls and counter-tests

1. The null model shall use the same geometry and passive coefficients, but all saturated and vesicular terms are set to zero. If the passive alone can reproduce the reference, H1 is not identified by the order of magnitude.
2. Each directed pathway is reported separately: luminal passive in, luminal passive out, transporter in, transporter out, abluminal in, abluminal out, receptor endocytosis and vesicle release.
3. Endothelial metabolite flow and endothelial retention are reported in the same time resolution. The test that can fail the model is that increased total brain flow without an independently directed flux or metabolite profile can be explained by the passive model; then causal identification must be marked `UNKNOWN`.
4. Test failures are counted as: incorrect unit, negative pool, integrated mass balance above tolerance, failed analytical limit case or displacement of reference/boundary after starting the first run.

## Data and limitations

No internal data or assumed measurement data are used. Published Hertz reference is used for external comparison only. Other parameters are expressly assumptions or general physiological ranges and may not be presented as new measurements. GLUT1/GLUT3 and the generic efflux-/transcytosparametern are mechanistic hypotheses for D-glucose; substrate-specific efflux and receptor identity are `UNKNOWN` without inhibitor or tracer data. The model execution is local and 1 threaded.
