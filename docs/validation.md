# Monash Myeloma Model — Validation

## Overview

The model is checked at three layers:

1.  **In-simulation invariants** — `core/validation.do`, run inside every simulation: lightweight ordering/bounds sanity checks (e.g. event times non-decreasing, durations non-negative). Catches gross engine errors on every run.
2.  **Engine verification** — `core/tests/`: unit tests for the Mata/survival primitives plus extreme-value (stress) testing. Verifies the *machinery* with no calibration data.
3.  **Calibration validation** — the **out-of-sample (70/30) analysis** in `analyses/default/`. This is the **mainstay validation, re-run for each model version** to confirm the model reproduces observed outcomes and nothing has regressed. The shared comparison engine is `analyses/default/validate_outcomes.do`.

**This document is the source of truth for how the model is validated.** The strategy is deliberately *holistic*: rather than unit-testing each of the 50 risk equations in isolation, the mainstay simulates the whole composite model and checks it against **real held-out patient outcomes** (Layer 3). For a predictive microsimulation this is a stronger test than per-component checks — a mis-specified equation would surface as an out-of-sample miss — so the OOS validation subsumes per-equation accuracy testing. Layers 1–2 guard the *machinery* (in-run invariants, Mata/survival primitives, stress tests); Layer 3 guards the *predictions*. (An earlier, aspirational per-component test-suite plan was superseded by this approach and retired.)

## Layer 2 — engine verification (`core/tests/`)

| File | Role |
|----|----|
| `core/tests/extreme_value.do` | Extreme-value / stress test: perturb one risk-equation intercept to an extreme and assert the engine responds in the right direction and stays bounded |
| `core/tests/test_reproducibility.do` | Determinism guard: a fixed seed + fixed cohort must give byte-identical output (CRN matrix identical across builds; `run_pipeline` output signature identical on re-run) |
| `core/tests/test_mata_functions.do` | Unit check: ordered-logit and survival helpers |
| `core/tests/test_survival_functions.do` | Unit check: exponential / Weibull / Gompertz inverse-CDF sampling |

`extreme_value.do` needs **no MRDR data** — it runs the engine on a 3k slice of the in-repo `synthetic_1995_2040` cohort with the base coefficients, perturbing one intercept at a time in Mata, and self-checks the result. The latest run passes **7/7**: OS hazard → ∞ gives median OS ≈ 0 months and → 0 gives ≈ 367 months (the age limit); ASCT probability → 1 gives ≈ 98% transplanted and → 0 gives 0%; TXD and TFI hazard → ∞ collapse to median ≈ 0; and a monotone sweep of the OS intercept gives a smoothly decreasing median OS. Run it from the repository root:

``` stata
do "core/tests/extreme_value.do"
```

> Known engine-robustness gap (surfaced by this test): at *intermediate-extreme* OS hazard (intercept shift ≈ +2, not the +20 boundary) the engine throws `r(3301) subscript invalid` — high early mortality empties a downstream line and an unguarded `selectindex` fails. The harness self-protects (`capture` per point); the engine guard (`rows(idx) > 0` before indexing) is an open item.

## Layer 3 — out-of-sample validation (`analyses/default/`, `$scenario outsample`)

The reference analysis's **out-of-sample scenario**. Trains the model on a random **70%** of MRDR patients, predicts the held-out **30%**, and compares the predictions to those patients' **observed** outcomes — testing generalisation, not just in-sample fit. See `analyses/default/README.md` for the full layout; the prep steps (split → train-70% imputation/risk-equations → 30% targets/cohort) run against the restricted registry data and are HPC-suited.

### How to run

The full pipeline — split, train-70% fit, 30% targets + cohort, simulate, and compare — is the out-of-sample track of `analyses/default/run.do` (its numbered steps O0–O6); run that top to bottom, or the individual steps it lists:

``` stata
// The whole deterministic out-of-sample pipeline (see the O-numbered steps inside):
do "analyses/default/run.do"

// ... or just the final compare (steps O0-O4 already done):
//   step O5  do "analyses/default/simulate.do" 0 . . outsample   // simulate the held-out 30% with the 70%-trained coefficients
//   step O6  do "analyses/default/validate_outsample.do"          // compare the simulated 30% to the observed targets
```

`validate_outsample.do` points the shared comparison engine `validate_outcomes.do` at the OOS targets (`$val_targets`) and the OOS simulation (`$val_simfile`); the engine imports the target CSVs inline, loads the simulated dataset, runs the checks below, and prints a pass/fail summary.

### Targets

`analyses/default/targets/` holds 13 CSVs (the held-out 30%'s observed outcomes), built by `prep/generate_benchmarks.do` from the test-fold imputed data: `os_l1_noasct.csv`, `os_asct.csv`, `os_l2.csv`, `os_l3.csv`, `bcr.csv`, `txd_l1_noasct.csv`, `txd_l1_asct.csv`, `txd_l2.csv`, `tfi_l1_noasct.csv`, `tfi_l1_asct.csv`, `tfi_l2.csv`, `tfi_l3.csv`, `pathways.csv`. Overall-survival files carry N, median and annual survival percentages; response carries N and the CR/VGPR/PR/MR/SD/PD percentages by line; treatment-duration and treatment-free-interval files carry N, mean, median, quartiles and the `M12`/`M24` horizon-survival columns; pathways carries the ASCT and subsequent-line reach rates. A further target, `os_wholepop_curve.csv`, holds the observed whole-population monthly Kaplan–Meier survivor (with Greenwood SE) over 0–120 months — the reference for the OS validation curve (below).

The survival, treatment-duration and treatment-free-interval targets are all estimated with censoring-aware survival methods (Kaplan–Meier / `stsum`), so they are comparable to the run-to-death simulation despite the registry's incomplete follow-up. The **subsequent-line reach rates** in `pathways.csv` are estimated **one transition at a time, conditionally** — the Aalen–Johansen probability of reaching line *L* given the patient reached line *L−1*, with *death before reaching L* as the competing event and the *previous line's reach date* as the origin. Validating each transition on its own keeps line-to-line errors from compounding down the pathway (an earlier cumulative-from-L1 estimate let a single early miss cascade into every later line); the validator and the simulation are compared on the same conditional quantity. Each conditional rate is still a competing-risks cumulative incidence rather than a crude "ever reached ÷ total" count, which would understate reach because recently-diagnosed registry patients still in an earlier line would usually progress given more follow-up.

**The ASCT rate uses a conditional denominator too, but a crude one.** Transplant is decided at the end of L1, so the target in `pathways.csv` is the share transplanted **among patients who reach the end of L1** (an L1-end event, `Event 11`), not among all diagnosed patients; patients who die during induction never reach the decision, and dividing by everyone understates the rate the simulation is actually being asked to reproduce. The model's ASCT logit is fit on the same conditional population, and `validate_outcomes.do` scores it on the matching denominator (patients alive at L1 end). Unlike the later transitions it is a crude proportion rather than an Aalen–Johansen estimate, because transplant is observed for essentially every patient and is not subject to the follow-up censoring that affects the deeper lines (`prep/generate_benchmarks.do` → Pathways).

### Checks and tolerances

`validate_outcomes.do` runs five families of checks and counts passes and failures:

| Family | Metric | Tolerance |
|----|----|----|
| Overall survival | 3-year and 5-year survival %, by response and line (L1 no-ASCT, ASCT, L2, L3) | within ±10 percentage points |
| Best clinical response | per-category % by line (L1, ASCT, L2, L3) | within ±5 percentage points |
| Treatment duration (TXD) | \% still on treatment at 12 & 24 months, by response and line (L1 no-ASCT, ASCT) | within ±10 percentage points |
| Treatment-free interval (TFI) | \% still treatment-free at 12 & 24 months, by response and line (L1 no-ASCT, ASCT, L2), plus 12-month response ordering (CR \> VGPR \> PR) at L2 | within ±10 percentage points |
| Pathways | ASCT rate and L2–L5 conditional reach rates (P(reach L \| reached L−1)) | within ±5 percentage points |
| Maintenance duration | KM median `MND_L1` in months, by regimen (lenalidomide / thalidomide / other) | within ±25%, **relative** |

**A minimum-N floor of 20 applies to every cell-based benchmark.** `prep/generate_benchmarks.do`
blanks the estimates in any cell below it and keeps N, so the thinness is visible and the cell is
skipped rather than scored. The 30% out-of-sample fold splits some cells to almost nothing, and a
Kaplan–Meier estimate on 17 patients is noise that presents as a model failure: L4 BCR=1 is N=17 out
of sample against N=54 in sample, the two folds disagreeing by 21 points about the same registry
quantity while the simulation sits near 70% in both. Where the registry cannot say what the truth
is, neither should the target.

**Maintenance is scored on duration, not on share of the gap.** An earlier version scored
`MND_L1 / TFI_L1` by regimen × gap band. Both quantities need an observed L2, so the registry side
could only be computed on patients who relapsed while the simulated side has a closed gap for
everybody — different populations, and systematically so, because at long gaps the relapsers are
exactly the patients whose maintenance ran close to relapse. The model scored 0/7 against that
target. Duration needs no closed gap (patients still on maintenance are censored, which is what KM
is for), uses the whole maintenance population rather than the ~37% who relapsed, and is the
quantity the cost engine actually bills. The cost of the change is that the regimen × gap
interaction is no longer testable on the registry side; that is recorded in `refractory.md` §5(8)
rather than papered over.

TXD and TFI are validated by **survival at fixed horizons** (% still on treatment / still treatment-free at 12 and 24 months), the same way overall survival is checked — not by the median. The median is unreliable for these outcomes because they are extremely right-skewed and heavily censored (TFI 20–54%): poor responders progress almost immediately (median ≈ 0, so a ratio test is unstable), while good responders have so much administrative censoring that the KM median is itself an unreliable, likely under-estimate. Survival at 12/24 months is well inside observed follow-up, so it is estimated robustly for both extremes and is directly comparable between the simulation and the registry.

### Interpreting results

The script ends with a summary of tests run, passed and failed. As a guide: a pass rate above 90% indicates the model reproduces the held-out outcomes well; 75–90% warrants reviewing the failed checks; below 75% indicates a structural issue to investigate before relying on the run.

The latest deterministic point-estimate run (July 2026) scores **83.3% (150/180)** out of sample; the *same* model checked in-sample scores **91.3% (168/184)**. The test count rose from 170 to 180 when the lenalidomide-refractory checks were added (prevalence by line, and OS from L2 by refractory status), so the rate is not directly comparable with the 86.5% (147/170) recorded before them: six of the ten added checks are the refractory ones, of which five now pass. The in-sample 91.3% is the best this model has recorded.

**The aggregate that matters passes comfortably.** Out-of-sample whole-population OS from diagnosis is within **1.6, 2.5 and 3.7 percentage points** at 3, 5 and 10 years. That is the arbiter for the specification, and it is the check a change should be judged on. The 10-year figure rose from 1.8pp when the duration ceilings moved from observed-ends-only to all records: a censored patient still on treatment is evidence that durations of at least that length occur, which is what a curtailment ceiling needs to know.

The remaining fails are data artefacts with known causes, not structural problems: the SMM-tail response misses, missing treatment-end dates inflating the TXD target, the older incidence cohort versus the registry, and the L1 no-ASCT CR under-prediction that vanishes in-sample. Twenty-six of the thirty out-of-sample fails are response-stratified checks of this kind.

**One is worth reading as real signal.** Three-year overall survival in the lenalidomide-refractory group is too favourable — 61.6% against an observed 35.1%. The mechanism is understood and is not a defect in the applied penalty: the fitted hazard ratio averages about 1.64 across the segments that carry it and the engine delivers roughly that, but the observed contrast is *marginal*, and registry refractory patients are sicker in ways the model does not observe. The most cytogenetically informative variable the registry holds, `FISHRisk`, was tested and adds nothing (AUC +0.002), so this is unmeasured biology rather than a missing covariate. Both refractory logits already discriminate at or above the level of a published early-relapse risk model. Recorded in `refractory.md` rather than resolved.

**Two quantities this validation does not score**, both introduced by the same rework and both recorded as open: the simulated mean maintenance duration sits well above its median, and cost bills the mean while this benchmark scores the median; and roughly a quarter of maintenance patients now receive a gap barely longer than their maintenance, hence no treatment-free interval at all.

### Bootstrap prediction intervals and the OS validation curve

The point estimate above uses one 70%-trained coefficient set. The **calibration** question — does each held-out observed value fall inside the *interval* the model predicts? — is answered by `analyses/default/bootstrap_validation.do` on the HPC, using **500 bootstrap simulations** of the held-out 30% (each a patient-cluster resample of the 70%, re-imputed, re-fit and re-simulated). By the **percentile method** it asks whether each observed target lies within the bootstrap 95% interval [p2.5, p97.5]. Latest coverage: **105/171 (61.4%)** — in line with the model's history, dominated by the tight parameter-only intervals plus the known data artefacts above; aggregate OS passes.

The same run also builds the **whole-population OS validation curve** (the 2024 PLOS ONE Fig 2 for this model): the held-out validation cohort's Kaplan–Meier **95% CI** against the simulated cohort's **95% CI** over 120 months, plus a **monthly two-sample z-test p-value** (observed Greenwood SE + simulated bootstrap SD) testing the difference. Result: **no significant difference in 118/120 months (98.3%)** — the simulated OS 95% CI sits within the observed CI across the full 10 years (p < 0.05 only in the first ~2–3 months), improving on the 2024 paper's 90%. It writes `results/os_wholepop_curve_validation.csv` and `os_wholepop_curve.png`; `analyses/default/plot_os_curve.do` redraws the figure locally from the CSV when HPC batch graphics are off.

## When to re-validate

Re-run the OOS validation **for each model version**, and whenever the training data are updated, the imputation strategy changes, variable definitions change, or the model structure changes (for example, the response categories). The targets are regenerated as part of the prep steps, so each version is validated against a freshly-derived held-out set.
