# Monash Myeloma Model — Technical Review

**Version:** 3.0 · **Updated:** 2026-08-11 · **Tooling:** Stata 15+ (Mata)

This document describes the architecture and implementation of the Monash Myeloma Model. Where a method underpins a specific analysis (notably Calibrated Transport), the authoritative specification lives in that analysis's own README; this review covers the shared engine.

## Executive summary

The Monash Myeloma Model is a discrete-event microsimulation of the multiple myeloma treatment journey, from diagnosis through up to nine lines of therapy and death. It predicts best clinical response, treatment duration, treatment-free interval and overall survival from 50 risk equations estimated on the Australia and New Zealand Myeloma and Related Diseases Registry (MRDR), and attaches costs and quality-adjusted life years for health-economic evaluation.

The engine is implemented in Stata's Mata language as a vectorised microsimulation: every patient is advanced through the same event simultaneously using matrix operations. All stochastic events draw from a pre-generated common-random-number (CRN) matrix, so two arms of a comparison see identical randomness for each patient — the basis for variance-reduced cost-effectiveness comparison.

## Architecture overview

### Design philosophy

The model holds the simulated population entirely in Mata memory and operates on it column-by-column across the patient pathway. The governing convention (set in `core/mata_setup.do`) is:

- **Vectors** hold fixed patient characteristics (age, sex, ECOG, R-ISS, comorbidity, ASCT intent, etc.), one element per patient.
- **Matrices** hold pathway-varying outcomes, one row per patient and one column per pathway point.

A run is driven by a small set of `global` macros set in a per-analysis dispatcher; there is no positional command-line interface.

### Pathway points (OMC)

The patient journey is discretised into 19 outcome-milestone checkpoints (OMC). Column 1 is diagnosis (DN); columns 2–19 are the start and end of each line, L1S/L1E through L9S/L9E. The engine walks these points in order, and outcome matrices are dimensioned `Obs × 19` to match.

## Run interface

Each analysis ships a dispatcher do-file (`analyses/<name>/simulate.do`) that owns a configuration block of globals, loads the core programs, loads the relevant coefficient set, and runs the pipeline; a per-analysis `run.do` runbook drives the dispatcher (looping the arms of a comparison, or chaining prep steps). Interactive runs are globals-only; `run.do` and the HPC arrays additionally pass a few optional positional overrides (`boot`, `min_bs`, `max_bs`, `scenario`). Dispatchers assume the **working directory is the repository root** — all paths are repo-root-relative.

The configuration block (canonical form in `analyses/default/simulate.do`):

| Global | Meaning | Default |
|----|----|----|
| `$analysis` | Analysis name (also sets the four `*_path` globals) | `default` |
| `$int` | Intervention label (two-arm analyses use `$int1`/`$int0`) | `all` |
| `$line` | Line of therapy assessed (`0` = full pathway from diagnosis; 1–9) | `0` |
| `$coeffs` | Coefficient set loaded via `mata matuse` (e.g. `full`, `train`) | `full` |
| `$data` | Patient cohort: `synthetic` (incidence population), `train`/`test` (registry folds), or `predicted` | `synthetic` |
| `$min_year` / `$max_year` | Diagnosis-year range | `1995` / `2040` |
| `$min_id` / `$max_id` | Patient ID range | `1` / `101212` |
| `$cost_year` | Price year for costs (AUD) | `2026` |
| `$drate` | Annual discount rate (PBAC = 5%) | `0.05` |
| `$report` | Generate PDF report (0/1) | `0` |
| `$scenario` | Scenario label (woven into output paths) | `""` |
| `$boot` | Bootstrap flag (0/1) | `0` |
| `$min_bs` / `$max_bs` | Bootstrap iteration range | `""` |

Every dispatcher loads the core programs with `run "core/…"` and runs the lean simulation pass via the shared program **`run_pipeline`** (`core/run_pipeline.do`), which performs `load_patients` → `mata_setup` → `simulation` → `process_data`, after sourcing `core/mata_functions.do` and `core/rng_slots.do` when the Mata functions aren't already compiled. `run_pipeline` deliberately excludes CSV export and saving, so each dispatcher composes those around it: save the simulated `.dta`, run `core/validation.do`, and optionally `core/generate_report.do` (`$report`) / `core/export_results.do`.

## Simulation flow

A typical run proceeds:

1.  **`core/load_patients.do`** (`load_patients`) — reads the patient `.dta` (a `synthetic_1995_2040_<n>.dta` cohort or a predicted `patients_<analysis>_<line>.dta`), filters on diagnosis year, disease state and ID range, and **resets `ID = _n`** so row order is canonical. This ordering is load-bearing for CRN alignment.
2.  **`core/mata_setup.do`** (`mata_setup`) — builds the Mata characteristic vectors and outcome matrices from the Stata data, and constructs the CRN matrix `mRN`. It asserts `ID == _n` (errors otherwise).
3.  **`core/simulation_engine.do`** (`simulation`) — the deterministic event loop. At each of the 19 OMC points it sets the current `OMC` and `Line` and executes the relevant `core/outcomes/sim_*.do` module, filling the outcome matrices.
4.  **`core/process_data.do`** (`process_data`) — drops `mRN`, stacks the outcome matrices into a summary matrix, writes them back to a flat Stata dataset with long variable names, and computes dates, costs, QALYs and discounting.

The dispatcher then saves the dataset, runs the in-run invariant checks in `core/validation.do`, and optionally writes CSV exports and/or a PDF report.

## Data structures

All built in `core/mata_setup.do`. Outcome matrices are `Obs × 19` unless noted, with paired row/column label matrices.

| Matrix | Shape | Contents |
|----|----|----|
| `mState` | `Obs×2` | Entry disease state, diagnosis date |
| `mAge` | `Obs×19` | Age at each pathway point |
| `mOS` | `Obs×19` | Simulated overall-survival time from diagnosis |
| `mTNE` | `Obs×19` | Time to next event (duration of the current OMC) |
| `mTSD` | `Obs×19` | Cumulative time since diagnosis (`TSD_DN = 0`) |
| `mMOR` | `Obs×19` | Mortality flag (1 dead, 0 alive, `.` not reached) |
| `mOC` | `Obs×2` | Final outcome: death/censor time, mortality flag |
| `mTXR` | `Obs×9` | Treatment regimen code per line |
| `mTXD` | `Obs×9` | Treatment duration (months) per line |
| `mBCR` | `Obs×10` | Best clinical response per line, plus post-ASCT response |
| `mTFI` | `Obs×9` | Treatment-free interval (TFI for line *L* sits in column *L+1*) |

Fixed characteristics are held as vectors: `vID`, `vAge`/`vAge2` (updated to current age each OMC), `vMale`, ECOG dummies, R-ISS and ISS dummies, four comorbidity indicators (chronic kidney disease, cardiac, pulmonary, diabetes), age-threshold indicators (`vAge70`, `vAge75`), the constant `vCons`, and ASCT/maintenance receipt vectors (`vSCT_DN`, `vSCT_L1`, `vMNT`). `process_data` reassembles these into a single summary matrix and exports them as long-named Stata variables (`OS_DN`, `BCR_L1`, `OC_TIME`, …).

## Risk-equation / outcomes architecture

The 50 risk equations are organised into self-contained modules under `core/outcomes/`. Each module filters to eligible patients, assembles a design matrix from the characteristic vectors, extracts the relevant coefficient block from an externally loaded Mata matrix, computes a linear predictor, draws a CRN uniform, maps it to an outcome, and writes back into the appropriate matrix. Most carry a design/coefficient dimension guard (`exit(459)` on mismatch).

| Module | Outcome | Method |
|----|----|----|
| `sim_asct_dn` | ASCT intent at diagnosis | Logistic |
| `sim_asct_l1` | ASCT receipt at L1 | Logistic (gated on regimen/response) |
| `sim_tfi_dn` | TFI diagnosis → L1 | Parametric survival |
| `sim_tfi_l1` | TFI L1 → L2 | Parametric survival, split by ASCT status |
| `sim_tfi` | TFI for L2+ | Parametric survival, line-specific |
| `sim_txr` | Treatment regimen per line | Multinomial logit (pooled "Other" fallback) |
| `sim_bcr` | Best clinical response (1=CR … 6=PD) | Ordered logit (6 categories, 5 cutpoints) |
| `sim_bcr_asct` | Post-ASCT response | Ordered logit |
| `sim_txd_l1` | L1 treatment duration | Three-segment spline survival (fixed+ASCT), plus continuous-therapy branch |
| `sim_txd` | L2+ treatment duration | Parametric survival, line-specific |
| `sim_os` | Overall survival | Per-line parametric survival: one Weibull per pathway stage, clocked from that stage's start |
| `sim_mort` | Death at this OMC | Deterministic (time crosses OS) |
| `sim_age` | Age update and age-limit deaths | Deterministic |
| `sim_mnt` | Maintenance receipt | Logistic, split by ASCT status |

The parametric survival families implemented in `core/mata_functions.do` are **exponential, Weibull, Gompertz, log-normal and log-logistic** (closed-form inverse-CDF sampling; `ereg` is an alias for exponential). Overall survival and treatment duration use Weibull; the treatment-free interval uses log-normal. The family for each equation is selected by a string global (`fbOS_L1S`, `fbDN_TFI`, …). Response is modelled by **ordered logit** (`calcOrdLogitProbs` + `assignOrdOutcome`) and regimen choice by **multinomial logit**; ASCT and maintenance receipt are logistic. The L1 treatment-duration equation uses three conditional spline segments for the fixed-duration plus ASCT group. The four comorbidity flags (renal, cardiac, pulmonary, diabetes) enter overall survival and both ASCT logits.

## Mortality and survival tracking

Survival is the master clock. `sim_os` fits a **separate parametric model to each pathway stage** — diagnosis, each line's start and end, with L1 end split by transplant status — each clocked from that stage's own entry event. At a line's first stage the elapsed time is zero, so the draw is an **unconditional** fresh survival on that stage's clock; the result is then stored back on the diagnosis clock (stage origin + residual) so mortality is compared like-for-like. (Lines 6+ share one model clocked from L6 start, so their later stages do condition on survival since L6.) `sim_mort` then declares death at an OMC when cumulative time-since-diagnosis crosses the drawn survival time, sets `mMOR`/`mOC`, and clips the realised duration (`mTXD` or `mTFI`) to the time actually experienced.

This is a single-cause time-to-event design rather than a formal competing-risks framework, with two additional absorbing mechanisms: an **age cap** (`Limit = 100`) that bounds age at death and can retroactively end a patient in the prior OMC, and a **terminal absorption** at L9E that forces any remaining survivors to die. `core/validation.do` enforces the survival invariants (no death-flag reversal, non-negative outcome times, monotone time-since-diagnosis, no outcomes recorded after death).

## Common random numbers (CRN)

The model uses common random numbers to reduce the variance of incremental (between-arm) cost-effectiveness estimates. Every stochastic event reads a fixed column of a pre-drawn uniform matrix for a fixed patient (row), so the same uniform feeds the same event for the same patient in every arm of a comparison; only the coefficients differ between arms, which isolates the treatment effect.

`core/rng_slots.do` is the single source of truth for the column layout of the CRN matrix. It defines a total width `rn_K()` = 74 columns, partitioned into named blocks with per-event accessor functions returning the absolute column for an (event, point) pair — for example `rn_os(omc)` (cols 1–19), `rn_bcr(line)` (cols 20–28), `rn_txr(line)` (cols 30–38), `rn_txd_l1(seg)` (cols 39–43), and a reserved override block `rn_override(1..8)` (cols 67–74). The accessor used throughout the outcome modules is:

``` mata
real colvector rnDraw(real colvector idx, real scalar slot) {
    external real matrix mRN
    return(mRN[idx, slot])
}
```

The matrix itself is built in `mata_setup`:

``` stata
if ("$crn_seed_base" == "") global crn_seed_base 20260615
local _crn_seed = $crn_seed_base + `_b'   // `_b' = bootstrap index (0 if not bootstrapping)
set seed `_crn_seed'
mata: mRN = runiform(Obs, rn_K())          // Obs × 74
```

Cross-arm alignment rests on three guarantees: an identical seed for both arms within a replication; identical cohort row order (enforced by `ID = _n` in `load_patients` and asserted in `mata_setup`); and fixed (patient × slot) addressing. Across bootstrap replications the seed is offset by the iteration index, so Monte-Carlo noise is independent between replications while arms stay aligned within each. CRN is unconditional — there is no `runiform()` fallback and no runtime toggle. The matrix is dropped in `process_data` to free peak memory.

The slot registry imposes two rules: sequential draws for the same patient must use distinct slots (e.g. the three L1 treatment-duration spline draws), and an override that *replaces* a core event reuses that event's slot, while a genuinely new stochastic event takes a column from the reserved override block.

## Bootstrap uncertainty

Parameter uncertainty is propagated through the coefficients. When `$boot == 1`, the dispatcher loops over iterations `$min_bs..$max_bs`; each iteration loads a resampled coefficient set (`coefficients_<set>_B<b>.mmat`), re-runs the full pipeline with the CRN seed offset by `b`, and saves a per-iteration dataset under `…/bootstrap/`. Each replication therefore combines a resampled coefficient vector (parameter uncertainty) with independent Monte-Carlo noise, while the two arms remain CRN-aligned within the replication.

`core/process_bootstraps.do` aggregates a two-arm comparison: for each saved replication it extracts per-arm means (costs, discounted QALYs, mean overall survival, treatment duration) and the response/subsequent-line distributions, computes incrementals and the ICER, and reports means with 2.5/97.5-percentile bootstrap confidence intervals.

## Analyses

- **`default`** — the reference analysis, run two ways via `$scenario`. **Projection** (`$scenario ""`): the full treatment landscape (all observed MRDR regimens), fit on the whole registry (`coeffs=full`) and applied to the synthetic incidence population (`data=synthetic`) — used for projections and as the baseline for the health-economic models; output `analyses/default/simulated/all_0_synthetic.dta`. **Out-of-sample validation** (`$scenario "outsample"`): the 70%-train fit (`coeffs=train`) applied to the held-out real 30% (`data=test`), compared to observed targets by `validate_outsample.do` — the mainstay calibration check.
- **`transport_dvd`** — a two-arm comparative cost-effectiveness analysis of DVd versus Vd at line 2, built on the Calibrated Transport method (below). Runs under three scenarios (`A_trial`, `B_transport`, `C_mrdr`), each saved to its own `simulated/<scenario>/` subtree.

## Calibrated Transport

Calibrated Transport predicts the real-world Australian second-line outcomes of a regimen — here DVd at line 2 — before funding, when no domestic real-world data on the new regimen exist. The method exploits a comparator observed in **both** the trial (CASTOR) and the registry (MRDR): Vd. This shared anchor calibrates the trial DVd-versus-Vd effect into the Australian setting, re-expressing it as an absolute real-world prediction. Implementation fits an ordered-logit response model on a stacked dataset of the registry Vd anchor and the resampled trial arms (with line-of-therapy indicators in the on-disk model), then overrides line-2 response in the simulation via `core/outcomes/sim_bcr_override.do`, drawing CRN-aligned uniforms so the DVd and Vd arms are paired patient-by-patient. The three scenarios contrast the traditional trial transfer (`A_trial`), Calibrated Transport (`B_transport`) and the observed MRDR benchmark (`C_mrdr`).

The definitive method specification, including the exact regression form and the validation results, lives in `analyses/transport_dvd/README.md` and the associated paper. A supporting Monte-Carlo precision workflow justifies the simulated cohort size:

- `cohort_pool.do` builds a reusable, arm-agnostic line-2 entry pool once.
- `ce_cohort.do` draws the fixed-size production cohort (50,000) from that pool.
- `ce_precision.do` estimates the per-patient SD of the incremental outcome using common random numbers, giving the Monte-Carlo SD of a size-N mean.
- `ce_sample_size.do` combines that per-patient SD with bootstrap parameter uncertainty (an O'Hagan/Stevenson/Madan ANOVA decomposition) to report the required cohort size.
- `bootstrap_summary.do` aggregates the cross-scenario bootstrap into response distributions, prediction error versus the observed benchmark, and ICERs.

## Health-economic integration

`process_data` attaches costs and utilities to the simulated pathway and discounts both at `$drate` (PBAC default 5%) to the `$cost_year` price base. Costs are accumulated by treatment and line; quality-adjusted life years are derived from response- and line-specific utilities applied over the simulated time. Two-arm analyses collapse mean discounted cost and QALY by arm to form the incremental cost-effectiveness ratio.

## Outputs

`core/export_results.do` (`export_results`) writes flat CSV outputs for downstream use — a response distribution (`bcr_*.csv`), an economic summary (`econ_*.csv`), and a per-patient export (`patients_*.csv`) — into `$simulated_path/$scenario/`, without modifying memory. It is point-estimate only: it exits early when `$boot == 1`, since bootstrap output is aggregated separately. It is wired into the pipeline by the orchestrators that call it (e.g. `transport_dvd/simulate.do`), not by `default/simulate.do`.

`core/generate_report.do` produces a `putpdf` report (titled "Monash Myeloma Model v3.0") when `$report == 1`, covering the patient sample, treatments, overall survival (with figures), lines of therapy, treatment duration, treatment-free interval and economic outcomes.

## Validation

Validation is described in `docs/validation.md`. In brief, three layers: `core/validation.do` runs lighter invariant checks within each simulation; `core/tests/` holds engine-verification tests (the Mata/survival unit checks plus `extreme_value.do` stress testing); and the **out-of-sample (70/30) analysis** in `analyses/default/` is the mainstay calibration validation — it trains on 70% of MRDR, predicts the held-out 30%, and compares to observed outcomes across five families (overall survival, response, treatment duration, treatment-free interval and pathways) via the shared engine `analyses/default/validate_outcomes.do` with documented tolerances. (The earlier in-sample registry acceptance test was retired; its debug diagnostics live in `scratch/`.)

## Performance

The vectorised Mata implementation processes all patients simultaneously through matrix operations, which allows large cohorts (tens of thousands of patients, with bootstrap replication) to be simulated in practical time. Memory is the main constraint at scale; the CRN matrix is the largest transient structure and is released in `process_data` once outcomes are finalised.

## Common pitfalls and debugging

- **Dimension mismatch (`exit(459)`).** Each outcome module checks that its design matrix width matches the loaded coefficient block. A mismatch usually means the wrong coefficient set was loaded for the analysis, or a characteristic vector is missing.
- **Broken CRN alignment.** If `ID != _n` after loading (for example, an extra filter or sort applied outside `load_patients`), `mata_setup` errors. Preserve the canonical row order, or arms will no longer be paired.
- **Response probabilities not summing to one.** Ordered-logit cutpoints must be loaded in order; a truncated or mis-ordered coefficient block produces invalid probabilities.
- **Extreme survival predictions.** Check the survival-family global (`fb*`) matches the fitted distribution; an exponential global applied to Weibull coefficients yields implausible tails.
- **Working directory.** Dispatchers assume the repository root. Run from there (the built-in error names the missing `core/…` file if you do not).
