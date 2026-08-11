# Changelog

All notable changes to the Monash Myeloma Model project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Planned
- Integration with R for post-processing analysis
- Extended documentation for health economic applications
- Daratumumab **SC** regimen option (a regimen-file switch; cheaper at this body weight and closer to current practice)

## [3.0] - 2026-08-11

### Added
- **Calibrated Transport methods** (`analyses/transport_dvd/`): out-of-trial outcome prediction (e.g. DVd at L2) with a common-random-number engine, sample-size workflow, and cohort pipeline.
- **Common Random Numbers (CRN)**: aligned RNG across treatment arms via an `mRN` slot registry (`core/rng_slots.do`) and `rnDraw` migration, for variance-reduced cost-effectiveness comparisons.
- **Per-line overall survival**: OS re-specified as a separate parametric model per line/stage of therapy, each clocked from that line's own start (`OS_DN`, `OS_L1`/`_NoASCT`/`_ASCT`, `OS_L2..L5` (+`_End`), `OS_L6plus`), replacing the single from-diagnosis survival curve — removes an accumulated-time bias that inflated survival for poor responders at later lines. Engine: `core/outcomes/sim_os.do` (per-stage firing map, diagnosis-clock storage).
- **Individual comorbidity covariates**: the OS and both ASCT-eligibility equations now carry four individual comorbidity flags — renal impairment (`CM_CKD`, derived from imputed eGFR), cardiac (`CM_CRD`), pulmonary (`CM_PLM`) and diabetes (`CM_DBT`) — replacing the earlier single ordinal comorbidity score (`CMc`). Engine plumbing in `core/mata_setup.do`, `sim_os.do`, `sim_asct_*.do`, `process_data.do`.
- **Consolidated pipeline**: dispatchers unified onto `core/run_pipeline.do`; added `analyses/transport_dvd/ce_sample_size.do`.
- **Standardised CSV result exports**: machine-readable results for downstream/programmatic access (R/Python post-processing, dashboards, assistant-driven manuscript drafting) instead of manual copy-to-Excel.
  - **`core/export_results.do`**: engine-level export of CSVs common to every analysis (per-patient summary, BCR distribution, mean cost/QALY/LY). Runs by default as part of the simulation pipeline (immediately after `process_data`, once per arm; skipped during bootstrap), reading `core/process_data.do` outputs into `simulated/<scenario>/`. First adopted by the `dvd_method` dispatcher.
  - **`analyses/<name>/results/` contract**: each analysis exposes a single `results/` folder of final (cross-scenario) CSVs plus a `results.md` narrating the key figures — the canonical downstream read surface. Analysis-specific and cross-scenario aggregation live under `analyses/<name>/`, not `core/`.
- **Refractory-status capture in the extraction** (`prep/data_extraction.do` + new `prep/sub/MRDR/build_refractory.do`). `MRDR Long` now carries, per chemotherapy line, an IMWG refractory flag (`refr_line`: non-response on the line, or a progression event within 60 days of it); per patient, `refr_len_mnt` and `refr_len_mnt_date`; and two carried-forward **lenalidomide-refractory** flags evaluated *as at entry to each line* (strictly prior lines only) — `refr_len_tx_in` (treatment-dose) and `refr_len_mnt_in` (maintenance-dose, line-resolved by date), kept separate because the form the drug sits in is the dose proxy. *(Names as shipped; see Changed above. `LineProg`, `LineProgDate`, `MNT_Refr`, `MNT_Len` and the patient-level maintenance drug were dropped before release.)* Progression uses a hybrid definition (Review structured ∪ treatment-record; the treatment record is essential — ~390 patients progress with no Review progression filed). Additive only: `MNT` and the event skeleton are unchanged, so `risk_equations.do` runs identically. Neither flag is consumed by the shipped engine: this entry is data capture only, and the generation/consumption wiring is built but held (see Known issues). Definitions, evidence and design decisions: `docs/refractory.md`.
- **`docs/refractory.md`** — the specification of record for the refractory subsystem. It also carries the specification for `prep/data_extraction.do` and `prep/sub/`, which are git-ignored under the data-governance rule and therefore have no other home in the repo.
- **Maintenance duration: capture and equations** (`prep/data_extraction.do`, `prep/risk_equations.do`, `analyses/default/outcomes/mnr_{full,train}.do`), addressing the maintenance over-costing under Known issues below. `MRDR Long` now carries, per patient, `MND_L1` (maintenance duration delivered in the billed L1 gap, days - the benchmark target) and `MNR_L1` (maintenance regimen), and the extraction now KEEPS the L1 maintenance start/end events (110/111) in the event skeleton instead of dropping them. `risk_equations.do` gains **`gen_mnr`**, mirroring `gen_txr`, and two equations beside the `MNT` logit complete whether / which / how long: **`L1_MNR`** (mlogit, lenalidomide vs thalidomide) and **`L1_MND`** (parametric survival on the duration, split by transplant like `L1_TFI`: `i.MNR_L1`, Age, `i.ECOGcc`, `i.RISS`, and `i.BCR_SCT` for the ASCT arm / `i.BCR_L1` for the no-ASCT arm). SIMPLE-FIRST scope: lenalidomide and thalidomide only, no 'other' regimen. `L1_MND` carries `ln(TFI_L1)` with a regimen interaction so lenalidomide duration scales with the gap (runs to progression) and thalidomide stays fixed; the engine caps the drawn duration at the realised TFI_L1. The restored maintenance events subdivide the L1-to-L2 span but do not move any other equation's origin/failure (which resolve by date), so OS/TFI/TXD are untouched. Consumed by the engine: see Fixed below. Design and evidence: `docs/refractory.md` sections 4.4 and 7.
- **Reproducible treatment-cost engine** (`prep/`): treatment costs rebuilt from first principles as the PBS **Dispensed Price for Maximum Quantity**, derived from a dated PBS Schedule extract rather than a manual spreadsheet — `build_cost_index.do` (ABS CPI deflator), `extract_pbs_costs.do` (Schedule → drug prices), `extract_pbs_restrictions.do` (eligibility reference) and `treatment_costs.do [year] [wholepack|prorata]` → `prep/inputs/treatment_costs_<year>.csv`. Non-treatment costs are phase-of-care (Yap 2025), with the initial phase netted of the transplant admission to avoid double-counting. Perspective is the Australian health system, so the full DPMQ is used with the co-payment **not** netted off; oral packs default to whole packs (dispensing wastage costed, per PBAC convention), with pro-rata as a sensitivity. Derivations in `docs/economic_inputs.md`.
  - **`$cost_year` falls back to the latest available file** when the requested year has no cost CSV. Check the cost year before re-running any near-submission analysis: 2026 generic price disclosure moved drug costs materially, so a re-run can silently reprice an older analysis.

### Changed
- **Rebranded** from *EpiMAP Myeloma* to **Monash Myeloma Model**; GitHub repository renamed `CHE-Monash/EpiMAP-Myeloma` → `CHE-Monash/Myeloma-Model` (old URLs auto-redirect). Published papers and DOIs retain the EpiMAP Myeloma name.
- **Version bump to 3.0**: consolidates the v2.1 vectorised engine with the Calibrated Transport/CRN methods, the July 2026 calibration work (per-line OS + individual comorbidities), the reproducible PBS treatment-cost engine and the refractory/maintenance subsystem into a single major release. From v3.0 onward, major versions increment with each published paper (see Version Naming Convention).
- **The refractory variables are one family, `refr_*`, and the maintenance outcome moved to the extraction.** Three changes in one pass, because they share a rebuild. (a) **Cull**: `MNT_TTM_L1`, `MNT_Refr`, `LineProg` and `LineProgDate` were built, kept and never read, and are removed; `EpRefr` (now `refr_ep`) stays, since it still builds the maintenance flag. (b) **Move**: `LENREFR_MNT` was fitted on an outcome derived at fit time in `risk_equations.do`, while `data_extraction.do` built a differently-defined `MNT_LenRefr_L1` that nothing consumed. The event-skeleton definition now lives in the extraction as `refr_len_mnt_l1` and the fit-time block is gone, taking `mi update`, `mi passive` and `mi xeq` out of `risk_equations.do` with it - those were writing Stata `__mitmpfile*.dta` into the working directory, and the `mi xeq` built a per-imputation extract on every bootstrap replicate to print a diagnostic table. The outcome is now **missing** where maintenance has not ended rather than 0, so `logit` drops the undetermined patients without a separate flag, and it is derived before the bootstrap date shift rather than after it. The two definitions disagreed on 84 of 545 determined rows in both directions, and the retired one coded 485 still-on-maintenance patients as not refractory (`scratch/refractory/mntrefr_defcheck.do`). (c) **Rename**: `refr_line`, `refr_len_mnt`, `refr_len_mnt_date`, `refr_len_mnt_l1`, `refr_len_tx_in`, `refr_len_mnt_in`, `refr_len_in`, `refr_len_grp`, `refr_len_oth`, `refr_ep`, `refr_len_l1..l9`, `bcr_grp_l1`. Mata objects keep UPPER (`LENREFR_TX`, `LENREFR_MNT`, `LENREFR_regimens`, `LENREFR_pother`) and the engine keeps `vLenRefr_in` / `mLenRefr_in`, so case now marks the layer rather than the concept. **`refr_len_l1..l9` are the columns the cohort pools carry, so every pool must be rebuilt.** Naming rules and the definition evidence: `docs/refractory.md`.
- **`run_pipeline.do` loads the four stage programs itself** (`load_patients`, `mata_setup`, `simulation`, `process_data`), so an orchestrator needs `run "core/run_pipeline.do"` alone instead of the five-line preamble that was duplicated in eight of them. They are run in the file, outside the program body, so they are parsed once rather than on every call. `export_results.do` deliberately stays with the callers.
- **`LENREFR_MNT` refitted on the rebuilt MI data is weak**: `Prob > F = 0.073` and the response term non-significant, against the monotone coefficients and AUC 0.708 previously recorded. The outcome is unchanged (the definition check reproduces 156/545 exactly), so the cause is the covariates - most likely the per-line rebuild of the `BCR_L1` imputation. The response term is **kept regardless**, per the pre-specified-baseline rule. `docs/refractory.md` 3.5.
- **Best clinical response is imputed PER LINE, congenially with the analysis** (`prep/multiple_imputation.do`). The single pooled model (`BCR = ... i.CLine`) conditioned on the line NUMBER, not the previous RESPONSE, while every risk equation conditions on the previous response - which attenuated the L1-to-L2 association by a third to a half. Now `impute_bcr` imputes L1, then `BCR_SCT`, then L2 each conditioning on what precedes it, with L3-L9 pooled (the paraprotein deltas nearly determine response and quasi-separate the ologit at per-line sample sizes). Each line imputes a dedicated working variable rather than the shared `BCR`, which otherwise appears on both sides of its own equation.
- **`L1_MND` simplified to two pooled arms** (lenalidomide, thalidomide), dropping the transplant split and the response covariate. Response was tested exhaustively once `BCR_SCT` was corrected - post-transplant response, L1 response collapsed and full 6-level, the two combined, and a transplant interaction - and none predicts maintenance duration (all LR p > 0.4, every AIC worse than baseline). The apparent signal that justified the split was an artefact of the `BCR_SCT` coding defect.
- **Duration ceilings are taken from ALL records, censored included** (`save_max`; `save_max_obs` retired). A patient censored while still on treatment is evidence that durations of at least that length occur, which is what a curtailment ceiling needs; taking the maximum over observed ends only understated it and truncated the simulated tail. Affects every TXD, TFI and MND ceiling.
- **Simulation cohorts migrated from `population_*` to `synthetic_*`.** The `population` cohort token is retired: `core/load_patients.do` now resolves `$data` to `synthetic` / `synthetic_<n>` (the incidence cohorts), `$cohort_file` (an explicit override) or a predicted patient file. `patients/population_1995_2040_*.dta` are superseded and deleted — they carried the retired ordinal comorbidity score (`CMc`) and the unused `CM_LVR` / `CM_PNR` / `CM_MLG` flags, and their covariates came from an imputation model that still included them. `patients/synthetic_1995_2040_1..10.dta` (from `prep/synthetic_1995_2040.do`) replace them.
  - The `patients/population_historical.csv` and `population_forecast.csv` **incidence inputs** (AIHW and Daffodil Centre) keep their names: they are inputs *to* the cohort build, not cohorts.
  - Line-entry cohort pools (`analyses/*/patients/cohort_pool.do`) now loop the synthetic cohorts. A pool built from the old `population_*` files is not reproducible against the current inputs and must be rebuilt.
- **Removed the dead `$data == "population"` early-exits** from `core/simulation_engine.do` (four line-truncation branches belonging to the archived `base_model` analysis; they never fired for the `population_<n>` / `synthetic_<n>` cohort-pool tokens).
- **Two-arm report shows patient counts in the treatment-pathways table.** The number receiving each line (n) now sits directly below the Reached % in each arm's column (`core/generate_report.do`), so the proportion and the absolute count are read together.
- **Two-arm report: cost-over-time figure.** New page showing mean undiscounted cost per patient by year since the decision line, split treatment vs non-treatment for each arm (`core/generate_report.do`). Treatment cost is allocated across each line's relative-time window; non-treatment across the survival window at the phase rate. A per-patient trajectory (no uptake) intended as the building block for a calendar-year budget-impact analysis.

### Fixed
- **`BCR_SCT == 0` conflated "no transplant" with "transplanted, response not recorded"** (`prep/multiple_imputation.do`). The zero-fill applied to every missing value despite its own comment saying it was for non-transplant patients, so 25.7% of the transplant arm sat in a category that meant two different things. Four equations handled it three different ways, and it masked a failing imputation: the chained block at `Event0 == 100` left 548 of 2,135 transplanted patients unimputed without erroring. `BCR_SCT` is now imputed in its own right, using `BCR_L1` as a predictor, and 0 means no transplant only - an invariant `prep/mi_diagnostics.do` checks directly.
- **The MND benchmark compared delivered exposure against a death-censored target.** `risk_equations.do` censors death when fitting the duration (cause-specific, with the engine imposing mortality separately) but the engine EXPORTS delivered exposure, so `generate_benchmarks.do` was scoring two different estimands. Death (`Event1 == 104`) is now a failure in the MND stsets. The fit deliberately still censors it: matching them would double-count death.
- **Maintenance cost was overstated by 69%** (`core/process_data.do`). `cost_tx_mnt` billed the blended `cMNT` across the entire `TFI_L1`, i.e. it assumed the patient was on maintenance for the whole L1-to-L2 gap. It now bills a modelled episode: `sim_mnr.do` draws the regimen (lenalidomide or thalidomide), `sim_mnd.do` draws the **duration** from a parametric survival model (`streg`, log-normal, split by transplant like `L1_TFI` with `i.MNR_L1`, response, `i.ECOGcc` and `i.RISS`), and `process_data.do` **caps** it at the realised `TFI_L1` before billing. The `L1_MND` fit is a normal `stset` on the maintenance start/end events (110/111) restored to the skeleton - `origin(Event1==110) failure(Event1==20 111)` - so `stset` derives the survival time and failure cleanly with no separate censoring variable. The fit adds `ln(TFI_L1)` (complete gaps) with a regimen interaction so lenalidomide duration scales with the gap and thalidomide stays fixed; the no-ASCT arm is restricted to responders (BCR CR/VGPR/PR/MR), since SD/PD do not receive maintenance and would otherwise empty a factor cell on the small sample - `sim_mnd.do` drops the same patients. The cap turns a draw that overshoots the gap into continuous-to-progression maintenance and inherits `sim_mort`'s death curtailment. **Known shortcoming**: the share validates only directionally. Fitted on complete gaps (median ~24 months) but simulated out to ~98-month gaps, the drawn duration undershoots at long gaps - the simulated lenalidomide share rises with the gap (right direction, OOS mid-bands within ~6-8pp of the registry) but too shallowly, so band 4 (42mo+, about half the maintenance cohort) under-shares (in-sample 0.44 against a registry 0.83) while very short gaps over-share via the cap. This under-bills long-gap maintenance and is accepted as the gap-extrapolation limitation (5(5)), not a bug. The whole model still validates out-of-sample at **83.9% (146/174)**, unchanged by the maintenance and refractory work. SIMPLE-FIRST: lenalidomide and thalidomide only, no window and no TTM start offset (an earlier share-of-the-window `betareg` design was withdrawn before release). Affects the **`default`** full-pathway analysis and budget-impact work; `transport_dvd` and `car_t` are `$line 2` and never costed maintenance. Two limitations remain, both pre-existing: pricing is still the blended `cMNT` (only lenalidomide has a PBS maintenance listing), and **later-line maintenance is costed nowhere**. The `rn_mnr()` / `rn_mnd()` slots re-lay-out `mRN` (`rn_K()` 74 -> 76), so patient-level results are not comparable with runs before this change. Evidence and the rejected full decomposition: `docs/refractory.md` section 7.
- **Best clinical response (BCR) was collapsing to a single imputation.** In `prep/multiple_imputation.do` the direct-column LOCF carry-forward filled BCR's m=0 master, so the subsequent `mi update` reset the per-imputation values to that single master value — only ~24 of ~5,000 imputed response rows varied across the 10 imputations (every covariate imputed correctly). Fixed by carrying forward only the imputation columns for BCR (`_cf … nomaster`), leaving the master missing as a registered imputed variable requires; the response now varies correctly (FMI ≈ 0.1–0.45). The two ASCT equations whose estimation sample conditions on imputed BCR (`if BCR != 6` at L1 end, `if BCR_L1 != 6` at ASCT; `prep/risk_equations.do`) now use `esampvaryok`, since the eligible sample legitimately varies across imputations. Out-of-sample validation is unchanged at **139/172 (80.8%)** — the pooled point estimates are stable, so the correction acts on the between-imputation variance (standard errors / prediction intervals), not the point predictions. Coefficients regenerated.

### Known issues
- **Lenalidomide-refractory status is now SIMULATED (both arms shipped).** Previously held out because the generation mix was wrong. Both causes were found and fixed. (a) The maintenance arm used a deterministic tail rule - refractory if the simulated gap minus the maintenance duration fell under 60 days - which was a correct definition but a broken mechanism: `sim_tfi_l1.do` draws the gap truncated just above the duration, so the tail piles up at the bound for LONG maintenance, selecting good-prognosis patients. Replaced by a fitted logit (`LENREFR_MNT`) drawn once per patient. (b) The treatment arm's gate counted only regimens in the analysis's modelled list, so everything in TXR code 0 ('other') was treated as non-lenalidomide - but that bucket is 50.6% lenalidomide at L1 and 33.7% at L4, and the L4/L5+ lists contain no lenalidomide regimen at all. Exposure is now a probability (1 for a modelled code, `P(len | other, line)` for 'other'), folded into the existing draw. **Result:** four of five prevalence checks pass (previously one), refractory-arm composition close to the registry's, and both folds improved. **Remaining:** 3-year OS in the refractory group is still too favourable (61.6% vs 35.1%); the applied hazard is roughly what the fit specifies, and the residual is confounding the model cannot observe - `FISHRisk`, the registry's only cytogenetic variable, was tested and adds nothing (AUC +0.002). `docs/refractory.md`; evidence in `scratch/refractory/_notes.md`.

### Licensing
- **Relicensed to a dual, source-available model from v3.0**: the software is now offered under the **PolyForm Noncommercial License 1.0.0** (free for academic/noncommercial use), with **commercial use — including industry-sponsored regulatory/reimbursement (e.g. PBAC) submissions — by separate licence from Monash University**. Registry-derived data and fitted parameters are additionally governed by the MRDR data agreement. See `LICENSING.md`.

### Incorporated from v2.1
- Vectorised Mata engine, modular `mata_setup.do`, and the comprehensive validation suite — see the [2.1] entry below for detail.

## [2.1] - 2025-10-27

### Added
- **Vectorised Implementation**: Complete rewrite of core simulation engine using Mata vector and matrix operations
- **`core/vector_setup.do`**: New modular vector setup module for efficient data preparation
- **Comprehensive Test Suite**: New validation framework in `tests/` directory
  - `validate_vectors.do`: Validates vectorised implementation against original
  - Additional outcome-specific validation tests
- **Improved Error Handling**: Enhanced validation of vector dimensions and data consistency
- **Performance Metrics**: Built-in validation summaries with detailed reporting

### Changed
- **Core Engine**: Replaced patient-level loops with vectorised Mata operations throughout
- **Data Processing**: All patient characteristics now processed as vectors for simultaneous operations
- **Matrix Operations**: Optimised matrix algebra for risk equation calculations
- **Code Structure**: Modular architecture with clear separation between setup, computation, and validation
- **Repository Organisation**: Modernised Git-based versioning (removed version folders)
- **File Naming**: 
  - `EpiMAP_Start.do` → `run.do` (clearer purpose)
  - Simplified main dispatcher naming
- **Documentation**: Updated README with vectorisation details and performance notes

### Performance
- **Execution Speed**: Significantly faster for large cohorts (10,000+ patients)
- **Memory Efficiency**: Reduced memory overhead through bulk operations
- **Scalability**: Better handling of bootstrap iterations and large-scale simulations

### Fixed
- Memory allocation issues in large simulation runs (now handled via vectorisation)
- Numerical precision edge cases through consistent vector operations
- Random number generation consistency across parallel operations

### Validation
- All vectorised outcomes validated to produce identical results to v2.0
- Comprehensive testing confirms bit-for-bit equivalence with loop-based implementation
- Extended validation suite for ongoing quality assurance

### Technical Details
- Implementation uses Mata's native vector and matrix operations
- Pre-allocated vectors for all patient characteristics (age, sex, ECOG, R-ISS, comorbidities)
- Matrix-based risk equation calculations with element-wise operations
- Consistent random number seeding for reproducibility

## [2.0] - 2025-01

### Added
- Reorganised repository structure with clear version folders (v1.0/, v2.0/)
- Enhanced documentation with detailed user guide
- Improved parameter validation and error checking
- New treatment pathway options for later lines of therapy
- Comprehensive testing framework for model validation
- Better handling of edge cases in survival calculations

### Changed
- **BREAKING**: Repository structure reorganised - models now in version-specific folders
- **BREAKING**: Updated file naming conventions for better clarity
- Improved simulation performance through code optimisation
- Enhanced random number generation for better reproducibility
- Updated coefficient matrix structure for additional parameters
- Refined treatment-free interval calculations

### Fixed
- Corrected edge case in ASCT eligibility determination
- Fixed rare numerical precision issues in survival probability calculations
- Resolved memory allocation issues in large simulation runs
- Improved handling of missing data in patient characteristics

### Model Updates
- Updated risk equations based on latest MRDR data analysis
- Refined treatment regimen proportions to reflect current clinical practice
- Enhanced Best Clinical Response prediction accuracy
- Improved Overall Survival curve fitting

### Documentation
- Complete rewrite of user documentation
- Added model specification document with technical details
- Created parameter reference guide
- Added validation report with benchmark results

## [1.0] - 2024-08

### Added
- Initial public release of EpiMAP Myeloma simulation model
- Discrete-event simulation framework for multiple myeloma outcomes
- 30 risk equations based on MRDR patient-level data
- Support for up to 9 Lines of Therapy (LoTs)
- Hypothetical patient dataset with 1,000 patients for testing
- Complete documentation and usage instructions

### Model Features
- **Patient Characteristics**: Age, sex, ECOG performance score, ISS stage
- **Treatment Pathways**: Comprehensive modelling of treatment sequences
- **Survival Analysis**: Parametric survival models for overall survival
- **Clinical Response**: Ordered logit models for Best Clinical Response
- **ASCT Support**: Separate analysis paths for transplant-eligible patients
- **Maintenance Therapy**: Modelling of post-induction maintenance treatment

### Risk Equations
1. Overall Survival
2. Planned ASCT eligibility
3. Diagnosis to treatment interval
4. LoT 1 chemotherapy regimen selection
5-7. LoT 1 chemotherapy duration (with manual splines for ASCT patients)
8. LoT 1 chemotherapy duration (non-ASCT patients)
9. LoT 1 Best Clinical Response
10. Receipt of ASCT
11. ASCT Best Clinical Response
12. Receipt of maintenance therapy
13-14. LoT 1 to LoT 2 treatment-free intervals (ASCT vs non-ASCT)
15-17. LoT 2 treatment pathways and outcomes
18-30. LoTs 3-6+ treatment pathways and outcomes

### Chemotherapy Regimens
- **LoT 1**: VCd (58%), VRd (15%), Other (26%)
- **LoT 2**: Rd (16%), DVd (11%), Other (73%)
- **LoT 3+**: Averaged survival benefit approach

### Technical Specifications
- **Platform**: Stata 15.0 or higher required
- **Programming**: Stata with Mata matrix operations
- **Data Format**: Stata .dta files for patient data
- **Coefficients**: Mata .mmat matrix files for risk equations
- **Output**: Comprehensive simulated patient outcomes dataset

### Validation
- Model calibrated against MRDR registry data (2009-2023)
- Survival curves validated against observed patient outcomes
- Treatment pathway distributions match registry patterns
- Out-of-sample validation with 70/30 split (2,884/1,237 patients)
- 100 bootstrap iterations for robustness testing
- No significant difference in 90% of 120 months post-diagnosis

---

## Migration Guides

### Upgrading from v2.0 to v2.1

**Good News**: v2.1 is fully backward compatible with v2.0 inputs!

**What's Different**:
- Internal implementation is vectorised, but all inputs/outputs remain the same
- No changes to data formats, coefficient files, or simulation parameters
- Results are identical (validated bit-for-bit equivalence)
- Significantly faster performance, especially for large cohorts

**Action Required**:
- None for basic usage - existing scripts will work unchanged
- Optional: Review new test suite in `tests/` for validation examples
- Optional: Update file references if using old naming (`EpiMAP_Start.do` → `run.do`)

**Performance Benefits**:
- ~2-5x faster for typical simulations (depends on cohort size)
- Better scaling for bootstrap analyses
- More efficient memory usage

### Upgrading from v1.0 to v2.1

**Breaking Changes**:
- Repository structure has changed (no more version folders)
- File paths for analyses updated
- Coefficient file organisation modified

**Action Required**:
1. Review new repository structure
2. Update file paths in custom scripts
3. Verify coefficient file locations in `analyses/*/data/coefficients/`
4. Test simulations with small cohort first

**New Features Available**:
- Extended treatment regimens
- Vectorised performance
- Comprehensive validation tools

## Release Notes

Each release includes:
- **Complete Model**: All necessary files to run simulations
- **Documentation**: User guides, technical specifications, and examples
- **Validation Results**: Benchmark tests and model performance metrics
- **Test Data**: Example datasets for testing and validation

## Contact

- **Model Questions**: adam.irving@monash.edu
- **Technical Issues**: Create an issue on GitHub
- **MRDR Data Access**: Visit mrdr.net.au
- **Collaboration**: Contact the research team via email

## Version Naming Convention

From v3.0 onward, **major versions increment with each published paper** (the model state used for that paper) or a significant methodological upgrade; minor versions cover interim features and fixes between papers:
- **Major version** (e.g. 2.0 → 3.0): the model as used for a new published paper, or a significant methodological/architectural upgrade.
- **Minor version** (e.g. 3.0 → 3.1): new features, improvements, or fixes between papers; backward compatible where possible.

Examples:
- v1.0: Initial release
- v2.0: Reorganised structure (breaking changes)
- v2.1: Vectorised implementation (backward compatible)
- v3.0: Calibrated Transport & CRN methods; PBS cost engine; refractory/maintenance subsystem; rebrand to Monash Myeloma Model
