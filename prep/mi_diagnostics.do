**********
* Monash Myeloma Model - MI diagnostics
*
* Purpose: Missing-data diagnostics for the imputed BASELINE COVARIATES that enter the risk
*          equations -- demographics (Male), performance status (ECOGcc), R-ISS (and its inputs),
*          and comorbidities (CM_CKD/CRD/PLM/DBT) -- plus the first-line response (BCR_L1), read off
*          the saved MRDR Long MI.dta at each variable's own record. Reports per variable: the
*          pre-imputation percent missing, then the relative increase in variance (RVI), the
*          fraction of missing information (FMI), the relative efficiency of the current number of
*          imputations (RE) and the imputation-adjusted degrees of freedom (DF) -- the figures that
*          justify the number of imputations and go in a manuscript's missing-data methods.
* Usage:   do "prep/mi_diagnostics.do"   (needs the MRDR drive; reads $data_path/MRDR Long MI.dta,
*          the main-model output of multiple_imputation.do run with $imp set to the reporting M).
* Notes:   Diagnostics only -- reads the imputed data, never rebuilds it. Deliberately kept out of
*          multiple_imputation.do so the production build (including the 500-iteration bootstrap on
*          the HPC) stays lean and its logs uncluttered. FMI and RE only settle as M grows, so read
*          them off the main M=10 model, not a 2-imputation run. These covariates also enter the
*          OS (streg) and ASCT (logit) equations; for the FMI/DF of those COEFFICIENTS as fitted,
*          add `vartable` / `dftable` to the matching `mi estimate` calls in prep/risk_equations.do.
**********

clear all
if "$repo_path" != "" cd "$repo_path"
capture run "config.do"
set linesize 200                  // keep the wide mi vartable/dftable on one line in the log

* Write a plain-text log to scratch/ (local + git-ignored via *.log; the MRDR drive is not readable
* off-Stata) so the results can be reviewed after a drive run. Named log -> coexists with any log
* already open in an interactive session. Convention documented in CLAUDE.md.
cap mkdir "scratch"
cap log close mi_diag
log using "scratch/mi_diagnostics.log", replace text name(mi_diag)

* Fail fast (and close the log) if the MRDR share is not mounted -- otherwise the first `use` throws
* a bare r(601) and leaves this log open to catch whatever runs next in the session.
capture confirm file "${data_path}/MRDR Long.dta"
if _rc {
    di as error "MRDR data not found at ${data_path}/ -- is the MRDR share mounted?"
    di as error "  Mount/authenticate the drive, confirm 'MRDR Long.dta' is visible, then re-run."
    cap log close mi_diag
    exit 601
}

* Baseline covariates entering the risk equations, at the diagnosis record (Event0 == 3, one row
* per patient), grouped by how the FMI is best summarised:
*   binary / continuous -> mean ;  ordinal (>2 levels) -> proportion (FMI per level).
local mean_vars  Male FISHRisk Albumin SerumB2Microglobulin LactateDehydrogenase CM_CRD CM_PLM CM_DBT eGFR CM_CKD   // 0/1 + continuous; CM_CKD passive (from eGFR). Labs (Albumin/B2M/LDH) feed R-ISS.
local prop_vars  ECOGcc RISS                                      // ordinal; ECOGcc imputed, RISS passive
* Directly-imputed inputs listed for the missingness table (labs feed ISS/R-ISS):
local miss_vars  Male FISHRisk ECOGcc eGFR CM_CRD CM_PLM CM_DBT Albumin SerumB2Microglobulin LactateDehydrogenase

* -- Pre-imputation missingness at diagnosis (context) --
use "${data_path}/MRDR Long.dta", clear
di as text _n(2) "{hline 74}"
di as text "Percent missing at diagnosis (Event0 == 3), before imputation"
di as text "{hline 74}"
misstable summarize `miss_vars' if Event0 == 3

* Response (BCR) missingness. BCR is stochastically imputed only at treatment starts (CStart==1 &
* Duration != ., per multiple_imputation.do); at the line-start records (Event0==L0) missing BCR is
* instead LOCF-carried. An earlier version of this comment said that carryforward is DETERMINISTIC
* and "returns FMI=0" - MEASURED 23 July 2026 AND IT IS NOT. BCR_L1 comes back with FMI up to 0.769
* and BCR_L2 up to 0.791 (see the per-line table further down). _cf fills the imputation columns
* separately with `nomaster', so imputation-to-imputation variation survives the carryforward and
* the equations that use response are NOT losing their uncertainty. Do not re-derive the FMI=0
* claim - a plan to restructure the whole BCR chain was written on it before anyone checked.
di as text _n "  BCR at treatment starts (CStart==1 & Duration != .) -- the imputation sample:"
misstable summarize BCR if CStart == 1 & Duration != .
foreach L in 1 2 {
    di as text _n "  BCR at L`L' start record (Event0 == `L'0) -- LOCF-carried, deterministic:"
    misstable summarize BCR if Event0 == `L'0
}

* THE POST-TRANSPLANT ROW, which this file did not look at until 23 July 2026 and which is the one
* place BCR is imputed by its OWN block (multiple_imputation.do ~line 183, "if Event0 == 100"),
* wrapped in `cap noi' so a failure prints once and is then lost. 548 of 2,135 transplanted patients
* had BCR still missing there at m = 1 - a quarter of the arm - and it was invisible from here for
* two reasons: this loop only covered Event0 == L0, and BCR_SCT reported ZERO missing in mi describe
* because the old `replace BCR_SCT = 0 if BCR_SCT == .' had already absorbed it.
* Downstream, four equations then treated that 0 as either a category or an exclusion. Watch BOTH
* lines: BCR missing here means the imputation block is not covering these rows; BCR_SCT missing
* here is the CORRECT presentation of that, not a new fault.
di as text _n "  BCR at the post-SCT response record (Event0 == 100) -- its own imputation block:"
misstable summarize BCR if Event0 == 100
di as text "  ... and among transplanted patients specifically:"
capture noisily misstable summarize BCR if Event0 == 100 & SCT == 1


* -- MI diagnostics at the current M --
use "${data_path}/MRDR Long MI.dta", clear
mi describe

* BCR_SCT invariant: 0 must mean "no transplant" and nothing else. Checked INSIDE an imputation -
* BCR_SCT is now registered imputed, so m = 0 legitimately retains the original gaps and reading the
* master would report every one of them as a failure. `mi xeq 1:' is style-agnostic; this file's
* data is mi set WIDE, where _mi_m does not exist at all.
di as text _n(2) "{hline 74}"
di as text "BCR_SCT -- 0 must count SCT == 0 patients exactly (checked in m = 1)"
di as text "{hline 74}"
capture noisily mi xeq 1: tab BCR_SCT SCT, missing

di as text _n(2) "{hline 74}"
di as text "Binary / continuous covariates -- variance information (RVI, FMI, rel. efficiency)"
di as text "{hline 74}"
mi estimate, vartable: mean `mean_vars' if Event0 == 3

di as text _n(2) "{hline 74}"
di as text "Binary / continuous covariates -- degrees of freedom (DF)"
di as text "{hline 74}"
mi estimate, dftable: mean `mean_vars' if Event0 == 3

di as text _n(2) "{hline 74}"
di as text "Ordinal covariates (ECOGcc, R-ISS) -- variance information, FMI per level"
di as text "{hline 74}"
mi estimate, vartable: proportion `prop_vars' if Event0 == 3

* Best clinical response, measured over the imputation sample where BCR carries genuine imputation
* uncertainty (CStart==1 & Duration != .): pooled, then by line via CLine (the current-line predictor
* the BCR imputation model uses). Measuring at Event0==10/20 gave FMI=0 (LOCF-carried there).
di as text _n(2) "{hline 74}"
di as text "Best clinical response -- imputation sample (CStart==1 & Duration != .), pooled"
di as text "{hline 74}"
mi estimate, vartable: proportion BCR if CStart == 1 & Duration != .

* THE PER-LINE VARIABLES - what the equations and the engine actually consume.
*
* BCR above carries real imputation uncertainty (FMI 0.06-0.81). The question this answers is
* whether BCR_L1..L9 and BCR_SCT INHERIT it. They are passive derivations of BCR followed by an LOCF
* carryforward, and this file's own header asserts they come back at FMI = 0 - deterministic, so no
* uncertainty about response reaches OS, TXD, TXR, MNT, TFI or MND. That assertion has never been
* printed here; the tables above deliberately measure at CStart instead, which is where FMI is
* non-zero. Measure it rather than assume it.
*
* If FMI is ~0 here, the case for imputing the per-line variables directly is strong and is the main
* argument in scratch/bcr/_notes.md. If it is materially non-zero, that argument largely falls away
* and the refactor is a tidiness change. BCR_SCT is now imputed in its own right, so it is the
* control: it SHOULD show non-zero FMI, and if it does not, the conversion did not take.
di as text _n(2) "{hline 74}"
di as text "Per-line response -- does the imputation uncertainty survive the LOCF carryforward?"
di as text "{hline 74}"
foreach L in 1 2 3 {
    capture confirm variable BCR_L`L'
    if !_rc {
        di as text _n "  BCR_L`L' at its own line-start record (Event0 == `L'0):"
        capture noisily mi estimate, vartable: proportion BCR_L`L' if Event0 == `L'0
    }
}
di as text _n "  BCR_SCT among transplanted patients -- the control, now imputed directly:"
capture noisily mi estimate, vartable: proportion BCR_SCT if Event0 == 100
foreach L in 1 2 {
    di as text _n(2) "{hline 74}"
    di as text "Best clinical response at line `L' (CStart==1 & Duration != . & CLine==`L')"
    di as text "{hline 74}"
    mi estimate, vartable: proportion BCR if CStart == 1 & Duration != . & CLine == `L'
}

* SANITY: does BCR actually vary across the imputations? Zero between-imputation variance everywhere,
* even where ~26% is imputed, usually means the draws collapsed. Count rows whose m-th imputed value
* differs from the first imputation, for BCR and for a known-varying imputed control (eGFR). If BCR is
* ~0 while eGFR is large, the response imputations are (near-)identical -- a real issue in the response
* imputation, not a diagnostic artefact. (M = 10 imputations in the main model.)
di as text _n(2) "{hline 74}"
di as text "SANITY: rows differing from imputation 1 (BCR vs the eGFR control)"
di as text "{hline 74}"
tempvar bd ed
qui gen byte `bd' = 0
qui gen byte `ed' = 0
forvalues m = 2/10 {
    capture confirm variable _`m'_BCR
    if !_rc  qui replace `bd' = 1 if _`m'_BCR  != _1_BCR
    capture confirm variable _`m'_eGFR
    if !_rc  qui replace `ed' = 1 if _`m'_eGFR != _1_eGFR
}
qui count if `bd'
di as text "  rows where BCR differs across imputations:  " %9.0fc r(N)
qui count if `ed'
di as text "  rows where eGFR differs across imputations: " %9.0fc r(N)

di as text _n "Diagnostics log saved to scratch/mi_diagnostics.log"
cap log close mi_diag
