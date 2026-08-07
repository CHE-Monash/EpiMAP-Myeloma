**********
* Monash Myeloma Model - Cohort Pool (template)
*
* Purpose: build-once line-entry pool for a LINE-SPECIFIC decision analysis. Simulate each independent
*          incident population from diagnosis, keep everyone alive at end of L(line-1) who reaches line L
*          inside the case-mix window, blank the line-L regimen (arm-agnostic pool), rebase the diagnosis
*          clock, and save. draw_cohort.do then samples a fixed-size cohort without re-simulating.
* Usage:   do cohort_pool.do -- edit the Configuration block first. Skip for whole-population analyses
*          (they point $data at a population cohort and go straight to simulate.do).
* Notes:   Case-mix is set by the window (who reaches L) and held fixed -- do NOT widen it to change size
*          (that shifts case-mix); size is draw_cohort.do's job. Needs fitted coefficients_<coeffs>.mmat
*          and patients/synthetic_<min>_<max>_<s>.dta (prep/synthetic_1995_2040.do). Output ->
*          analyses/template/patients/cohort_pool_<line>.dta. Example: analyses/transport_dvd/patients/sample_size/.
**********

clear all
set more off

if "$repo_path" != "" cd "$repo_path"   // cd to repo root only if config.do set it
capture run "config.do"     // machine-specific paths (git-ignored; see config.example.do)
macro drop _all

**********
* Configuration  -- EDIT THESE
**********

local n_samples  = 10          // independent population files to pool (max 10)
local start_year = 2020        // case-mix window: who reaches line L in these diagnosis-cohort years
local end_year   = 2025

global analysis  "template"    // <- your analysis
global line      "2"           // <- the DECISION line L (patients must reach it)
global coeffs    "template"    // <- coefficient set (coefficients_<coeffs>.mmat)
global cost_year "2020"
global drate     "0.05"
global boot      "0"

global coefficients_path  "analyses/$analysis/coefficients"
global outcomes_path      "analyses/$analysis/outcomes"
global patients_path      "analyses/$analysis/patients"
global simulated_path     "analyses/$analysis/simulated"

local pool_file "$patients_path/cohort_pool_${line}.dta"

**********
* Load core programs (incl. the shared engine pass)
**********

run "core/run_pipeline.do"

**********
* Build pool: simulate each independent population -> line-L entry -> filter
**********

di as text _n "=== Building cohort pool from `n_samples' independent populations ==="

forval s = 1/`n_samples' {

	mata: mata clear

	global int      "all"            // natural-history setting only; the line-L regimen is blanked below
	global data     "synthetic_`s'" // INDEPENDENT synthetic cohort file s (NOT "synthetic", which is file 1)
	global min_year "1995"
	global max_year "2020"
	global min_id   "1"
	global max_id   "101212"

	qui mata: mata matuse "$coefficients_path/coefficients_$coeffs"

	run_pipeline

	// Keep patients alive at end of the prior line who reach line L in the case-mix window
	keep if MOR_L`= ${line} - 1'E == 0
	keep if YearL${line} >= `start_year' & YearL${line} <= `end_year'

	// Prepare each row as a valid line-L simulation-entry cohort:
	replace State = ${line} * 2                              // enter the pathway at the start of line L
	replace Age_L${line}S = .                                // let the engine recompute age at entry
	replace TXR_L${line} = .                                 // blank the line-L regimen -> arm-agnostic pool
	replace ID = _n
	replace DateDN = td(1jan2020) - (TSD_L${line}S * 12)     // rebase the diagnosis clock to the entry time
	replace YearDN = yofd(DateDN)
	drop DateL* DateMOR YearL* YearMOR c* q* OC_TIME_L TSD_*_ref   // drop realised-pathway/cost columns
	cap drop DateSCT YearSCT

	gen Sample = `s'
	save "$patients_path/pool_part_`s'.dta", replace
}

// Combine the parts into one pool
use "$patients_path/pool_part_1.dta", clear
erase "$patients_path/pool_part_1.dta"
forval s = 2/`n_samples' {
	append using "$patients_path/pool_part_`s'.dta", nolabel
	erase "$patients_path/pool_part_`s'.dta"
}

replace ID = _n
save "`pool_file'", replace

di as text "Pool built: " as result _N as text " line-${line} entrants -> `pool_file'"
tab Sample, missing
summarize TSD_L${line}S, detail
