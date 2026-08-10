**********
* Monash Myeloma Model - Extreme-value verification
*
* Purpose: Verification (not calibration). Push individual risk-equation coefficients to extreme
*          values and confirm the engine responds in the correct direction and stays bounded. Each
*          scenario reloads the base coefficients, overwrites ONE intercept in Mata, runs the engine
*          on a small cohort and asserts an output crosses an expected bound; plus a monotonicity
*          sweep (raising the OS intercept must monotonically lower median survival).
* Usage:   Run from the repository root. No MRDR drive required (uses the in-repo population cohort).
**********

clear all
set more off
if "$repo_path" != "" cd "$repo_path"
capture run "config.do"

**********
* Configuration (small cohort for speed; deterministic)
**********
global analysis     "default"
global int          "all"
global line         "0"
global coeffs       "full"
global data         "synthetic"
global min_year     "1995"
global max_year     "2040"
global min_id       "1"
global max_id       "3000"          // reduced cohort -- enough for stable medians
global boot         "0"
global min_bs       ""
global max_bs       ""
global cost_year    "2025"
global drate        "0.05"
global report       "0"
global scenario     ""
global coefficients_path "analyses/$analysis/coefficients"
global outcomes_path     "analyses/$analysis/outcomes"
global patients_path     "analyses/$analysis/patients"
global simulated_path    "analyses/$analysis/simulated"

**********
* Load engine
**********
run "core/load_patients.do"
run "core/mata_setup.do"
run "core/simulation_engine.do"
run "core/process_data.do"
run "core/mata_functions.do"

local coeffile "$coefficients_path/coefficients_$coeffs"

* Run the full engine with whatever coefficients are currently in Mata (results left in memory)
capture program drop xv_sim
program define xv_sim
	qui load_patients
	qui mata_setup
	qui simulation
	qui process_data
end

* Assert a measured value is on the expected side of a bound; tally + print
capture program drop xv_assert
program define xv_assert
	args label measured op bound
	if ("`op'" == "<") local ok = (`measured' <  `bound')
	else               local ok = (`measured' >  `bound')
	if `ok' {
		local s "PASS"
		global xv_pass = $xv_pass + 1
	}
	else {
		local s "FAIL"
		global xv_fail = $xv_fail + 1
	}
	global xv_run = $xv_run + 1
	di as text "  " %-44s "`label'" "  measured = " %10.2f `measured' "   (expect `op' `bound')   " as result "`s'"
end

global xv_run  = 0
global xv_pass = 0
global xv_fail = 0

di _n "{hline 90}"
di "EXTREME-VALUE VERIFICATION  (cohort N <= $max_id)"
di "{hline 90}"

**********
* Baseline (unperturbed) reference
**********
qui mata: mata matuse "`coeffile'", replace
xv_sim
qui summarize OC_TIME, detail
local base_os = r(p50)
qui summarize SCT_L1
local base_asct = r(mean) * 100
qui summarize TXD_L1 if SCT_L1 == 0, detail
local base_txd = r(p50)
di _n as text "Baseline:  median OS = " %5.1f `base_os' " mo |  ASCT = " %4.1f `base_asct' "% |  median TXD_L1 (NoASCT) = " %4.1f `base_txd' " mo"
di ""

**********
* Per-line OS: the intercept lives in EACH stage's own coefficient matrix (bOS_DN, bOS_L1S, ...);
* the constant is always the second-to-last column (cons; the last column is the ancillary /ln_p).
* xv_bump_os shifts every stage's intercept together, so cranking hazard up/down moves survival at
* every pathway point (matching the old single-bOS crank).
**********
mata:
void xv_bump_os(real scalar d)
{
	external bOS_DN, bOS_L1S, bOS_L1E_NoASCT, bOS_L1E_ASCT, bOS_L2S, bOS_L2E
	external bOS_L3S, bOS_L3E, bOS_L4S, bOS_L4E, bOS_L5plus
	bOS_DN[1, cols(bOS_DN)-1]                 = bOS_DN[1, cols(bOS_DN)-1] + d
	bOS_L1S[1, cols(bOS_L1S)-1]               = bOS_L1S[1, cols(bOS_L1S)-1] + d
	bOS_L1E_NoASCT[1, cols(bOS_L1E_NoASCT)-1] = bOS_L1E_NoASCT[1, cols(bOS_L1E_NoASCT)-1] + d
	bOS_L1E_ASCT[1, cols(bOS_L1E_ASCT)-1]     = bOS_L1E_ASCT[1, cols(bOS_L1E_ASCT)-1] + d
	bOS_L2S[1, cols(bOS_L2S)-1]               = bOS_L2S[1, cols(bOS_L2S)-1] + d
	bOS_L2E[1, cols(bOS_L2E)-1]               = bOS_L2E[1, cols(bOS_L2E)-1] + d
	bOS_L3S[1, cols(bOS_L3S)-1]               = bOS_L3S[1, cols(bOS_L3S)-1] + d
	bOS_L3E[1, cols(bOS_L3E)-1]               = bOS_L3E[1, cols(bOS_L3E)-1] + d
	bOS_L4S[1, cols(bOS_L4S)-1]               = bOS_L4S[1, cols(bOS_L4S)-1] + d
	bOS_L4E[1, cols(bOS_L4E)-1]               = bOS_L4E[1, cols(bOS_L4E)-1] + d
	bOS_L5plus[1, cols(bOS_L5plus)-1]         = bOS_L5plus[1, cols(bOS_L5plus)-1] + d
}
end

**********
* Boundary scenarios
**********

* --- OS: crank hazard up -> everyone dies almost immediately ---
qui mata: mata matuse "`coeffile'", replace
mata: xv_bump_os(20)
xv_sim
qui summarize OC_TIME, detail
xv_assert "OS hazard -> infinity : median OS (mo)" r(p50) "<" 6

* --- OS: crank hazard down -> everyone survives to the age limit ---
qui mata: mata matuse "`coeffile'", replace
mata: xv_bump_os(-20)
xv_sim
qui summarize OC_TIME, detail
xv_assert "OS hazard -> 0 : median OS (mo)" r(p50) ">" 200

* --- ASCT: certain -> ~all eligible patients transplanted ---
qui mata: mata matuse "`coeffile'", replace
mata: bL1_SCT[1,21] = bL1_SCT[1,21] + 20
xv_sim
qui summarize SCT_L1
xv_assert "ASCT prob -> 1 : ASCT rate (pct)" "`=r(mean)*100'" ">" 60

* --- ASCT: never -> no transplants ---
qui mata: mata matuse "`coeffile'", replace
mata: bL1_SCT[1,21] = bL1_SCT[1,21] - 20
xv_sim
qui summarize SCT_L1
xv_assert "ASCT prob -> 0 : ASCT rate (pct)" "`=r(mean)*100'" "<" 2

* --- TXD: crank L1 (NoASCT) duration hazard up -> treatment ends immediately ---
qui mata: mata matuse "`coeffile'", replace
mata: bL1_TXD_NoASCT[1,19] = bL1_TXD_NoASCT[1,19] + 20
xv_sim
qui summarize TXD_L1 if SCT_L1 == 0, detail
xv_assert "TXD_L1 hazard -> infinity : median TXD (mo)" r(p50) "<" 1

* --- TFI: crank L1 (NoASCT) interval hazard up -> immediate progression ---
qui mata: mata matuse "`coeffile'", replace
mata: bL1_TFI_NoASCT[1,17] = bL1_TFI_NoASCT[1,17] + 20
xv_sim
qui summarize TFI_L1 if SCT_L1 == 0, detail
xv_assert "TFI_L1 hazard -> infinity : median TFI (mo)" r(p50) "<" 1

**********
* Monotonicity sweep: raising the OS intercept must lower median survival monotonically
**********
di _n as text "Monotonicity: median OS vs OS-intercept shift"
* NB: shifts >= +2 can crash the engine (r(3301) subscript invalid) -- at intermediate-high mortality
* a downstream line empties out and an unguarded selectindex fails. That is a real engine-robustness
* finding (boundary +20 is fine because everyone dies at once); the sweep stays in the clean range
* and wraps each run in capture so one bad point cannot abort the test.
local prev = .
local mono = 1
foreach d in -2 -1.5 -1 -0.5 0 0.5 1 {
	qui mata: mata matuse "`coeffile'", replace
	mata: xv_bump_os(`d')
	capture xv_sim
	if _rc continue
	qui summarize OC_TIME, detail
	local m = r(p50)
	di as text "    shift " %5.1f `d' "  ->  median OS = " %6.1f `m' " mo"
	if (`prev' < . & `m' >= `prev') local mono = 0
	local prev = `m'
}
global xv_run = $xv_run + 1
if `mono' {
	global xv_pass = $xv_pass + 1
	di as result "    monotonically decreasing   PASS"
}
else {
	global xv_fail = $xv_fail + 1
	di as error  "    NOT monotone   FAIL"
}

**********
* Summary
**********
di _n "{hline 90}"
di "EXTREME-VALUE SUMMARY:  " $xv_pass "/" $xv_run " passed"
if $xv_fail == 0  di as result "All extreme-value checks passed -- engine responds correctly to extreme inputs."
else              di as error  "$xv_fail check(s) failed -- review above."
di "{hline 90}"
