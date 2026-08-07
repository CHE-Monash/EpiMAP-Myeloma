**********
* Monash Myeloma Model - Validate Outcomes
*
* Purpose: Shared comparison engine. Compares a simulated dataset to observed-outcome target CSVs (OS,
*          BCR, TXD/TFI by horizon survival, pathways by competing-risks CIF, ASCT among L1-end
*          reachers) and prints a pass/fail summary. Driven by $val_targets (target CSV directory) and
*          $val_simfile (simulated .dta).
* Notes:   Normally invoked via analyses/default/validate_outsample.do, which sets both globals and calls it.
**********

if "$repo_path" != "" cd "$repo_path"   // cd to repo root only if config.do set it; a bare cd "" goes to home on Mac/Unix
capture run "config.do"     // machine-specific paths (git-ignored; see config.example.do)

* The caller MUST set both the target-CSV directory and the simulated dataset -- there is no in-sample
* default (the standalone registry acceptance test was retired; see _notes -> Validation). Normally
* invoked via analyses/default/validate_outsample.do, not directly.
if "$val_targets" == "" | "$val_simfile" == "" {
	di as error "Set globals \$val_targets (target CSV directory) and \$val_simfile (simulated .dta) first."
	di as error "This comparison engine is normally run via analyses/default/validate_outsample.do."
	exit 198
}


**********
// Load Benchmarks
**********

di _n "{hline 80}"
di "LOADING BENCHMARKS"
di "{hline 80}"

// NOTE: the Censored column is carried through mkmat only so the column positions below line up
// with the target CSVs -- NOTHING in this file scores it. The OS tests read Y3/Y5 (cols 5 and 7);
// the TXD/TFI tests read M12/M24 (cols 7 and 8). No pass/fail verdict depends on Censored, so
// don't read one into it: it is registry follow-up bookkeeping, not a model prediction.

// OS benchmarks
import delimited "${val_targets}/os_l1_noasct.csv", clear case(preserve)
mkmat N Median Y1 Y2 Y3 Y4 Y5 Y6 Y7 Y8 Y10 Censored, matrix(OS_L1_NoASCT_bench)

import delimited "${val_targets}/os_asct.csv", clear case(preserve)
mkmat N Median Y1 Y2 Y3 Y4 Y5 Y6 Y7 Y8 Y10 Censored, matrix(OS_ASCT_bench)

import delimited "${val_targets}/os_l2.csv", clear case(preserve)
mkmat N Median Y1 Y2 Y3 Y4 Y5 Y6 Y7 Y8 Y10 Censored, matrix(OS_L2_bench)

import delimited "${val_targets}/os_l3.csv", clear case(preserve)
mkmat N Median Y1 Y2 Y3 Y4 Y5 Y6 Y7 Y8 Y10 Censored, matrix(OS_L3_bench)

// GUARD: a benchmark whose entire N column is 0/blank means its target CSV is EMPTY (failed to
// regenerate) -- NOT that the model is wrong. Flag it loudly so an empty target can't masquerade as a
// model regression (this once turned a fine model into a "55.8%" all-FAIL OS-by-BCR run).
global empty_targets ""
foreach bm in OS_L1_NoASCT OS_ASCT OS_L2 OS_L3 {
	mata: st_numscalar("nbench", colsum(editmissing(st_matrix("`bm'_bench")[., 1], 0)))
	if nbench == 0 {
		global empty_targets "$empty_targets `bm'"
		di as error "  WARNING: benchmark `bm' has N=0 in every stratum -- its target CSV is EMPTY."
		di as error "           This is a TARGET-regeneration failure, not a model regression; re-run the targets."
	}
}

// Whole-population OS from diagnosis (optional; older target sets without it still run)
capture confirm file "${val_targets}/os_wholepop.csv"
if _rc == 0 {
	import delimited "${val_targets}/os_wholepop.csv", clear case(preserve)
	mkmat N Median Y1 Y2 Y3 Y4 Y5 Y6 Y7 Y8 Y10 Censored, matrix(OS_All_bench)
}

// Whole-population OS by comorbidity burden (0/1/2+); optional
capture confirm file "${val_targets}/os_wholepop_cm.csv"
if _rc == 0 {
	import delimited "${val_targets}/os_wholepop_cm.csv", clear case(preserve)
	mkmat N Median Y1 Y2 Y3 Y4 Y5 Y6 Y7 Y8 Y10 Censored, matrix(OS_CM_bench)
}

// Lenalidomide-refractory (treatment lines); optional. Prevalence by line and OS from L2 by status.
capture confirm file "${val_targets}/lenrefr.csv"
if _rc == 0 {
	import delimited "${val_targets}/lenrefr.csv", clear case(preserve)
	mkmat N PctRefr, matrix(LENREFR_bench)
}
capture confirm file "${val_targets}/os_lenrefr.csv"
if _rc == 0 {
	import delimited "${val_targets}/os_lenrefr.csv", clear case(preserve)
	mkmat N Median Y1 Y2 Y3 Y4 Y5 Y6 Y7 Y8 Y10 Censored, matrix(OS_LENREFR_bench)
}

// BCR distributions
import delimited "${val_targets}/bcr.csv", clear case(preserve)
mkmat N CR VG PR MR SD PD, matrix(BCR_bench)

// TXD benchmarks
import delimited "${val_targets}/txd_l1_noasct.csv", clear case(preserve)
mkmat N Mean Median P25 P75 Censored M12 M24, matrix(TXD_L1_NoASCT_bench)

import delimited "${val_targets}/txd_l1_asct.csv", clear case(preserve)
mkmat N Mean Median P25 P75 Censored M12 M24, matrix(TXD_L1_ASCT_bench)

// MND benchmark - maintenance duration (KM median months) by regimen. Optional, so target sets
// generated before the maintenance work still validate. NOTE the column set changed when the
// share-of-gap metric was retired: an older mnd_l1.csv carries MNR GapBand N Mean Median P25 P75
// and will not mkmat here, which is deliberate - it would otherwise be scored as if it were the
// new metric. Regenerate the targets.
capture confirm file "${val_targets}/mnd_l1.csv"
if _rc == 0 {
	import delimited "${val_targets}/mnd_l1.csv", clear case(preserve)
	capture mkmat MNR N Failures Median P25 P75, matrix(MND_L1_bench)
	if _rc {
		di as error "  WARNING: mnd_l1.csv is in the OLD share-by-gap-band format and was not loaded."
		di as error "           Re-run prep/generate_benchmarks.do to rebuild it as duration."
	}
}

* L2-L4 TXD (now with M12/M24 on-treatment cols); optional so older target sets still run
foreach L in 2 3 4 {
	capture confirm file "${val_targets}/txd_l`L'.csv"
	if _rc == 0 {
		import delimited "${val_targets}/txd_l`L'.csv", clear case(preserve)
		mkmat N Mean Median P25 P75 Censored M12 M24, matrix(TXD_L`L'_bench)
	}
}

// TFI benchmarks
import delimited "${val_targets}/tfi_l1_noasct.csv", clear case(preserve)
mkmat N Mean Median P25 P75 Censored M12 M24, matrix(TFI_L1_NoASCT_bench)

import delimited "${val_targets}/tfi_l1_asct.csv", clear case(preserve)
mkmat N Mean Median P25 P75 Censored M12 M24, matrix(TFI_L1_ASCT_bench)

import delimited "${val_targets}/tfi_l2.csv", clear case(preserve)
mkmat N Mean Median P25 P75 Censored M12 M24, matrix(TFI_L2_bench)

import delimited "${val_targets}/tfi_l3.csv", clear case(preserve)
mkmat N Mean Median P25 P75 Censored, matrix(TFI_L3_bench)

// Pathways
import delimited "${val_targets}/pathways.csv", clear case(preserve)
mkmat N ASCT L2 L3 L4 L5 L6 L7 L8, matrix(Pathways_bench)

di "Benchmarks loaded successfully"

**********
// Load Simulation
**********

use "$val_simfile", clear

local tests_run = 0
local tests_passed = 0
local tests_failed = 0


**********
// Overall Survival
**********

local tolerance = 0.10

// OS whole-population (from diagnosis, all patients) at 3/5/10yr -- the aggregate check the
// per-line/BCR breakdown lacks. Only runs if the os_wholepop.csv target is present.
capture confirm matrix OS_All_bench
if _rc == 0 {
	n di "OS WHOLE-POPULATION (from diagnosis)  | Benchmark | Simulated | Diff   | Pass?"
	qui stset OC_TIME, failure(OC_MORT==1) id(ID)
	qui sts generate surv_temp = s
	foreach yr in 3 5 10 {
		local col = cond(`yr'==3, 5, cond(`yr'==5, 7, 11))
		local mo  = `yr' * 12
		local bench = OS_All_bench[1, `col'] * 100
		qui summarize surv_temp if _t <= `mo'
		if r(N) > 0 & !missing(`bench') {
			local sim = r(min) * 100
			local diff = `sim' - `bench'
			if abs(`diff') <= `tolerance' * 100 {
				local status "PASS"
				local tests_passed = `tests_passed' + 1
			}
			else {
				local status "FAIL"
				local tests_failed = `tests_failed' + 1
			}
			local tests_run = `tests_run' + 1
			n di "  " %2.0f `yr' "-year (all)                      | " %8.1f `bench' "% | " %8.1f `sim' "% | " %5.1f `diff' "% | `status'"
		}
	}
	drop surv_temp
	n di _n
}

// OS whole-population by comorbidity burden (CKD+CRD+PLM+DBT: 0/1/2+) -- tests whether the model's
// comorbidity differentiation is calibrated (the aggregate check above cannot see this).
capture confirm matrix OS_CM_bench
if _rc == 0 {
	local _cmok = 1
	foreach v in CM_CKD CM_CRD CM_PLM CM_DBT {
		capture confirm variable `v'
		if _rc local _cmok = 0
	}
	if `_cmok' {
		n di "OS BY COMORBIDITY BURDEN (from diagnosis)  | Benchmark | Simulated | Diff   | Pass?"
		qui capture drop _cmn _cmg
		qui gen byte _cmn = CM_CKD + CM_CRD + CM_PLM + CM_DBT
		qui gen byte _cmg = cond(_cmn==0, 0, cond(_cmn==1, 1, 2))
		qui stset OC_TIME, failure(OC_MORT==1) id(ID)
		forvalues g = 0/2 {
			local lbl : word `=`g'+1' of "none" "one " "two+"
			qui capture drop surv_temp
			qui sts generate surv_temp = s if _cmg == `g'
			foreach yr in 3 5 10 {
				local col = cond(`yr'==3, 5, cond(`yr'==5, 7, 11))
				local mo  = `yr' * 12
				local bench = OS_CM_bench[`=`g'+1', `col'] * 100
				qui summarize surv_temp if _t <= `mo' & _cmg == `g'
				if r(N) > 0 & !missing(`bench') {
					local sim = r(min) * 100
					local diff = `sim' - `bench'
					local status = cond(abs(`diff') <= `tolerance' * 100, "PASS", "FAIL")
					if "`status'" == "PASS" local tests_passed = `tests_passed' + 1
					else                    local tests_failed = `tests_failed' + 1
					local tests_run = `tests_run' + 1
					n di "  " %2.0f `yr' "-year (`lbl')                    | " %8.1f `bench' "% | " %8.1f `sim' "% | " %5.1f `diff' "% | `status'"
				}
			}
		}
		qui drop surv_temp _cmn _cmg
		n di _n
	}
}

// Lenalidomide-refractory prevalence by line (entry state).
// LenRefr_L* IS the union: the engine keeps ONE latched flag that both sim_lenrefr.do (treatment)
// and sim_mnt_refr.do (maintenance) write into, since the two carry the same conditional OS hazard
// (test Tx = Mnt p = 0.97). An earlier version assembled the union here from two exported flags. -- scores the GENERATION model:
// does the engine make the right share refractory as it accrues across lines (docs/refractory.md 4.7)?
capture confirm matrix LENREFR_bench
if _rc == 0 {
	local _lrok = 1
	forvalues l = 1/6 {
		capture confirm variable LenRefr_L`l'
		if _rc local _lrok = 0
	}
	if _rc local _lrok = 0
	if `_lrok' {
		n di "LEN-REFRACTORY PREVALENCE BY LINE (true=Tx|Mnt) | Benchmark | Simulated | Diff   | Pass?"
		forvalues l = 1/6 {
			local bench = LENREFR_bench[`l', 2]
			qui count if !missing(LenRefr_L`l')
			local denom = r(N)
			if `denom' > 0 & !missing(`bench') {
				// TRUE len-refractory = treatment OR maintenance; maintenance applies only from L2
				qui count if !missing(LenRefr_L`l') & LenRefr_L`l' == 1
				local sim = 100 * r(N) / `denom'
				local diff = `sim' - `bench'
				local status = cond(abs(`diff') <= `tolerance' * 100, "PASS", "FAIL")
				if "`status'" == "PASS" local tests_passed = `tests_passed' + 1
				else                    local tests_failed = `tests_failed' + 1
				local tests_run = `tests_run' + 1
				n di "  L`l'                                       | " %8.1f `bench' "% | " %8.1f `sim' "% | " %5.1f `diff' "% | `status'"
			}
		}
		n di _n
	}
	else n di "LEN-REFRACTORY prevalence: SKIPPED - LenRefr_L* not in the simulated data (re-run simulate)." _n
}
else n di "LEN-REFRACTORY prevalence: SKIPPED - no lenrefr.csv target (re-run test_targets / generate_benchmarks)." _n

// OS from L2 start, split by len-refractory-at-L2 (refr_len_l2) -- scores the CONSUMPTION /
// redistribution: does the simulated subgroup OS gap match observed (docs/refractory.md 5.6)?
capture confirm matrix OS_LENREFR_bench
if _rc == 0 {
	capture confirm variable refr_len_l2
	if _rc == 0 {
		n di "OS FROM L2 BY LEN-REFRACTORY STATUS (true=Tx|Mnt) | Benchmark | Simulated | Diff   | Pass?"
		// TRUE len-refractory at L2 entry = treatment OR maintenance; missing for non-L2-reachers
		qui cap drop _lrU2
		qui gen byte _lrU2 = (refr_len_l2 == 1) if !missing(refr_len_l2)
		qui cap drop OC_TIME_L2S
		qui gen OC_TIME_L2S = OC_TIME - TSD_L2S
		qui stset OC_TIME_L2S, failure(OC_MORT==1) id(ID)
		forvalues r = 0/1 {
			local lbl : word `=`r'+1' of "not-refr" "refr    "
			qui capture drop surv_temp
			qui sts generate surv_temp = s if _lrU2 == `r'
			foreach yr in 3 5 10 {
				local col = cond(`yr'==3, 5, cond(`yr'==5, 7, 11))
				local mo  = `yr' * 12
				local bench = OS_LENREFR_bench[`=`r'+1', `col'] * 100
				qui summarize surv_temp if _t <= `mo' & _lrU2 == `r'
				if r(N) > 0 & !missing(`bench') {
					local sim = r(min) * 100
					local diff = `sim' - `bench'
					local status = cond(abs(`diff') <= `tolerance' * 100, "PASS", "FAIL")
					if "`status'" == "PASS" local tests_passed = `tests_passed' + 1
					else                    local tests_failed = `tests_failed' + 1
					local tests_run = `tests_run' + 1
					n di "  " %2.0f `yr' "-year (`lbl')                  | " %8.1f `bench' "% | " %8.1f `sim' "% | " %5.1f `diff' "% | `status'"
				}
			}
			qui drop surv_temp
		}
		n di _n
	}
	else n di "OS by len-refractory status: SKIPPED - refr_len_l2 not in the simulated data (re-run simulate)." _n
}
else n di "OS by len-refractory status: SKIPPED - no os_lenrefr.csv target (re-run test_targets / generate_benchmarks)." _n

// OS by BCR_L1 (NoASCT)
qui cap gen OC_TIME_L1S = OC_TIME - TSD_L1S
qui stset OC_TIME_L1S, failure(OC_MORT==1) id(ID)

* BCR | Timepoint | Benchmark | Simulated | Diff   | Pass?
qui forvalues bcr = 1/6 {
	// 3-year survival
	local bench_3yr = OS_L1_NoASCT_bench[`bcr', 5] * 100
	
	qui sts generate surv_temp = s if BCR_L1 == `bcr' & SCT_L1 == 0
	qui summarize surv_temp if BCR_L1 == `bcr' & SCT_L1 == 0 & _t <= 36
	if r(N) > 0 {
		local sim_3yr = r(min) * 100
		local diff = `sim_3yr' - `bench_3yr'
		
		if abs(`diff') <= `tolerance' * 100 {
			local status "PASS"
			local tests_passed = `tests_passed' + 1
		}
		else {
			local status "FAIL"
			local tests_failed = `tests_failed' + 1
		}
		
		local tests_run = `tests_run' + 1
		
		n di "  BCR=`bcr' | 3-year    | " %8.1f `bench_3yr' "% | " %8.1f `sim_3yr' "% | " %5.1f `diff' "% | `status'"
	}
	drop surv_temp

	// 5-year survival
	local bench_5yr = OS_L1_NoASCT_bench[`bcr', 7] * 100
	
	qui sts generate surv_temp = s if BCR_L1 == `bcr' & SCT_L1 == 0
	qui summarize surv_temp if BCR_L1 == `bcr' & SCT_L1 == 0 & _t <= 60
	qui if r(N) > 0 {
		local sim_5yr = r(min) * 100
		local diff = `sim_5yr' - `bench_5yr'
		
		if abs(`diff') <= `tolerance' * 100 {
			local status "PASS"
			local tests_passed = `tests_passed' + 1
		}
		else {
			local status "FAIL"
			local tests_failed = `tests_failed' + 1
		}
		
		local tests_run = `tests_run' + 1
		
		n di "  BCR=`bcr' | 5-year    | " %8.1f `bench_5yr' "% | " %8.1f `sim_5yr' "% | " %5.1f `diff' "% | `status'"
		n di _n
	}
	drop surv_temp
}

// OS by BCR_SCT
qui cap gen OC_TIME_L1E = OC_TIME - TSD_L1E
qui stset OC_TIME_L1E, failure(OC_MORT==1) id(ID)

* BCR | Timepoint | Benchmark | Simulated | Diff   | Pass?
qui forvalues bcr = 1/4 {
	// 3-year survival
	local bench_3yr = OS_ASCT_bench[`bcr', 5] * 100
	
	qui sts generate surv_temp = s if BCR_SCT == `bcr'
	qui summarize surv_temp if BCR_SCT == `bcr' & _t <= 36
	if r(N) > 0 {
		local sim_3yr = r(min) * 100
		local diff = `sim_3yr' - `bench_3yr'
		
		if abs(`diff') <= `tolerance' * 100 {
			local status "PASS"
			local tests_passed = `tests_passed' + 1
		}
		else {
			local status "FAIL"
			local tests_failed = `tests_failed' + 1
		}
		
		local tests_run = `tests_run' + 1
		
		n di "  BCR=`bcr' | 3-year    | " %8.1f `bench_3yr' "% | " %8.1f `sim_3yr' "% | " %5.1f `diff' "% | `status'"
	}
	drop surv_temp
	
	// 5-year survival
	local bench_5yr = OS_ASCT_bench[`bcr', 7] * 100
	
	qui sts generate surv_temp = s if BCR_SCT == `bcr'
	qui summarize surv_temp if BCR_SCT == `bcr' & _t <= 60
	qui if r(N) > 0 {
		local sim_5yr = r(min) * 100
		local diff = `sim_5yr' - `bench_5yr'
		
		if abs(`diff') <= `tolerance' * 100 {
			local status "PASS"
			local tests_passed = `tests_passed' + 1
		}
		else {
			local status "FAIL"
			local tests_failed = `tests_failed' + 1
		}
		
		local tests_run = `tests_run' + 1
		
		n di "  BCR=`bcr' | 5-year    | " %8.1f `bench_5yr' "% | " %8.1f `sim_5yr' "% | " %5.1f `diff' "% | `status'"
		n di _n
	}
	drop surv_temp
}

// OS by BCR_L2
qui cap gen OC_TIME_L2S = OC_TIME - TSD_L2S
qui stset OC_TIME_L2S, failure(OC_MORT==1) id(ID)

* BCR | Timepoint | Benchmark | Simulated | Diff   | Pass?
qui forvalues bcr = 1/6 {
	// 3-year survival
	local bench_3yr = OS_L2_bench[`bcr', 5] * 100
	
	qui sts generate surv_temp = s if BCR_L2 == `bcr'
	qui summarize surv_temp if BCR_L2 == `bcr' & _t <= 36
	if r(N) > 0 {
		local sim_3yr = r(min) * 100
		local diff = `sim_3yr' - `bench_3yr'
		
		if abs(`diff') <= `tolerance' * 100 {
			local status "PASS"
			local tests_passed = `tests_passed' + 1
		}
		else {
			local status "FAIL"
			local tests_failed = `tests_failed' + 1
		}
		
		local tests_run = `tests_run' + 1
		
		n di "  BCR=`bcr' | 3-year    | " %8.1f `bench_3yr' "% | " %8.1f `sim_3yr' "% | " %5.1f `diff' "% | `status'"
	}
	drop surv_temp
	
	// 5-year survival
	local bench_5yr = OS_L2_bench[`bcr', 7] * 100
	
	qui sts generate surv_temp = s if BCR_L2 == `bcr'
	qui summarize surv_temp if BCR_L2 == `bcr' & _t <= 60
	qui if r(N) > 0 {
		local sim_5yr = r(min) * 100
		local diff = `sim_5yr' - `bench_5yr'
		
		if abs(`diff') <= `tolerance' * 100 {
			local status "PASS"
			local tests_passed = `tests_passed' + 1
		}
		else {
			local status "FAIL"
			local tests_failed = `tests_failed' + 1
		}
		
		local tests_run = `tests_run' + 1
		
		n di "  BCR=`bcr' | 5-year    | " %8.1f `bench_5yr' "% | " %8.1f `sim_5yr' "% | " %5.1f `diff' "% | `status'"
		n di _n
	}
	drop surv_temp
}

// OS by BCR_L3
qui cap gen OC_TIME_L3S = OC_TIME - TSD_L3S
qui stset OC_TIME_L3S, failure(OC_MORT==1) id(ID)

* BCR | Timepoint | Benchmark | Simulated | Diff   | Pass?
qui forvalues bcr = 1/6 {
	// 3-year survival
	local bench_3yr = OS_L3_bench[`bcr', 5] * 100
	
	qui sts generate surv_temp = s if BCR_L3 == `bcr'
	qui summarize surv_temp if BCR_L3 == `bcr' & _t <= 36
	if r(N) > 0 {
		local sim_3yr = r(min) * 100
		local diff = `sim_3yr' - `bench_3yr'
		
		if abs(`diff') <= `tolerance' * 100 {
			local status "PASS"
			local tests_passed = `tests_passed' + 1
		}
		else {
			local status "FAIL"
			local tests_failed = `tests_failed' + 1
		}
		
		local tests_run = `tests_run' + 1
		
		n di "  BCR=`bcr' | 3-year    | " %8.1f `bench_3yr' "% | " %8.1f `sim_3yr' "% | " %5.1f `diff' "% | `status'"
	}
	drop surv_temp
	
	// 5-year survival
	local bench_5yr = OS_L3_bench[`bcr', 7] * 100
	
	qui sts generate surv_temp = s if BCR_L3 == `bcr'
	qui summarize surv_temp if BCR_L3 == `bcr' & _t <= 60
	qui if r(N) > 0 {
		local sim_5yr = r(min) * 100
		local diff = `sim_5yr' - `bench_5yr'
		
		if abs(`diff') <= `tolerance' * 100 {
			local status "PASS"
			local tests_passed = `tests_passed' + 1
		}
		else {
			local status "FAIL"
			local tests_failed = `tests_failed' + 1
		}
		
		local tests_run = `tests_run' + 1
		
		n di "  BCR=`bcr' | 5-year    | " %8.1f `bench_5yr' "% | " %8.1f `sim_5yr' "% | " %5.1f `diff' "% | `status'"
		n di _n
	}
	drop surv_temp
}

**********
// BCR
**********

local tolerance = 5

// L1
qui count if !missing(BCR_L1)
local n_sim = r(N)

* Line | Metric | Benchmark | Simulated | Diff | Pass?

qui forvalues bcr = 1/6 {
	qui count if BCR_L1 == `bcr'
	local pct_sim = r(N) / `n_sim' * 100
	local pct_bench = BCR_bench[1, `bcr'+1]
	local diff = `pct_sim' - `pct_bench'
	
	if abs(`diff') <= `tolerance' {
		local status "PASS"
		local tests_passed = `tests_passed' + 1
	}
	else {
		local status "FAIL"
		local tests_failed = `tests_failed' + 1
	}
	
	local tests_run = `tests_run' + 1
	
	n di "  L1   | BCR=`bcr' | " %8.1f `pct_bench' "% | " %8.1f `pct_sim' "% | " %5.1f `diff' "% | `status'"
}

//ASCT
qui count if !missing(BCR_SCT)
local n_sim = r(N)

qui forvalues bcr = 1/4 {
	qui count if BCR_SCT == `bcr'
	local pct_sim = r(N) / `n_sim' * 100
	local pct_bench = BCR_bench[2, `bcr'+1]
	local diff = `pct_sim' - `pct_bench'
	
	if abs(`diff') <= `tolerance' {
		local status "PASS"
		local tests_passed = `tests_passed' + 1
	}
	else {
		local status "FAIL"
		local tests_failed = `tests_failed' + 1
	}
	
	local tests_run = `tests_run' + 1
	
	n di "  ASCT | BCR=`bcr' | " %8.1f `pct_bench' "% | " %8.1f `pct_sim' "% | " %5.1f `diff' "% | `status'"
}

// L2
qui count if !missing(BCR_L2)
local n_sim = r(N)

qui forvalues bcr = 1/6 {
	qui count if BCR_L2 == `bcr'
	local pct_sim = r(N) / `n_sim' * 100
	local pct_bench = BCR_bench[3, `bcr'+1]
	local diff = `pct_sim' - `pct_bench'
	
	if abs(`diff') <= `tolerance' {
		local status "PASS"
		local tests_passed = `tests_passed' + 1
	}
	else {
		local status "FAIL"
		local tests_failed = `tests_failed' + 1
	}
	
	local tests_run = `tests_run' + 1
	
	n di "  L2   | BCR=`bcr' | " %8.1f `pct_bench' "% | " %8.1f `pct_sim' "% | " %5.1f `diff' "% | `status'"
}

// L3
qui count if !missing(BCR_L3)
local n_sim = r(N)

qui forvalues bcr = 1/6 {
	qui count if BCR_L3 == `bcr'
	local pct_sim = r(N) / `n_sim' * 100
	local pct_bench = BCR_bench[4, `bcr'+1]
	local diff = `pct_sim' - `pct_bench'
	
	if abs(`diff') <= `tolerance' {
		local status "PASS"
		local tests_passed = `tests_passed' + 1
	}
	else {
		local status "FAIL"
		local tests_failed = `tests_failed' + 1
	}
	
	local tests_run = `tests_run' + 1
	
	n di " L3   | BCR=`bcr' | " %8.1f `pct_bench' "% | " %8.1f `pct_sim' "% | " %5.1f `diff' "% | `status'"
}

**********
// TXD  --  % still on treatment at 12 / 24 months vs benchmark (+/- 10 pp)
**********

di _n _n "{hline 80}"
di "TXD PREDICTIONS (% still on treatment at 12 / 24 months)"
di "{hline 80}"
di _n "Line | BCR | Horizon | Benchmark | Simulated |  Diff | Pass?"
di    "-----|-----|---------|-----------|-----------|-------|------"

local tol_pp = 10

// TXD_L1 NoASCT
qui stset TXD_L1, failure(MOR_L1S==0) id(ID)
qui forvalues bcr = 1/6 {
    capture drop surv_h
    qui sts generate surv_h = s if BCR_L1 == `bcr' & SCT_L1 == 0
    local c = 7
    foreach h in 12 24 {
        local bench = TXD_L1_NoASCT_bench[`bcr', `c']
        qui summarize surv_h if BCR_L1 == `bcr' & SCT_L1 == 0 & _t <= `h'
        if r(N) > 0 & !missing(`bench') {
            local sim = r(min) * 100
            local diff = `sim' - `bench'
            if abs(`diff') <= `tol_pp' {
                local status "PASS"
                local tests_passed = `tests_passed' + 1
            }
            else {
                local status "FAIL"
                local tests_failed = `tests_failed' + 1
            }
            local tests_run = `tests_run' + 1
            n di "L1NA | " %3.0f `bcr' " |  " %4.0f `h' "mo | " %8.1f `bench' "% | " %8.1f `sim' "% | " %5.1f `diff' " | `status'"
        }
        local ++c
    }
    qui drop surv_h
}

// TXD_L1 ASCT
qui forvalues bcr = 1/4 {
    capture drop surv_h
    qui sts generate surv_h = s if BCR_SCT == `bcr' & SCT_L1 == 1
    local c = 7
    foreach h in 12 24 {
        local bench = TXD_L1_ASCT_bench[`bcr', `c']
        qui summarize surv_h if BCR_SCT == `bcr' & SCT_L1 == 1 & _t <= `h'
        if r(N) > 0 & !missing(`bench') {
            local sim = r(min) * 100
            local diff = `sim' - `bench'
            if abs(`diff') <= `tol_pp' {
                local status "PASS"
                local tests_passed = `tests_passed' + 1
            }
            else {
                local status "FAIL"
                local tests_failed = `tests_failed' + 1
            }
            local tests_run = `tests_run' + 1
            n di "L1AS | " %3.0f `bcr' " |  " %4.0f `h' "mo | " %8.1f `bench' "% | " %8.1f `sim' "% | " %5.1f `diff' " | `status'"
        }
        local ++c
    }
    qui drop surv_h
}

// MND - L1 maintenance DURATION by regimen, in months.
//
// This replaced a share-of-the-gap benchmark (MND_L1 / TFI_L1 by regimen x gap band) that was not
// comparable between the two sides: both quantities need an observed L2, so the registry side was
// computed on relapsers only while the simulated side had a closed gap for everybody. See the long
// note in prep/generate_benchmarks.do and docs/refractory.md 5(8).
//
// The registry side is now a Kaplan-Meier median, which uses the whole maintenance population
// including patients still on the drug at the cut. The simulated side has no censoring - every
// drawn maintenance episode is complete - so a plain median is the right comparator for it.
//
// THALIDOMIDE IS ON AN 18-MONTH CLOCK on both sides, but by two different mechanisms: the registry
// target is CENSORED at 18 (prep/generate_benchmarks.do) and the simulated draw is CAPPED at 18
// (core/outcomes/sim_mnd.do). Those are not the same operation - a cap puts mass exactly at 18, a
// censoring puts none there - but they agree on the MEDIAN as long as it falls below 18, which it
// does on both sides (registry 11.14, simulated ~11.4). If either median ever drifts above 18 the
// comparison stops being like-for-like and this check needs rethinking, not just re-reading.
//
// Tolerance is RELATIVE, not absolute, because the two regimens differ by an order of magnitude
// (thalidomide is a fixed course, lenalidomide runs to progression). A fixed +/- N months would be
// trivial for one and impossible for the other.
capture confirm matrix MND_L1_bench
if _rc == 0 {
	local tol_mnd = 0.25
	n di ""
	n di "MND_L1 maintenance duration by regimen (median months; tol +/- " %3.0f `tol_mnd'*100 "%)"
	n di "grp  |    bench |      sim |   diff | status"

	qui forvalues r = 1/3 {
		local g     = MND_L1_bench[`r', 1]
		local bn    = MND_L1_bench[`r', 2]
		local bench = MND_L1_bench[`r', 4]

		// Groups follow $MNR_L1 "1 5": lenalidomide, thalidomide, everything else pooled to 0.
		capture drop mnd_cell
		qui gen byte mnd_cell = 0
		if `g' == 0 qui replace mnd_cell = 1 if MNT == 1 & !mi(MND_L1) & !inlist(MNR_L1, 1, 5)
		else        qui replace mnd_cell = 1 if MNT == 1 & !mi(MND_L1) & MNR_L1 == `g'

		qui count if mnd_cell == 1
		local simn = r(N)

		// A missing bench is either a cell below the N floor or a KM median the registry's
		// follow-up never reached. Both mean the registry cannot say, so neither do we.
		if !missing(`bench') & `simn' > 0 {
			qui summarize MND_L1 if mnd_cell == 1, detail
			local sim  = r(p50)
			local diff = `sim' - `bench'
			if abs(`diff') <= `tol_mnd' * `bench' {
				local status "PASS"
				local tests_passed = `tests_passed' + 1
			}
			else {
				local status "FAIL"
				local tests_failed = `tests_failed' + 1
			}
			local tests_run = `tests_run' + 1
			n di %4.0f `g' " | " %8.2f `bench' " | " %8.2f `sim' " | " %6.2f `diff' " | `status'"
		}
		else if `simn' == 0 {
			// Expected for group 0: sim_mnr only ever assigns the levels the analysis declared
			// ($MNR_L1 "1 5"), so the engine produces no 'other' maintenance at all (docs 7.4).
			// The registry median can be perfectly good and still have nothing to score against.
			n di %4.0f `g' " | " %8.2f `bench' " |        . |      . | skipped (engine produces no such regimen)"
		}
		else {
			n di %4.0f `g' " | " %8.2f `bench' " |        . |      . | skipped (no scoreable registry median, N = " %4.0f `bn' ")"
		}
	}
	capture drop mnd_cell
}

// TXD_L2 / L3 / L4  (no ASCT split at later lines); optional -- only if the target was generated
foreach L in 2 3 4 {
    capture confirm matrix TXD_L`L'_bench
    if _rc continue
    qui stset TXD_L`L', failure(MOR_L`L'S==0) id(ID)
    qui forvalues bcr = 1/6 {
        capture drop surv_h
        qui sts generate surv_h = s if BCR_L`L' == `bcr'
        local c = 7
        foreach h in 12 24 {
            local bench = TXD_L`L'_bench[`bcr', `c']
            qui summarize surv_h if BCR_L`L' == `bcr' & _t <= `h'
            if r(N) > 0 & !missing(`bench') {
                local sim = r(min) * 100
                local diff = `sim' - `bench'
                if abs(`diff') <= `tol_pp' {
                    local status "PASS"
                    local tests_passed = `tests_passed' + 1
                }
                else {
                    local status "FAIL"
                    local tests_failed = `tests_failed' + 1
                }
                local tests_run = `tests_run' + 1
                n di "L`L'   | " %3.0f `bcr' " |  " %4.0f `h' "mo | " %8.1f `bench' "% | " %8.1f `sim' "% | " %5.1f `diff' " | `status'"
            }
            local ++c
        }
        qui drop surv_h
    }
}


**********
// TFI  --  % still treatment-free at 12 / 24 months vs benchmark (+/- 10 pp)
**********

di _n _n "{hline 80}"
di "TFI PREDICTIONS (% still treatment-free at 12 / 24 months)"
di "{hline 80}"
di _n "Line | BCR | Horizon | Benchmark | Simulated |  Diff | Pass?"
di    "-----|-----|---------|-----------|-----------|-------|------"

// TFI_L1 NoASCT
qui stset TFI_L1, failure(MOR_L1E==0) id(ID)
qui forvalues bcr = 1/6 {
    capture drop surv_h
    qui sts generate surv_h = s if BCR_L1 == `bcr' & SCT_L1 == 0
    local c = 7
    foreach h in 12 24 {
        local bench = TFI_L1_NoASCT_bench[`bcr', `c']
        qui summarize surv_h if BCR_L1 == `bcr' & SCT_L1 == 0 & _t <= `h'
        if r(N) > 0 & !missing(`bench') {
            local sim = r(min) * 100
            local diff = `sim' - `bench'
            if abs(`diff') <= `tol_pp' {
                local status "PASS"
                local tests_passed = `tests_passed' + 1
            }
            else {
                local status "FAIL"
                local tests_failed = `tests_failed' + 1
            }
            local tests_run = `tests_run' + 1
            n di "L1NA | " %3.0f `bcr' " |  " %4.0f `h' "mo | " %8.1f `bench' "% | " %8.1f `sim' "% | " %5.1f `diff' " | `status'"
        }
        local ++c
    }
    qui drop surv_h
}

// TFI_L1 ASCT
qui forvalues bcr = 1/4 {
    capture drop surv_h
    qui sts generate surv_h = s if BCR_SCT == `bcr' & SCT_L1 == 1
    local c = 7
    foreach h in 12 24 {
        local bench = TFI_L1_ASCT_bench[`bcr', `c']
        qui summarize surv_h if BCR_SCT == `bcr' & SCT_L1 == 1 & _t <= `h'
        if r(N) > 0 & !missing(`bench') {
            local sim = r(min) * 100
            local diff = `sim' - `bench'
            if abs(`diff') <= `tol_pp' {
                local status "PASS"
                local tests_passed = `tests_passed' + 1
            }
            else {
                local status "FAIL"
                local tests_failed = `tests_failed' + 1
            }
            local tests_run = `tests_run' + 1
            n di "L1AS | " %3.0f `bcr' " |  " %4.0f `h' "mo | " %8.1f `bench' "% | " %8.1f `sim' "% | " %5.1f `diff' " | `status'"
        }
        local ++c
    }
    qui drop surv_h
}

// TFI_L2
qui stset TFI_L2, failure(MOR_L2E==0) id(ID)
qui forvalues bcr = 1/6 {
    capture drop surv_h
    qui sts generate surv_h = s if BCR_L2 == `bcr'
    local c = 7
    foreach h in 12 24 {
        local bench = TFI_L2_bench[`bcr', `c']
        qui summarize surv_h if BCR_L2 == `bcr' & _t <= `h'
        if r(N) > 0 & !missing(`bench') {
            local sim = r(min) * 100
            local diff = `sim' - `bench'
            if abs(`diff') <= `tol_pp' {
                local status "PASS"
                local tests_passed = `tests_passed' + 1
            }
            else {
                local status "FAIL"
                local tests_failed = `tests_failed' + 1
            }
            local tests_run = `tests_run' + 1
            n di "L2   | " %3.0f `bcr' " |  " %4.0f `h' "mo | " %8.1f `bench' "% | " %8.1f `sim' "% | " %5.1f `diff' " | `status'"
        }
        local ++c
    }
    qui drop surv_h
}

// TFI_L2 response ordering (better response -> longer treatment-free interval), tested as
// % treatment-free at 12 months for the clearly-separated categories CR > VGPR > PR.
// (MR/SD progress almost immediately and are not reliably ordered -- even the registry medians
//  do not separate them, so the strict 5-way ordering is not asserted.)
di _n "Checking TFI_L2 ordering at 12 months (CR > VG > PR)..."
foreach bcr in 1 2 3 {
    capture drop surv_h
    qui sts generate surv_h = s if BCR_L2 == `bcr'
    qui summarize surv_h if BCR_L2 == `bcr' & _t <= 12
    local s12_`bcr' = cond(r(N) > 0, r(min) * 100, .)
    qui drop surv_h
}
if (`s12_1' > `s12_2') & (`s12_2' > `s12_3') {
    local status "PASS"
    local tests_passed = `tests_passed' + 1
}
else {
    local status "FAIL"
    local tests_failed = `tests_failed' + 1
}
local tests_run = `tests_run' + 1
di "  CR=" %4.1f `s12_1' "% > VG=" %4.1f `s12_2' "% > PR=" %4.1f `s12_3' "%   Result: `status'"


**********
// Pathways
**********

di _n _n "{hline 80}"
di "TEST 5: TREATMENT PATHWAYS"
di "{hline 80}"

di _n "Line | Benchmark | Simulated | Diff | Pass?"
di "-----|-----------|-----------|------|-------"

qui count
local n_total = r(N)
local tolerance = 5

// ASCT -- transplant is a decision taken AT THE END OF L1, so the rate is among patients who
// reach L1 end (alive past L1 induction), matching the benchmark denominator in generate_benchmarks.
// (Dividing by all patients would dilute it by the ~8% who die during induction and never reach
//  the transplant decision.)
qui count if OC_TIME > TSD_L1E & !missing(TSD_L1E)
local n_l1e = r(N)
qui count if SCT_L1 == 1
local pct_sim = r(N) / `n_l1e' * 100
local pct_bench = Pathways_bench[1, 2]
local diff = `pct_sim' - `pct_bench'

if abs(`diff') <= `tolerance' {
	local status "PASS"
	local tests_passed = `tests_passed' + 1
}
else {
	local status "FAIL"
	local tests_failed = `tests_failed' + 1
}

local tests_run = `tests_run' + 1

di "ASCT | " %8.1f `pct_bench' "% | " %8.1f `pct_sim' "% | " %5.1f `diff' "% | `status'"

// Lines 2-5: conditional per-transition progression P(reach L | reached L-1); denominator = patients
// who reached the previous line (BCR_L{prev}). Matches the benchmark's per-transition AJ CIF.
qui forvalues line = 2/5 {
	local prev = `line' - 1
	qui count if !missing(BCR_L`prev')
	local denom = r(N)
	qui count if !missing(BCR_L`line')
	local pct_sim = cond(`denom' > 0, r(N) / `denom' * 100, .)
	local pct_bench = Pathways_bench[1, `line'+1]
	local diff = `pct_sim' - `pct_bench'

	if abs(`diff') <= `tolerance' {
		local status "PASS"
		local tests_passed = `tests_passed' + 1
	}
	else {
		local status "FAIL"
		local tests_failed = `tests_failed' + 1
	}

	local tests_run = `tests_run' + 1

	n di "L`prev'->L`line'  | " %8.1f `pct_bench' "% | " %8.1f `pct_sim' "% | " %5.1f `diff' "% | `status'"
}

**********
// Summary
**********

di _n _n "{hline 80}"
di "VALIDATION SUMMARY"
di "{hline 80}"
di "Tests run:    " %3.0f `tests_run'
di "Tests passed: " %3.0f `tests_passed' " (" %5.1f `tests_passed'/`tests_run'*100 "%)"
di "Tests failed: " %3.0f `tests_failed' " (" %5.1f `tests_failed'/`tests_run'*100 "%)"
di "{hline 80}"

if "$empty_targets" != "" {
	di as error _n "⚠ EMPTY TARGET BENCHMARK(S): $empty_targets"
	di as error "  These scored as FAILs because their target CSV is empty, not because the model is wrong."
	di as error "  Regenerate the targets before reading the score as a model result."
}

if `tests_failed' == 0 {
	di _n "✓ ALL TESTS PASSED"
	di "Simulation accurately reproduces MRDR training data"
}
else {
	di _n "⚠ SOME TESTS FAILED"
	di "Review failed tests above for model issues"
}
