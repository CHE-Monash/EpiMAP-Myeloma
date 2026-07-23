**********
* Monash Myeloma Model - Generate Benchmarks
*
* Purpose: Extract validation benchmarks from the MRDR MI data (in-sample by default, or held-out
*          OOS targets when passed input/output args).
* Usage:   do generate_benchmarks.do [input_MI.dta] [output_dir]   (no args -> in-sample defaults)
**********

* Optional positional args (read into locals, which survive clear all) let this script build either
* the in-sample benchmarks (no args -> defaults below) or the held-out OOS targets, when called from
* analyses/default/prep/test_targets.do:   arg 1 = input MI dataset, arg 2 = output directory.
local bench_in  `"`1'"'
local bench_out `"`2'"'

clear all

* Guard: output paths below are relative to the repo root (see prep/README.md).
* config.example.do is a tracked file that exists only at the repo root.
capture confirm file "config.example.do"
if _rc {
	di as error "Run generate_benchmarks.do from the repository root."
	di as error "Current directory: " c(pwd)
	exit 601
}

if "$repo_path" != "" cd "$repo_path"   // cd to repo root only if config.do set it; a bare cd "" goes to home on Mac/Unix

* Log. Adam runs this and Claude reads the result, per CLAUDE.md - and without one a failure here is
* INVISIBLE: the targets simply keep their previous contents and the next validation scores against
* a stale file. That happened on 23 July, twice, and cost two full validation cycles before the
* timestamps gave it away. Named so it coexists with an interactive session's own log.
cap mkdir "scratch"
cap log close gbm
log using "scratch/generate_benchmarks.log", replace text name(gbm)

capture run "config.do"     // machine-specific paths: $data_path (git-ignored)

local Date "$data_cut"

* Defaults (no args) reproduce the main in-sample benchmarks.
if `"`bench_in'"'  == "" local bench_in  "${data_path}/MRDR Long MI.dta"
if `"`bench_out'"' == "" local bench_out "scratch/benchmarks"

use "`bench_in'", replace

mi extract 1, clear

global timepoints "12 24 36 48 60 72 84 96 120"
global years "1 2 3 4 5 6 7 8 10"

* Helper: write the registry KM survival (%) at 12 and 24 months into columns 7-8 of a TXD/TFI
* benchmark matrix, for the active stset and the given subgroup condition. Used so the validator
* can test TXD/TFI by survival-at-horizon (censoring-robust) rather than the median, which is
* unreliable where the distribution is heavily censored / near-zero (see docs/validation.md).
capture program drop bench_horizons
program define bench_horizons
	args M row ifc
	* Build the KM survivor for this cell. Small / event-free cells (common in the 30% OOS test
	* fold) can make sts generate fail to create the variable -- if so, leave M12/M24 missing
	* (the validator skips missing horizon benchmarks) rather than crashing.
	capture drop _bh_surv
	capture quietly sts generate _bh_surv = s if `ifc'
	capture confirm variable _bh_surv
	if _rc exit
	quietly summarize _bh_surv if (`ifc') & _t <= 12
	if r(N) > 0  matrix `M'[`row', 7] = r(min) * 100
	quietly summarize _bh_surv if (`ifc') & _t <= 24
	if r(N) > 0  matrix `M'[`row', 8] = r(min) * 100
	capture drop _bh_surv
end

* Helper: write the censored % into column 6 of a TXD/TFI benchmark matrix, for the active stset and
* the given subgroup condition. Column 1 counts one row per patient (the record at the stset origin,
* _t0 == 0), so the numerator has to be patients too. Counting _d == 0 counts in-spell ROWS, and a
* patient with several records inside the spell contributes one per record -- that is what pushed the
* TFI_L1_ASCT censored % over 100, ASCT patients carrying the most episode records. A spell fails at
* most once, so _d == 1 counts failed PATIENTS; the censored ones are column 1 less those.
// Helper: blank every estimate in a benchmark row whose N is below the floor, leaving N itself so
// the thinness stays visible. The validator skips missing benchmarks, so a blanked row is silently
// not scored rather than scored badly.
//
// WHY THIS EXISTS. The 30% out-of-sample fold splits some cells to almost nothing, and a KM
// estimate on 17 patients is noise that scores as a model failure. L4 BCR=1 is the worked example:
// N=17 out of sample with 12-month on-treatment 35.3%, against N=54 and 56.4% in sample - the two
// folds disagree by 21 points about the same registry quantity, while the simulation sits near 70%
// in both. That produced six L4 FAILs which were read as a model defect and cost a day (see
// scratch/maintenance/_notes.md). Below the floor the registry cannot say what the truth is, so
// neither should the target.
//
// The floor is deliberately low. It removes cells that are indefensible, not cells that are merely
// small; borderline cells stay in and are meant to be argued about.
capture program drop bench_floor
program define bench_floor
	args M floorN
	local nr = rowsof(`M')
	local nc = colsof(`M')
	forvalues r = 1/`nr' {
		local n = `M'[`r', 1]
		if !missing(`n') & `n' < `floorN' {
			forvalues c = 2/`nc' {
				matrix `M'[`r', `c'] = .
			}
		}
	}
end

capture program drop bench_censored
program define bench_censored
	args M row ifc
	quietly count if (`ifc') & _d == 1
	local fail = r(N)
	local n = `M'[`row', 1]
	if `n' > 0  matrix `M'[`row', 6] = (`n' - `fail') / `n' * 100
end

**********
// Overall Survival
**********

bysort ID_BS: gen first_record = (_n == 1)
// last_record was `_n == _N-1' - the SECOND-to-last row, not the last - and was used to count
// censored patients as "in-spell rows that did not fail". That double-counts patients with several
// records in the spell and silently depends on a record-position assumption nobody documented; it
// is what put the transplant TFI censored figure above 100%. Censored patients are now derived the
// way bench_censored derives them (a spell fails at most once, so censored = N - failures), so the
// variable is gone rather than fixed. Do not reintroduce it.

// OS by BCR_L1 (NoASCT patients only)
stset Date1 if(F_OS != 1 & SCT == 0), id(ID_BS) origin(Event1 == 10) failure(Event1 == 104) scale(30.4375)

matrix OS_L1_NoASCT = J(6, 12, .)
matrix colnames OS_L1_NoASCT = "N" "Median" "Y1" "Y2" "Y3" "Y4" "Y5" "Y6" "Y7" "Y8" "Y10" "Censored"
matrix rownames OS_L1_NoASCT = "CR" "VG" "PR" "MR" "SD" "PD"

forvalues bcr = 1/6 {
	quietly count if BCR_L1 == `bcr' & SCT == 0 & _t0 == 0   // anchor at the L1-start origin: BCR_L1..L9 are LOCF forward from each line's start (not on the diagnosis record)
	local n = r(N)
	matrix OS_L1_NoASCT[`bcr', 1] = `n'
	
	if `n' > 0 {
		quietly stsum if BCR_L1 == `bcr' & SCT == 0
		matrix OS_L1_NoASCT[`bcr', 2] = r(p50)
		
		quietly sts generate surv_temp = s if BCR_L1 == `bcr' & SCT == 0
		quietly summarize _t if BCR_L1 == `bcr' & SCT == 0
		local tmax = r(max)

		local col = 3
		foreach tp of global timepoints {
			quietly summarize surv_temp if BCR_L1 == `bcr' & SCT == 0 & _t <= `tp'
			if r(N) > 0 & `tp' <= `tmax' {
				matrix OS_L1_NoASCT[`bcr', `col'] = r(min)
			}
			local ++col
		}
		
		drop surv_temp
		
		quietly count if BCR_L1 == `bcr' & SCT == 0 & _d == 1
		matrix OS_L1_NoASCT[`bcr', 12] = (`n' - r(N)) / `n' * 100
	}
}

// OS by BCR_ASCT (ASCT patients only)
stset Date1 if(F_OS != 1 & SCT == 1), id(ID_BS) origin(Event1 == 11) failure(Event1 == 104) scale(30.4375)

matrix OS_ASCT = J(4, 12, .)
matrix colnames OS_ASCT = "N" "Median" "Y1" "Y2" "Y3" "Y4" "Y5" "Y6" "Y7" "Y8" "Y10" "Censored"
matrix rownames OS_ASCT = "CR" "VG" "PR" "MR"

forvalues bcr = 1/4 {
	quietly count if BCR_SCT == `bcr' & !missing(BCR_SCT) & first_record == 1
	local n = r(N)
	matrix OS_ASCT[`bcr', 1] = `n'
	
	if `n' > 0 {
		quietly stsum if BCR_SCT == `bcr'
		matrix OS_ASCT[`bcr', 2] = r(p50)
		
		quietly sts generate surv_temp = s if BCR_SCT == `bcr'
		quietly summarize _t if BCR_SCT == `bcr'
		local tmax = r(max)

		local col = 3
		foreach tp of global timepoints {
			quietly summarize surv_temp if BCR_SCT == `bcr' & _t <= `tp'
			if r(N) > 0 & `tp' <= `tmax' {
				matrix OS_ASCT[`bcr', `col'] = r(min)
			}
			local ++col
		}
		
		drop surv_temp
		
		quietly count if BCR_SCT == `bcr' & _d == 1
		matrix OS_ASCT[`bcr', 12] = (`n' - r(N)) / `n' * 100
	}
}

// OS by BCR_L2
stset Date1 if(F_OS != 1), id(ID_BS) origin(Event1 == 20) failure(Event1 == 104) scale(30.4375)

matrix OS_L2 = J(6, 12, .)
matrix colnames OS_L2 = "N" "Median" "Y1" "Y2" "Y3" "Y4" "Y5" "Y6" "Y7" "Y8" "Y10" "Censored"
matrix rownames OS_L2 = "CR" "VG" "PR" "MR" "SD" "PD"

forvalues bcr = 1/6 {
	quietly count if BCR_L2 == `bcr' & !missing(BCR_L2) & _t0 == 0   // anchor at the L2-start origin (BCR_L2 is LOCF forward from L2 start)
	local n = r(N)
	matrix OS_L2[`bcr', 1] = `n'
	
	if `n' > 0 {
		quietly stsum if BCR_L2 == `bcr'
		matrix OS_L2[`bcr', 2] = r(p50)
		
		quietly sts generate surv_temp = s if BCR_L2 == `bcr'
		quietly summarize _t if BCR_L2 == `bcr'
		local tmax = r(max)

		local col = 3
		foreach tp of global timepoints {
			quietly summarize surv_temp if BCR_L2 == `bcr' & _t <= `tp'
			if r(N) > 0 & `tp' <= `tmax' {
				matrix OS_L2[`bcr', `col'] = r(min)
			}
			local ++col
		}
		
		drop surv_temp
		
		quietly count if BCR_L2 == `bcr' & _d == 1
		matrix OS_L2[`bcr', 12] = (`n' - r(N)) / `n' * 100
	}
}

// OS by BCR_L3
stset Date1 if(F_OS != 1), id(ID_BS) origin(Event1 == 30) failure(Event1 == 104) scale(30.4375)

matrix OS_L3 = J(6, 12, .)
matrix colnames OS_L3 = "N" "Median" "Y1" "Y2" "Y3" "Y4" "Y5" "Y6" "Y7" "Y8" "Y10" "Censored"
matrix rownames OS_L3 = "CR" "VG" "PR" "MR" "SD" "PD"

forvalues bcr = 1/6 {
	quietly count if BCR_L3 == `bcr' & !missing(BCR_L3) & _t0 == 0   // anchor at the L3-start origin (BCR_L3 is LOCF forward from L3 start)
	local n = r(N)
	matrix OS_L3[`bcr', 1] = `n'
	
	if `n' > 0 {
		quietly stsum if BCR_L3 == `bcr'
		matrix OS_L3[`bcr', 2] = r(p50)
		
		quietly sts generate surv_temp = s if BCR_L3 == `bcr'
		quietly summarize _t if BCR_L3 == `bcr'
		local tmax = r(max)

		local col = 3
		foreach tp of global timepoints {
			quietly summarize surv_temp if BCR_L3 == `bcr' & _t <= `tp'
			if r(N) > 0 & `tp' <= `tmax' {
				matrix OS_L3[`bcr', `col'] = r(min)
			}
			local ++col
		}
		
		drop surv_temp
		
		quietly count if BCR_L3 == `bcr' & _d == 1
		matrix OS_L3[`bcr', 12] = (`n' - r(N)) / `n' * 100
	}
}

**********
// Lenalidomide-refractory (treatment lines)
**********
// Validates the len-refractory wiring (docs/refractory.md 4.7, 4.4). Target is TRUE len-refractory =
// treatment OR maintenance (LenRefr_Tx_in | LenRefr_Mnt_in), so the checks register BOTH generation
// models together (the maintenance half is ~as large as treatment - 5(6)). Prevalence-by-line scores
// generation; OS-by-status scores the consumption/redistribution. Guarded on both flags so a fold
// built before either existed skips these rather than erroring.
capture confirm variable LenRefr_Tx_in
local have_lenrefr = (_rc == 0)
capture confirm variable LenRefr_Mnt_in
if _rc local have_lenrefr = 0

if `have_lenrefr' {

	// Prevalence of TRUE len-refractory (treatment OR maintenance) AS AT ENTRY to each line (= the
	// sim's LenRefr_L`l', which IS the union - the engine keeps one latched flag that both arms write
	// into), one value per patient per line. Both registry flags are held within a
	// line, so their value on a patient's line-l rows is the entry-to-l state; egen max over those rows
	// recovers it (missing if the line is unreached). L1 is 0 by construction; L2+ carry the accrual.
	matrix LENREFR = J(6, 2, .)
	matrix colnames LENREFR = "N" "PctRefr"
	matrix rownames LENREFR = "L1" "L2" "L3" "L4" "L5" "L6"
	forvalues l = 1/6 {
		capture drop lr_ln lr_pt
		gen byte lr_ln = (LenRefr_Tx_in == 1 | LenRefr_Mnt_in == 1) if Line == `l'
		bysort ID_BS: egen lr_pt = max(lr_ln)
		quietly count if !missing(lr_pt) & first_record == 1
		matrix LENREFR[`l', 1] = r(N)
		if r(N) > 0 {
			quietly count if lr_pt == 1 & first_record == 1
			matrix LENREFR[`l', 2] = r(N) / LENREFR[`l', 1] * 100
		}
		drop lr_ln lr_pt
	}

	// OS from L2 start, split by TRUE len-refractory as at L2 entry (treatment OR maintenance, vs
	// neither). Mirrors the OS-by-BCR_L2 benchmark; the sim scores OS_L2S by LenRefr_L2 (the union)
	// the same way. The direct check on the subgroup OS split (5.6) the whole-population OS cannot see.
	capture drop lr2 lr2_pt
	gen byte lr2 = (LenRefr_Tx_in == 1 | LenRefr_Mnt_in == 1) if Line == 2
	bysort ID_BS: egen lr2_pt = max(lr2)

	stset Date1 if(F_OS != 1), id(ID_BS) origin(Event1 == 20) failure(Event1 == 104) scale(30.4375)

	matrix OS_LENREFR = J(2, 12, .)
	matrix colnames OS_LENREFR = "N" "Median" "Y1" "Y2" "Y3" "Y4" "Y5" "Y6" "Y7" "Y8" "Y10" "Censored"
	matrix rownames OS_LENREFR = "NotRefr" "Refr"

	forvalues r = 0/1 {
		local row = `r' + 1
		quietly count if lr2_pt == `r' & _t0 == 0
		local n = r(N)
		matrix OS_LENREFR[`row', 1] = `n'
		if `n' > 0 {
			quietly stsum if lr2_pt == `r'
			matrix OS_LENREFR[`row', 2] = r(p50)

			quietly sts generate surv_temp = s if lr2_pt == `r'
			quietly summarize _t if lr2_pt == `r'
			local tmax = r(max)

			local col = 3
			foreach tp of global timepoints {
				quietly summarize surv_temp if lr2_pt == `r' & _t <= `tp'
				if r(N) > 0 & `tp' <= `tmax' {
					matrix OS_LENREFR[`row', `col'] = r(min)
				}
				local ++col
			}
			drop surv_temp

			// Censored PATIENTS = N less failures. A spell fails at most once, so counting _d == 1
			// counts patients; the earlier idiom counted in-spell ROWS at a fixed record position
			// (`last_record'), which double-counted patients with several records and was retired
			// from this file. This block came back from a branch that predates that change.
			quietly count if lr2_pt == `r' & _d == 1
			matrix OS_LENREFR[`row', 12] = (`n' - r(N)) / `n' * 100
		}
	}
	drop lr2 lr2_pt
}

**********
// BCR
**********

matrix BCR = J(9, 7, .)
matrix colnames BCR = "N" "CR" "VG" "PR" "MR" "SD" "PD"
matrix rownames BCR = "L1" "ASCT" "L2" "L3" "L4" "L5" "L6" "L7" "L8"

// L1 (count patients, not records). BCR_L1..L9 are LOCF-forward from each line's start (not on the
// diagnosis record), so reduce to one value per patient (egen max over the LOCF'd column) then count
// once per patient at first_record.
capture drop bcr_pt
bysort ID_BS: egen bcr_pt = max(BCR_L1)
quietly count if !missing(bcr_pt) & first_record == 1
matrix BCR[1,1] = r(N)
forvalues bcr = 1/6 {
	quietly count if bcr_pt == `bcr' & first_record == 1
	matrix BCR[1, `bcr'+1] = r(N) / BCR[1,1] * 100
}
drop bcr_pt

// ASCT (only 4 categories)
quietly count if !missing(BCR_SCT) & BCR_SCT != 0 & first_record == 1
matrix BCR[2,1] = r(N)
forvalues bcr = 1/4 {
	quietly count if BCR_SCT == `bcr' & first_record == 1
	matrix BCR[2, `bcr'+1] = r(N) / BCR[2,1] * 100
}

// L2-L8 (same one-value-per-patient reduction as L1, since BCR_L`line' is LOCF-forward)
forvalues line = 2/8 {
	local row = `line' + 1
	capture drop bcr_pt
	bysort ID_BS: egen bcr_pt = max(BCR_L`line')
	quietly count if !missing(bcr_pt) & first_record == 1
	matrix BCR[`row',1] = r(N)

	if r(N) > 0 {
		forvalues bcr = 1/6 {
			quietly count if bcr_pt == `bcr' & first_record == 1
			matrix BCR[`row', `bcr'+1] = r(N) / BCR[`row',1] * 100
		}
	}
	drop bcr_pt
}

**********
// TXD
**********

// TXD_L1 by BCR_L1 - NoASCT patients
stset Date1, id(ID_BS) origin(Event1 == 10) failure(Event1 == 11) scale(30.4375)

matrix TXD_L1_NoASCT = J(6, 8, .)
matrix colnames TXD_L1_NoASCT = "N" "Mean" "Median" "P25" "P75" "Censored" "M12" "M24"
matrix rownames TXD_L1_NoASCT = "CR" "VG" "PR" "MR" "SD" "PD"

forvalues bcr = 1/6 {
	quietly count if BCR_L1 == `bcr' & SCT == 0 & _t0 == 0
	matrix TXD_L1_NoASCT[`bcr', 1] = r(N)

	if r(N) > 0 {
		quietly stci if BCR_L1 == `bcr' & SCT == 0, rmean
		matrix TXD_L1_NoASCT[`bcr', 2] = r(rmean)

		quietly stsum if BCR_L1 == `bcr' & SCT == 0
		matrix TXD_L1_NoASCT[`bcr', 3] = r(p50)
		matrix TXD_L1_NoASCT[`bcr', 4] = r(p25)
		matrix TXD_L1_NoASCT[`bcr', 5] = r(p75)

		bench_censored TXD_L1_NoASCT `bcr' "BCR_L1 == `bcr' & SCT == 0"
		bench_horizons TXD_L1_NoASCT `bcr' "BCR_L1 == `bcr' & SCT == 0"
	}
}

// TXD_L1 by BCR_SCT - ASCT patients (ASCT response has 4 categories: BCR 5/6 are folded into 4 in
// multiple_imputation.do. Stratify by BCR_SCT (1-4) to match OS_ASCT / TFI_L1_ASCT and the validator.)
matrix TXD_L1_ASCT = J(4, 8, .)
matrix colnames TXD_L1_ASCT = "N" "Mean" "Median" "P25" "P75" "Censored" "M12" "M24"
matrix rownames TXD_L1_ASCT = "CR" "VG" "PR" "MR"

forvalues bcr = 1/4 {
	quietly count if BCR_SCT == `bcr' & !missing(BCR_SCT) & _t0 == 0
	matrix TXD_L1_ASCT[`bcr', 1] = r(N)

	if r(N) > 0 {
		quietly stci if BCR_SCT == `bcr', rmean
		matrix TXD_L1_ASCT[`bcr', 2] = r(rmean)

		quietly stsum if BCR_SCT == `bcr'
		matrix TXD_L1_ASCT[`bcr', 3] = r(p50)
		matrix TXD_L1_ASCT[`bcr', 4] = r(p25)
		matrix TXD_L1_ASCT[`bcr', 5] = r(p75)

		bench_censored TXD_L1_ASCT `bcr' "BCR_SCT == `bcr'"
		bench_horizons TXD_L1_ASCT `bcr' "BCR_SCT == `bcr'"
	}
}

// TXD_L2 by BCR_L2
stset Date1, id(ID_BS) origin(Event1 == 20) failure(Event1 == 21) scale(30.4375)

matrix TXD_L2 = J(6, 8, .)
matrix colnames TXD_L2 = "N" "Mean" "Median" "P25" "P75" "Censored" "M12" "M24"
matrix rownames TXD_L2 = "CR" "VG" "PR" "MR" "SD" "PD"

forvalues bcr = 1/6 {
	quietly count if BCR_L2 == `bcr' & !missing(BCR_L2) & _t0 == 0
	matrix TXD_L2[`bcr', 1] = r(N)

	if r(N) > 0 {
		quietly stci if BCR_L2 == `bcr', rmean
		matrix TXD_L2[`bcr', 2] = r(rmean)

		quietly stsum if BCR_L2 == `bcr'
		matrix TXD_L2[`bcr', 3] = r(p50)
		matrix TXD_L2[`bcr', 4] = r(p25)
		matrix TXD_L2[`bcr', 5] = r(p75)

		bench_censored TXD_L2 `bcr' "BCR_L2 == `bcr'"
		bench_horizons TXD_L2 `bcr' "BCR_L2 == `bcr'"
	}
}

// TXD_L3 by BCR_L3
stset Date1, id(ID_BS) origin(Event1 == 30) failure(Event1 == 31) scale(30.4375)

matrix TXD_L3 = J(6, 8, .)
matrix colnames TXD_L3 = "N" "Mean" "Median" "P25" "P75" "Censored" "M12" "M24"
matrix rownames TXD_L3 = "CR" "VG" "PR" "MR" "SD" "PD"

forvalues bcr = 1/6 {
	quietly count if BCR_L3 == `bcr' & !missing(BCR_L3) & _t0 == 0
	matrix TXD_L3[`bcr', 1] = r(N)

	if r(N) > 0 {
		quietly stci if BCR_L3 == `bcr', rmean
		matrix TXD_L3[`bcr', 2] = r(rmean)

		quietly stsum if BCR_L3 == `bcr'
		matrix TXD_L3[`bcr', 3] = r(p50)
		matrix TXD_L3[`bcr', 4] = r(p25)
		matrix TXD_L3[`bcr', 5] = r(p75)

		bench_censored TXD_L3 `bcr' "BCR_L3 == `bcr'"
		bench_horizons TXD_L3 `bcr' "BCR_L3 == `bcr'"
	}
}

// TXD_L4 by BCR_L4
stset Date1, id(ID_BS) origin(Event1 == 40) failure(Event1 == 41) scale(30.4375)

matrix TXD_L4 = J(6, 8, .)
matrix colnames TXD_L4 = "N" "Mean" "Median" "P25" "P75" "Censored" "M12" "M24"
matrix rownames TXD_L4 = "CR" "VG" "PR" "MR" "SD" "PD"

forvalues bcr = 1/6 {
	quietly count if BCR_L4 == `bcr' & !missing(BCR_L4) & _t0 == 0
	matrix TXD_L4[`bcr', 1] = r(N)

	if r(N) > 0 {
		quietly stci if BCR_L4 == `bcr', rmean
		matrix TXD_L4[`bcr', 2] = r(rmean)

		quietly stsum if BCR_L4 == `bcr'
		matrix TXD_L4[`bcr', 3] = r(p50)
		matrix TXD_L4[`bcr', 4] = r(p25)
		matrix TXD_L4[`bcr', 5] = r(p75)

		bench_censored TXD_L4 `bcr' "BCR_L4 == `bcr'"
		bench_horizons TXD_L4 `bcr' "BCR_L4 == `bcr'"
	}
}
**********
// TFI
**********

// TFI_L1 by BCR_L1 - NoASCT patients
stset Date1, id(ID_BS) origin(Event1 == 11) failure(Event1 == 20) scale(30.4375)

matrix TFI_L1_NoASCT = J(6, 8, .)
matrix colnames TFI_L1_NoASCT = "N" "Mean" "Median" "P25" "P75" "Censored" "M12" "M24"
matrix rownames TFI_L1_NoASCT = "CR" "VG" "PR" "MR" "SD" "PD"

forvalues bcr = 1/6 {
	quietly count if BCR_L1 == `bcr' & SCT == 0 & _t0 == 0
	matrix TFI_L1_NoASCT[`bcr', 1] = r(N)

	if r(N) > 0 {
		quietly stci if BCR_L1 == `bcr' & SCT == 0, rmean
		matrix TFI_L1_NoASCT[`bcr', 2] = r(rmean)

		quietly stsum if BCR_L1 == `bcr' & SCT == 0
		matrix TFI_L1_NoASCT[`bcr', 3] = r(p50)
		matrix TFI_L1_NoASCT[`bcr', 4] = r(p25)
		matrix TFI_L1_NoASCT[`bcr', 5] = r(p75)

		bench_censored TFI_L1_NoASCT `bcr' "BCR_L1 == `bcr' & SCT == 0"
		bench_horizons TFI_L1_NoASCT `bcr' "BCR_L1 == `bcr' & SCT == 0"
	}
}

// TFI_L1 by BCR_ASCT - ASCT patients
matrix TFI_L1_ASCT = J(4, 8, .)
matrix colnames TFI_L1_ASCT = "N" "Mean" "Median" "P25" "P75" "Censored" "M12" "M24"
matrix rownames TFI_L1_ASCT = "CR" "VG" "PR" "MR"

forvalues bcr = 1/4 {
	quietly count if BCR_SCT == `bcr' & !missing(BCR_SCT) & _t0 == 0
	matrix TFI_L1_ASCT[`bcr', 1] = r(N)

	if r(N) > 0 {
		quietly stci if BCR_SCT == `bcr', rmean
		matrix TFI_L1_ASCT[`bcr', 2] = r(rmean)

		quietly stsum if BCR_SCT == `bcr'
		matrix TFI_L1_ASCT[`bcr', 3] = r(p50)
		matrix TFI_L1_ASCT[`bcr', 4] = r(p25)
		matrix TFI_L1_ASCT[`bcr', 5] = r(p75)

		bench_censored TFI_L1_ASCT `bcr' "BCR_SCT == `bcr'"
		bench_horizons TFI_L1_ASCT `bcr' "BCR_SCT == `bcr'"
	}
}

// TFI_L2 by BCR_L2
stset Date1, id(ID_BS) origin(Event1 == 21) failure(Event1 == 30) scale(30.4375)

matrix TFI_L2 = J(6, 8, .)
matrix colnames TFI_L2 = "N" "Mean" "Median" "P25" "P75" "Censored" "M12" "M24"
matrix rownames TFI_L2 = "CR" "VG" "PR" "MR" "SD" "PD"

forvalues bcr = 1/6 {
	quietly count if BCR_L2 == `bcr' & !missing(BCR_L2) & _t0 == 0
	matrix TFI_L2[`bcr', 1] = r(N)

	if r(N) > 0 {
		quietly stci if BCR_L2 == `bcr', rmean
		matrix TFI_L2[`bcr', 2] = r(rmean)

		quietly stsum if BCR_L2 == `bcr'
		matrix TFI_L2[`bcr', 3] = r(p50)
		matrix TFI_L2[`bcr', 4] = r(p25)
		matrix TFI_L2[`bcr', 5] = r(p75)

		bench_censored TFI_L2 `bcr' "BCR_L2 == `bcr'"
		bench_horizons TFI_L2 `bcr' "BCR_L2 == `bcr'"
	}
}

// TFI_L3 by BCR_L3
stset Date1, id(ID_BS) origin(Event1 == 31) failure(Event1 == 40) scale(30.4375)

matrix TFI_L3 = J(6, 6, .)
matrix colnames TFI_L3 = "N" "Mean" "Median" "P25" "P75" "Censored"
matrix rownames TFI_L3 = "CR" "VG" "PR" "MR" "SD" "PD"

forvalues bcr = 1/6 {
	quietly count if BCR_L3 == `bcr' & !missing(BCR_L3) & _t0 == 0
	matrix TFI_L3[`bcr', 1] = r(N)
	
	if r(N) > 0 {
		quietly stci if BCR_L3 == `bcr', rmean
		matrix TFI_L3[`bcr', 2] = r(rmean)
		
		quietly stsum if BCR_L3 == `bcr'
		matrix TFI_L3[`bcr', 3] = r(p50)
		matrix TFI_L3[`bcr', 4] = r(p25)
		matrix TFI_L3[`bcr', 5] = r(p75)
		
		bench_censored TFI_L3 `bcr' "BCR_L3 == `bcr'"
	}
}

**********
// Pathways  (censoring-aware: competing-risks cumulative incidence)
**********
//
// The registry has incomplete follow-up, so a crude "ever reached / total" count understates the
// CONDITIONAL per-transition progression: P(reach line L | reached line L-1). For each line we
// estimate the Aalen-Johansen cumulative incidence of reaching it with *death before reaching it* as
// a competing event, but with origin = the date the PREVIOUS line was reached and the risk set
// restricted to patients who reached that previous line. This isolates each L-1 -> L step instead of
// the cumulative-from-L1 reach, so errors don't compound line by line and the follow-up truncation is
// far milder. The simulation posts the matching ratio (reached L / reached L-1) per resample.
//
// ASCT is left as a crude proportion: transplant is decided during L1, so it is observed for
// essentially every patient and is not subject to the cumulative follow-up censoring above.

capture mata: mata drop aj_cif()
mata:
// Aalen-Johansen cumulative incidence for cause 1 (reach) at end of follow-up, with cause 2 (death)
// competing. time = sojourn from origin (any monotone scale); status: 1=reach, 2=death, 0=censored.
real scalar aj_cif(real colvector time, real colvector status)
{
	real colvector keep, t, s, ord
	real scalar n, atrisk, S, cif, i, j, ti, d1, d2, dc
	keep = (time :< .) :& (status :< .)
	t = select(time, keep); s = select(status, keep)
	n = rows(t)
	if (n == 0) return(.)
	ord = order(t, 1); t = t[ord]; s = s[ord]
	atrisk = n; S = 1; cif = 0; i = 1
	while (i <= n) {
		ti = t[i]; d1 = 0; d2 = 0; dc = 0; j = i
		while (j <= n) {                          // group ties at time ti (j-index guarded first;
			if (t[j] != ti) break                 // Mata & is not short-circuit, so don't test t[j] in the while)
			if      (s[j] == 1) d1++
			else if (s[j] == 2) d2++
			else                dc++
			j++
		}
		if (atrisk > 0) {
			cif = cif + S * (d1 / atrisk)        // increment of CIF for cause 1
			S   = S   * (1 - (d1 + d2) / atrisk)  // overall "still in earlier line & alive" survival
		}
		atrisk = atrisk - (d1 + d2 + dc)
		i = j
	}
	return(cif)
}
end

matrix Pathways = J(1, 9, .)
matrix colnames Pathways = "N" "ASCT" "L2" "L3" "L4" "L5" "L6" "L7" "L8"
matrix rownames Pathways = "Percent"

preserve

// Total patients (denominator)
quietly count if first_record == 1
matrix Pathways[1,1] = r(N)

// ASCT -- transplant is decided AT THE END OF L1, so the rate is among patients who reach L1 end
// (have an L1-end event, Event 11), NOT all diagnosed patients. The model's ASCT logit is fit on
// the same conditional population (Event1 == 11), and validate_simulation.do uses the matching
// denominator (patients alive at L1 end).
bysort ID_BS: egen temp_asct = max(SCT == 1)
bysort ID_BS: egen temp_l1e  = max(Event0 == 11 | Event1 == 11)
quietly count if temp_l1e == 1 & first_record == 1
local n_l1e = r(N)
quietly count if temp_asct == 1 & first_record == 1
matrix Pathways[1,2] = r(N) / `n_l1e' * 100
drop temp_asct temp_l1e

// Per-patient competing-risks fields (constant within ID_BS)
tempvar rowv
//  origin = L1 start date
gen double `rowv' = Date0 if Event0 == 10
bysort ID_BS: egen double cr_origin = min(`rowv')
drop `rowv'
//  death date (recorded as Event0==104, or as the next event Event1==104 on the terminal row)
gen double `rowv' = .
replace `rowv' = Date0 if Event0 == 104
replace `rowv' = Date1 if Event1 == 104 & missing(`rowv')
bysort ID_BS: egen double cr_death = min(`rowv')
drop `rowv'
//  last follow-up date (latest date observed for the patient)
gen double `rowv' = max(Date0, Date1)
bysort ID_BS: egen double cr_lastfu = max(`rowv')
drop `rowv'
//  first date each line is reached (Event0 == line*10)
forvalues line = 2/8 {
	gen double `rowv' = Date0 if Event0 == `line' * 10
	bysort ID_BS: egen double cr_reach`line' = min(`rowv')
	drop `rowv'
}

// Collapse to one record per patient
keep if first_record == 1

// Lines 2-8: CONDITIONAL per-transition progression -- P(reach L | reached L-1), Aalen-Johansen
// cumulative incidence with death competing. Origin = the date the PREVIOUS line was reached, and the
// risk set is patients who reached that previous line (prevdate non-missing). L2's previous line is L1
// (origin cr_origin), so L2 is unchanged; L3+ are now conditioned on the prior line rather than
// cumulative from L1 -- isolates each transition and removes the cumulative-cascade artefact.
forvalues line = 2/8 {
	local prev = `line' - 1
	if (`prev' == 1) local prevdate cr_origin
	else             local prevdate cr_reach`prev'
	tempvar status time
	gen byte   `status' = cond(!missing(cr_reach`line'), 1, cond(!missing(cr_death), 2, 0)) if !missing(`prevdate')
	gen double `time'   = (cond(`status'==1, cr_reach`line', cond(`status'==2, cr_death, cr_lastfu)) - `prevdate') if !missing(`prevdate')
	replace    `time'   = 0 if `time' < 0     // guard rare same-day / ordering artefacts
	mata: st_numscalar("cr_cif", aj_cif(st_data(., "`time'"), st_data(., "`status'")))
	matrix Pathways[1, `line'+1] = cr_cif * 100
	drop `status' `time'
}

restore

// Whole-population OS from diagnosis (all patients; the aggregate check the per-line/BCR OS lacks).
// Same estimators as the stratified OS above, but clocked from diagnosis (origin = Event1==3) with
// no line/BCR filter -> the overall survival curve at 3/5/10 yr for the whole cohort.
stset Date1 if(F_OS != 1), id(ID_BS) origin(Event1 == 3) failure(Event1 == 104) scale(30.4375)
matrix OS_All = J(1, 12, .)
matrix colnames OS_All = "N" "Median" "Y1" "Y2" "Y3" "Y4" "Y5" "Y6" "Y7" "Y8" "Y10" "Censored"
matrix rownames OS_All = "All"
quietly count if first_record == 1
local n = r(N)
matrix OS_All[1, 1] = `n'
if `n' > 0 {
	quietly stsum
	matrix OS_All[1, 2] = r(p50)
	quietly sts generate surv_temp = s
	quietly summarize _t
	local tmax = r(max)
	local col = 3
	foreach tp of global timepoints {
		quietly summarize surv_temp if _t <= `tp'
		if r(N) > 0 & `tp' <= `tmax' matrix OS_All[1, `col'] = r(min)
		local ++col
	}
	drop surv_temp
	quietly count if _d == 1
	matrix OS_All[1, 12] = (`n' - r(N)) / `n' * 100
}

// Whole-population OS monthly KM curve (survivor + Greenwood SE) at 0..120 months from diagnosis.
// These are the inputs for the 2024-style validation figure (observed KM 95% CI vs simulated 95% CI
// + monthly p-value). Exported as a target CSV so bootstrap_validation.do -- which runs on the HPC
// without the drive -- can read the observed curve. Uses the same whole-pop stset as OS_All above.
preserve
	keep if _st == 1
	bysort ID_BS (_t): keep if _n == _N          // one row per patient: _t = OS exit (mo), _d = died
	keep _t _d
	gen double _tot = 1
	collapse (sum) d = _d (sum) tot = _tot, by(_t)   // per distinct exit time: deaths, # exiting
	sort _t
	egen double _Ntot = total(tot)
	gen double _cum = sum(tot)
	gen double n = _Ntot - _cum + tot            // number at risk just before _t (no delayed entry)
	gen double _fac = 1 - d/n
	gen byte   _zero = sum(_fac <= 0) > 0         // running flag: has the survivor hit 0?
	gen double _s = exp(sum(ln(cond(_fac > 0, _fac, 1))))   // product-limit survivor
	replace _s = 0 if _zero
	gen double _gt = cond(n > d, d/(n*(n - d)), 0)          // Greenwood increment
	gen double _gcum = sum(_gt)
	tempfile _curvecsv
	tempname _cvf
	postfile `_cvf' int month double s_obs double se_obs using "`_curvecsv'", replace
	forvalues m = 0/120 {
		qui summarize _s    if _t <= `m', meanonly
		local sm = cond(r(N) > 0, r(min), 1)     // S(m): survivor is non-increasing -> min over t<=m
		qui summarize _gcum if _t <= `m', meanonly
		local gm = cond(r(N) > 0, r(max), 0)     // Greenwood cum: non-decreasing -> max over t<=m
		post `_cvf' (`m') (`sm') (`sm' * sqrt(`gm'))
	}
	postclose `_cvf'
	use "`_curvecsv'", clear
	export delimited using "`bench_out'/os_wholepop_curve.csv", replace
restore

// Whole-population OS stratified by baseline comorbidity burden (CKD+CRD+PLM+DBT: 0 / 1 / 2+).
// Mirrors OS_All. Adding the four comorbidities to OS leaves the aggregate ~unchanged (it only
// redistributes hazard), so THIS stratified target is what tests whether the model's comorbidity
// differentiation is calibrated out-of-sample.
local have_cm = 0
capture confirm variable CM_CKD
if _rc == 0 {
	local have_cm = 1
	capture drop _cmn _cmnp _cmg
	quietly gen byte _cmn = CM_CKD + CM_CRD + CM_PLM + CM_DBT
	// Patient-level burden. The flags are carried FORWARD from diagnosis, so pre-diagnosis rows
	// (e.g. the DOB record = first_record) are missing; a row-wise group would then dump every
	// patient into "2+". CM_* are constant per patient from diagnosis, so max over the patient
	// recovers that value; patients with no non-missing flag stay missing (excluded).
	bysort ID_BS: egen byte _cmnp = max(_cmn)
	quietly gen byte _cmg = cond(_cmnp == 0, 0, cond(_cmnp == 1, 1, 2)) if !missing(_cmnp)
	matrix OS_CM = J(3, 12, .)
	matrix colnames OS_CM = "N" "Median" "Y1" "Y2" "Y3" "Y4" "Y5" "Y6" "Y7" "Y8" "Y10" "Censored"
	matrix rownames OS_CM = "CM0" "CM1" "CM2plus"
	forvalues g = 0/2 {
		local r = `g' + 1
		stset Date1 if(F_OS != 1 & _cmg == `g'), id(ID_BS) origin(Event1 == 3) failure(Event1 == 104) scale(30.4375)
		quietly count if first_record == 1 & _cmg == `g'
		local n = r(N)
		matrix OS_CM[`r', 1] = `n'
		if `n' > 0 {
			quietly stsum
			matrix OS_CM[`r', 2] = r(p50)
			quietly sts generate surv_temp = s
			quietly summarize _t
			local tmax = r(max)
			local col = 3
			foreach tp of global timepoints {
				quietly summarize surv_temp if _t <= `tp'
				if r(N) > 0 & `tp' <= `tmax' matrix OS_CM[`r', `col'] = r(min)
				local ++col
			}
			drop surv_temp
			quietly count if _d == 1 & _cmg == `g'
			matrix OS_CM[`r', 12] = (`n' - r(N)) / `n' * 100
		}
	}
	drop _cmn _cmnp _cmg
}



// MND benchmark - L1 maintenance DURATION by regimen, censoring-aware.
//
// REPLACED a share-of-the-gap metric that was not comparable between the two sides: both quantities
// need an observed L2, so the registry side was computed on relapsers only while the simulated side
// has a closed gap for everybody. The model scored 0/7 against it. Duration needs no closed gap, uses
// the whole maintenance population, and is what the cost engine bills. What is lost is the regimen x
// gap interaction, which cannot be recovered on the registry side - see docs/refractory.md 5(8).
//
// Regimen groups match analyses/default/outcomes/mnr_full.do ($MNR_L1 "1 5"): lenalidomide,
// thalidomide, everything else pooled to 0. If that list changes this must change with it.
//
// DEATH IS A FAILURE HERE AND CENSORING IN THE FIT. That asymmetry is deliberate - do NOT
// "correct" prep/risk_equations.do to match, it would double-count death.
//
// The fit is cause-specific: censoring death estimates how long maintenance runs for MAINTENANCE
// reasons, and the engine then imposes death itself (sim_mort curtails TFI_L1, process_data clips
// MND_L1 to it). Draw the cause-specific duration, kill patients separately. That is the right
// division of labour and both halves are correct.
//
// But it means the two sides measure different things. The engine EXPORTS delivered exposure,
// while a death-censored KM estimates the duration that would be seen if nobody died. This target
// benchmarked the first against the second, so a gap was guaranteed and it was large: with death
// censored the registry RMST(158) is 51.0 months against a simulated 35.7. Adding death as a
// failure makes the target delivered exposure, which is what the engine produces and what cost
// bills.
//
// The model itself was never wrong. The DRAWN duration (MND_raw, exported by process_data.do for
// exactly this check) reproduces the death-censored registry at every horizon - E[min(MND_raw,c)]
// of 35.6 / 42.7 / 48.5 against RMST of 35.23 / 42.91 / 51.01 at c = 84 / 120 / 158. Curtailment
// then costs 12.8 months, and that is real exposure lost, not an error. See ../scratch/maintenance
// /_notes.md and docs/refractory.md.
//
// One caveat to state rather than fix: censoring death in the fit assumes death is uninformative
// about maintenance duration, and it is not - patients who die are sicker and would have come off
// maintenance sooner anyway. So the fitted duration is biased slightly long.
preserve
	// stset the maintenance episode as risk_equations.do fits L1_MND: origin at the maintenance
	// start, failure at the recorded end or the next line. Patients still on maintenance at the cut
	// never fail and are censored.
	// THALIDOMIDE IS CENSORED AT 18 MONTHS, matching the fit - both sides must measure the same
	// quantity or the miss is purely definitional. TRAP: exit() is in the TIME VARIABLE's scale and
	// `origin' is unavailable inside it when origin() is an event condition, so it references the
	// maintenance-start date in days.
	capture drop _mxMNDo
	capture drop MND_origin
	quietly gen double _mxMNDo = Date1 if Event1 == 110
	quietly egen double MND_origin = min(_mxMNDo), by(ID_BS)

	// A stset covering everyone, purely so the _t guard below has something to confirm. Each group
	// re-stsets itself inside the loop.
	capture stset Date1 if(MNT == 1), ///
		id(ID_BS) failure(Event1 == 20 111 104) origin(Event1 == 110) scale(30.4375)

	matrix MND_L1 = J(3, 5, .)
	matrix colnames MND_L1 = "N" "Failures" "Median" "P25" "P75"
	matrix rownames MND_L1 = "len" "thal" "other"

	// GUARD: if the maintenance events are absent from this cut (an older extraction), stset fails
	// or yields no subjects. Leave the matrix missing rather than exporting zeros - the validator
	// skips missing benchmarks, and an empty target that scores as 0 is worse than no target.
	capture confirm variable _t
	if _rc == 0 {
		local mrow = 0
		foreach g in 1 5 0 {
			local ++mrow
			local gcond "MNR_L1 == `g'"
			if `g' == 0 local gcond "!inlist(MNR_L1, 1, 5)"

			// EACH GROUP SETS ITS OWN CLOCK, explicitly. Thalidomide is scored 18-month-censored
			// and the others are not, so a single stset outside the loop plus a special case for
			// thalidomide would leave group 0 - which comes AFTER thalidomide in this loop -
			// running on the thalidomide-only risk set and returning nothing.
			if `g' == 5 {
				capture stset Date1 if(MNT == 1 & MNR_L1 == 5), ///
					id(ID_BS) failure(Event1 == 20 111 104) origin(Event1 == 110) scale(30.4375) ///
					exit(time MND_origin + `=18 * 30.4375')
			}
			else {
				capture stset Date1 if(MNT == 1), ///
					id(ID_BS) failure(Event1 == 20 111 104) origin(Event1 == 110) scale(30.4375)
			}

			quietly count if (`gcond') & _t0 == 0
			matrix MND_L1[`mrow', 1] = r(N)
			quietly count if (`gcond') & _d == 1
			matrix MND_L1[`mrow', 2] = r(N)

			// stci returns the KM median and its quartiles; it leaves them missing when the
			// survivor never reaches the relevant probability, which is the honest answer for a
			// heavily-censored cell and is skipped downstream rather than guessed at.
			capture quietly stci if `gcond'
			if _rc == 0 & r(N) > 0 {
				matrix MND_L1[`mrow', 3] = r(p50)
				// Two traps here, both of which shipped once. stci's p(#) is the #th percentile
				// of SURVIVAL TIME, so p(25) is the SHORTER duration and belongs in the P25
				// column. And it returns the result in r(p25) / r(p75), NOT r(p50) - reading
				// r(p50) after a p(#) call silently exports an empty column.
				capture quietly stci if `gcond', p(25)
				if _rc == 0 matrix MND_L1[`mrow', 4] = r(p25)
				capture quietly stci if `gcond', p(75)
				if _rc == 0 matrix MND_L1[`mrow', 5] = r(p75)
			}
		}
	}
restore

**********
// 5b. MINIMUM-N FLOOR
**********

// Applied to every benchmark AFTER it is built and BEFORE it is exported, so the floor is one
// decision in one place rather than a condition repeated in twenty loops. N stays; the estimates
// go. Pathways is excluded: its rows are whole-cohort reach rates, not per-response cells, so they
// are never thin, and OS_All / the KM curve are single whole-population rows for the same reason.
local bench_min = 20

foreach M in OS_L1_NoASCT OS_ASCT OS_L2 OS_L3 OS_CM ///
             TXD_L1_NoASCT TXD_L1_ASCT TXD_L2 TXD_L3 TXD_L4 ///
             TFI_L1_NoASCT TFI_L1_ASCT TFI_L2 TFI_L3 ///
             BCR MND_L1 {
	capture confirm matrix `M'
	if _rc == 0 bench_floor `M' `bench_min'
}

**********
// 6. EXPORT TO CSV
**********

preserve

// Lenalidomide-refractory benchmarks: prevalence by line, and OS from L2 by refractory status
if `have_lenrefr' {
	clear
	svmat LENREFR, names(col)
	gen Line = _n
	order Line
	export delimited using "`bench_out'/lenrefr.csv", replace

	clear
	svmat OS_LENREFR, names(col)
	gen Refr = _n - 1
	order Refr
	export delimited using "`bench_out'/os_lenrefr.csv", replace
}

// OS benchmarks
clear
svmat OS_L1_NoASCT, names(col)
gen BCR = _n
order BCR
export delimited using "`bench_out'/os_l1_noasct.csv", replace

clear
svmat OS_ASCT, names(col)
gen BCR = _n
order BCR
export delimited using "`bench_out'/os_asct.csv", replace

clear
svmat OS_L2, names(col)
gen BCR = _n
order BCR
export delimited using "`bench_out'/os_l2.csv", replace

clear
svmat OS_L3, names(col)
gen BCR = _n
order BCR
export delimited using "`bench_out'/os_l3.csv", replace

clear
svmat OS_All, names(col)
gen BCR = _n
order BCR
export delimited using "`bench_out'/os_wholepop.csv", replace

// Whole-population OS by comorbidity burden (0 / 1 / 2+); rows CM0/CM1/CM2plus
if `have_cm' {
	clear
	svmat OS_CM, names(col)
	gen CM = _n - 1
	order CM
	export delimited using "`bench_out'/os_wholepop_cm.csv", replace
}


// TFI benchmarks
clear
svmat TFI_L1_NoASCT, names(col)
gen BCR = _n
order BCR
export delimited using "`bench_out'/tfi_l1_noasct.csv", replace

clear
svmat TFI_L1_ASCT, names(col)
gen BCR = _n
order BCR
export delimited using "`bench_out'/tfi_l1_asct.csv", replace

clear
svmat TFI_L2, names(col)
gen BCR = _n
order BCR
export delimited using "`bench_out'/tfi_l2.csv", replace

clear
svmat TFI_L3, names(col)
gen BCR = _n
order BCR
export delimited using "`bench_out'/tfi_l3.csv", replace

// TXD benchmarks
clear
svmat TXD_L1_NoASCT, names(col)
gen BCR = _n
order BCR
export delimited using "`bench_out'/txd_l1_noasct.csv", replace

clear
svmat TXD_L1_ASCT, names(col)
gen BCR = _n
order BCR
export delimited using "`bench_out'/txd_l1_asct.csv", replace

clear
svmat TXD_L2, names(col)
gen BCR = _n
order BCR
export delimited using "`bench_out'/txd_l2.csv", replace

clear
svmat TXD_L3, names(col)
gen BCR = _n
order BCR
export delimited using "`bench_out'/txd_l3.csv", replace

clear
svmat TXD_L4, names(col)
gen BCR = _n
order BCR
export delimited using "`bench_out'/txd_l4.csv", replace

clear
svmat MND_L1, names(col)
gen MNR = .
local mrow = 0
foreach g in 1 5 0 {
	local ++mrow
	quietly replace MNR = `g' in `mrow'
}
order MNR
export delimited using "`bench_out'/mnd_l1.csv", replace

// BCR distributions
clear
svmat BCR, names(col)
gen Line = _n
order Line
export delimited using "`bench_out'/bcr.csv", replace

// Pathways
clear
svmat Pathways, names(col)
export delimited using "`bench_out'/pathways.csv", replace

restore

cap log close gbm
