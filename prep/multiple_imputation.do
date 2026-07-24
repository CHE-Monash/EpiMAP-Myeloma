**********
* Monash Myeloma Model - Multiple Imputation
*
* Purpose: Multiple imputation of MRDR Long. Performance-optimised variant intended to produce
*          IDENTICAL output to the original for a given seed.
* Notes:   Optimisation rationale (carryforward/broadcast helpers, and a deliberately-deferred
*          structural win) documented below.
**********

* -------------------------------------------------------------------------------------------------
* PERFORMANCE. The post-imputation carryforward was the bottleneck (~27 `mi xeq 0/$imp: bysort'
* calls, each re-sorting the full long dataset M+1 times). The _cf / _bcast_idbs helpers below work
* directly on the wide per-imputation columns instead, sorting ONCE: `mi xeq' does not preserve sort
* order across its passes, and replacing a non-sort-key never clears the sort. Output is identical
* for a given seed.
*
* DEFERRED: imputing baseline covariates on a one-row-per-patient `keep if Event0==3' extract and
* mi merge-ing back would shrink the mi object ~10x, but the current carryforward is forward-LOCF
* rather than a pure broadcast, so it needs its own equivalence check first.
* -------------------------------------------------------------------------------------------------

clear
clear mata

if "$repo_path" != "" cd "$repo_path"   // cd to repo root only if config.do set it; a bare cd "" goes to home on Mac/Unix
capture run "config.do"     // machine-specific paths: $data_path (git-ignored)

local Data "$data_cut"

**********
// Settings
global imp `1'
global boot `2'
global min_bs `3'
global max_bs `4'
global sample `5'   // "" = full cohort (main model); "train"/"test" = OOS fold (analyses/default/)

* OOS routing: when $sample is set, restrict to that fold (split crosswalk written by
* analyses/default/prep/split.do) and write outputs under ${data_path}/oos/. Empty = main model.
if "$sample" == "" {
	global mi_outdir ""
	global mi_outtag ""
}
else {
	global mi_outdir "oos/"
	global mi_outtag "_$sample"
	capture mkdir "${data_path}/oos"
	capture mkdir "${data_path}/oos/bootstrap"
}

**********
// Helpers: carry values across a patient's rows by operating directly on the wide per-imputation
// columns (`_m_var') plus the master (`var'), sorting once instead of looping `mi xeq 0/$imp'.
//
// _cf           temporal LOCF within ID_BS, across m=0..M. Caller sorts first. Optional extra `if'.
//               Third arg "nomaster" skips the m=0 fill: for a REGISTERED IMPUTED variable, filling
//               the master at imputed rows makes `mi update' treat them as observed and collapse the
//               imputations to that one value. Pass it for imputed vars, omit for derived ones.
// _bcast_idbs   _cf forward then backward - one value broadcast to all of a patient's rows.
cap program drop _cf
program define _cf
	args v extra nomaster
	if "`extra'" != "" local extra "& `extra'"
	if "`nomaster'" == "" qui by ID_BS: replace `v' = `v'[_n-1] if `v' == . `extra'
	forvalues m = 1/$imp {
		qui by ID_BS: replace _`m'_`v' = _`m'_`v'[_n-1] if _`m'_`v' == . `extra'
	}
end

// _bcast_idbs: broadcast a single per-patient value onto ALL of that patient's rows (fill forward
//      then backward in Date0). Identical to the old `bysort ID_BS (VAR)` value-sort broadcast when
//      there is one non-missing value per patient (true for the diagnosis-/line-/SCT-anchored vars
//      it is used on). Two sorts total for the whole variable list, regardless of M.
cap program drop _bcast_idbs
program define _bcast_idbs
	local vlist "`0'"
	sort ID_BS Date0
	foreach v of local vlist {
		_cf `v'
	}
	// backward pass: sort on a negated date (plain sort -> `by ID_BS:` guaranteed valid; gsort's
	// descending key can leave the by-list unusable). Temp var is dropped before any mi command.
	tempvar nd
	gen double `nd' = -Date0
	sort ID_BS `nd'
	foreach v of local vlist {
		_cf `v'
	}
	drop `nd'
	sort ID_BS Date0
end


**********
// MI Settings
cap program drop mi_settings
program define mi_settings
	qui cap drop ISS
	qui cap drop RISS
	qui cap drop CM_CKD   // created in data_extraction from observed eGFR; re-derived below from imputed eGFR

	mi set wide
	mi register imputed Albumin AlkalinePhosphatase BMPlasmaCells LactateDehydrogenase SerumB2Microglobulin SerumCalcium SerumCreatinine eGFR EQ5D_Diagnosis LTHaemoglobinGL WhiteCellCount NeutrophillCount PlateletCount CRABScore CM_CRD CM_PLM CM_DBT Male FISHRisk ExtraMedullaryD LyticLesion ECOGcc Para Lambda Kappa FLC dPara dLambda dKappa dFLC BCR
	mi register regular Age CLine
	mi describe
end

// Impute diagnosis
cap program drop impute_diagnosis
program define impute_diagnosis
	args RN_diag

	cap noi mi impute chained (regress) Albumin AlkalinePhosphatase BMPlasmaCells LactateDehydrogenase SerumB2Microglobulin SerumCalcium SerumCreatinine eGFR EQ5D_Diagnosis LTHaemoglobinGL WhiteCellCount NeutrophillCount PlateletCount CRABScore Para Lambda Kappa FLC (logit, augment) Male CM_CRD CM_PLM CM_DBT FISHRisk ExtraMedullaryD LyticLesion (ologit, augment) ECOGcc = Age if Event0 == 3, add($imp) rseed(`RN_diag')
	if _rc {
		exit _rc
	}
		// ISS
		qui mi passive: gen ISS = 1 if SerumB2Microglobulin < 3.5 & Albumin >= 35 & Event0 == 3
		qui mi passive: replace ISS = 3 if SerumB2Microglobulin >= 5.5 & Event0 == 3
		qui mi passive: replace ISS = 2 if ISS == . & Event0 == 3
		label variable ISS "ISS (MI)"

		// LDHRisk
		qui mi passive: gen LDHRisk = 0 if LactateDehydrogenase <= LDHUpperLimit & Event0 == 3
		qui mi passive: replace LDHRisk = 1 if LactateDehydrogenase > LDHUpperLimit & Event0 == 3
		label variable LDHRisk "LDH Risk (MI)"

		// CM_CKD
		qui mi passive: gen CM_CKD = 1 if eGFR <= 59 & eGFR != . & Event0 == 3
		qui mi passive: replace CM_CKD = 0 if eGFR > 59 & Event0 == 3
		qui mi passive: replace CM_CKD = 0 if CM_CKD == . & Event0 == 3
		label variable CM_CKD "Chronic Kidney Disease"

		// Carryforward (LOCF within patient by Date0). Sort once; direct-column fills (no mi xeq).
		local vars "Male ECOGcc ISS CM_CKD CM_CRD CM_PLM CM_DBT"
		sort ID_BS Date0
		foreach v of local vars {
			_cf `v'
		}

		// RISS
		qui mi passive: gen RISS = 1 if ISS == 1 & LDHRisk == 0 & FISHRisk == 0 & Event0 == 3
		qui mi passive: replace RISS = 3 if ISS == 3 & (LDHRisk == 1 | FISHRisk == 1) & Event0 == 3
		qui mi passive: replace RISS = 2 if RISS == . & Event0 == 3
		_bcast_idbs RISS
		label variable RISS "Revised ISS (MI)"
end

// Impute BCR TXR
cap program drop impute_bcr_txr
program define impute_bcr_txr
	args RN_txr

	// TXR imputation
		cap noi mi impute chained (regress) dPara dLambda dKappa dFLC (ologit, augment) BCR = Age i.ECOGcc i.RISS i.CLine if CStart == 1 & Duration != ., replace rseed(`RN_txr')
		if _rc {
			exit _rc
		}
		
		// BCR_L1..L9: BCR at each line's start, copied FORWARD only
		forvalues l = 1/9 {
			qui mi passive: gen BCR_L`l' = BCR if Event0 == `l'0
		}
		sort ID_BS Date0
		forvalues l = 1/9 {
			_cf BCR_L`l'
		}
		forvalues l = 1/9 {
			label values BCR_L`l' BCR_label
		}

		// pBCR: previous line's response, needed only for the pooled L6+ BCR regressions, where it
		// is exactly BCR_L{line-1}. Left missing at L2 - risk_equations.do uses i.BCR_L1 there.
		qui mi passive: gen pBCR = .
		qui mi passive: replace pBCR = BCR_L5 if Event0 == 60
		qui mi passive: replace pBCR = BCR_L6 if Event0 == 70
		qui mi passive: replace pBCR = BCR_L7 if Event0 == 80
		qui mi passive: replace pBCR = BCR_L8 if Event0 == 90
		label values pBCR BCR_label
		
end

// Impute BCR SCT
cap program drop impute_bcr_sct
program define impute_bcr_sct
	args RN_sct

		// BCR_SCT - post-transplant response. ONE imputation, on BCR_SCT itself.
		//
		// A chained model used to impute BCR at Event0 == 100 with BCR_SCT derived from it, but it
		// left 548 of 2,135 transplanted patients (25.7%) missing without erroring, and the old
		// "replace BCR_SCT = 0 if BCR_SCT == ." absorbed the hole into the no-transplant code - so a
		// quarter of the transplant arm looked like a category, and four equations built different
		// workarounds on it (scratch/maintenance/_notes.md).
		//
		// Imputed at Event0 == 100 only - one row per patient, then broadcast; imputing the long
		// form would give one patient different responses on different rows. Uses BCR_L1, which is
		// why this must run after impute_bcr_txr.
		//
		// Assumes a missing assessment is MAR. If it is informative (died or progressed before
		// assessment) this biases toward better responses.
		qui mi passive: gen BCR_SCT = BCR if Event0 == 100
		mi unregister BCR_SCT
		mi register imputed BCR_SCT

		// Collapse BEFORE imputing, so imputed values can only be levels the engine can draw
		// (sim_bcr_asct.do categoryValues = 1,2,3,4). Small n at 5/6.
		qui mi xeq 0/$imp: replace BCR_SCT = 4 if inlist(BCR_SCT, 5, 6) & Event0 == 100

		// CHAINED, with the paraprotein deltas. They have no downstream consumer - they are
		// unregistered and dropped from the keep list below - and exist purely to inform the
		// response imputation, which is what they are clinically: the change in paraprotein and
		// light chains IS what determines best response. Chained also means their own missingness
		// does not listwise-delete rows from the BCR_SCT equation, and it keeps them imputed at
		// Event0 == 100 as the model this replaced did.
		cap noi mi impute chained (regress) dPara dLambda dKappa dFLC (ologit, augment) BCR_SCT ///
			= Age i.ECOGcc i.RISS i.BCR_L1 if Event0 == 100, replace rseed(`RN_sct')
		if _rc {
			di as error "  BCR_SCT imputation failed (rc = " _rc "): transplanted patients with no"
			di as error "  recorded response will drop from the SCT == 1 equations."
		}

		// Write back to BCR, then carry forward, so rows after the transplant see the POST-transplant
		// response rather than the pre-transplant one. At m = 0 BCR_SCT keeps its gaps and BCR is
		// already equal where it does not, so the master is untouched; nomaster on _cf for the same
		// reason - filling an imputed variable's master makes mi update treat it as observed.
		qui mi xeq 0/$imp: replace BCR = BCR_SCT if Event0 == 100 & !mi(BCR_SCT)
		sort ID_BS Date0
		_cf BCR "Duration != ." nomaster

		_bcast_idbs BCR_SCT
		qui mi xeq 0/$imp: replace BCR_SCT = 0 if BCR_SCT == . & SCT == 0
		qui mi xeq 0/$imp: label values BCR_SCT BCR_label

		// The invariant: BCR_SCT == 0 must count SCT == 0 patients EXACTLY. Anything else means the
		// zero-fill has caught missing data again, which is the fault this whole block exists to fix.
		//
		// CHECK INSIDE AN IMPUTATION, not the master. BCR_SCT is now a REGISTERED IMPUTED variable,
		// so m = 0 correctly keeps the original gaps and only m = 1..M are complete - that is how mi
		// works, not a fault. An earlier version of this check read the master and reported all 4,330
		// records as "still missing after imputation" when the imputation had in fact filled every
		// one of the 548 patients. `mi xeq 1:' is style-agnostic (wide or flong) where reading
		// _mi_m or _1_BCR_SCT directly is not.
		mi xeq 1: qui count if SCT == 1 & BCR_SCT == 0
		local _nbad = r(N)
		mi xeq 1: qui count if SCT == 1 & mi(BCR_SCT)
		local _nmiss = r(N)
		if `_nbad' > 0 {
			di as error "  BCR_SCT: " `_nbad' " TRANSPLANTED records still coded 0 - the invariant is broken."
		}
		if `_nmiss' > 0 {
			di as error "  BCR_SCT: " `_nmiss' " transplanted records still missing IN m = 1 after" ///
				" imputation - the model did not cover them; they will drop from the SCT == 1 equations."
		}
		if `_nbad' == 0 & `_nmiss' == 0 {
			di as txt "  BCR_SCT: complete in m = 1 for all transplanted records; 0 means no transplant only."
			di as txt "           (m = 0 retains the original gaps, as an imputed variable should.)"
		}
end

// Finalise. Was left at top level by the split into separate programs, where it would have run at
// file-read time with no data loaded - and the bootstrap branch never reached it at all, so that
// output would have kept every auxiliary variable.
cap program drop finalise_mi
program define finalise_mi

	// AFTER every direct-column write. _cf and _bcast_idbs bypass mi to write the `_m_var' columns,
	// so _mi_miss is stale until this runs; and BEFORE the unregister, which needs a consistent mi
	// object. This is the one mi update in the file, in the same position the original had it.
	mi update

	// Unregister, keep, sort & order
	mi unregister AlkalinePhosphatase BMPlasmaCells SerumCalcium SerumCreatinine EQ5D_Diagnosis LTHaemoglobinGL WhiteCellCount NeutrophillCount PlateletCount CRABScore ExtraMedullaryD LyticLesion Para Lambda Kappa FLC dPara dLambda dKappa dFLC
	keep ID ID_BS Event* Date* Age* Male ECOGcc ISS RISS SCT MNT MND_L1 MNR_L1 LineRefr LenRefr_Tx_in LenRefr_Mnt_in MNT_LenRefr_L1 CM* BCR* pBCR* Reg* OS Line Duration CID CLine CStart CEnd Country F_* CN_* Year Albumin SerumB2Microglobulin LactateDehydrogenase LDHUpperLimit LDHRisk FISHRisk eGFR _* Bortezomib Carfilzomib Cisplatin Cyclophosphamide Daratumamab Dexamethasone Doxorubicin Elotuzamab Etoposide Lenalidomide Melphalan Methylprednisolone Panobinostat Prednisolone Thalidomide Pomalidomide Ixazomib TXD* TFI*
	sort ID_BS Date0
	order $core Age Male ECOGcc RISS BCR Reg Regimen Line Duration
end

**********
// Execute based on arguments

if "$boot" == "0" {

	// Create temp folder on the shared drive (the locally-synced repo confuses Google Drive;
	// $data_path is mounted whenever this script runs).
	local repo = c(pwd)
	cap mkdir "~/temp"
	cd "~/temp"

	// Open MRDR Long Data
	use "${data_path}/MRDR Long.dta"
	cap drop CM_LVR CM_PNR CM_MLG   // unused comorbidities; dropped before mi set
	gen ID_BS = ID

	// OOS: restrict to the requested fold (train/test) before imputing
	if "$sample" != "" {
		merge m:1 ID using "${data_path}/oos/oos_split.dta", keep(match) keepusing(fold) nogen
		keep if fold == "$sample"
		drop fold
	}

	// Seeds. Each program takes only the ones it uses, so a miscount is a syntax error rather than
	// a silently empty rseed().
	local RN_diag = 3949
	local RN_txr  = 6192
	local RN_sct  = 5117

	mi_settings
	impute_diagnosis `RN_diag'
	impute_bcr_txr   `RN_txr'
	impute_bcr_sct   `RN_sct'
	finalise_mi

	// Save Long MI
	save "${data_path}/${mi_outdir}MRDR Long MI${mi_outtag}.dta", replace

	// Create Wide MI to create synthetic patient datasets
	keep if Event0 == 3

		// Convert _1_var to var_1
		foreach v in Male ECOGcc RISS ISS LDHRisk FISHRisk CM_CKD CM_CRD CM_PLM CM_DBT {
			forvalues i = 1/$imp {
				gen `v'_`i' = _`i'_`v'
				drop _`i'_`v'
			}
		}

		// Drop non-MI variables
		keep ID Event0 Date0 Age Male* ECOGc* RISS* ISS* LDHRisk* FISHRisk* CM_* _mi_miss
		drop Male ECOGcc RISS ISS LDHRisk FISHRisk CM_CKD CM_CRD CM_PLM CM_DBT

		// Turn MI imps into rows
		mi unset
		reshape long Male_ ECOGcc_ RISS_ ISS_ LDHRisk_ FISHRisk_ CM_CKD_ CM_CRD_ CM_PLM_ CM_DBT_, i(ID) j(Imp)
		rename Male_ Male
		rename ECOGcc_ ECOGcc
		rename RISS_ RISS
		rename ISS_ ISS
		rename LDHRisk_ LDHRisk
		rename FISHRisk_ FISHRisk
		rename CM_CKD_ CM_CKD
		rename CM_CRD_ CM_CRD
		rename CM_PLM_ CM_PLM
		rename CM_DBT_ CM_DBT
		drop mi_miss

		// Save Wide MI
		gen MRDR = 1
		order MRDR ID Imp Event0 Date0 Age Male ECOGcc RISS ISS LDHRisk FISHRisk CM_CKD CM_CRD CM_PLM CM_DBT
		save "${data_path}/${mi_outdir}MRDR Wide MI${mi_outtag}.dta", replace

	// Return to the repo root, then delete temp folder
	// (must cd out first, or Stata is left sitting in the just-removed directory)
	cd "`repo'"
	cap rmdir "~/temp"
}
else if "$boot" == "1" {
	forval b = $min_bs / $max_bs {

		// Open MRDR Long Data
		use "${data_path}/MRDR Long.dta"
		cap drop CM_LVR CM_PNR CM_MLG   // unused comorbidities (engine uses only CM_CKD/CRD/PLM/DBT)

		// OOS: restrict to the requested fold before resampling/imputing
		if "$sample" != "" {
			merge m:1 ID using "${data_path}/oos/oos_split.dta", keep(match) keepusing(fold) nogen
			keep if fold == "$sample"
			drop fold
		}

		// Create per-iteration temp folder on the shared drive - needed for M3 array jobs
		// (keeps scratch off the locally-synced repo; $data_path is mounted at run time)
		local curdir = c(pwd)
		capture mkdir "${data_path}/temp"
		capture mkdir "${data_path}/temp/job_`b'"
		cd "${data_path}/temp/job_`b'"

		// Base random numbers
		local RN_diag = 7394
		local RN_txr  = 1392
		local RN_sct  = 8856

		// Retry settings
		local maxtries = 50
		local success = 0
		local try = 0

		// Snapshot the clean (pre-resample) data for restarts
		preserve

		while `success' == 0 & `try' < `maxtries' {
			local ++try

			// Restore clean data before each attempt
			restore, preserve

			// Attempt-specific, deterministic seeds
			local RS  = 5839 * `b' + 100003 * `try'
			local a_diag = `RN_diag' + 911 * `try'
			local a_txr  = `RN_txr'  + 911 * `try'
			local a_sct  = `RN_sct'  + 911 * `try'

			// Resample with this attempt's seed
			set seed `RS'
			bsample, cluster(ID) idcluster(ID_BS)

			// Try imputation; trap perfect-predictor / convergence failures.
			// The `multiple_imputation' wrapper was removed when the stages were split into
			// separate programs, so this called a program that does not exist. Bootstrap always
			// runs all three stages - the stage flags are a development convenience and
			// checkpointing a resampled dataset would be meaningless.
			capture noisily {
				mi_settings
				impute_diagnosis `a_diag'
				impute_bcr_txr   `a_txr'
				impute_bcr_sct   `a_sct'
				finalise_mi
			}

			if _rc == 0 {
				local success = 1
			}
			else {
				di as txt "Bootstrap `b': attempt `try' failed (rc=" _rc "), retrying with new seeds"
			}
		}

		// Drop the snapshot
		restore, not

		if `success' == 0 {
			di as error "Bootstrap `b' failed after `maxtries' attempts"
			exit 459
		}

		// Save Long MI
		cd "`curdir'"
		save "${data_path}/${mi_outdir}bootstrap/MRDR Long MI${mi_outtag} B`b'.dta", replace

		// Delete temporary folder
		cap rmdir "${data_path}/temp/job_`b'"
	}

	// Remove the now-empty parent temp dir
	cap rmdir "${data_path}/temp"
}
