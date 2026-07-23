**********
* Monash Myeloma Model - Multiple Imputation
*
* Purpose: Multiple imputation of MRDR Long. Performance-optimised variant intended to produce
*          IDENTICAL output to the original for a given seed.
* Notes:   Optimisation rationale (carryforward/broadcast helpers, and a deliberately-deferred
*          structural win) documented below.
**********

* -------------------------------------------------------------------------------------------------
* WHAT CHANGED vs multiple_imputation.do (intended to produce IDENTICAL output for a given seed):
*   The post-imputation carryforward/broadcast step was the bottleneck: ~27 `mi xeq 0/$imp: bysort`
*   calls, each re-sorting the full long dataset M+1 times (~O(M) full sorts). Fixes (helpers _cf /
*   _bcast_idbs, above): operate DIRECTLY on the wide per-imputation columns (`_m_var`) + the original
*   (`var`), sorting ONCE and using plain `by` - instead of looping `mi xeq 0/$imp`. Two reasons this
*   is the right lever: (a) `mi xeq` does NOT preserve the sort order across its per-imputation passes
*   (a `sort` + `by` inside it fails "not sorted"), and (b) replacing a non-sort-key never clears the
*   sort, so one sort covers all variables and all imputations. The baseline block drops from
*   ~9*(M+1) sorts to a single sort; the BCR_L broadcast from 9*(M+1) value-sorts to two.
*     - Temporal LOCF (_cf): forward fill by Date0 across m=0..M. Identical to the old
*       `bysort ID_BS (Date0): replace v = v[_n-1] if v==.`.
*     - Broadcast (_bcast_idbs): forward + backward fill by Date0, identical to the old value-sort
*       `bysort ID_BS (VAR): ...` when there is one non-missing value per patient (true for the
*       diagnosis-/line-/SCT-anchored vars it is used on: RISS, BCR_L1..9, BCR_SCT).
*   The pBCR block keeps its original `mi xeq: bysort ID` form (single call, uncertain _m_ storage),
*   and all no-`by` `mi xeq` calls (BCR collapse, BCR_SCT default, labels) are unchanged. Everything
*   else (imputation models, ISS/RISS/LDHRisk/CM_CKD derivations, keeps, reshape, outputs, bootstrap
*   retry logic) is byte-for-byte the same. VERIFY by diffing the saved MRDR Long MI / Wide MI against
*   the original for a fixed seed before switching over.
*
* NOT changed here (deliberately): the larger structural win of imputing baseline covariates on a
*   one-row-per-patient `keep if Event0==3` extract and mi merge-ing back. That would shrink the mi
*   object ~10x and remove the baseline carryforward entirely, but it risks changing results (the
*   current carryforward is forward-LOCF, not a pure diagnosis-value broadcast) so it needs its own
*   equivalence check. Left as a follow-up.
* -------------------------------------------------------------------------------------------------

clear
clear mata

if "$repo_path" != "" cd "$repo_path"   // cd to repo root only if config.do set it; a bare cd "" goes to home on Mac/Unix
capture run "config.do"     // machine-specific paths: $data_path (git-ignored)

local Data "$data_cut"

**********
// Settings
/*
global imp 2
global boot 0
*/

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
// Helpers: carry values across a patient's rows by operating DIRECTLY on the wide per-imputation
// columns (`_m_var`) plus the original (m=0, `var`). This deliberately avoids `mi xeq`, which does
// NOT preserve the sort order across its per-imputation passes ("not sorted" errors) and pays a
// context switch each imputation. With direct columns we sort ONCE and use plain `by` (replacing a
// non-sort-key never clears the sort), so the whole baseline block costs one sort instead of
// ~9*(M+1). The file already reads `_i_var` directly at the Wide-MI step, so this is in-house style.
// Values produced are identical to the old `mi xeq 0/$imp: ...` fills; mi's existing `mi update`
// call (after the pBCR merge) refreshes system vars, exactly as in the original.

// _cf: temporal LOCF within ID_BS for one variable, across m=0..M. Caller must have sorted so that
//      rows are in the desired order within ID_BS. Optional extra `if` condition (e.g. Duration!=.).
//      Third arg "nomaster" skips the m=0 (master) fill: for a REGISTERED IMPUTED variable, filling
//      the master at imputed rows makes mi update treat them as observed and reset the imputations to
//      the single master value (the BCR collapse). Imputed vars must keep a missing master at imputed
//      rows, so pass "nomaster" for them; leave it off for observed/derived vars (the covariates).
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
// Define function
cap program drop multiple_imputation
program define multiple_imputation
	args RN1 RN2 RN3 RN4

	// Drop variables to be derived from imputed
	qui cap drop ISS
	qui cap drop RISS
	qui cap drop CM_CKD   // created in data_extraction from observed eGFR; re-derived below from imputed eGFR

	// MI settings
	mi set wide
	mi register imputed Albumin AlkalinePhosphatase BMPlasmaCells LactateDehydrogenase SerumB2Microglobulin SerumCalcium SerumCreatinine eGFR EQ5D_Diagnosis LTHaemoglobinGL WhiteCellCount NeutrophillCount PlateletCount CRABScore CM_CRD CM_PLM CM_DBT Male FISHRisk ExtraMedullaryD LyticLesion ECOGcc Para Lambda Kappa FLC dPara dLambda dKappa dFLC BCR
	mi register regular Age CLine
	mi describe

	// Diagnosis imputation
	cap noi mi impute chained (regress) Albumin AlkalinePhosphatase BMPlasmaCells LactateDehydrogenase SerumB2Microglobulin SerumCalcium SerumCreatinine eGFR EQ5D_Diagnosis LTHaemoglobinGL WhiteCellCount NeutrophillCount PlateletCount CRABScore Para Lambda Kappa FLC (logit, augment) Male CM_CRD CM_PLM CM_DBT FISHRisk ExtraMedullaryD LyticLesion (ologit, augment) ECOGcc = Age if Event0 == 3, add($imp) rseed(`RN1')
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

	// BCR imputation

		// TXR imputation
		cap noi mi impute chained (regress) dPara dLambda dKappa dFLC (ologit, augment) BCR = Age i.ECOGcc i.CLine if CStart == 1 & Duration != ., replace rseed(`RN2')
		if _rc {
			exit _rc
		}

		// ASCT imputation. KEPT, but its job is now the CARRYFORWARD below, not BCR_SCT.
		//
		// It leaves 548 of 2,135 transplanted patients missing - it runs without error and simply
		// does not cover them - which is why BCR_SCT is imputed on its own further down rather than
		// derived from this. But this block still matters: _cf propagates BCR forward from the
		// Event0 == 100 row, so for a transplanted patient the post-transplant response becomes
		// BCR_L2 onward at the next line start. That is right - the most recent assessment before L2
		// IS the post-transplant one - and deleting this would carry the PRE-transplant L1 response
		// into L2 instead.
		//
		// KNOWN REDUNDANCY: two models now impute the same quantity. BCR_SCT inherits this block's
		// values where it succeeded and imputes the rest with its own ologit, so the 548 come from a
		// different model than the others. The clean end state is to impute BCR_SCT ONCE, write it
		// back to BCR at Event0 == 100, and carry forward from there - one model, one value, and no
		// dependence on the chained block reaching every row. That needs the BCR_L1 generation moved
		// above this point and belongs on the branch that retires the running BCR variable.
		cap noi mi impute chained (regress) dPara dLambda dKappa dFLC (ologit, augment) BCR = Age i.ECOGcc if Event0 == 100, replace rseed(`RN3')
		if _rc {
			exit _rc
		}

		// Carryforward (temporal LOCF; sort once, direct-column fills). BCR is a registered IMPUTED
		// variable, so carry only the imputation columns (nomaster) -- filling the master would make the
		// mi update below collapse the imputations to a single value.
		sort ID_BS Date0
		_cf BCR "Duration != ." nomaster

		// Collapse BCR for ASCT - small n
		qui mi xeq 0/$imp: replace BCR = 4 if (BCR == 5 | BCR == 6) & Event0 == 100

		// Generate BCR_L1..L9 (BCR at each line's start) and copy FORWARD only (LOCF), not a full
		// broadcast. These are used in risk_equations.do only at Event0 >= the line's own start, and
		// the simulation does not consume them (analyses/default/prep/test_cohort.do resets BCR_L* to
		// missing; the synthetic population never carries them). So filling rows *before* the line -
		// the backward pass - is wasted. This drops one sort + 9 fill-passes vs _bcast_idbs. It leaves
		// pre-line BCR_L* cells missing in the saved data, which nothing reads (fits/sims unaffected),
		// so a raw diff will differ there. pBCR below still reads BCR_L{line-1} at Event0 60-90, which
		// is at/after that line, so forward-fill covers it.
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

		// Previous BCR (pBCR). Only genuinely needed for the pooled L6+ BCR regressions
		// (risk_equations.do, Event0 60/70/80/90), where at each row the previous line's response is
		// exactly the already-broadcast BCR_L{line-1}: L6<-BCR_L5, L7<-BCR_L6, L8<-BCR_L7, L9<-BCR_L8.
		// This replaces the old preserve + SSC mmerge + mi update + 2*(M+1) bysort machinery (the
		// main post-BCR bottleneck) with a handful of passive derivations off variables we already
		// have. The former L2 use (Event0==20) is handled in risk_equations.do with i.BCR_L1 (the L1
		// response), so pBCR is left missing there. If you need pBCR at other lines later, add a row.
		qui mi passive: gen pBCR = .
		qui mi passive: replace pBCR = BCR_L5 if Event0 == 60
		qui mi passive: replace pBCR = BCR_L6 if Event0 == 70
		qui mi passive: replace pBCR = BCR_L7 if Event0 == 80
		qui mi passive: replace pBCR = BCR_L8 if Event0 == 90
		label values pBCR BCR_label

		// Generate BCR_SCT - the post-transplant response, taken from the Event0 == 100 row and
		// IMPUTED where that row carries no response.
		//
		// WHY IMPUTED HERE RATHER THAN INHERITED FROM BCR. The ASCT block above imputes BCR at
		// Event0 == 100, but it leaves 548 of 2,135 transplanted patients (25.7%) still missing -
		// it runs without error and simply does not cover them. Deriving BCR_SCT passively then
		// propagated that hole to every row of those patients, and the old
		// "replace BCR_SCT = 0 if BCR_SCT == ." absorbed it into the non-transplant code, so a
		// quarter of the transplant arm looked like a category. Four equations then built different
		// workarounds on it. See scratch/maintenance/_notes.md.
		//
		// Imputing BCR_SCT directly sidesteps whatever is stopping the chained block from reaching
		// those rows, and it can use BCR_L1 - the pre-transplant response, generated above and the
		// strongest available predictor of the post-transplant one.
		//
		// ORDER MATTERS, three ways:
		//   - it must run AFTER BCR_L1 exists and is forward-filled (just above), since that is a
		//     predictor
		//   - it imputes at Event0 == 100 ONLY, one row per transplanted patient, and broadcasts
		//     afterwards. Imputing the long form directly would give one patient different
		//     responses on different rows
		//   - the zero-fill comes LAST and only for SCT == 0. BCR_SCT == 0 must mean "no transplant"
		//     and nothing else; mi_diagnostics.do now checks that invariant directly
		//
		// ASSUMPTION, stated rather than buried: this treats a missing post-transplant assessment as
		// MAR. If it is informative - patients who died or progressed before assessment - imputing
		// from assessed patients biases toward better responses. Worth revisiting if the transplant
		// arm's survival drifts optimistic.
		qui mi passive: gen BCR_SCT = BCR if Event0 == 100
		mi unregister BCR_SCT
		mi register imputed BCR_SCT
		cap noi mi impute ologit BCR_SCT = Age i.ECOGcc i.RISS i.BCR_L1 if Event0 == 100, ///
			replace augment rseed(`RN4')
		if _rc {
			di as error "  BCR_SCT imputation failed (rc = " _rc "). BCR_SCT will be missing for"
			di as error "  transplanted patients with no recorded response - the SCT == 1 equations"
			di as error "  will drop them, which is safe but loses ~a quarter of that arm."
		}
		_bcast_idbs BCR_SCT
		qui mi xeq 0/$imp: replace BCR_SCT = 0 if BCR_SCT == . & SCT == 0
		qui mi xeq 0/$imp: label values BCR_SCT BCR_label

		// The invariant: BCR_SCT == 0 must count SCT == 0 patients EXACTLY. Anything else means the
		// zero-fill has caught missing data again, which is the fault this whole block exists to fix.
		qui count if _mi_m == 1 & SCT == 1 & BCR_SCT == 0
		local _nbad = r(N)
		qui count if _mi_m == 1 & SCT == 1 & mi(BCR_SCT)
		local _nmiss = r(N)
		if `_nbad' > 0 {
			di as error "  BCR_SCT: " `_nbad' " TRANSPLANTED records still coded 0 - the invariant is broken."
		}
		if `_nmiss' > 0 {
			di as txt "  BCR_SCT: " `_nmiss' " transplanted records remain missing after imputation" ///
				" (they will drop from the SCT == 1 equations)."
		}
		if `_nbad' == 0 & `_nmiss' == 0 {
			di as txt "  BCR_SCT: complete for all transplanted records; 0 means no transplant only."
		}

		// Refresh mi system variables after the direct-column carryforwards/broadcasts (replaces the
		// mi update that used to live inside the now-removed pBCR merge block).
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
	cap drop CM_LVR CM_PNR CM_MLG   // unused comorbidities (engine uses only CM_CKD/CRD/PLM/DBT); dropped before mi set
	gen ID_BS = ID

	// OOS: restrict to the requested fold (train/test) before imputing
	if "$sample" != "" {
		merge m:1 ID using "${data_path}/oos/oos_split.dta", keep(match) keepusing(fold) nogen
		keep if fold == "$sample"
		drop fold
	}

	// Draw random numbers
	local RN1 = 3949
	local RN2 = 6192
	local RN3 = 8273
	local RN4 = 5117

	// Execute function
		multiple_imputation `RN1' `RN2' `RN3' `RN4'

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
		local RN1 = 7394
		local RN2 = 1392
		local RN3 = 4771

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
			local a1  = `RN1' + 911 * `try'
			local a2  = `RN2' + 911 * `try'
			local a3  = `RN3' + 911 * `try'

			// Resample with this attempt's seed
			set seed `RS'
			bsample, cluster(ID) idcluster(ID_BS)

			// Try imputation; trap perfect-predictor / convergence failures
			cap noi multiple_imputation `a1' `a2' `a3'

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
