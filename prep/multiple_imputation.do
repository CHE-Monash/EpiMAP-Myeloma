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
cap program drop impute_bcr
program define impute_bcr
	args RN_l1 RN_sct RN_l2 RN_l3 RN_l4 RN_l5 RN_l6

	// PER-LINE response imputation, each line conditioning on the previous line's response so the
	// imputation mirrors the risk-equation model at that line (congeniality). The pooled model this
	// replaced (BCR = ... i.CLine, no previous response) attenuated the L1->L2 association by a
	// third to a half (scratch/bcr_congeniality.log).
	//
	// THE OUTCOME IS A DEDICATED PER-LINE VARIABLE (iL{l}), NOT the shared BCR. Imputing BCR while
	// conditioning on i.BCR_L{prev} - a passive function of BCR - puts BCR on both sides of its own
	// equation, and mi refuses it ("missing imputed values produced", r(498)). Imputing iL{l}
	// instead, then writing it back into BCR, breaks the self-reference - the same pattern the
	// BCR_SCT block already uses.
	//
	// L1-L5 SEPARATE (cells support it), L6-L9 POOLED on the L5 response (later cells too thin and
	// fold-dependent to split). BCR_SCT is INTERLEAVED after L1, because L2 conditions on it.
	local aux "dPara dLambda dKappa dFLC"
	local base "Age i.ECOGcc i.RISS"

	// ---- L1 (Event0 == 10): baseline ----
	qui gen iL1 = BCR if Event0 == 10
	mi register imputed iL1
	cap noi mi impute chained (regress) `aux' (ologit, augment) iL1 = `base' ///
		if Event0 == 10 & CStart == 1 & Duration != ., replace rseed(`RN_l1')
	if _rc {
		exit _rc
	}
	qui mi xeq 0/$imp: replace BCR = iL1 if Event0 == 10 & !mi(iL1)
	cap mi unregister iL1
	forvalues m = 1/$imp {
		cap drop _`m'_iL1
	}
	cap drop iL1
	qui mi passive: gen BCR_L1 = BCR if Event0 == 10
	sort ID_BS Date0
	_cf BCR_L1
	label values BCR_L1 BCR_label
	// ---- BCR_SCT (Event0 == 100): baseline + BCR_L1. Interleaved so L2 can condition on it. ----
	qui mi passive: gen BCR_SCT = BCR if Event0 == 100
	mi unregister BCR_SCT
	mi register imputed BCR_SCT

	// Collapse BEFORE imputing, so imputed values can only be levels the engine draws
	// (sim_bcr_asct.do categoryValues = 1,2,3,4).
	qui mi xeq 0/$imp: replace BCR_SCT = 4 if inlist(BCR_SCT, 5, 6) & Event0 == 100

	// Chained on the paraprotein deltas: they have no downstream consumer and exist purely to
	// inform this imputation - the change in paraprotein and light chains is what determines
	// response. Chained also keeps their missingness from listwise-deleting BCR_SCT rows.
	cap noi mi impute chained (regress) dPara dLambda dKappa dFLC (ologit, augment) BCR_SCT ///
		= Age i.ECOGcc i.RISS i.BCR_L1 if Event0 == 100, replace rseed(`RN_sct')
	if _rc {
		di as error "  BCR_SCT imputation failed (rc = " _rc "): transplanted patients with no"
		di as error "  recorded response will drop from the SCT == 1 equations."
	}

	// Write back to BCR before carrying forward, so later rows see the POST-transplant response.
	// m = 0 is untouched: BCR_SCT keeps its gaps there and equals BCR elsewhere. nomaster for the
	// usual reason - filling an imputed variable's master makes mi update treat it as observed.
	qui mi xeq 0/$imp: replace BCR = BCR_SCT if Event0 == 100 & !mi(BCR_SCT)
	sort ID_BS Date0
	_cf BCR "Duration != ." nomaster

	_bcast_idbs BCR_SCT
	qui mi xeq 0/$imp: replace BCR_SCT = 0 if BCR_SCT == . & SCT == 0
	qui mi xeq 0/$imp: label values BCR_SCT BCR_label

	// Invariant: BCR_SCT == 0 must count SCT == 0 patients exactly, or the zero-fill has caught
	// missing data again. Checked in m = 1 - the master legitimately keeps its gaps.
	mi xeq 1: qui count if SCT == 1 & BCR_SCT == 0
	local _nbad = r(N)
	mi xeq 1: qui count if SCT == 1 & mi(BCR_SCT)
	local _nmiss = r(N)
	if `_nbad' > 0 {
		di as error "  BCR_SCT: " `_nbad' " TRANSPLANTED records coded 0 - invariant broken."
	}
	if `_nmiss' > 0 {
		di as error "  BCR_SCT: " `_nmiss' " transplanted records still missing in m = 1 - the model" ///
			" did not cover them; they will drop from the SCT == 1 equations."
	}
	if `_nbad' == 0 & `_nmiss' == 0 {
		di as txt "  BCR_SCT: complete in m = 1 for all transplanted records; 0 means no transplant only."
	}
	// ---- L2 (Event0 == 20): baseline + BCR_L1 + BCR_SCT ----
	qui gen iL2 = BCR if Event0 == 20
	mi register imputed iL2
	cap noi mi impute chained (regress) `aux' (ologit, augment) iL2 = `base' i.BCR_L1 i.BCR_SCT ///
		if Event0 == 20 & CStart == 1 & Duration != ., replace rseed(`RN_l2')
	if _rc {
		exit _rc
	}
	qui mi xeq 0/$imp: replace BCR = iL2 if Event0 == 20 & !mi(iL2)
	cap mi unregister iL2
	forvalues m = 1/$imp {
		cap drop _`m'_iL2
	}
	cap drop iL2
	qui mi passive: gen BCR_L2 = BCR if Event0 == 20
	sort ID_BS Date0
	_cf BCR_L2
	label values BCR_L2 BCR_label
	// ---- L3 (Event0 == 30): baseline + previous response (LOCF) ----
	// Previous response, carried forward (LOCF), COMPLETE. Strict i.BCR_L{l-1} has gaps on this
	// line's sample - patients who reach line l whose earlier line fell outside the imputation
	// restriction (CStart == 1 & Duration != .) never had it filled - and mi rejects a predictor
	// with any missing. Carrying the last assembled response forward fills the gaps (L1 is imputed
	// first, so there is always a fallback). _pr at a line row = the most recent response BEFORE it,
	// since this line's own response is not assembled yet.
	cap drop _pr
	qui mi passive: gen _pr = BCR
	sort ID_BS Date0
	_cf _pr
	label values _pr BCR_label
	qui gen iL3 = BCR if Event0 == 30
	mi register imputed iL3
	cap noi mi impute chained (regress) `aux' (ologit, augment) iL3 = `base' i._pr ///
		if Event0 == 30 & CStart == 1 & Duration != ., replace rseed(`RN_l3')
	if _rc {
		exit _rc
	}
	qui mi xeq 0/$imp: replace BCR = iL3 if Event0 == 30 & !mi(iL3)
	cap mi unregister iL3
	forvalues m = 1/$imp {
		cap drop _`m'_iL3
	}
	cap drop iL3
	qui mi passive: gen BCR_L3 = BCR if Event0 == 30
	sort ID_BS Date0
	_cf BCR_L3
	label values BCR_L3 BCR_label
	// ---- L4 (Event0 == 40): baseline + previous response (LOCF) ----
	// Previous response, carried forward (LOCF), COMPLETE. Strict i.BCR_L{l-1} has gaps on this
	// line's sample - patients who reach line l whose earlier line fell outside the imputation
	// restriction (CStart == 1 & Duration != .) never had it filled - and mi rejects a predictor
	// with any missing. Carrying the last assembled response forward fills the gaps (L1 is imputed
	// first, so there is always a fallback). _pr at a line row = the most recent response BEFORE it,
	// since this line's own response is not assembled yet.
	cap drop _pr
	qui mi passive: gen _pr = BCR
	sort ID_BS Date0
	_cf _pr
	label values _pr BCR_label
	qui gen iL4 = BCR if Event0 == 40
	mi register imputed iL4
	cap noi mi impute chained (regress) `aux' (ologit, augment) iL4 = `base' i._pr ///
		if Event0 == 40 & CStart == 1 & Duration != ., replace rseed(`RN_l4')
	if _rc {
		exit _rc
	}
	qui mi xeq 0/$imp: replace BCR = iL4 if Event0 == 40 & !mi(iL4)
	cap mi unregister iL4
	forvalues m = 1/$imp {
		cap drop _`m'_iL4
	}
	cap drop iL4
	qui mi passive: gen BCR_L4 = BCR if Event0 == 40
	sort ID_BS Date0
	_cf BCR_L4
	label values BCR_L4 BCR_label
	// ---- L5 (Event0 == 50): baseline + previous response (LOCF) ----
	// Previous response, carried forward (LOCF), COMPLETE. Strict i.BCR_L{l-1} has gaps on this
	// line's sample - patients who reach line l whose earlier line fell outside the imputation
	// restriction (CStart == 1 & Duration != .) never had it filled - and mi rejects a predictor
	// with any missing. Carrying the last assembled response forward fills the gaps (L1 is imputed
	// first, so there is always a fallback). _pr at a line row = the most recent response BEFORE it,
	// since this line's own response is not assembled yet.
	cap drop _pr
	qui mi passive: gen _pr = BCR
	sort ID_BS Date0
	_cf _pr
	label values _pr BCR_label
	qui gen iL5 = BCR if Event0 == 50
	mi register imputed iL5
	cap noi mi impute chained (regress) `aux' (ologit, augment) iL5 = `base' i._pr ///
		if Event0 == 50 & CStart == 1 & Duration != ., replace rseed(`RN_l5')
	if _rc {
		exit _rc
	}
	qui mi xeq 0/$imp: replace BCR = iL5 if Event0 == 50 & !mi(iL5)
	cap mi unregister iL5
	forvalues m = 1/$imp {
		cap drop _`m'_iL5
	}
	cap drop iL5
	qui mi passive: gen BCR_L5 = BCR if Event0 == 50
	sort ID_BS Date0
	_cf BCR_L5
	label values BCR_L5 BCR_label
	// ---- L6-L9 POOLED (Event0 60-90): baseline + previous response (LOCF) ----
	// Previous response, carried forward (LOCF), COMPLETE. Strict i.BCR_L{l-1} has gaps on this
	// line's sample - patients who reach line l whose earlier line fell outside the imputation
	// restriction (CStart == 1 & Duration != .) never had it filled - and mi rejects a predictor
	// with any missing. Carrying the last assembled response forward fills the gaps (L1 is imputed
	// first, so there is always a fallback). _pr at a line row = the most recent response BEFORE it,
	// since this line's own response is not assembled yet.
	cap drop _pr
	qui mi passive: gen _pr = BCR
	sort ID_BS Date0
	_cf _pr
	label values _pr BCR_label
	qui gen iL69 = BCR if inlist(Event0, 60, 70, 80, 90)
	mi register imputed iL69
	cap noi mi impute chained (regress) `aux' (ologit, augment) iL69 = `base' i._pr ///
		if inlist(Event0, 60, 70, 80, 90) & CStart == 1 & Duration != ., replace rseed(`RN_l6')
	if _rc {
		exit _rc
	}
	qui mi xeq 0/$imp: replace BCR = iL69 if inlist(Event0, 60, 70, 80, 90) & !mi(iL69)
	cap mi unregister iL69
	forvalues m = 1/$imp {
		cap drop _`m'_iL69
	}
	cap drop iL69
	forvalues l = 6/9 {
		qui mi passive: gen BCR_L`l' = BCR if Event0 == `l'0
	}
	sort ID_BS Date0
	forvalues l = 6/9 {
		_cf BCR_L`l'
		label values BCR_L`l' BCR_label
	}

	// pBCR for the risk equations' pooled L6+ regressions: the IMPUTED previous response
	// (BCR_L{line-1}), replacing the OBSERVED pBCR that data_extraction.do builds as a predictor.
	cap drop pBCR
	qui mi passive: gen pBCR = .
	qui mi passive: replace pBCR = BCR_L5 if Event0 == 60
	qui mi passive: replace pBCR = BCR_L6 if Event0 == 70
	qui mi passive: replace pBCR = BCR_L7 if Event0 == 80
	qui mi passive: replace pBCR = BCR_L8 if Event0 == 90
	label values pBCR BCR_label

	cap drop _pr

	// Completeness: every line-start response must be imputed. Fires if a per-line impute missed a line.
	mi xeq 1: qui count if inlist(Event0, 10, 20, 30, 40, 50, 60, 70, 80, 90) & mi(BCR) & CStart == 1 & Duration != .
	if r(N) > 0 {
		di as error "  impute_bcr: " r(N) " line-start responses still missing in m = 1 - a per-line"
		di as error "  impute did not cover its line. Do NOT trust this imputation."
	}
	else {
		di as txt "  impute_bcr: all line-start responses imputed (per-line L1-L5, pooled L6+)."
	}
end

cap program drop finalise_mi
program define finalise_mi

	// AFTER every direct-column write. _cf and _bcast_idbs bypass mi to write the `_m_var' columns,
	// so _mi_miss is stale until this runs; and BEFORE the unregister, which needs a consistent mi object.
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
	local RN_l1   = 6192
	local RN_sct  = 5117
	local RN_l2   = 2731
	local RN_l3   = 4409
	local RN_l4   = 8102
	local RN_l5   = 1567
	local RN_l6   = 9284

	mi_settings
	impute_diagnosis `RN_diag'
	impute_bcr       `RN_l1' `RN_sct' `RN_l2' `RN_l3' `RN_l4' `RN_l5' `RN_l6'
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
	cd "`repo'"
	cap rmdir "~/temp"
}
else if "$boot" == "1" {
	forval b = $min_bs / $max_bs {

		// Open MRDR Long Data
		use "${data_path}/MRDR Long.dta"
		cap drop CM_LVR CM_PNR CM_MLG   // unused comorbidities

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
		local RN_l1   = 1392
		local RN_sct  = 8856
		local RN_l2   = 3517
		local RN_l3   = 6640
		local RN_l4   = 2298
		local RN_l5   = 7051
		local RN_l6   = 4483

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
			local a_l1   = `RN_l1'   + 911 * `try'
			local a_sct  = `RN_sct'  + 911 * `try'
			local a_l2   = `RN_l2'   + 911 * `try'
			local a_l3   = `RN_l3'   + 911 * `try'
			local a_l4   = `RN_l4'   + 911 * `try'
			local a_l5   = `RN_l5'   + 911 * `try'
			local a_l6   = `RN_l6'   + 911 * `try'

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
				impute_bcr       `a_l1' `a_sct' `a_l2' `a_l3' `a_l4' `a_l5' `a_l6'
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
