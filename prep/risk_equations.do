**********
* Monash Myeloma Model - Risk Equations
*
* Purpose: Fit every model risk equation (grouped by outcome type) from the MI data and matsave
*          the coefficients to analyses/$analysis/coefficients/. Called by the prep pipeline.
* Usage:   do risk_equations.do <analysis> <coeffs> <minyr> <maxyr> <boot> <minbs> <maxbs> [sample]
**********

clear
clear mata

local Data "$data_cut"

**********
// Capture arguments
global analysis `1'
global coeffs `2'
global min_year `3'
global max_year `4'
global boot `5'
global min_bs `6'
global max_bs `7'
global sample `8'   // "" = standard MI data; "train"/"test" = OOS fold under ${data_path}/oos/

cap mkdir "scratch"
cap log close req
log using "scratch/risk_equations_$coeffs.log", replace text name(req)

if "$repo_path" != "" cd "$repo_path"   // cd to repo root only if config.do set it; a bare cd "" goes to home on Mac/Unix
capture run "config.do"     // machine-specific paths: $data_path (git-ignored)

* OOS routing: read the OOS-fold imputed data when $sample is set (output is already keyed by
* $analysis, so OOS coefficients land in analyses/default/coefficients/). Empty = main model.
if "$sample" == "" {
	global re_indir ""
	global re_intag ""
}
else {
	global re_indir "oos/"
	global re_intag "_$sample"
}

**********
// Define functions

cap program drop save_coefs
program define save_coefs
	args mat

	mata: b`mat' = st_matrix("e(b_mi)")
	mata: rb`mat' = st_matrixrowstripe("e(b_mi)")
	mata: cb`mat' = st_matrixcolstripe("e(b_mi)")

	scalar ncol = colsof(e(b_mi))
	forval i = 1/`=ncol'{
		mata: cb`mat'[`i',1] = "`i'"
	}
	scalar drop ncol
	global Coeffs $Coeffs b`mat' rb`mat' cb`mat'

	// For parametric survival analysis - save functional form
	if ("e(ecmd_mi)" != "") {
		mata: fb`mat' = st_global("e(ecmd_mi)")
		global Coeffs $Coeffs fb`mat'
	}

	// For multinomial/ordered logt: save outcome levels
	capture mata: o`mat' = st_matrix("e(out)")
	if _rc == 0 {
		global Coeffs $Coeffs o`mat'
	}

	ereturn clear
end

cap program drop save_max
program define save_max
	// Ceiling over ALL records, censored included. Correct for a MAXIMUM: a patient censored while
	// still on treatment proves durations of at least that length occur, which is exactly what a
	// curtailment ceiling needs. Taking the max over observed ends only (the retired save_max_obs)
	// understated it - most for the heavily-censored quantities - and truncated the simulated tail.
	args mat

	qui summarize _t if _st == 1, meanonly
	di as txt "  ceiling `mat': " %7.1f r(max) " months over " %6.0f r(N) " records (censored included)"
	mata: max`mat' = st_numscalar("r(max)")
	global Coeffs $Coeffs max`mat'
end


cap program drop gen_mnr
program define gen_mnr
	// MNR_L1 arrives from the extraction as a RAW 6-level maintenance drug code (see
	// docs/refractory.md 2). Keep that as MNR_drug and rebuild MNR_L1 as the levels this
	// analysis models, per $MNR_L1 from outcomes/mnr_<coeffs>.do; anything unlisted falls to
	// 0 = 'other'. Mirrors gen_txr, and is what lets one extraction serve both the historical
	// window the OOS validation scores against and a current-paradigm window: the maintenance
	// mix stepped in 2020, so no single fixed list serves both (docs/refractory.md 7.4).
	// Re-entrant: rebuilds from MNR_drug if already called.
	cap confirm variable MNR_drug
	if _rc {
		cap confirm variable MNR_L1
		if _rc {
			di as err "gen_mnr: MNR_L1 is not in the MI data."
			di as err "         Rebuild MRDR Long (prep/data_extraction.do), then re-run"
			di as err "         prep/multiple_imputation.do - its keep list must name MNR_L1."
			exit 111
		}
		rename MNR_L1 MNR_drug
	}
	cap drop MNR_L1
	gen byte MNR_L1 = 0 if !mi(MNR_drug)
	foreach r of global MNR_L1 {
		replace MNR_L1 = `r' if MNR_drug == `r'
	}
	cap label values MNR_L1 MNR_label
end

cap program drop gen_txr
program define gen_txr
	// Build TXR_L1..L9 from the per-line regimen code lists in $TXR_L1..$TXR_L9 (declared by the
	// analysis's outcomes/txr_<coeffs>.do). Any regimen not listed for a line falls into 0 ('other').
	forval l = 1/9 {
		cap drop TXR_L`l'
		gen TXR_L`l' = 0 if Event0 == `l'0
		foreach r of global TXR_L`l' {
			replace TXR_L`l' = `r' if Event0 == `l'0 & Regimen == `r'
		}
		bysort ID (TXR_L`l'): replace TXR_L`l' = TXR_L`l'[_n-1] if TXR_L`l' == .
	}
end

cap program drop risk_equations
program define risk_equations

	// OS distribution
	global dOS "w"

	// TXD distribution
	global dTXD "w"

	// TFI distribution
	global dTFI = "lognormal"

	// TXR variables
	forval l = 1/9 {
		global TXR_L`l' ""
	}
	qui do "analyses/$analysis/outcomes/txr_$coeffs.do"
	gen_txr

	// MNR variable. Not every analysis declares a maintenance regimen list; without one every
	// regimen falls to 'other' and L1_MNR is skipped below, which is right for the $line 2
	// analyses where maintenance is never costed (docs/refractory.md 7.3).
	global MNR_L1 ""
	cap qui do "analyses/$analysis/outcomes/mnr_$coeffs.do"
	gen_mnr

	// Reset global
	global Coeffs

	// Risk equations are grouped BY OUTCOME TYPE (not by line). Each block below is a
	// self-contained estimation (its own stset/if), so grouping is purely organisational -
	// $Coeffs/matsave key by name, so order does not affect the saved coefficients. Within a
	// group the equations run in pathway order (diagnosis -> L1 -> ... -> LX). See the engine's
	// core/outcomes/sim_*.do for how each is consumed.

	***** LENALIDOMIDE-REFRACTORY, TREATMENT LINES (LENREFR_TX) *****
	di "Lenalidomide-refractory (treatment lines)"
	// LenRefr_Tx is a LATCHED state: once refractory to a lenalidomide treatment line the patient
	// stays refractory (cumulative-any, and downstream TXR/OS read only that latched flag). So the
	// only event to model is the 0 -> 1 FLIP, which can only happen to a NOT-YET-REFRACTORY patient.
	// The engine (sim_lenrefr.do) therefore draws only where LenRefr_Tx_in == 0:
	//     BCR in {5,6}  -> refractory        (definitional, no equation)
	//     BCR in {1-4}  -> Bernoulli(this logit)
	// and leaves already-refractory patients latched at 1. This logit is fitted on exactly that
	// population - not-yet-refractory (LenRefr_Tx_in == 0), responding (BCR 1-4), len line-start
	// rows - so prior state is NOT a covariate (it is 0 by construction here). Fitting it on the
	// full sample with LenRefr_Tx_in as a covariate, as an earlier draft did, mismatches the engine's
	// application population (the 4.4 "test with what the engine will have" rule).
	//
	// POOLED across lines with the line collapsed to L1 / L2 / L3+ (LENREFR_line): docs/refractory.md
	// 3.5 settles that shape - full line resolution and the L4+ split each add nothing (LR p > 0.95),
	// and line is the dominant predictor. Carries the pre-specified baseline set REGARDLESS of
	// significance (selecting on in-sample p-values does not replicate out of sample; the OOS
	// validation is the arbiter), plus i.BCR.
	//
	// FIT vs ENGINE gate: the fit conditions on Lenalidomide == 1 (the true drug binary, the clean
	// residual-arm definition). The engine has no drug binary, so sim_lenrefr applies this where the
	// DRAWN regimen is a lenalidomide code (analysis-declared, 'other' treated as non-len). That
	// asymmetry is deliberate and noted in docs/refractory.md 3.5 / 4.
	cap drop LENREFR_line
	gen byte LENREFR_line = min(Line, 3)
	// esampvaryok: the residual arm (BCR 1-4) is defined on imputed BCR, so its membership varies
	// across imputations - legitimate, and pooled the same way the ASCT equations are (docs note).
	mi estimate, esampvaryok: logit LineRefr Age Age2 Male i.ECOGcc i.RISS CM_CKD CM_CRD CM_PLM CM_DBT ///
		i.LENREFR_line i.BCR ///
		if(CStart == 1 & inlist(Event0, 10, 20, 30, 40, 50, 60, 70, 80, 90) & Lenalidomide == 1 & inlist(BCR, 1, 2, 3, 4) & LenRefr_Tx_in == 0)
	save_coefs LENREFR_TX
	mata: _matrix_list(bLENREFR_TX, rbLENREFR_TX, cbLENREFR_TX)

	// The len-regimen gate travels WITH the coefficients, not as a global, so a simulation cannot
	// silently run the refractory logit against a different regimen list from the one it was fitted
	// beside. sim_lenrefr.do reads it via get_lenrefr_regimens(); absent, the treatment arm is a
	// no-op rather than an error.
	if ("$LENREFR_regimens" != "") {
		mata: LENREFR_regimens = strtoreal(tokens(st_global("LENREFR_regimens")))
		global Coeffs $Coeffs LENREFR_regimens
	}

	***** LENALIDOMIDE-REFRACTORY, L1 MAINTENANCE (LENREFR_MNT) *****
	di "Lenalidomide-refractory (L1 maintenance)"
	// Replaces the deterministic TAIL RULE that sim_mnt_refr.do used to apply. That rule flagged a
	// patient when the SIMULATED tail (TFI_L1 - MND_L1) fell under 60 days - a correct DEFINITION but
	// a broken MECHANISM, because sim_tfi_l1 draws the gap truncated just above the maintenance
	// duration, so the tail piles up at the bound for LONG maintenance. It therefore selected
	// long-maintenance (good-prognosis) patients: refractory prevalence came out at 3.7% against a
	// registry 15.5% at L2, and the flagged group SURVIVED BETTER than the registry's despite the OS
	// penalty being applied correctly. Evidence: scratch/refractory/mnt_refr_logit.log.
	//
	// OUTCOME: the tail rule applied HONESTLY to observed data - progressed to L2 within 60 days of
	// maintenance ending, among patients whose maintenance ENDED. Still-on-maintenance patients are
	// undetermined and excluded; the engine runs everyone to completion, so the ended-maintenance
	// group is the congenial fitting sample. Base rate 28.6% (156/545).
	//   Note this is DATED progression only. mntrefr_reconcile.do's ~36% additionally counts patients
	//   flagged by the registry's cessation-reason column with no dated progression - which the engine
	//   has no field for and no mechanism could reproduce. 28.6% is the simulation-reachable target.
	//
	// SPECIFICATION, chosen by test not assumption (scratch/refractory/mnt_refr_logit.log):
	//   i.BCR_L1 beats the post-transplant response and the combined variable on BOTH AIC and AUC
	//     (327.3/0.708 vs 332.3/0.691), and within transplant patients the L1 response out-predicts
	//     BCR_SCT (AUC 0.681 vs 0.659). The induction response is simply the better predictor.
	//   POOLED across transplant with an SCT main effect: the response x transplant interaction is
	//     chi2(2) = 1.47, p = 0.48, so the effects do not differ - SCT shifts the level (-1.38,
	//     strongly protective), not the slopes. Same shape the MND split turned out to have.
	//   Response COLLAPSED to CR/VG/PR/poor. Forced: SD was a perfect predictor on 5 observations and
	//     was dropped from the uncollapsed fit, yet sim_bcr.do CAN draw it - so an uncollapsed model
	//     would leave simulated SD/PD patients falling through to the CR base level.
	// Coefficients are right-signed and monotone: worse response -> higher odds (CR base, VG +0.66,
	// PR +1.14, poor +1.17), older age higher, transplant protective. AUC 0.708.
	// ONE statement, not gen-then-replace: in mi-set data the multi-statement form left the variable
	// unregistered and mi estimate could not find it on m = 1. min(BCR_L1, 4) IS the CR/VG/PR/poor
	// collapse (4,5,6 all map to 4), and it matches the idiom used elsewhere in this file.
	cap drop MNTREFR_bcr
	qui gen byte MNTREFR_bcr = min(BCR_L1, 4) if !mi(BCR_L1)

	// Build the outcome from event dates: maintenance end = the earlier of a recorded cessation
	// (Event1 == 111) and L2 start (Event1 == 20); refractory = L2 within 60 days of that end.
	cap drop MNTREFR_r111
	cap drop MNTREFR_r20
	cap drop MNTREFR_cess
	cap drop MNTREFR_l2
	cap drop MNTREFR_end
	cap drop MNTREFR_ended
	cap drop MNTREFR
	qui gen double MNTREFR_r111 = Date1 if Event1 == 111        // this row's cessation date, if any
	qui gen double MNTREFR_r20  = Date1 if Event1 == 20         // this row's L2 start, if any
	qui egen double MNTREFR_cess = min(MNTREFR_r111), by(ID_BS) // patient's cessation date
	qui egen double MNTREFR_l2   = min(MNTREFR_r20),  by(ID_BS) // patient's L2 start
	qui gen double MNTREFR_end = min(MNTREFR_cess, MNTREFR_l2)  // maintenance ended at the earlier
	qui gen byte MNTREFR_ended = (!mi(MNTREFR_cess) | !mi(MNTREFR_l2))
	qui gen byte MNTREFR = (!mi(MNTREFR_l2) & (MNTREFR_l2 - MNTREFR_end) < 60) ///
		if MNTREFR_ended == 1

	// One row per patient: the fit is patient-level, so restrict to the maintenance-start row rather
	// than letting a patient contribute every row of their history.
	mi estimate, esampvaryok: logit MNTREFR Age Age2 Male i.ECOGcc i.RISS ///
		CM_CKD CM_CRD CM_PLM CM_DBT SCT i.MNTREFR_bcr ///
		if(Event1 == 110 & MNT == 1 & MNR_L1 == 1 & MNTREFR_ended == 1)
	save_coefs LENREFR_MNT
	mata: _matrix_list(bLENREFR_MNT, rbLENREFR_MNT, cbLENREFR_MNT)

	// COLLAPSED len-refractory covariate. Treatment- and maintenance-dose refractoriness carry the
	// same conditional OS hazard (HR ~1.64, test Tx = Mnt p = 0.97;
	// scratch/refractory/os_lenrefr_check.do), so ONE flag is fitted and the engine keeps one latched
	// state rather than two. Simpler, and the two-flag wiring did not apply the penalty correctly.
	cap drop LenRefr_any
	gen byte LenRefr_any = (LenRefr_Tx_in == 1 | LenRefr_Mnt_in == 1)

	***** OVERALL SURVIVAL *****
	di "Overall Survival"
	// Per-line OS family: one Weibull per pathway stage, window-censored via origin()/exit().
	// Same covariates throughout (Age Age2 Male i.ECOGcc i.RISS CM_* + stage BCR); stages differ
	// only by the stset window and BCR source. See core/outcomes/sim_os.do for the firing map.

	// OS_DN: DN -> L1S
	mi stset Date1 if(F_OS != 1), id(ID_BS) failure(Event1 == 104) origin(Event1 == 3) exit(Event1 == 10) scale(30.4375)
	mi estimate: streg Age Age2 Male i.ECOGcc i.RISS CM_CKD CM_CRD CM_PLM CM_DBT, d($dOS)
	save_coefs OS_DN
	mata: _matrix_list(bOS_DN, rbOS_DN, cbOS_DN)

	// OS_L1S: L1S -> L1E
	mi stset Date1 if(F_OS != 1), id(ID_BS) failure(Event1 == 104) origin(Event1 == 10) exit(Event1 == 11) scale(30.4375)
	mi estimate: streg Age Age2 Male i.ECOGcc i.RISS CM_CKD CM_CRD CM_PLM CM_DBT i.BCR_L1, d($dOS)
	save_coefs OS_L1S
	mata: _matrix_list(bOS_L1S, rbOS_L1S, cbOS_L1S)

	// OS_L1E_NoASCT: L1E -> L2S, No ASCT
	mi stset Date1 if(F_OS != 1 & SCT == 0), id(ID_BS) failure(Event1 == 104) origin(Event1 == 11) exit(Event1 == 20) scale(30.4375)
	mi estimate: streg Age Age2 Male i.ECOGcc i.RISS CM_CKD CM_CRD CM_PLM CM_DBT i.BCR_L1, d($dOS)
	save_coefs OS_L1E_NoASCT
	mata: _matrix_list(bOS_L1E_NoASCT, rbOS_L1E_NoASCT, cbOS_L1E_NoASCT)

	// OS_L1E_ASCT: L1E -> L2S, ASCT (BCR_SCT, 4 levels)
	mi stset Date1 if(F_OS != 1 & SCT == 1), id(ID_BS) failure(Event1 == 104) origin(Event1 == 11) exit(Event1 == 20) scale(30.4375)
	mi estimate: streg Age Age2 Male i.ECOGcc i.RISS CM_CKD CM_CRD CM_PLM CM_DBT i.BCR_SCT, d($dOS)
	save_coefs OS_L1E_ASCT
	mata: _matrix_list(bOS_L1E_ASCT, rbOS_L1E_ASCT, cbOS_L1E_ASCT)

	// OS_L2S: L2S -> L2E
	mi stset Date1 if(F_OS != 1), id(ID_BS) failure(Event1 == 104) origin(Event1 == 20) exit(Event1 == 21) scale(30.4375)
	mi estimate: streg Age Age2 Male i.ECOGcc i.RISS CM_CKD CM_CRD CM_PLM CM_DBT i.BCR_L2 LenRefr_any, d($dOS)
	save_coefs OS_L2S
	mata: _matrix_list(bOS_L2S, rbOS_L2S, cbOS_L2S)

	// OS_L2E: L2E -> L3S
	mi stset Date1 if(F_OS != 1), id(ID_BS) failure(Event1 == 104) origin(Event1 == 21) exit(Event1 == 30) scale(30.4375)
	mi estimate: streg Age Age2 Male i.ECOGcc i.RISS CM_CKD CM_CRD CM_PLM CM_DBT i.BCR_L2 LenRefr_any, d($dOS)
	save_coefs OS_L2E
	mata: _matrix_list(bOS_L2E, rbOS_L2E, cbOS_L2E)

	// OS_L3S: L3S -> L3E
	mi stset Date1 if(F_OS != 1), id(ID_BS) failure(Event1 == 104) origin(Event1 == 30) exit(Event1 == 31) scale(30.4375)
	mi estimate: streg Age Age2 Male i.ECOGcc i.RISS CM_CKD CM_CRD CM_PLM CM_DBT i.BCR_L3 LenRefr_any, d($dOS)
	save_coefs OS_L3S
	mata: _matrix_list(bOS_L3S, rbOS_L3S, cbOS_L3S)

	// OS_L3E: L3E -> L4S
	mi stset Date1 if(F_OS != 1), id(ID_BS) failure(Event1 == 104) origin(Event1 == 31) exit(Event1 == 40) scale(30.4375)
	mi estimate: streg Age Age2 Male i.ECOGcc i.RISS CM_CKD CM_CRD CM_PLM CM_DBT i.BCR_L3 LenRefr_any, d($dOS)
	save_coefs OS_L3E
	mata: _matrix_list(bOS_L3E, rbOS_L3E, cbOS_L3E)

	// OS_L4S: L4S -> L4E
	mi stset Date1 if(F_OS != 1), id(ID_BS) failure(Event1 == 104) origin(Event1 == 40) exit(Event1 == 41) scale(30.4375)
	mi estimate: streg Age Age2 Male i.ECOGcc i.RISS CM_CKD CM_CRD CM_PLM CM_DBT i.BCR_L4 LenRefr_any, d($dOS)
	save_coefs OS_L4S
	mata: _matrix_list(bOS_L4S, rbOS_L4S, cbOS_L4S)

	// OS_L4E: L4E -> L5S
	mi stset Date1 if(F_OS != 1), id(ID_BS) failure(Event1 == 104) origin(Event1 == 41) exit(Event1 == 50) scale(30.4375)
	mi estimate: streg Age Age2 Male i.ECOGcc i.RISS CM_CKD CM_CRD CM_PLM CM_DBT i.BCR_L4 LenRefr_any, d($dOS)
	save_coefs OS_L4E
	mata: _matrix_list(bOS_L4E, rbOS_L4E, cbOS_L4E)

	// OS_L5S: L5S -> L5E
	mi stset Date1 if(F_OS != 1), id(ID_BS) failure(Event1 == 104) origin(Event1 == 50) exit(Event1 == 51) scale(30.4375)
	mi estimate: streg Age Age2 Male i.ECOGcc i.RISS CM_CKD CM_CRD CM_PLM CM_DBT i.BCR_L5, d($dOS)
	save_coefs OS_L5S
	mata: _matrix_list(bOS_L5S, rbOS_L5S, cbOS_L5S)

	// OS_L5E: L5E -> L6S
	mi stset Date1 if(F_OS != 1), id(ID_BS) failure(Event1 == 104) origin(Event1 == 51) exit(Event1 == 60) scale(30.4375)
	mi estimate: streg Age Age2 Male i.ECOGcc i.RISS CM_CKD CM_CRD CM_PLM CM_DBT i.BCR_L5, d($dOS)
	save_coefs OS_L5E
	mata: _matrix_list(bOS_L5E, rbOS_L5E, cbOS_L5E)

	// OS_L6plus: L6S onward (single conditional model for the sparse deep tail; running BCR)
	mi stset Date1 if(F_OS != 1), id(ID_BS) failure(Event1 == 104) origin(Event1 == 60) scale(30.4375)
	mi estimate: streg Age Age2 Male i.ECOGcc i.RISS CM_CKD CM_CRD CM_PLM CM_DBT i.BCR, d($dOS)
	save_coefs OS_L6plus
	mata: _matrix_list(bOS_L6plus, rbOS_L6plus, cbOS_L6plus)

	***** ASCT (TRANSPLANT DECISION) *****
	di "ASCT"
	// Transplant-eligibility logits: at diagnosis (DN_SCT) and at L1 end (L1_SCT).
	// Comorbidities: 4 individual flags (CKD/cardiac/lung/diabetes).

	// DN
	mi estimate: logit SCT Age Age2 Male i.ECOGcc i.RISS Age70 Age75 CM_CKD CM_CRD CM_PLM CM_DBT if(Event0 == 3)
	save_coefs DN_SCT
	mata: _matrix_list(bDN_SCT, rbDN_SCT, cbDN_SCT)

	// L1
	mi estimate, esampvaryok: logit SCT Age Age2 Male i.ECOGcc i.RISS i.BCR Age70 Age75 CM_CKD CM_CRD CM_PLM CM_DBT if(Event1 == 11 & TXR_L1 != 7 & BCR != 6)
	save_coefs L1_SCT
	mata: _matrix_list(bL1_SCT, rbL1_SCT, cbL1_SCT)

	***** TREATMENT REGIMEN (TXR) *****
	di "Treatment Regimen"
	// Regimen choice is availability-driven; L1 carries Age Age2 SCT (transplant gates some
	// regimens), L2+ carry Age only - Male/ECOG/RISS/prior-BCR carried no signal (see
	// scratch/txr_predictor_check.do). Engine design mPat in core/outcomes/sim_txr.do must match.
	// If a line has only 'other' (r(r) == 1) an empty o<L>_TXR = 0 is stored instead.

	// L1
	qui tab TXR_L1 // Check if modelling specific regimens
	if `r(r)' > 1 {
		mi estimate: mlogit TXR_L1 Age Age2 SCT if(Event0 == 10 & yofd(Date0) >= $min_year & yofd(Date0) <= $max_year), baseoutcome(0)
		save_coefs L1_TXR
		mata: _matrix_list(bL1_TXR, rbL1_TXR, cbL1_TXR)
	}
	else { // If all patients receive 'other'
		mata: oL1_TXR = 0
		global Coeffs $Coeffs oL1_TXR
	}

	// L2
	qui tab TXR_L2 // Check if modelling specific regimens
	if `r(r)' > 1 {
		mi estimate: mlogit TXR_L2 Age Age2 LenRefr_any if(Event0 == 20 & yofd(Date0) >= $min_year & yofd(Date0) <= $max_year), baseoutcome(0)
		save_coefs L2_TXR
		mata: _matrix_list(bL2_TXR, rbL2_TXR, cbL2_TXR)
	}
	else { // If all patients receive 'other'
		mata: oL2_TXR = 0
		global Coeffs $Coeffs oL2_TXR
	}

	// L3
	qui tab TXR_L3 // Check if modelling specific regimens
	if `r(r)' > 1 {
		mi estimate: mlogit TXR_L3 Age Age2 LenRefr_any if(Event0 == 30 & yofd(Date0) >= $min_year & yofd(Date0) <= $max_year), baseoutcome(0)
		save_coefs L3_TXR
		mata: _matrix_list(bL3_TXR, rbL3_TXR, cbL3_TXR)
	}
	else { // If all patients receive 'other'
		mata: oL3_TXR = 0
		global Coeffs $Coeffs oL3_TXR
	}

	// L4
	qui tab TXR_L4 // Check if modelling specific regimens
	if `r(r)' > 1 {
		mi estimate: mlogit TXR_L4 Age Age2 LenRefr_any if(Event0 == 40 & yofd(Date0) >= $min_year & yofd(Date0) <= $max_year), baseoutcome(0)
		save_coefs L4_TXR
		mata: _matrix_list(bL4_TXR, rbL4_TXR, cbL4_TXR)
	}
	else { // If all patients receive 'other'
		mata: oL4_TXR = 0
		global Coeffs $Coeffs oL4_TXR
	}

	// L5
	qui tab TXR_L5 // Check if modelling specific regimens
	if `r(r)' > 1 {
		mi estimate: mlogit TXR_L5 Age Age2 if(Event0 == 50 & yofd(Date0) >= $min_year & yofd(Date0) <= $max_year), baseoutcome(0)
		save_coefs L5_TXR
		mata: _matrix_list(bL5_TXR, rbL5_TXR, cbL5_TXR)
	}
	else { // If all patients receive 'other'
		mata: oL5_TXR = 0
		global Coeffs $Coeffs oL5_TXR
	}

	***** BEST CLINICAL RESPONSE (BCR) *****
	di "Best Clinical Response"
	// Ordered logit for response at each line; each line conditions on the prior line's response
	// (and regimen). SCT_BCR is the post-transplant response (conditional on BCR_L1).

	// L1
	mi estimate: ologit BCR Age Age2 Male i.ECOGcc i.RISS SCT i.TXR_L1 if(Event0 == 10)
	save_coefs L1_BCR
	mata: _matrix_list(bL1_BCR, rbL1_BCR, cbL1_BCR)

	// ASCT (post-transplant response)
	mi estimate, esampvaryok: ologit BCR Age Age2 Male i.ECOGcc i.RISS i.BCR_L1 if(Event0 == 100 & BCR_L1 != 6)
	save_coefs SCT_BCR
	mata: _matrix_list(bSCT_BCR, rbSCT_BCR, cbSCT_BCR)

	// L2
	mi estimate: ologit BCR Age Age2 Male i.ECOGcc i.RISS i.BCR_L1 i.BCR_SCT i.TXR_L2 if(Event0 == 20)
	save_coefs L2_BCR
	mata: _matrix_list(bL2_BCR, rbL2_BCR, cbL2_BCR)

	// L3
	mi estimate: ologit BCR Age Age2 Male i.ECOGcc i.RISS i.BCR_L2 i.TXR_L3 if(Event0 == 30)
	save_coefs L3_BCR
	mata: _matrix_list(bL3_BCR, rbL3_BCR, cbL3_BCR)

	// L4
	mi estimate: ologit BCR Age Age2 Male i.ECOGcc i.RISS i.BCR_L3 i.TXR_L4 if(Event0 == 40)
	save_coefs L4_BCR
	mata: _matrix_list(bL4_BCR, rbL4_BCR, cbL4_BCR)

	// L5
	mi estimate: ologit BCR Age Age2 Male i.ECOGcc i.RISS i.BCR_L4 if(Event0 == 50)
	save_coefs L5_BCR
	mata: _matrix_list(bL5_BCR, rbL5_BCR, cbL5_BCR)

	// LX (lines 6-9)
	mi estimate: ologit BCR Age Age2 Male i.ECOGcc i.RISS i.pBCR if(Event0 == 60 | Event0 == 70 | Event0 == 80 | Event0 == 90)
	save_coefs LX_BCR
	mata: _matrix_list(bLX_BCR, rbLX_BCR, cbLX_BCR)

	// L6 (superseded by LX; kept for reference)
/*
	mi estimate: ologit BCR Age Age2 Male i.ECOGcc i.RISS i.pBCR if(Event0 == 60)
	save_coefs L6_BCR
	mata: _matrix_list(bL6_BCR, rbL6_BCR, cbL6_BCR)
*/

	***** TREATMENT DURATION (TXD) *****
	di "Treatment Duration"
	// Time on treatment within a line (Weibull). L1 splits by transplant: ASCT fixed-duration is a
	// 3-spline fit (cut-offs 60/120 mo), NoASCT is a single fit, and Continuous (Rd) is separate.
	// L2+ are single fits conditional on that line's response and regimen.

	// L1 - Fixed w/ ASCT (3 splines)
	// Cut off 1 - 60 for optimal fit with 3 splines, 90 for 2 splines
	scalar L1_TXD_ASCT_C1 = 60
	mata: L1_TXD_ASCT_C1 = st_numscalar("L1_TXD_ASCT_C1")

	// Cut off 2 - 120 for optimal fit with 3 splines
	scalar L1_TXD_ASCT_C2 = 120
	mata: L1_TXD_ASCT_C2 = st_numscalar("L1_TXD_ASCT_C2")

	global Coeffs $Coeffs L1_TXD_ASCT_C1 L1_TXD_ASCT_C2

	// Spline 1
	gen S1_exit = Date1 + `=L1_TXD_ASCT_C1'
	mi stset Date1 if(TXR_L1 != 7 & SCT == 1), id(ID_BS) failure(Event1 == 11) origin(Event1 == 10) exit(S1_exit) scale(30.4375)
	mi estimate: streg Age Age2 Male i.ECOGcc i.RISS i.TXR_L1, d($dTXD)
	save_coefs L1_TXD_ASCT_S1
	mata: _matrix_list(bL1_TXD_ASCT_S1, rbL1_TXD_ASCT_S1, cbL1_TXD_ASCT_S1)

	// Spline 2
	gen S2_enter = S1_exit + 1
	gen S2_exit = Date1 + `=L1_TXD_ASCT_C2'
	mi stset Date1 if(TXR_L1 != 7 & SCT == 1), id(ID_BS) failure(Event1 == 11) origin(Event1 == 10) enter(S2_enter) exit(S2_exit) scale(30.4375)
	mi estimate: streg Age Age2 Male i.ECOGcc i.RISS i.TXR_L1, d($dTXD)
	save_coefs L1_TXD_ASCT_S2
	mata: _matrix_list(bL1_TXD_ASCT_S2, rbL1_TXD_ASCT_S2, cbL1_TXD_ASCT_S2)

	// Spline 3
	gen S3_enter = S2_exit + 1
	mi stset Date1 if(TXR_L1 != 7 & SCT == 1), id(ID_BS) failure(Event1 == 11) origin(Event1 == 10) enter(S3_enter) scale(30.4375)
	mi estimate: streg Age Age2 Male i.ECOGcc i.RISS i.TXR_L1, d($dTXD)
	save_coefs L1_TXD_ASCT_S3
	drop S1_exit S2_enter S2_exit S3_enter
	mata: _matrix_list(bL1_TXD_ASCT_S3, rbL1_TXD_ASCT_S3, cbL1_TXD_ASCT_S3)

	// L1 - Fixed w/o ASCT
	mi stset Date1 if(TXR_L1 != 7 & SCT == 0), id(ID_BS) failure(Event1 == 11) origin(Event1 == 10) scale(30.4375)
	mi estimate: streg Age Age2 Male i.ECOGcc i.RISS i.BCR i.TXR_L1, d($dTXD)
	save_coefs L1_TXD_NoASCT
	mata: _matrix_list(bL1_TXD_NoASCT, rbL1_TXD_NoASCT, cbL1_TXD_NoASCT)

	// L1 - Continuous (Rd)
	count if(TXR_L1 == 7)
	if r(N) > 0 { // If Rd is included
		mi stset Date1 if(TXR_L1 == 7), id(ID_BS) failure(Event1 == 11) origin(Event1 == 10) scale(30.4375)
		save_max L1_TXD_Cont
		mi estimate: streg Age Age2 Male i.ECOGcc i.RISS i.BCR, d($dTXD)
		save_coefs L1_TXD_Cont
	}
	else { // If Rd is not included (need empty matrices)
		mata: bL1_TXD_Cont = .
		mata: rbL1_TXD_Cont = .
		mata: cbL1_TXD_Cont = .
		mata: fbL1_TXD_Cont = .
		mata: maxL1_TXD_Cont = .
		global Coeffs $Coeffs bL1_TXD_Cont rbL1_TXD_Cont cbL1_TXD_Cont fbL1_TXD_Cont maxL1_TXD_Cont
	}
	*mata: _matrix_list(bL1_TXD_Cont, rbL1_TXD_Cont, cbL1_TXD_Cont)

	// L2
	mi stset Date1, id(ID_BS) failure(Event1 == 21) origin(Event1 == 20) scale(30.4375)
	save_max L2_TXD
	mi estimate: streg Age Age2 Male i.ECOGcc i.RISS i.BCR_L2 i.TXR_L2, d($dTXD)
	save_coefs L2_TXD
	mata: _matrix_list(bL2_TXD, rbL2_TXD, cbL2_TXD)

	// L3
	mi stset Date1, id(ID_BS) failure(Event1 == 31) origin(Event1 == 30) scale(30.4375)
	save_max L3_TXD
	mi estimate: streg Age Age2 Male i.ECOGcc i.RISS i.BCR_L3 i.TXR_L3, d($dTXD)
	save_coefs L3_TXD
	mata: _matrix_list(bL3_TXD, rbL3_TXD, cbL3_TXD)

	// L4
	mi stset Date1, id(ID_BS) failure(Event1 == 41) origin(Event1 == 40) scale(30.4375)
	save_max L4_TXD
	mi estimate: streg Age Age2 Male i.ECOGcc i.RISS i.BCR_L4 i.TXR_L4, d($dTXD)
	save_coefs L4_TXD
	mata: _matrix_list(bL4_TXD, rbL4_TXD, cbL4_TXD)

	// L5
	mi stset Date1, id(ID_BS) failure(Event1 == 51) origin(Event1 == 50) scale(30.4375)
	save_max L5_TXD
	mi estimate: streg Age Age2 Male i.ECOGcc i.RISS i.BCR_L5, d($dTXD)
	save_coefs L5_TXD
	mata: _matrix_list(bL5_TXD, rbL5_TXD, cbL5_TXD)

	// LX (lines 6-9)
	qui cap drop fail
	qui gen fail = 1 if(Event1 == 61 | Event1 == 71 | Event1 == 81 | Event1 == 91)
	qui replace fail = 0 if(Event1 == 60 | Event1 == 70 | Event1 == 80 | Event1 == 90)
	qui bysort ID_BS (Date0): replace fail = fail[_n-1] if(fail == .)
	mi stset Date1, id(ID_BS) failure(fail) origin(Event1 == 60) exit(time .) scale(30.4375)
	save_max LX_TXD
	mi estimate: streg Age Age2 Male i.ECOGcc i.RISS i.BCR, d($dTXD)
	save_coefs LX_TXD
	mata: _matrix_list(bLX_TXD, rbLX_TXD, cbLX_TXD)

	// L6 (superseded by LX; kept for reference)
/*
	mi stset Date1, id(ID_BS) failure(Event1 == 61) origin(Event1 == 60) scale(30.4375)
	save_max L6_TXD
	mi estimate: streg Age Age2 Male i.ECOGcc i.RISS i.BCR, d($dTXD)
	save_coefs L6_TXD
	mata: _matrix_list(bL6_TXD, rbL6_TXD, cbL6_TXD)
*/

	***** MAINTENANCE (MNT) *****
	di "Maintenance"
	// Maintenance-therapy logit after L1 induction, split by transplant status.

	// ASCT
	mi estimate: logit MNT Age Age2 Male i.ECOGcc i.RISS i.TXR_L1 i.BCR_SCT if(Event1 == 11 & SCT == 1)
	save_coefs MNT_ASCT
	mata: _matrix_list(bMNT_ASCT, rbMNT_ASCT, cbMNT_ASCT)

	// No ASCT
	mi estimate: logit MNT Age Age2 Male i.ECOGcc i.RISS i.TXR_L1 i.BCR_L1 if(Event1 == 11 & SCT == 0)
	save_coefs MNT_NoASCT
	mata: _matrix_list(bMNT_NoASCT, rbMNT_NoASCT, cbMNT_NoASCT)

	***** MAINTENANCE REGIMEN AND DURATION (MNR / MND) *****
	di "Maintenance regimen and duration"
	// The maintenance analogue of TXR_L1 / TXD_L1, sat here because both condition on the MNT
	// logit above: MNT decides WHETHER, L1_MNR WHICH, L1_MND HOW LONG. Both fit at the same
	// Event1 == 11 row as MNT, and only among MNT == 1. Consumed only by cost_tx_mnt, so the
	// out-of-sample validation fits them but never uses them. See docs/refractory.md 7.4.

	// Which regimen. Lenalidomide and thalidomide only (docs/refractory.md 7.4) - lenalidomide
	// is the base, thalidomide the alternative, so this is effectively a logit and the engine
	// never produces an 'other' maintenance regimen. Year-windowed exactly as TXR_L1 is, so the
	// mix reflects the era; the list comes from outcomes/mnr_$coeffs.do ("1 5"). If a window
	// leaves only one regimen (r(r) == 1, e.g. thalidomide empty in a modern window) then
	// oL1_MNR = 1 is stored and sim_mnr assigns every maintenance patient lenalidomide.
	qui tab MNR_L1 if(Event1 == 11 & MNT == 1 & inlist(MNR_L1, 1, 5))
	if `r(r)' > 1 {
		mi estimate: mlogit MNR_L1 Age Age2 Male i.ECOGcc i.RISS SCT ///
			if(Event1 == 11 & MNT == 1 & inlist(MNR_L1, 1, 5) & yofd(Date0) >= $min_year & yofd(Date0) <= $max_year), baseoutcome(1)
		save_coefs L1_MNR
		mata: _matrix_list(bL1_MNR, rbL1_MNR, cbL1_MNR)
	}
	else {
		mata: oL1_MNR = 1
		global Coeffs $Coeffs oL1_MNR
	}

	// How long: maintenance DURATION via parametric survival on the maintenance EVENTS in the
	// skeleton, exactly as L1_TFI is fitted (docs/refractory.md 4.4). The extraction keeps the L1
	// maintenance start (110) and end (111) events, so this is a normal stset:
	//     origin(Event1 == 110)          maintenance start
	//     failure(Event1 == 20 111)      maintenance end - the recorded end (111), or L2 start (20)
	//
	// SPLIT BY REGIMEN, POOLED ACROSS TRANSPLANT, NO GAP TERM. The two drugs are different processes
	// with non-overlapping ancillaries (thalidomide a fixed course, lenalidomide running to
	// progression); the two transplant arms are near-identical, so SCT is a covariate instead. BCR
	// drops out as a consequence - the arms key on different response variables.
	//
	// The gap covariate is gone because ln(TFI_L1) is missing without an observed L2, which
	// restricted this fit to a quarter of the maintenance population. Dropping it alone is NOT
	// enough: it leaves duration and gap independent and the billing cap does the damage instead.
	// sim_tfi_l1.do supplies the dependence by drawing the gap truncated below at the maintenance.
	//
	// THALIDOMIDE IS CENSORED AT 18 MONTHS - a documented judgement that recorded ends beyond that
	// are stale records, not a measurement. Censored rather than excluded: excluding leaves a
	// truncated sample and an untruncated family fitted to it undershoots. generate_benchmarks.do
	// carries the same exit() so both sides measure the same quantity.
	//
	// Evidence for all of the above, and the five specifications rejected on the way, are in
	// scratch/maintenance/_notes.md. Conclusions in docs/refractory.md 5(7) and 7.
	//
	// TRAP: exit() is in the TIME VARIABLE's scale, not analysis time, and `origin' is unavailable
	// inside it when origin() is an event condition. It must reference the origin date in days.
	capture drop MND_origin
	qui gen double _mxMNDo = Date1 if Event1 == 110
	qui egen double MND_origin = min(_mxMNDo), by(ID_BS)
	qui drop _mxMNDo

	// Lenalidomide, POOLED across transplant, baseline covariates only.
	//
	// NO response term and NO transplant split. Response depth was tested exhaustively once BCR_SCT
	// was correctly imputed (it had been mis-coded, which is what created an apparent signal): the
	// post-transplant response BCR_SCT, the L1 response BCR_L1 collapsed AND full 6-level, the two
	// combined, and a transplant interaction - none predicts maintenance duration (all LR p > 0.4,
	// every AIC worse than baseline). Comorbidities likewise (p = 0.23). The p = 0.0021 that once
	// justified splitting the arm to carry BCR_SCT was purely the coding artefact.
	// Evidence: scratch/mnd_bcr_pooled.log. Response still drives TFI and OS at full 6 levels; it
	// simply does not add to the maintenance-DURATION draw within the selected maintenance cohort.
	//
	// CEILING and FAMILY, both settled: $dTFI is lognormal, the 157.96-month ceiling is a genuine
	// recorded cessation (RMST(158) is the only unflagged KM estimate), and the drawn duration
	// reproduces the registry restricted mean at every horizon. Do not lower the ceiling.
	// Evidence: scratch/mnd_dist.log, scratch/mnd_failtype.log, ../scratch/maintenance/_notes.md.
	mi stset Date1 if(MNT == 1 & MNR_L1 == 1), ///
		id(ID_BS) failure(Event1 == 20 111) origin(Event1 == 110) scale(30.4375)
	save_max L1_MND_LEN
	mi estimate: streg Age Age2 Male i.ECOGcc i.RISS SCT, d($dTFI)
	save_coefs L1_MND_LEN
	mata: _matrix_list(bL1_MND_LEN, rbL1_MND_LEN, cbL1_MND_LEN)

	// Thalidomide, censored at 18 months. Pooled across transplant, no BCR: n is roughly 311 against
	// lenalidomide's 1,020, too small to add a factor. Thalidomide maintenance also ended with the
	// 2020 PBS listing, so it matters to historical validation rather than any current analysis.
	mi stset Date1 if(MNT == 1 & MNR_L1 == 5), ///
		id(ID_BS) failure(Event1 == 20 111) origin(Event1 == 110) scale(30.4375) ///
		exit(time MND_origin + `=18 * 30.4375')
	// The engine caps the thalidomide draw at 18 months, so the next line overwrites whatever the
	// ceiling computes; the save_max call is kept only so the refit log reports it like the others.
	save_max L1_MND_THAL
	mata: maxL1_MND_THAL = 18
	mi estimate: streg Age Age2 Male i.ECOGcc i.RISS SCT, d($dTFI)
	save_coefs L1_MND_THAL
	mata: _matrix_list(bL1_MND_THAL, rbL1_MND_THAL, cbL1_MND_THAL)

	***** TREATMENT-FREE INTERVAL (TFI) *****
	di "Treatment-free Interval"
	// Gap from the end of one line's treatment to the start of the next (log-normal). DN_TFI is the
	// diagnosis-to-L1 interval; L1 splits by transplant (with maintenance MNT); L2+ are single fits.

	// DN
	mi stset Date1, id(ID_BS) failure(Event1 == 10) origin(Event1 == 3) scale(30.4375)
	mi estimate: streg Age Age2 Male i.ECOGcc i.RISS SCT, d($dTFI)
	save_coefs DN_TFI
	mata: _matrix_list(bDN_TFI, rbDN_TFI, cbDN_TFI)

	// L1 - ASCT
	// inrange(BCR_SCT, 1, 4), NOT "BCR_SCT != 0". Since multiple_imputation.do stopped coding
	// transplanted-but-unrecorded patients as 0, they are MISSING - and in Stata `. != 0' is TRUE,
	// so the old condition would have quietly let them into the risk set. They would still drop out
	// of the estimation listwise, but not before affecting _st == 1, and therefore the ceiling.
	mi stset Date1 if(SCT == 1 & inrange(BCR_SCT, 1, 4)), id(ID_BS) failure(Event1 == 20) origin(Event1 == 11) scale(30.4375)
	save_max L1_TFI_ASCT
	mi estimate: streg Age Age2 Male i.ECOGcc i.RISS MNT i.BCR_SCT, d($dTFI)
	save_coefs L1_TFI_ASCT
	mata: _matrix_list(bL1_TFI_ASCT, rbL1_TFI_ASCT, cbL1_TFI_ASCT)

	// L1 - No ASCT
	mi stset Date1 if(SCT == 0), id(ID_BS) failure(Event1 == 20) origin(Event1 == 11) scale(30.4375)
	save_max L1_TFI_NoASCT
	mi estimate: streg Age Age2 Male i.ECOGcc i.RISS MNT i.BCR_L1, d($dTFI)
	save_coefs L1_TFI_NoASCT
	mata: _matrix_list(bL1_TFI_NoASCT, rbL1_TFI_NoASCT, cbL1_TFI_NoASCT)

	// L2
	mi stset Date1, id(ID_BS) failure(Event1 == 30) origin(Event1 == 21) scale(30.4375)
	save_max L2_TFI
	mi estimate: streg Age Age2 Male i.ECOGcc i.RISS i.BCR_L2, d($dTFI)
	save_coefs L2_TFI
	mata: _matrix_list(bL2_TFI, rbL2_TFI, cbL2_TFI)

	// L3
	mi stset Date1, id(ID_BS) failure(Event1 == 40) origin(Event1 == 31) scale(30.4375)
	save_max L3_TFI
	mi estimate: streg Age Age2 Male i.ECOGcc i.RISS i.BCR_L3, d($dTFI)
	save_coefs L3_TFI
	mata: _matrix_list(bL3_TFI, rbL3_TFI, cbL3_TFI)

	// L4
	mi stset Date1, id(ID_BS) failure(Event1 == 50) origin(Event1 == 41) scale(30.4375)
	save_max L4_TFI
	mi estimate: streg Age Age2 Male i.ECOGcc i.RISS i.BCR_L4, d($dTFI)
	save_coefs L4_TFI
	mata: _matrix_list(bL4_TFI, rbL4_TFI, cbL4_TFI)

	// L5
	mi stset Date1, id(ID_BS) failure(Event1 == 60) origin(Event1 == 51) scale(30.4375)
	save_max L5_TFI
	mi estimate: streg Age Age2 Male i.ECOGcc i.RISS i.BCR_L5, d($dTFI)
	save_coefs L5_TFI
	mata: _matrix_list(bL5_TFI, rbL5_TFI, cbL5_TFI)

	// LX (lines 6-9)
	qui cap drop fail
	qui gen fail = 1 if(Event1 == 70 | Event1 == 80 | Event1 == 90)
	qui replace fail = 0 if(Event1 == 61 | Event1 == 71 | Event1 == 81)
	qui bysort ID_BS (Date0): replace fail = fail[_n-1] if(fail == .)
	mi stset Date1, id(ID_BS) failure(fail) origin(Event1 == 61) exit(time .) scale(30.4375)
	save_max LX_TFI
	mi estimate: streg Age Age2 Male i.ECOGcc i.RISS i.BCR, d($dTFI)
	save_coefs LX_TFI
	mata: _matrix_list(bLX_TFI, rbLX_TFI, cbLX_TFI)

	// L6 (superseded by LX; kept for reference)
/*
	mi stset Date1, id(ID_BS) failure(Event1 == 70) origin(Event1 == 61) scale(30.4375)
	save_max L6_TFI
	mi estimate: streg Age Age2 Male i.ECOGcc i.RISS i.BCR, d($dTFI)
	save_coefs L6_TFI
	mata: _matrix_list(bL6_TFI, rbL6_TFI, cbL6_TFI)
*/

end

**********
// Execute based on arguments

if("$boot" == "0") {

	// Load data
	use "${data_path}/${re_indir}MRDR Long MI${re_intag}.dta", replace

	// Execute function
	risk_equations

	// Save matrices
	mata: mata matsave "analyses/$analysis/coefficients/coefficients_$coeffs" $Coeffs, replace
}
else if("$boot" == "1") {

	forval b = $min_bs/$max_bs {

		// Load data
		use "${data_path}/${re_indir}bootstrap/MRDR Long MI${re_intag} B`b'.dta", replace

		// Shift dates so bootstrap samples of the same ID_BS do not have the same Date0/Date1
		qui bysort ID Event0 (ID_BS): gen tag0 = _n if Date0[_n-1] == Date0 & Date0 != .
		qui bysort ID Event1 (ID_BS): gen tag1 = _n if Date1[_n-1] == Date1 & Date1 != .
		qui replace Date0 = Date0 + tag0 - 1 if tag0 != .
		qui replace Date1 = Date1 + tag1 - 1 if tag1 != .
		qui drop tag0 tag1

		// Execute function
		risk_equations

		// Save matrices
		mata: mata matsave "analyses/$analysis/coefficients/bootstrap/coefficients_${coeffs}_B`b'" $Coeffs, replace
	}
}

cap log close req
