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
//log using "scratch/risk_equations_$coeffs.log", replace text name(req)

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

	// HEALTH CHECK. An unidentified block does not error - OS_L5S shipped SEs of 462 and hazard
	// ratios of 50-140 in every generation until four bootstrap replicates happened to die. This
	// flags it at the point of fit. Silent on pass; one line per problem.
	capture {
		tempname V
		capture matrix `V' = e(V_mi)
		if _rc matrix `V' = e(V)
		local mx = 0
		forvalues i = 1/`=colsof(`V')' {
			local v = `V'[`i', `i']
			if `v' < . & `v' > 0 & sqrt(`v') > `mx' local mx = sqrt(`v')
		}
		if `mx' > 50 ///
			di as error "  CHECK `mat': max SE " %10.1f `mx' " - block likely unidentified (see docs/refractory.md 3.5)."
		local npar = colsof(`V')
		if !missing(e(N_fail)) & `npar' > 0 {
			if e(N_fail) / `npar' < 5 ///
				di as error "  CHECK `mat': " e(N_fail) " events for " `npar' " parameters (" %4.1f e(N_fail)/`npar' " each)."
		}
	}

	ereturn clear
end

cap program drop save_max
program define save_max
	// Ceiling over ALL records, censored included: a patient censored while still on treatment proves
	// durations of at least that length occur. Observed ends only understates it.
	args mat

	qui summarize _t if _st == 1, meanonly
	mata: max`mat' = st_numscalar("r(max)")
	global Coeffs $Coeffs max`mat'
end

cap program drop gen_mnr
program define gen_mnr
	// Keep the raw 6-level drug code as MNR_drug and rebuild MNR_L1 as the levels this analysis
	// models per $MNR_L1; anything unlisted falls to 0 = 'other'. Mirrors gen_txr. Re-entrant.
	// docs/refractory.md 2 / 7.4.
	cap confirm variable MNR_drug
	if _rc {
		cap confirm variable MNR_L1
		if _rc {
			di as err "gen_mnr: MNR_L1 not in the MI data - rebuild the extraction and re-run the MI."
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

	// Regimen lists: one file per analysis declares the per-line treatment lists, the L1 maintenance
	// list and the len-refractory gate. Reset first so a stale value cannot leak in from an earlier
	// run in the same session. An empty $MNR_L1 is not supported and errors below.
	forval l = 1/9 {
		global TXR_L`l' ""
	}
	global MNR_L1 ""
	qui do "analyses/$analysis/outcomes/txr_$coeffs.do"
	gen_txr
	gen_mnr

	// Reset global
	global Coeffs

	// Grouped by outcome type, each block a self-contained estimation. $Coeffs/matsave key by name,
	// so order does not affect the saved coefficients. See core/outcomes/sim_*.do for consumption.

	***** LENALIDOMIDE-REFRACTORY, TREATMENT LINES (LENREFR_TX) *****
	// refr_len_tx_in is a latched 0 -> 1 state, so only the flip is modelled, and only where the patient is
	// not yet refractory: BCR in {5,6} is definitional, BCR 1-4 draws this logit. Fitted on exactly that
	// population, so prior state is not a covariate. Line collapsed to L1/L2/L3+; baseline set carried
	// regardless of significance. The fit gates on Lenalidomide == 1, the engine on the drawn regimen -
	// deliberate. docs/refractory.md 3.5 / 4.
	cap drop refr_len_grp
	gen byte refr_len_grp = min(Line, 3)
	// esampvaryok: the residual arm is defined on imputed BCR, so its membership varies across m.
	mi estimate, esampvaryok: logit refr_line Age Age2 Male i.ECOGcc i.RISS CM_CKD CM_CRD CM_PLM CM_DBT ///
		i.refr_len_grp i.BCR ///
		if(CStart == 1 & inlist(Event0, 10, 20, 30, 40, 50, 60, 70, 80, 90) & Lenalidomide == 1 & inlist(BCR, 1, 2, 3, 4) & refr_len_tx_in == 0)
	save_coefs LENREFR_TX
	mata: _matrix_list(bLENREFR_TX, rbLENREFR_TX, cbLENREFR_TX)

	// The gate travels WITH the coefficients so a run cannot use a different regimen list from the one
	// it was fitted beside. sim_lenrefr.do reads it via get_lenrefr_regimens(); absent, it is a no-op.
	if ("$LENREFR_regimens" != "") {
		mata: LENREFR_regimens = strtoreal(tokens(st_global("LENREFR_regimens")))
		global Coeffs $Coeffs LENREFR_regimens
	}

	// P(lenalidomide | regimen NOT modelled), per line - the 'other'-bucket gate probability. Unlisted
	// regimens fall to TXR code 0, which the gate would otherwise treat as non-len; that bucket carries
	// a large lenalidomide share, skewed older and more refractory. Estimated here rather than
	// hard-coded because it depends on $TXR_L{line}, and travels with the coefficients. Does not touch
	// which regimen is drawn or for how long. Evidence: docs/refractory.md 5(6).
	cap drop refr_len_oth
	// Year-windowed as the TXR mlogits are: the bucket is an era composition, not a fixed quantity.
	// Below `minoth' records in the window the all-years value is used, rather than a p built on a
	// handful of patients.
	local minoth = 50
	mata: LENREFR_pother = J(1, 9, 0)
	forvalues l = 1/9 {
		qui gen byte refr_len_oth = 1 if Event0 == `l'0 & CStart == 1
		foreach r of global TXR_L`l' {
			qui replace refr_len_oth = 0 if Event0 == `l'0 & CStart == 1 & Regimen == `r'
		}

		// All-years estimate, kept as the fallback
		qui count if refr_len_oth == 1
		local nall = r(N)
		local pall = 0
		if `nall' > 0 {
			qui count if refr_len_oth == 1 & Lenalidomide == 1
			local pall = r(N)/`nall'
		}

		// Windowed estimate
		qui count if refr_len_oth == 1 & yofd(Date0) >= $min_year & yofd(Date0) <= $max_year
		local noth = r(N)
		local p = 0
		if `noth' > 0 {
			qui count if refr_len_oth == 1 & Lenalidomide == 1 & yofd(Date0) >= $min_year & yofd(Date0) <= $max_year
			local p = r(N)/`noth'
		}

		local note ""
		if `noth' < `minoth' {
			local p = `pall'
			local note "  <- only `noth' in window, using all-years"
		}
		mata: LENREFR_pother[1, `l'] = `p'
		qui drop refr_len_oth
	}
	global Coeffs $Coeffs LENREFR_pother

	***** LENALIDOMIDE-REFRACTORY, L1 MAINTENANCE (LENREFR_MNT) *****
	// Replaces the withdrawn tail rule, which selected long-maintenance good-prognosis patients because
	// sim_tfi_l1 draws the gap truncated just above the duration. Outcome refr_len_mnt_l1 is that rule
	// applied to OBSERVED data, built in data_extraction.do: progressed to L2 within 60 days of
	// maintenance ending. It is MISSING where maintenance has not ended, so logit drops the
	// undetermined patients without a separate flag. Pooled across transplant with an SCT main effect;
	// response collapsed to CR/VG/PR/poor because SD perfectly predicts on a handful of rows yet
	// sim_bcr.do can draw it. docs/refractory.md 4.4.
	// esampvaryok: bcr_grp_l1 derives from imputed BCR_L1, so its levels can vary across m.
	mi estimate, esampvaryok: logit refr_len_mnt_l1 Age Age2 Male i.ECOGcc i.RISS ///
		CM_CKD CM_CRD CM_PLM CM_DBT SCT i.bcr_grp_l1 ///
		if(Event1 == 110 & MNT == 1 & MNR_L1 == 1)
	save_coefs LENREFR_MNT
	mata: _matrix_list(bLENREFR_MNT, rbLENREFR_MNT, cbLENREFR_MNT)

	// Treatment- and maintenance-dose refractoriness carry the same conditional OS hazard, so one
	// collapsed flag is fitted and the engine keeps one latched state. docs/refractory.md 5.
	cap drop refr_len_in
	gen byte refr_len_in = (refr_len_tx_in == 1 | refr_len_mnt_in == 1)

	***** OVERALL SURVIVAL *****
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
	mi estimate: streg Age Age2 Male i.ECOGcc i.RISS CM_CKD CM_CRD CM_PLM CM_DBT i.BCR_L2 refr_len_in, d($dOS)
	save_coefs OS_L2S
	mata: _matrix_list(bOS_L2S, rbOS_L2S, cbOS_L2S)

	// OS_L2E: L2E -> L3S
	mi stset Date1 if(F_OS != 1), id(ID_BS) failure(Event1 == 104) origin(Event1 == 21) exit(Event1 == 30) scale(30.4375)
	mi estimate: streg Age Age2 Male i.ECOGcc i.RISS CM_CKD CM_CRD CM_PLM CM_DBT i.BCR_L2 refr_len_in, d($dOS)
	save_coefs OS_L2E
	mata: _matrix_list(bOS_L2E, rbOS_L2E, cbOS_L2E)

	// OS_L3S: L3S -> L3E
	mi stset Date1 if(F_OS != 1), id(ID_BS) failure(Event1 == 104) origin(Event1 == 30) exit(Event1 == 31) scale(30.4375)
	mi estimate: streg Age Age2 Male i.ECOGcc i.RISS CM_CKD CM_CRD CM_PLM CM_DBT i.BCR_L3 refr_len_in, d($dOS)
	save_coefs OS_L3S
	mata: _matrix_list(bOS_L3S, rbOS_L3S, cbOS_L3S)

	// OS_L3E: L3E -> L4S
	mi stset Date1 if(F_OS != 1), id(ID_BS) failure(Event1 == 104) origin(Event1 == 31) exit(Event1 == 40) scale(30.4375)
	mi estimate: streg Age Age2 Male i.ECOGcc i.RISS CM_CKD CM_CRD CM_PLM CM_DBT i.BCR_L3 refr_len_in, d($dOS)
	save_coefs OS_L3E
	mata: _matrix_list(bOS_L3E, rbOS_L3E, cbOS_L3E)

	// OS_L4S: L4S -> L4E
	mi stset Date1 if(F_OS != 1), id(ID_BS) failure(Event1 == 104) origin(Event1 == 40) exit(Event1 == 41) scale(30.4375)
	mi estimate: streg Age Age2 Male i.ECOGcc i.RISS CM_CKD CM_CRD CM_PLM CM_DBT i.BCR_L4 refr_len_in, d($dOS)
	save_coefs OS_L4S
	mata: _matrix_list(bOS_L4S, rbOS_L4S, cbOS_L4S)

	// OS_L4E: L4E -> L5S
	mi stset Date1 if(F_OS != 1), id(ID_BS) failure(Event1 == 104) origin(Event1 == 41) exit(Event1 == 50) scale(30.4375)
	mi estimate: streg Age Age2 Male i.ECOGcc i.RISS CM_CKD CM_CRD CM_PLM CM_DBT i.BCR_L4 refr_len_in, d($dOS)
	save_coefs OS_L4E
	mata: _matrix_list(bOS_L4E, rbOS_L4E, cbOS_L4E)

	// OS_L5plus: L5S onward, pooled, with an L7+ stage term. Replaces OS_L5S, OS_L5E and OS_L6plus.
	// The boundary was in the wrong place: L5 and L6 are indistinguishable (stage 0.07, p = 0.70)
	// while L7+ is not (0.78, p < 0.001, HR 2.2), so one indicator carries what three equations
	// carried badly. 291 deaths against OS_L5S's 49 - which is what fixes the separation, since
	// OS_L5S's CR cell drew ZERO deaths in 2 imputations of 10 and its BCR block was unidentified
	// (SE 462, HRs of 50-140) in every generation shipped to date. The BCR gradient is common
	// across stages (all interactions p > 0.17), so no response collapse is needed here.
	// Evidence: scratch/l5_pooling_test.log, scratch/os_l5s_check.log.
	mi stset Date1 if(F_OS != 1), id(ID_BS) failure(Event1 == 104) origin(Event1 == 50) scale(30.4375)
	mi estimate: streg Age Age2 Male i.ECOGcc i.RISS CM_CKD CM_CRD CM_PLM CM_DBT i.BCR L7plus, d($dOS)
	save_coefs OS_L5plus
	mata: _matrix_list(bOS_L5plus, rbOS_L5plus, cbOS_L5plus)

	***** ASCT (TRANSPLANT DECISION) *****
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
		mi estimate: mlogit TXR_L2 Age Age2 refr_len_in if(Event0 == 20 & yofd(Date0) >= $min_year & yofd(Date0) <= $max_year), baseoutcome(0)
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
		mi estimate: mlogit TXR_L3 Age Age2 refr_len_in if(Event0 == 30 & yofd(Date0) >= $min_year & yofd(Date0) <= $max_year), baseoutcome(0)
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
		mi estimate: mlogit TXR_L4 Age Age2 refr_len_in if(Event0 == 40 & yofd(Date0) >= $min_year & yofd(Date0) <= $max_year), baseoutcome(0)
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
	// Time on treatment within a line (Weibull). L1 splits by transplant: ASCT fixed-duration is a
	// 3-spline fit, NoASCT is a single fit, and Continuous (Rd) is separate.
	// L2+ are single fits conditional on that line's response and regimen.

	// L1 - Fixed w/ ASCT (3 splines)
	// Cut off 1 - 60 days for optimal fit but too few event, 90 days selected
	scalar L1_TXD_ASCT_C1 = 90
	mata: L1_TXD_ASCT_C1 = st_numscalar("L1_TXD_ASCT_C1")

	// Cut off 2 - 120 days for optimal fit with 3 splines
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
	// The maintenance analogue of TXR_L1 / TXD_L1, sat here because both condition on the MNT
	// logit above: MNT decides WHETHER, L1_MNR WHICH, L1_MND HOW LONG. Both fit at the same
	// Event1 == 11 row as MNT, and only among MNT == 1. Consumed only by cost_tx_mnt, so the
	// out-of-sample validation fits them but never uses them. See docs/refractory.md 7.4.

	// An analysis may declare lenalidomide alone: the thalidomide fits below then have no sample and
	// are skipped, oL1_MNR = 1, and sim_mnd's thalidomide branch is a no-op. Lenalidomide itself is
	// not optional - it is the mlogit base outcome and sim_mnr's single-regimen fallback.
	qui count if(MNT == 1 & MNR_L1 == 1)
	local nLen = r(N)
	qui count if(MNT == 1 & MNR_L1 == 5)
	local nThal = r(N)
	if `nLen' == 0 {
		di as err "risk_equations: \$MNR_L1 is '$MNR_L1' but no lenalidomide records remain. Declare 1 in analyses/$analysis/outcomes/txr_$coeffs.do."
		exit 2000
	}

	// Which regimen: lenalidomide (base) vs thalidomide, so effectively a logit. Year-windowed as
	// TXR_L1 is. The guard counts the estimation sample, so it catches a list that omits thalidomide
	// and one that keeps it but lands in a window where it is empty. docs/refractory.md 7.4.
	qui count if(Event1 == 11 & MNT == 1 & MNR_L1 == 5 & yofd(Date0) >= $min_year & yofd(Date0) <= $max_year)
	local nThalFit = r(N)
	qui count if(Event1 == 11 & MNT == 1 & MNR_L1 == 1 & yofd(Date0) >= $min_year & yofd(Date0) <= $max_year)
	if `nThalFit' > 0 & r(N) > 0 {
		mi estimate: mlogit MNR_L1 Age Age2 Male i.ECOGcc i.RISS SCT ///
			if(Event1 == 11 & MNT == 1 & inlist(MNR_L1, 1, 5) & yofd(Date0) >= $min_year & yofd(Date0) <= $max_year), baseoutcome(1)
		save_coefs L1_MNR
		mata: _matrix_list(bL1_MNR, rbL1_MNR, cbL1_MNR)
	}
	else {
		mata: oL1_MNR = 1
		global Coeffs $Coeffs oL1_MNR
	}

	// Maintenance DURATION via parametric survival on the 110/111 maintenance events, as L1_TFI is
	// fitted. Split by regimen, pooled across transplant (SCT a covariate), no gap term - ln(TFI_L1)
	// is missing without an observed L2; sim_tfi_l1.do supplies the dependence instead by drawing the
	// gap truncated below at the duration. Thalidomide is censored at 18 months, a judgement that
	// later recorded ends are stale; generate_benchmarks.do carries the same exit().
	// docs/refractory.md 4.4 / 5(7) / 7; rejected specifications in scratch/maintenance/notes.md.
	//
	// TRAP: exit() is in the TIME VARIABLE's scale, and `origin' is unavailable inside it when
	// origin() is an event condition, so it must reference the origin date in days.
	capture drop MND_origin
	qui gen double _mxMNDo = Date1 if Event1 == 110
	qui egen double MND_origin = min(_mxMNDo), by(ID_BS)
	qui drop _mxMNDo

	// Lenalidomide, pooled across transplant, baseline covariates only. No response term and no
	// transplant split: once BCR_SCT was correctly imputed, no response specification predicted
	// duration. Ceiling and family are settled - do not lower the ceiling.
	// Evidence: scratch/maintenance/notes.md.
	mi stset Date1 if(MNT == 1 & MNR_L1 == 1), ///
		id(ID_BS) failure(Event1 == 20 111) origin(Event1 == 110) scale(30.4375)
	save_max L1_MND_LEN
	mi estimate: streg Age Age2 Male i.ECOGcc i.RISS SCT, d($dTFI)
	save_coefs L1_MND_LEN
	mata: _matrix_list(bL1_MND_LEN, rbL1_MND_LEN, cbL1_MND_LEN)

	// Thalidomide, censored at 18 months. Pooled across transplant, no BCR: n is roughly 311 against
	// lenalidomide's 1,020, too small to add a factor. Thalidomide maintenance also ended with the
	// 2020 PBS listing, so it matters to historical validation rather than any current analysis -
	// which is why an analysis may legitimately leave it out of $MNR_L1 and skip this block.
	if `nThal' > 0 {
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
	}
	else {
	}

	***** TREATMENT-FREE INTERVAL (TFI) *****
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
	// Response COLLAPSED to CR vs not-CR. On n = 302 the six-level block was carrying a single
	// binary contrast - every non-CR coefficient fell between -1.23 and -1.43 against SEs near 0.85 -
	// and a level empty in one imputation but not another failed `mi estimate' ("omitted terms vary")
	// on 5 of 500 bootstrap replicates. sim_tfi.do reads the block width from the coefficients.
	mi estimate: streg Age Age2 Male i.ECOGcc i.RISS i.bcr_grp_l5, d($dTFI)
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
