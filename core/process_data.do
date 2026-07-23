**********
* Monash Myeloma Model - Process Data
*
* Purpose: Assemble the Mata simulation matrices into a per-patient dataset and derive
*          costs and QALYs starting at line $line ($line == 1: full pathway from diagnosis
*          with calendar-date discounting; $line > 1: from L$line onwards with relative-time
*          discounting, L$line start = time 0).
* Notes:   Variable naming: cost_{component}[_L#][_d] (tx, nt, total), qaly_{component}[_d];
*          _d suffix = discounted. Frees mRN (peak-memory stage).
**********

capture program drop process_data
program process_data

local L = $line
if `L' == 0 local L = 1   // Line 0 = full pathway from diagnosis
local maxL = 9

di as text "Processing Simulated Data (Starting Line: `L')"

* Free the CRN matrix - all draws are done by end of simulation; not needed here,
*   and process_data is the peak-memory stage (mSum assembly + getmata).
cap mata: mata drop mRN

* Create mSum in Mata
	mata: mSum = vID , vMale , vECOG , vRISS , vISS , vCKD , vCRD , vPLM , vDBT , vAge70 , vAge75 , vSCT_DN , vSCT_L1 , vMNT , vMNR , vMND , mLenRefr_in , ///
			mAge , mOS , mTNE , mTSD , mMOR , mOC , mTXR , mTXD , mBCR , mTFI , mState
	
* Column names for mSum, in assembly order below.
* (getmata errors on a name/column count mismatch, which guards this alignment.)
	local varnames ID Male ECOGcc RISS ISS CM_CKD CM_CRD CM_PLM CM_DBT Age70 Age75 SCT_DN SCT_L1 MNT MNR_L1 MND_L1 ///
		LenRefr_L1 LenRefr_L2 LenRefr_L3 LenRefr_L4 LenRefr_L5 LenRefr_L6 LenRefr_L7 LenRefr_L8 LenRefr_L9 ///
		Age_DN Age_L1S Age_L1E Age_L2S Age_L2E Age_L3S Age_L3E Age_L4S Age_L4E Age_L5S Age_L5E Age_L6S Age_L6E Age_L7S Age_L7E Age_L8S Age_L8E Age_L9S Age_L9E ///
		OS_DN OS_L1S OS_L1E OS_L2S OS_L2E OS_L3S OS_L3E OS_L4S OS_L4E OS_L5S OS_L5E OS_L6S OS_L6E OS_L7S OS_L7E OS_L8S OS_L8E OS_L9S OS_L9E ///
		TNE_DN TNE_L1S TNE_L1E TNE_L2S TNE_L2E TNE_L3S TNE_L3E TNE_L4S TNE_L4E TNE_L5S TNE_L5E TNE_L6S TNE_L6E TNE_L7S TNE_L7E TNE_L8S TNE_L8E TNE_L9S TNE_L9E ///
		TSD_DN TSD_L1S TSD_L1E TSD_L2S TSD_L2E TSD_L3S TSD_L3E TSD_L4S TSD_L4E TSD_L5S TSD_L5E TSD_L6S TSD_L6E TSD_L7S TSD_L7E TSD_L8S TSD_L8E TSD_L9S TSD_L9E ///
		MOR_DN MOR_L1S MOR_L1E MOR_L2S MOR_L2E MOR_L3S MOR_L3E MOR_L4S MOR_L4E MOR_L5S MOR_L5E MOR_L6S MOR_L6E MOR_L7S MOR_L7E MOR_L8S MOR_L8E MOR_L9S MOR_L9E ///
		OC_TIME OC_MORT ///
		TXR_L1 TXR_L2 TXR_L3 TXR_L4 TXR_L5 TXR_L6 TXR_L7 TXR_L8 TXR_L9 ///
		TXD_L1 TXD_L2 TXD_L3 TXD_L4 TXD_L5 TXD_L6 TXD_L7 TXD_L8 TXD_L9 ///
		BCR_L1 BCR_L2 BCR_L3 BCR_L4 BCR_L5 BCR_L6 BCR_L7 BCR_L8 BCR_L9 BCR_SCT ///
		TFI_DN TFI_L1 TFI_L2 TFI_L3 TFI_L4 TFI_L5 TFI_L6 TFI_L7 TFI_L8 ///
		State DateDN ///
		

* Write the Mata matrix straight to named variables. getmata reads mSum
* directly, bypassing the st_matrix()/svmat round-trip that dominated runtime
* (st_matrix on the full ~101k x 148 matrix was ~60s; getmata is < 1s).
	drop _all
	getmata (`varnames') = mSum, double
	cap mata: mata drop mSum
	
	format DateDN %td
	order ID Male ECOGcc RISS SCT_L1 MNT CM_CKD
		
* Label
	label values State State_lbl

**********
* Time Variables
**********

* Survival time from starting line to death (months)
	if `L' == 1 {
		qui gen OC_TIME_L = OC_TIME
	}
	else {
		qui gen OC_TIME_L = OC_TIME - TSD_L`L'S
	}

* Relative time markers (starting line = time 0)
	qui {
		if `L' == 1 {
			gen TSD_L1S_ref = TFI_DN if TFI_DN != .
			replace TSD_L1S_ref = 0 if TFI_DN == .
			gen TSD_L1E_ref = TSD_L1S_ref + TXD_L1 if TXD_L1 != .
			local first_l = 2
		}
		else if `L' > 1 {
			gen TSD_L`L'S_ref = 0
			gen TSD_L`L'E_ref = TXD_L`L' if TXD_L`L' != .
			local first_l = `=`L'+1'
		}
		local prev = `L'
		forval l = `first_l'/`maxL' {
			local tfi_idx = `=`l'-1'
			gen TSD_L`l'S_ref = TSD_L`prev'E_ref + TFI_L`tfi_idx' if TFI_L`tfi_idx' != .
			gen TSD_L`l'E_ref = TSD_L`l'S_ref + TXD_L`l' if TXD_L`l' != .
			local prev = `l'
		}
	}

* Generate Dates & Years from L onwards
	qui {
		forval l = `L'/`maxL' {
			gen DateL`l'S = DateDN + (TSD_L`l'S * 30.4375)
			gen DateL`l'E = DateDN + (TSD_L`l'E * 30.4375)
		}
		
		if `L' == 1 {
			gen DateSCT = DateL1E + 1 if SCT_L1 == 1
			gen YearSCT = yofd(DateSCT)
		}
		gen DateMOR = DateDN + (OC_TIME * 30.4375)
		format Date* %td
		
		gen YearDN = yofd(DateDN)
		forval l = `L'/`maxL' {
			gen YearL`l' = yofd(DateL`l'S)
		}
		gen YearMOR = yofd(DateMOR)
	}

**********
* Costs
**********

* ---- Per-cycle regimen + non-drug costs, from the derived PBS DPMQ table ----
* prep/treatment_costs.do builds prep/inputs/treatment_costs_<year>.csv (PBS first-principles DPMQ)
* from the committed PBS inputs. Read it into a FRAME so the in-memory (peak-size) simulation data is
* untouched. Uses $cost_year; if that year's file is absent, falls back to the latest available with a
* note. (Values were hardcoded here until Jul 2026 - see docs/economic_inputs.md.)
	local cyear = "$cost_year"
	if "`cyear'" == "" local cyear 2026
	local costfile "prep/inputs/treatment_costs_`cyear'.csv"
	capture confirm file "`costfile'"
	if _rc {
		local avail : dir "prep/inputs" files "treatment_costs_*.csv"
		local latest ""
		foreach f of local avail {
			if "`f'" > "`latest'" local latest "`f'"
		}
		if "`latest'" == "" {
			di as error "process_data: no prep/inputs/treatment_costs_*.csv - run prep/treatment_costs.do first"
			exit 601
		}
		if "$cost_fallback_note" != "`latest'" {   // note once per session (process_data runs many times)
			di as text "  (cost year `cyear' not found; using `latest')"
			global cost_fallback_note "`latest'"
		}
		local costfile "prep/inputs/`latest'"
	}
	capture frame drop _costs
	frame create _costs
	frame _costs: quietly import delimited "`costfile'", varnames(1) case(preserve) stringcols(_all) clear
	frame _costs: quietly destring value, replace     // to double (avoids float import rounding)
	foreach p in cVCd cVRd cRd cPd cVd cOther cMNT cASCT ///
	             cOther_L1 cOther_L2 cOther_L3 cOther_L4 cOther_L5 cOther_L6 ///
	             cKd_p1 cKd_p2 cDVd_p1 cDVd_p2 cDVd_p3 ///
	             cHosp_initial cHosp_continuing cHosp_terminal ///
	             cMBS_initial  cMBS_continuing  cMBS_terminal ///
	             cEmer_initial cEmer_continuing cEmer_terminal {
		frame _costs: quietly summarize value if parameter == "`p'", meanonly
		local `p' = r(mean)
	}
	frame drop _costs

	local ln_r = ln(1 + $drate)

* DVd and Kd are phase-based (front-loaded): the per-cycle cost changes over the course, so cost is
* allocated to phase windows (months from regimen start) at each phase's per-month rate rather than a
* single flat rate. DVd: load = cycles 1-3 (21-day, dara weekly), mid = cycles 4-8 (21-day, dara q3w),
* tail = cycles 9+ (28-day, dara-only). Kd: load = cycle 1 (28-day step-up), maint = cycles 2+ (28-day).
	local dvd_le = 3*21/30.4375                 // DVd phase-1 (load) ends (months)
	local dvd_me = 168/30.4375                  // DVd phase-2 (mid) ends (months) = 8 x 21-day cycles
	local dvd_lr = `cDVd_p1' * 30.4375/21       // DVd per-month rate: p1 / p2 / p3
	local dvd_mr = `cDVd_p2' * 30.4375/21
	local dvd_tr = `cDVd_p3' * 30.4375/28
	local kd_le  = 28/30.4375                   // Kd phase-1 (load) ends (months) = 1 x 28-day cycle
	local kd_lr  = `cKd_p1' * 30.4375/28
	local kd_mr  = `cKd_p2' * 30.4375/28

* Treatment costs by line (undiscounted)
	forval l = `L'/`maxL' {
		qui gen cost_tx_L`l' = 0
		
		qui replace cost_tx_L`l' = `cVCd' * min(4, TXD_L`l' * 30.4375 / 21) if TXR_L`l' == 4
		qui replace cost_tx_L`l' = `cVRd' * min(5, TXD_L`l' * 30.4375 / 21) if TXR_L`l' == 31
		qui replace cost_tx_L`l' = `cRd' * (TXD_L`l' * 30.4375 / 28) if TXR_L`l' == 7
		qui replace cost_tx_L`l' = `kd_lr'*min(TXD_L`l',`kd_le') + `kd_mr'*max(0,TXD_L`l'-`kd_le') if TXR_L`l' == 49
		qui replace cost_tx_L`l' = `dvd_lr'*min(TXD_L`l',`dvd_le') + `dvd_mr'*max(0,min(TXD_L`l',`dvd_me')-`dvd_le') + `dvd_tr'*max(0,TXD_L`l'-`dvd_me') if TXR_L`l' == 80
		qui replace cost_tx_L`l' = `cPd' * (TXD_L`l' * 30.4375 / 28) if TXR_L`l' == 56
		qui replace cost_tx_L`l' = `cVd' * min(8, TXD_L`l' * 30.4375 / 21) if TXR_L`l' == 5
		* Pooled "Other" (code 0) is costed at a USAGE-WEIGHTED, LINE-SPECIFIC rate where one is
		* available (prep/inputs/other_mix.csv via treatment_costs.do), because what is prescribed
		* outside the modelled regimens differs sharply by line. L6+ share line 6's rate. Falls back
		* to the single fixed c_Other blend when the per-line rates are absent.
		local _og = min(`l', 6)
		local _orate "`cOther_L`_og''"
		if ("`_orate'" == "" | "`_orate'" == ".") local _orate "`cOther'"
		qui replace cost_tx_L`l' = `_orate' * (TXD_L`l' * 30.4375 / 28) if TXR_L`l' == 0
	}

* L1-specific costs (ASCT, MNT)
	qui gen cost_tx_asct = 0
	qui gen cost_tx_mnt = 0
	if `L' == 1 {
		qui replace cost_tx_asct = `cASCT' if SCT_L1 == 1
		// Bill the maintenance actually delivered, not the whole L1-to-L2 gap. MND_L1 is the
		// maintenance DURATION (months) drawn by sim_mnd.do and carried out of the engine. Cap
		// it at the REALISED TFI_L1 - sim_mort has curtailed TFI_L1 at death by now, so a patient
		// who dies mid-gap is billed only for the maintenance they lived to receive, and a draw
		// that overshoots the gap becomes continuous-to-progression maintenance (the patient
		// stayed on maintenance until relapse; docs/refractory.md 4.4). Billing the whole TFI_L1
		// overstated maintenance by 69% population-wide (docs/refractory.md 7).
		// The !mi() guard matters: MND_L1 is missing for MNT == 0 and for every patient if
		// L1_MND was never fitted, and an unguarded bill would propagate missing through cost_tx
		// into every downstream total. Leaving cost_tx_mnt at its 0 initialisation is the safe
		// direction - no maintenance cost beats a silently missing total.
		// INTERIM: still the blended cMNT (roughly the all-years lenalidomide/thalidomide mix,
		// so wrong for any single year). MNR_L1 is simulated and carried out of the engine, but
		// regimen-specific pricing needs maintenance DPMQs that mostly do not exist: per the MSAG
		// guideline only lenalidomide has a PBS maintenance listing. See docs/refractory.md 7.4
		// and the Costs todo in the programme note.
		// SAFETY NET, no longer the mechanism. sim_tfi_l1.do draws the gap truncated below at the
		// maintenance duration, so MND_L1 > TFI_L1 should now be impossible. It is kept because the
		// alternative to a stray overshoot is billing maintenance the patient could not have had,
		// and because it still catches the death curtailment sim_mort applies after both draws.
		// If it fires on any material share of patients the truncation is not working - count it.
		qui count if MNT == 1 & !mi(MND_L1) & !mi(TFI_L1) & MND_L1 > TFI_L1 + 0.01
		if r(N) > 0 {
			di as error "  NOTE: MND_L1 exceeded TFI_L1 for " r(N) " patients despite the truncated"
			di as error "        TFI draw. Expected ~0 (death curtailment aside); investigate if large."
		}
		qui replace MND_L1 = TFI_L1 if MNT == 1 & !mi(MND_L1) & MND_L1 > TFI_L1
		qui replace cost_tx_mnt = `cMNT' * (MND_L1 * 30.4375 / 28) if MNT == 1 & !mi(MND_L1)
	}

* Total undiscounted treatment cost
	qui gen cost_tx = cost_tx_asct + cost_tx_mnt
	forval l = `L'/`maxL' {
		qui replace cost_tx = cost_tx + cost_tx_L`l'
	}

* Non-treatment costs (undiscounted) - Yap 2025 phase-based (initial / continuing / terminal).
* Phases are defined on the diagnosis clock: initial = first 12 months, terminal = last 12 months of
* life, continuing = between (Yap short-survivor rules: <2 yr -> last year terminal + remainder initial;
* <1 yr -> all terminal). We cost only the survival window from the starting line, [w0, OC_TIME], where
* w0 = OC_TIME - OC_TIME_L (=0 at L1). Each phase contributes rate x (months of overlap)/12.
	qui gen double _w0 = OC_TIME - OC_TIME_L
	qui gen double _i1 = cond(OC_TIME<12, 0, cond(OC_TIME<24, OC_TIME-12, 12))   // initial ends
	qui gen double _c1 = cond(OC_TIME>=24, OC_TIME-12, 12)                       // continuing ends
	qui gen double _t0 = cond(OC_TIME<12, 0, OC_TIME-12)                         // terminal starts
	qui gen double _os_i = max(0,   _w0)
	qui gen double _oe_i = min(_i1, OC_TIME)
	qui gen double _os_c = max(12,  _w0)
	qui gen double _oe_c = min(_c1, OC_TIME)
	qui gen double _os_t = max(_t0, _w0)
	qui gen double _oe_t = OC_TIME
	* months of overlap in each phase (computed once, applied to every component)
	qui gen double _mo_i = max(0, _oe_i-_os_i)/12
	qui gen double _mo_c = max(0, _oe_c-_os_c)/12
	qui gen double _mo_t = max(0, _oe_t-_os_t)/12
	* kept split by component (Yap 2025): admitted-hospital + out-of-hospital Medicare + emergency
	qui gen cost_nt_hosp = `cHosp_initial'*_mo_i + `cHosp_continuing'*_mo_c + `cHosp_terminal'*_mo_t
	qui gen cost_nt_mbs  = `cMBS_initial' *_mo_i + `cMBS_continuing' *_mo_c + `cMBS_terminal' *_mo_t
	qui gen cost_nt_emer = `cEmer_initial'*_mo_i + `cEmer_continuing'*_mo_c + `cEmer_terminal'*_mo_t
	qui gen cost_nt = cost_nt_hosp + cost_nt_mbs + cost_nt_emer

* Total undiscounted cost
	qui gen cost_total = cost_tx + cost_nt

* Discounted treatment costs
	qui gen cost_tx_d = 0
	
	* First line (starts at time 0 or after TFI_DN for L1)
	if `L' == 1 {
		qui gen cost_tx_L1_d = cost_tx_L1 * ((1 + $drate)^(-TSD_L1S_ref/12) - (1 + $drate)^(-TSD_L1E_ref/12)) / (`ln_r' * TXD_L1/12) if cost_tx_L1 > 0 & TXD_L1 > 0
		local first_l = 2
	}
	else {
		qui gen cost_tx_L`L'_d = cost_tx_L`L' * (1 - (1 + $drate)^(-TSD_L`L'E_ref/12)) / (`ln_r' * TXD_L`L'/12) if cost_tx_L`L' > 0 & TXD_L`L' > 0
		local first_l = `=`L'+1'
	}
	qui replace cost_tx_d = cost_tx_d + cost_tx_L`L'_d if cost_tx_L`L'_d != .
	
	* Subsequent lines
	forval l = `first_l'/`maxL' {
		qui gen cost_tx_L`l'_d = cost_tx_L`l' * ((1 + $drate)^(-TSD_L`l'S_ref/12) - (1 + $drate)^(-TSD_L`l'E_ref/12)) / (`ln_r' * TXD_L`l'/12) if cost_tx_L`l' > 0 & TXD_L`l' > 0
		qui replace cost_tx_d = cost_tx_d + cost_tx_L`l'_d if cost_tx_L`l'_d != .
	}
	
	* ASCT and MNT (L1 only)
	if `L' == 1 {
		qui gen cost_tx_asct_d = cost_tx_asct * (1 + $drate)^(-TSD_L1E_ref/12) if SCT_L1 == 1
		qui replace cost_tx_d = cost_tx_d + cost_tx_asct_d if cost_tx_asct_d != .
		
		qui gen cost_tx_mnt_d = cost_tx_mnt * ((1 + $drate)^(-TSD_L1E_ref/12) - (1 + $drate)^(-TSD_L2S_ref/12)) / (`ln_r' * TFI_L1/12) if cost_tx_mnt > 0 & TFI_L1 > 0
		qui replace cost_tx_d = cost_tx_d + cost_tx_mnt_d if cost_tx_mnt_d != .
	}
	
	* Non-treatment costs (discounted) - each phase discounted over its own sub-interval, in relative
	* time from the starting line (a = os - w0, b = oe - w0). The uniform-over-interval factor
	* ((1+r)^(-a/12) - (1+r)^(-b/12))/ln(1+r) tends to (b-a)/12 as r->0, matching the undiscounted amount.
	* per-phase discount factor (uniform-over-interval), computed once and applied to every component
	qui gen double _df_i = cond(_oe_i > _os_i, ((1+$drate)^(-(_os_i-_w0)/12) - (1+$drate)^(-(_oe_i-_w0)/12))/`ln_r', 0)
	qui gen double _df_c = cond(_oe_c > _os_c, ((1+$drate)^(-(_os_c-_w0)/12) - (1+$drate)^(-(_oe_c-_w0)/12))/`ln_r', 0)
	qui gen double _df_t = cond(_oe_t > _os_t, ((1+$drate)^(-(_os_t-_w0)/12) - (1+$drate)^(-(_oe_t-_w0)/12))/`ln_r', 0)
	qui gen cost_nt_hosp_d = `cHosp_initial'*_df_i + `cHosp_continuing'*_df_c + `cHosp_terminal'*_df_t
	qui gen cost_nt_mbs_d  = `cMBS_initial' *_df_i + `cMBS_continuing' *_df_c + `cMBS_terminal' *_df_t
	qui gen cost_nt_emer_d = `cEmer_initial'*_df_i + `cEmer_continuing'*_df_c + `cEmer_terminal'*_df_t
	qui gen cost_nt_d = cost_nt_hosp_d + cost_nt_mbs_d + cost_nt_emer_d

	* Undiscounted cost accrued from diagnosis to 5 years (60 months) - comparable to the Yap 2025
	* diagnosis-to-5-year excess cost. Only meaningful for a from-diagnosis run ($line == 1); each
	* stream is truncated to its overlap with [0, 60] months on the diagnosis clock (missing otherwise).
	if `L' == 1 {
		qui gen double cost_5yr = 0
		forval l = 1/`maxL' {
			qui replace cost_5yr = cost_5yr + cost_tx_L`l' * ///
				max(0, min(TSD_L`l'E_ref,60) - max(TSD_L`l'S_ref,0)) / (TSD_L`l'E_ref - TSD_L`l'S_ref) ///
				if cost_tx_L`l' > 0 & TSD_L`l'E_ref > TSD_L`l'S_ref
		}
		qui replace cost_5yr = cost_5yr + cost_tx_asct if SCT_L1 == 1 & TSD_L1E_ref <= 60
		qui replace cost_5yr = cost_5yr + cost_tx_mnt * ///
			max(0, min(TSD_L2S_ref,60) - max(TSD_L1E_ref,0)) / TFI_L1 if cost_tx_mnt > 0 & TFI_L1 > 0
		* non-treatment: phase rates over the window capped at min(60, OC_TIME) (reuse _i1/_c1/_t0)
		qui replace cost_5yr = cost_5yr ///
			+ (`cHosp_initial'   +`cMBS_initial'   +`cEmer_initial')    * max(0, min(_i1,60)-0 )/12 ///
			+ (`cHosp_continuing'+`cMBS_continuing'+`cEmer_continuing') * max(0, min(_c1,60)-12)/12 ///
			+ (`cHosp_terminal'  +`cMBS_terminal'  +`cEmer_terminal')   * max(0, min(OC_TIME,60)-_t0)/12
		* death within the window => whole lifetime falls in [0,60], so 5-yr cost = lifetime total
		qui replace cost_5yr = cost_total if OC_TIME <= 60
	}
	else qui gen double cost_5yr = .
	qui drop _w0 _i1 _c1 _t0 _os_i _oe_i _os_c _oe_c _os_t _oe_t _mo_i _mo_c _mo_t _df_i _df_c _df_t

* Total discounted costs
	qui gen cost_total_d = cost_tx_d + cost_nt_d

**********
* QALYs
**********

* L1 analysis
	if `L' == 1 {
		
		local uTFI = 0.72
		local uTXD_L1 = 0.63
		local uTXD_L2 = 0.67
		local uPostL2 = 0.63
		
		qui gen PreL2 = 0
		
		* TFI_DN
		qui gen qaly_tfi_DN = 0
		qui replace qaly_tfi_DN = (TFI_DN / 12) * `uTFI' if TFI_DN != .
		qui replace PreL2 = PreL2 + TFI_DN if TFI_DN != .
		
		* TXD_L1
		qui gen qaly_txd_L1 = 0
		qui replace qaly_txd_L1 = (TXD_L1 / 12) * `uTXD_L1' if TXD_L1 != .
		qui replace PreL2 = PreL2 + TXD_L1 if TXD_L1 != .
		
		* TFI_L1
		qui gen qaly_tfi_L1 = 0
		qui replace qaly_tfi_L1 = (TFI_L1 / 12) * `uTFI' if TFI_L1 != .
		qui replace PreL2 = PreL2 + TFI_L1 if TFI_L1 != .
		
		* TXD_L2
		qui gen qaly_txd_L2 = 0
		qui replace qaly_txd_L2 = (TXD_L2 / 12) * `uTXD_L2' if TXD_L2 != .
		qui replace PreL2 = PreL2 + TXD_L2 if TXD_L2 != .
		
		* PostL2
		qui gen qaly_post_L2 = 0
		qui replace qaly_post_L2 = ((OC_TIME_L - PreL2) / 12) * `uPostL2' if OC_TIME_L > PreL2 & PreL2 < .
		drop PreL2
		
		* Total undiscounted QALYs
		qui gen qaly_total = qaly_tfi_DN + qaly_txd_L1 + qaly_tfi_L1 + qaly_txd_L2 + qaly_post_L2
		
		* Discounting (using _ref time markers)
		
		* TFI_DN: from time 0 to TSD_L1S_ref
		qui gen qaly_tfi_DN_d = 0
		qui replace qaly_tfi_DN_d = qaly_tfi_DN * (1 - (1 + $drate)^(-TSD_L1S_ref/12)) / (`ln_r' * TSD_L1S_ref/12) if TSD_L1S_ref > 0 & qaly_tfi_DN > 0
		
		* TXD_L1: from TSD_L1S_ref to TSD_L1E_ref
		qui gen qaly_txd_L1_d = 0
		qui replace qaly_txd_L1_d = qaly_txd_L1 * ((1 + $drate)^(-TSD_L1S_ref/12) - (1 + $drate)^(-TSD_L1E_ref/12)) / (`ln_r' * TXD_L1/12) if TXD_L1 > 0 & qaly_txd_L1 > 0
		
		* TFI_L1: from TSD_L1E_ref to TSD_L2S_ref
		qui gen qaly_tfi_L1_d = 0
		qui replace qaly_tfi_L1_d = qaly_tfi_L1 * ((1 + $drate)^(-TSD_L1E_ref/12) - (1 + $drate)^(-TSD_L2S_ref/12)) / (`ln_r' * TFI_L1/12) if TFI_L1 > 0 & qaly_tfi_L1 > 0
		
		* TXD_L2: from TSD_L2S_ref to TSD_L2E_ref
		qui gen qaly_txd_L2_d = 0
		qui replace qaly_txd_L2_d = qaly_txd_L2 * ((1 + $drate)^(-TSD_L2S_ref/12) - (1 + $drate)^(-TSD_L2E_ref/12)) / (`ln_r' * TXD_L2/12) if TXD_L2 > 0 & qaly_txd_L2 > 0
		
		* PostL2: from TSD_L2E_ref to OC_TIME_L
		qui gen qaly_post_L2_d = 0
		qui replace qaly_post_L2_d = qaly_post_L2 * ((1 + $drate)^(-TSD_L2E_ref/12) - (1 + $drate)^(-OC_TIME_L/12)) / (`ln_r' * (OC_TIME_L - TSD_L2E_ref)/12) if qaly_post_L2 > 0 & OC_TIME_L > TSD_L2E_ref
		
		* Total discounted QALYs
		qui gen qaly_total_d = qaly_tfi_DN_d + qaly_txd_L1_d + qaly_tfi_L1_d + qaly_txd_L2_d + qaly_post_L2_d
	}
* L2 analysis
	else if `L' == 2 {
		
		local uTXD_L2 = 0.67
		local uPostL2 = 0.63
		
		* TXD_L2
		qui gen qaly_txd_L2 = (TXD_L2 / 12) * `uTXD_L2' if TXD_L2 != .
		qui replace qaly_txd_L2 = 0 if qaly_txd_L2 == .
		
		* Post L2
		qui gen qaly_post_L2 = ((OC_TIME_L - TSD_L2E_ref) / 12) * `uPostL2' if OC_TIME_L > TSD_L2E_ref
		qui replace qaly_post_L2 = 0 if qaly_post_L2 == .
		
		* Total undiscounted QALYs
		qui gen qaly_total = qaly_txd_L2 + qaly_post_L2
		
		* Discounting (TSD_L2S_ref = 0 by construction)
		
		* TXD_L2: from TSD_L2S_ref (=0) to TSD_L2E_ref
		qui gen qaly_txd_L2_d = qaly_txd_L2 * ((1 + $drate)^(-TSD_L2S_ref/12) - (1 + $drate)^(-TSD_L2E_ref/12)) / (`ln_r' * TXD_L2/12) if TXD_L2 > 0 & qaly_txd_L2 > 0
		qui replace qaly_txd_L2_d = 0 if qaly_txd_L2_d == .
		
		* Post: from TSD_L2E_ref to OC_TIME_L
		qui gen qaly_post_L2_d = qaly_post_L2 * ((1 + $drate)^(-TSD_L2E_ref/12) - (1 + $drate)^(-OC_TIME_L/12)) / (`ln_r' * (OC_TIME_L - TSD_L2E_ref)/12) if qaly_post_L2 > 0 & OC_TIME_L > TSD_L2E_ref
		qui replace qaly_post_L2_d = 0 if qaly_post_L2_d == .
		
		* Total discounted QALYs
		qui gen qaly_total_d = qaly_txd_L2_d + qaly_post_L2_d
	}
* L3+ analysis
	else {
		* L >= 3: Everything from line L start uses uPostL2
		local uPostL2 = 0.63
		
		* Total time from start of line L (TSD_L`L'S_ref = 0) to outcome
		qui gen qaly_total = (OC_TIME_L / 12) * `uPostL2' if OC_TIME_L > 0 & OC_TIME_L != .
		qui replace qaly_total = 0 if qaly_total == .
		
		* Discounting: from TSD_L`L'S_ref (=0) to OC_TIME_L
		qui gen qaly_total_d = qaly_total * ((1 + $drate)^(-TSD_L`L'S_ref/12) - (1 + $drate)^(-OC_TIME_L/12)) / (`ln_r' * OC_TIME_L/12) if qaly_total > 0 & OC_TIME_L > 0
		qui replace qaly_total_d = 0 if qaly_total_d == .
	}

end
