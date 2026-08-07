**********
* Monash Myeloma Model - Sim maintenance len-refractory (L1)
*
* Purpose: set the LATCHED lenalidomide-refractory state for patients who progress on, or within 60
*          days of, their L1 lenalidomide MAINTENANCE. Writes into the same vLenRefr_in the
*          treatment arm uses - one flag, not two (os_lenrefr_check.do settles that at p = 0.97).
*
* Notes:   FITTED LOGIT, drawn once per patient. Replaces a deterministic tail rule that flagged a
*          patient when the simulated tail (TFI_L1 - MND_L1) fell under 60 days. That rule was a
*          correct DEFINITION but a broken MECHANISM: sim_tfi_l1.do draws the gap truncated just
*          above the maintenance duration, so the tail piles up at the bound for LONG maintenance -
*          i.e. it selected long-maintenance, GOOD-prognosis patients. Prevalence came out at 3.7%
*          against a registry 15.5% at L2, and the flagged group survived BETTER than the registry's
*          despite the OS penalty being applied correctly (scratch/refractory/mnt_refr_logit.log).
*
*          The logit decouples the flag from that artefact: covariates set WHO is flagged, the
*          intercept sets HOW MANY. Fitted on the registry's determinate maintenance episodes
*          (28.6% base rate) with i.BCR_L1 collapsed to CR/VG/PR/poor, pooled across transplant with
*          an SCT main effect. Response is monotone in the right direction - worse L1 response means
*          higher odds - and transplant is protective. AUC 0.708.
*
*          ONE DRAW, not per line: maintenance refractoriness is determined once, at the end of L1
*          maintenance, and the flag is latched thereafter. Consumes rn_mntrefr() (one CRN column).
*
* ORDER:   AFTER sim_mnr.do (needs vMNR) and sim_bcr_asct.do; at L1E, after that line's OS, for the
*          same reason the treatment arm runs there - the flag applies from L2 onward, which is
*          satisfied because L2's equations run at a later OMC.
* Reads:   vMNT, vMNR, mBCR[.,1] (BCR_L1), vSCT_L1, the baseline and comorbidity vectors,
*          bLENREFR_MNT, rn_mntrefr().
* Writes:  vLenRefr_in (latched 0 -> 1).
*
*          Design MUST match the fit in prep/risk_equations.do:
*              logit MNTREFR Age Age2 Male i.ECOGcc i.RISS CM_CKD CM_CRD CM_PLM CM_DBT SCT
*                            i.bcr_grp_l1
*          e(b) order: Age Age2 Male, ECOG(0,1,2), RISS(1,2,3), CM x4, SCT, BCR(1,2,3,4), _cons.
*          Base-level dummies carry a 0 coefficient and are included as columns, as elsewhere.
**********

mata {
	vCoefM = get_mntrefr_coef()

	if (cols(vCoefM) > 0) {

		// Alive, reached this point, on LENALIDOMIDE maintenance, not already refractory.
		// sim_mnr never draws 'other', so this is len against thal.
		idx = selectindex((mMOR[., OMC-1] :== 0) :& (mState[., 1] :<= OMC) ///
		                  :& (vMNT :== 1) :& (vMNR :== 1) :& (vLenRefr_in :== 0))

		if (rows(idx) > 0) {

			// Response collapsed to CR / VG / PR / poor, matching bcr_grp_l1 in the fit. The
			// collapse is not cosmetic: SD was a perfect predictor on 5 observations and was dropped
			// from the uncollapsed fit, while sim_bcr.do CAN draw SD and PD - so without it those
			// patients would fall through every dummy to the CR base level, the BEST response.
			vB  = mBCR[idx, 1]
			vB1 = (vB :== 1)
			vB2 = (vB :== 2)
			vB3 = (vB :== 3)
			vB4 = (vB :>= 4)

			mPat = (vAge[idx], vAge2[idx], vMale[idx],
			        vECOG0[idx], vECOG1[idx], vECOG2[idx],
			        vRISS1[idx], vRISS2[idx], vRISS3[idx],
			        vCKD[idx], vCRD[idx], vPLM[idx], vDBT[idx],
			        vSCT_L1[idx],
			        vB1, vB2, vB3, vB4,
			        vCons[idx])

			// Guard: design columns must equal the coefficient count (plain logit, no ancillary).
			if (cols(mPat) != cols(vCoefM)) {
				errprintf("sim_mnt_refr: design/coefficient mismatch - mPat has %g columns but bLENREFR_MNT has %g. The fit implies %g BCR levels against the %g assumed here (design is Age Age2 Male + ECOG 3 + RISS 3 + CM 4 + SCT + BCR + _cons).\n",
					cols(mPat), cols(vCoefM), cols(vCoefM) - 15, 4)
				exit(459)
			}

			vXB = mPat * vCoefM'
			vPR = 1 :/ (1 :+ exp(-vXB))

			// One CRN draw per patient; flip where p > u. Latched: 1 stays 1.
			vRN = rnDraw(idx, rn_mntrefr())
			vLenRefr_in[idx] = vLenRefr_in[idx] :| (vPR :> vRN)
		}
	}
	else {
		// No model - the maintenance arm is a no-op rather than an error, matching sim_lenrefr.do's
		// behaviour when its own logit is absent. The treatment arm still contributes.
		errprintf("sim_mnt_refr: bLENREFR_MNT not found - maintenance refractoriness will not be generated. Re-run prep/risk_equations.do.\n")
	}
}

// Check for override file, execute if it exists
local override_file "${outcomes_path}/sim_mnt_refr_override.do"
capture confirm file "`override_file'"
if _rc == 0 {
	qui do `override_file'
}
