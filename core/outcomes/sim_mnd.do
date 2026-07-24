**********
* Monash Myeloma Model - Sim MND (L1 maintenance duration)
*
* Purpose: Draw L1 maintenance DURATION by parametric survival, among patients on maintenance
*          (MNT == 1). Continuous time in months.
*
* Notes:   Two arms, split by REGIMEN and pooled across transplant, with SCT as a covariate and no
*          BCR. No ln(TFI) covariate: it restricted the fit to patients with an observed L2. The
*          ordering (maintenance must fit inside the gap) is enforced downstream instead -
*          sim_tfi_l1.do draws the gap truncated below at the duration drawn here, so nothing is
*          clipped. Thalidomide is capped at 18 months, matching the censoring in its fit.
*          Reasoning and the rejected alternatives: prep/risk_equations.do and
*          scratch/maintenance/_notes.md.
*
* ORDER:   AFTER sim_mnr.do (needs vMNR) and sim_bcr_asct.do (needs mBCR).
*          BEFORE sim_tfi_l1.do, which now depends on vMND. This is a REVERSAL of the previous order
*          and the whole point of the design: maintenance first, then the gap that must contain it.
*
*          Must match the fits in prep/risk_equations.do:
*              len   streg Age Age2 Male i.ECOGcc i.RISS SCT
*              thal  streg Age Age2 Male i.ECOGcc i.RISS SCT   (exit at 18 months)
*          e(b) order: Age Age2 Male, ECOG(0,1,2), RISS(1,2,3), SCT, _cons, aux. The base-level
*          dummies are 0-coef in e(b), as sim_tfi_l1.do does.
**********

mata {
	vCoefL = get_mnd_coef_len()
	vCoefT = get_mnd_coef_thal()

	if (cols(vCoefL) > 0 | cols(vCoefT) > 0) {

		// Alive, eligible AND receiving maintenance. Same population sim_mnr drew for.
		idx = selectindex((mMOR[., OMC-1] :== 0) :& (mState[., 1] :<= OMC) :& (vMNT :== 1))
		if (rows(idx) > 0) {

			// ---- Lenalidomide ----
			if (cols(vCoefL) > 0) {
				iL = idx[selectindex(vMNR[idx] :== 1)]
				if (rows(iL) > 0) {
					mPatL = (vAge[iL], vAge2[iL], vMale[iL],
							 vECOG0[iL], vECOG1[iL], vECOG2[iL],
							 vRISS1[iL], vRISS2[iL], vRISS3[iL],
							 vSCT_L1[iL],
							 vCons[iL])
					nPredL = cols(mPatL)

					if (cols(vCoefL) != nPredL + 1) {
						errprintf("sim_mnd (len): design/coefficient mismatch - mPat has %g columns so %g were expected (mean + ancillary), but bL1_MND_LEN has %g. An ECOG/RISS level was likely empty in the fit.\n",
							nPredL, nPredL + 1, cols(vCoefL))
						exit(459)
					}

					vBetaL = vCoefL[1, 1..nPredL]'
					auxL   = vCoefL[1, cols(vCoefL)]
					vXBl   = mPatL * vBetaL
					vRNl   = rnDraw(iL, rn_mnd())
					vOCl   = calcSurvTime(vXBl, vRNl, fbL1_MND_LEN, auxL)
					vMND[iL] = rowmin((vOCl, J(rows(iL), 1, maxL1_MND_LEN)))
				}
			}

			// ---- Thalidomide ----
			if (cols(vCoefT) > 0) {
				iT = idx[selectindex(vMNR[idx] :== 5)]
				if (rows(iT) > 0) {
					mPatT = (vAge[iT], vAge2[iT], vMale[iT],
							 vECOG0[iT], vECOG1[iT], vECOG2[iT],
							 vRISS1[iT], vRISS2[iT], vRISS3[iT],
							 vSCT_L1[iT],
							 vCons[iT])
					nPredT = cols(mPatT)

					if (cols(vCoefT) != nPredT + 1) {
						errprintf("sim_mnd (thal): design/coefficient mismatch - mPat has %g columns so %g were expected (mean + ancillary), but bL1_MND_THAL has %g. An ECOG/RISS level was likely empty in the fit.\n",
							nPredT, nPredT + 1, cols(vCoefT))
						exit(459)
					}

					vBetaT = vCoefT[1, 1..nPredT]'
					auxT   = vCoefT[1, cols(vCoefT)]
					vXBt   = mPatT * vBetaT
					vRNt   = rnDraw(iT, rn_mnd())
					vOCt   = calcSurvTime(vXBt, vRNt, fbL1_MND_THAL, auxT)
					// maxL1_MND_THAL is set to 18 in risk_equations.do, overriding the observed
					// maximum, because the records beyond that point are the ones the fit has just
					// declared untrustworthy.
					vMND[iT] = rowmin((vOCt, J(rows(iT), 1, maxL1_MND_THAL)))
				}
			}
		}
	}
	else {
		// No model - leave vMND missing. process_data.do bills nothing where the duration is
		// missing, which is the safe direction: no maintenance cost beats silently reverting to
		// the old whole-gap bill. sim_tfi_l1.do also falls back to an untruncated draw.
		errprintf("sim_mnd: bL1_MND_LEN / bL1_MND_THAL not found - maintenance duration will not be costed. Re-run prep/risk_equations.do.\n")
	}
}

// Check for override file, execute if it exists
local override_file "${outcomes_path}/sim_mnd_override.do"
capture confirm file "`override_file'"
if _rc == 0 {
	qui do `override_file'
}
