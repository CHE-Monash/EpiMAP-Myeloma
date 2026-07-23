**********
* Monash Myeloma Model - Sim MND (L1 maintenance duration)
*
* Purpose: Draw L1 maintenance DURATION by parametric survival, among patients on maintenance
*          (MNT == 1). Continuous time in months.
*
* Notes:   THREE arms: lenalidomide with ASCT, lenalidomide without, and thalidomide pooled across
*          transplant. Lenalidomide is split so each arm can carry the response depth that exists
*          for it - i.BCR_SCT after transplant, i.BCR_L1 otherwise. The signal is entirely in the
*          transplant response (LR p = 0.0021 against 0.50 for BCR_L1), and a pooled equation
*          cannot carry it: BCR_SCT does not exist for never-transplanted patients, and the
*          registry's 0 code cannot stand in for them because the engine never draws 0 either.
*          Mirrors sim_tfi_l1.do, which splits the same way for the same reason.
*
*          No ln(TFI) covariate: it restricted the fit to patients with an observed L2. The
*          ordering (maintenance must fit inside the gap) is enforced downstream instead -
*          sim_tfi_l1.do draws the gap truncated below at the duration drawn here. Thalidomide is
*          capped at 18 months, matching the censoring in its fit; the lenalidomide ceilings are
*          the observed maxima and are NOT to be lowered (scratch/maintenance/_notes.md).
*
* ORDER:   AFTER sim_mnr.do (needs vMNR) and sim_bcr_asct.do (needs mBCR[.,10] = BCR_SCT).
*          BEFORE sim_tfi_l1.do, which depends on vMND. This is a REVERSAL of the previous order
*          and the whole point of the design: maintenance first, then the gap that must contain it.
*
*          Must match the fits in prep/risk_equations.do:
*              len ASCT    streg Age Age2 Male i.ECOGcc i.RISS i.BCR_SCT
*              len NoASCT  streg Age Age2 Male i.ECOGcc i.RISS i.BCR_L1
*              thal        streg Age Age2 Male i.ECOGcc i.RISS SCT        (exit at 18 months)
*          e(b) order follows the command: Age Age2 Male, ECOG(0,1,2), RISS(1,2,3), then the BCR
*          levels, then _cons, then aux. Base-level dummies carry a 0 coefficient in e(b) and are
*          included as columns, exactly as sim_tfi_l1.do does.
*
*          BCR LEVELS, set by what the ENGINE can draw, not by what the registry codes:
*              ASCT arm    BCR_SCT 1-4   sim_bcr_asct.do categoryValues = (1,2,3,4)
*              NoASCT arm  BCR_L1  1-6   sim_bcr.do      categoryValues = (1,2,3,4,5,6)
*          The registry also codes BCR_SCT == 0 ("transplanted, response not recorded", 22.7% of
*          transplanted maintenance records) but the engine never assigns it, so the fit excludes
*          it - otherwise it would become the base category and every simulated patient would be
*          measured against a group that does not exist. If a level is empty in the fit, e(b) is
*          short and the guard below fires rather than silently misaligning the design.
**********

mata {
	vCoefA = get_mnd_coef_len_asct()
	vCoefN = get_mnd_coef_len_noasct()
	vCoefT = get_mnd_coef_thal()

	if (cols(vCoefA) > 0 | cols(vCoefN) > 0 | cols(vCoefT) > 0) {

		// Alive, eligible AND receiving maintenance. Same population sim_mnr drew for.
		idx = selectindex((mMOR[., OMC-1] :== 0) :& (mState[., 1] :<= OMC) :& (vMNT :== 1))
		if (rows(idx) > 0) {

			// ---- Lenalidomide, ASCT ----
			if (cols(vCoefA) > 0) {
				iA = idx[selectindex((vMNR[idx] :== 1) :& (vSCT_L1[idx] :== 1))]
				if (rows(iA) > 0) {
					// BCR_SCT levels 1-4 only. sim_bcr_asct.do draws categoryValues = (1,2,3,4),
					// so 0 ("transplanted, response not recorded" in the registry) is unreachable
					// here and is excluded from the fit to match. Same as sim_tfi_l1.do's ASCT arm.
					vB1 = (mBCR[iA, 10] :== 1)
					vB2 = (mBCR[iA, 10] :== 2)
					vB3 = (mBCR[iA, 10] :== 3)
					vB4 = (mBCR[iA, 10] :== 4)

					mPatA = (vAge[iA], vAge2[iA], vMale[iA],
							 vECOG0[iA], vECOG1[iA], vECOG2[iA],
							 vRISS1[iA], vRISS2[iA], vRISS3[iA],
							 vB1, vB2, vB3, vB4,
							 vCons[iA])
					nPredA = cols(mPatA)

					if (cols(vCoefA) != nPredA + 1) {
						errprintf("sim_mnd (len ASCT): design/coefficient mismatch - mPat has %g columns so %g were expected (mean + ancillary), but bL1_MND_LEN_ASCT has %g. The fit implies %g BCR_SCT levels against the %g assumed here (design is Age Age2 Male + ECOG 3 + RISS 3 + BCR + _cons).\n",
							nPredA, nPredA + 1, cols(vCoefA), cols(vCoefA) - 11, 4)
						exit(459)
					}

					vBetaA = vCoefA[1, 1..nPredA]'
					auxA   = vCoefA[1, cols(vCoefA)]
					vXBa   = mPatA * vBetaA
					vRNa   = rnDraw(iA, rn_mnd())
					vOCa   = calcSurvTime(vXBa, vRNa, fbL1_MND_LEN_ASCT, auxA)
					vMND[iA] = rowmin((vOCa, J(rows(iA), 1, maxL1_MND_LEN_ASCT)))
				}
			}

			// ---- Lenalidomide, no ASCT ----
			if (cols(vCoefN) > 0) {
				iN = idx[selectindex((vMNR[idx] :== 1) :& (vSCT_L1[idx] :!= 1))]
				if (rows(iN) > 0) {
					vC1 = (mBCR[iN, 1] :== 1)
					vC2 = (mBCR[iN, 1] :== 2)
					vC3 = (mBCR[iN, 1] :== 3)
					vC4 = (mBCR[iN, 1] :== 4)
					vC5 = (mBCR[iN, 1] :== 5)
					vC6 = (mBCR[iN, 1] :== 6)

					mPatN = (vAge[iN], vAge2[iN], vMale[iN],
							 vECOG0[iN], vECOG1[iN], vECOG2[iN],
							 vRISS1[iN], vRISS2[iN], vRISS3[iN],
							 vC1, vC2, vC3, vC4, vC5, vC6,
							 vCons[iN])
					nPredN = cols(mPatN)

					if (cols(vCoefN) != nPredN + 1) {
						errprintf("sim_mnd (len NoASCT): design/coefficient mismatch - mPat has %g columns so %g were expected (mean + ancillary), but bL1_MND_LEN_NoASCT has %g. The fit implies %g BCR_L1 levels against the %g assumed here (design is Age Age2 Male + ECOG 3 + RISS 3 + BCR + _cons).\n",
							nPredN, nPredN + 1, cols(vCoefN), cols(vCoefN) - 11, 6)
						exit(459)
					}

					vBetaN = vCoefN[1, 1..nPredN]'
					auxN   = vCoefN[1, cols(vCoefN)]
					vXBn   = mPatN * vBetaN
					vRNn   = rnDraw(iN, rn_mnd())
					vOCn   = calcSurvTime(vXBn, vRNn, fbL1_MND_LEN_NoASCT, auxN)
					vMND[iN] = rowmin((vOCn, J(rows(iN), 1, maxL1_MND_LEN_NoASCT)))
				}
			}

			// ---- Thalidomide, pooled across transplant ----
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
		errprintf("sim_mnd: no L1_MND coefficients found - maintenance duration will not be costed. Re-run prep/risk_equations.do.\n")
	}
}

// Check for override file, execute if it exists
local override_file "${outcomes_path}/sim_mnd_override.do"
capture confirm file "`override_file'"
if _rc == 0 {
	qui do `override_file'
}
