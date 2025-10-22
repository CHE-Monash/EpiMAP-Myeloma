**********
*SIM TFI LX - Vectorised Implementation
*
* Purpose: Treatment-free Interval at Line 2+ End (time from LXE to L(X+1)S)
**********

mata {
	// Filters
	vEligible = (mMOR[.,OMC-1] :== 0) :& (mState[.,1] :<= OMC+1)
	idxEligible = selectindex(vEligible)
	
	if (rows(idxEligible) > 0) {
		// Generate random numbers for all patients
		vRN = runiform(rows(mMOR), 1)
		
		// Initialize outcome vector
		vOutcome = J(rows(mMOR), 1, .)
		
		// Select coefficient matrix based on Line
		if (Line == 2) {
			bCoef = bL2_CI
			fbCoef = fbL2_CI
			maxTFI = maxL2_CI
		}
		else if (Line == 3) {
			bCoef = bL3_CI
			fbCoef = fbL3_CI
			maxTFI = maxL3_CI
		}
		else if (Line == 4) {
			bCoef = bL4_CI
			fbCoef = fbL4_CI
			maxTFI = maxL4_CI
		}
		else if (Line >= 5) {
			bCoef = bLX_CI
			fbCoef = fbLX_CI
			maxTFI = maxLX_CI
		}
		
		// Patient vectors
		vAge_e = mAge[idxEligible, OMC]
		vAge2_e = vAge_e :^ 2
		vMale_e = vMale[idxEligible]
		vECOG0_e = (vECOG[idxEligible] :== 0)
		vECOG1_e = (vECOG[idxEligible] :== 1)
		vECOG2_e = (vECOG[idxEligible] :== 2)
		vRISS1_e = (vRISS[idxEligible] :== 1)
		vRISS2_e = (vRISS[idxEligible] :== 2)
		vRISS3_e = (vRISS[idxEligible] :== 3)
		vBCR1_e = (mBCR[idxEligible, Line] :== 1)
		vBCR2_e = (mBCR[idxEligible, Line] :== 2)
		vBCR3_e = (mBCR[idxEligible, Line] :== 3)
		vBCR4_e = (mBCR[idxEligible, Line] :== 4)
		vBCR5_e = (mBCR[idxEligible, Line] :== 5)
		vBCR6_e = (mBCR[idxEligible, Line] :== 6)
		vCons_e = vCons[idxEligible]
		
		// Patient matrix
		pMatrix = (vAge_e, vAge2_e, vMale_e, 
				   vECOG0_e, vECOG1_e, vECOG2_e,
		           vRISS1_e, vRISS2_e, vRISS3_e, 
				   vBCR1_e, vBCR2_e, vBCR3_e, vBCR4_e, vBCR5_e, vBCR6_e, 
				   vCons_e)
		
		// Extract coefficients
		nPredictors = cols(pMatrix)
		coef = bCoef[1, 1..nPredictors]'
		
		// Calculate XB
		vXB = pMatrix * coef
		
		// Calculate survival time
		vRN_e = vRN[idxEligible]
		vOutcome[idxEligible] = calcSurvTime(vXB, vRN_e, fbCoef, bCoef[1, cols(bCoef)])
		
		// Curtail if beyond maximum observed
		vOutcome[idxEligible] = rowmin((vOutcome[idxEligible], J(rows(idxEligible), 1, maxTFI)))
		
		// Handle prevalent patients
		prevalent = selectindex(mState[., 1] :> OMC + 1)
		if (rows(prevalent) > 0) {
			vOutcome[prevalent] = mTNE[prevalent, OMC] :* 365.25
		}
		
		// Update matrices
		mTNE[., OMC] = vOutcome :/ 365.25
		mTSD[., OMC+1] = mTSD[., OMC] + mTNE[., OMC]
		mTFI[., LX+1] = vOutcome
	}
}
