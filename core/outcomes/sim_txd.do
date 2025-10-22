**********
*SIM TXD LX Vector

*Purpose: Calculate treatment duration for L2 onwards
*Method: Parametric survival analysis
*Outcome: Continuous time (days)
**********

mata {
	vEligible = (mMOR[.,OMC-1] :== 0) :& (mState[.,1] :<= OMC+1)
	idxEligible = selectindex(vEligible)
	
	if (rows(idxEligible) > 0) {
		
		// Determine model stage
		if (Line == 1) {
			mCoef = bL2_CD
			dist = fbL2_CD
		}
		else if (Line == 2) {
			mCoef = bL3_CD
			dist = fbL3_CD
		}
		else if (Line == 3) {
			mCoef = bL4_CD
			dist = fbL4_CD
		}
		else if (Line >= 4) {
			mCoef = bLX_CD
			dist = fbLX_CD
		}
		
		vRN = runiform(rows(mMOR), 1)
		vOutcome = J(rows(mMOR), 1, .)
		
		// Patient characteristics
		vAge_e = mAge[idxEligible, OMC]
		vAge2_e = vAge_e :^ 2
		vMale_e = vMale[idxEligible]
		vECOG0_e = (vECOG[idxEligible] :== 0)
		vECOG1_e = (vECOG[idxEligible] :== 1)
		vECOG2_e = (vECOG[idxEligible] :== 2)
		vRISS1_e = (vRISS[idxEligible] :== 1)		
		vRISS2_e = (vRISS[idxEligible] :== 2)
		vRISS3_e = (vRISS[idxEligible] :== 3)
		vSCT_L1_e = vSCT_L1[idxEligible]
		vCons_e = vCons[idxEligible]
		
		// Build patient matrix
		mPat = (vAge_e, vAge2_e, vMale_e, 
				   vECOG0_e, vECOG1_e, vECOG2_e, 
				   vRISS1_e, vRISS2_e, vRISS3_e, 
				   vSCT_L1_e)
		
		// Previous BCR
		if (Line == 1 | Line == 2) { // Calculating for Line 2 or Line 3
			prevBCR = mBCR[idxEligible, 1]
			vBCR_CR_e = (prevBCR :== 1)
			vBCR_VG_e = (prevBCR :== 2)
			vBCR_PR_e = (prevBCR :== 3)
			vBCR_MR_e = (prevBCR :== 4)
			vBCR_SD_e = (prevBCR :== 5)
			vBCR_PD_e = (prevBCR :== 6)
			mPat = (mPat, vBCR_CR_e, vBCR_VG_e, vBCR_PR_e, vBCR_MR_e, vBCR_SD_e, vBCR_PD_e)
		}
		if (Line >= 3) {
			prevBCR = mBCR[idxEligible, Line-1]
			vBCR_CR_e = (prevBCR :== 1)
			vBCR_PR_e = (prevBCR :== 3)
			vBCR_SD_e = (prevBCR :== 5)
			mPat = (mPat, vBCR_CR_e, vBCR_PR_e, vBCR_SD_e)
		}
		
		// Add Treatment regimen dummies
		currentTX = mTXR[idxEligible, Line]
		
		if (Line == 1) oVector = oL2_CR
		else if (Line == 2) oVector = oL3_CR
		else if (Line == 3) oVector = oL4_CR
		else oVector = J(1, 0, .)
		
		if (cols(oVector) >= 1) {
			vTXR1_e = (currentTX :== oVector[1,1])
			mPat = (mPat, vTXR1_e)
		}
		if (cols(oVector) >= 2) {
			vTXR2_e = (currentTX :== oVector[1,2])
			mPat = (mPat, vTXR2_e)
		}
		if (cols(oVector) >= 3) {
			vTXR3_e = (currentTX :== oVector[1,3])
			mPat = (mPat, vTXR3_e)
		}
		
		mPat = (mPat, vCons_e)
		
		// Extract coefficients
		nPredictors = cols(mPat)
		vCoef = mCoef[1, 1..nPredictors]'
		aux = mCoef[1, cols(mCoef)]
		
		// Calculate XB and duration
		vXB = mPat * vCoef
		vRN_e = vRN[idxEligible]
		vOutcome[idxEligible] = calcSurvTime(vXB, vRN_e, dist, aux)
		
		// Store outcomes
		mTXD[idxEligible, Line] = vOutcome[idxEligible]
		mTNE[idxEligible, OMC] = vOutcome[idxEligible]
		mTSD[idxEligible, OMC] = mTSD[idxEligible, OMC-1] :+ vOutcome[idxEligible]
	}
}
