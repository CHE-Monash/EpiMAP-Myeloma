**********
*SIM TXD (L2+)
* Purpose: Calculate treatment duration for L2 onwards
* Method: Parametric survival analysis
* Outcome: Continuous time (days)
**********
mata {
	vEligible = (mMOR[.,OMC-1] :== 0) :& (mState[.,1] :<= OMC+1)
	idxEligible = selectindex(vEligible)
	
	if (rows(idxEligible) > 0) {
		
		// Determine model stage and BCR structure
		if (Line == 1) {
			vCoef = bL2_TXD
			dist = fbL2_TXD
			vTXR = oL2_TXR
			BCR_cat = 6 
		}
		else if (Line == 2) {
			vCoef = bL3_TXD
			dist = fbL3_TXD
			vTXR = oL3_TXR
			BCR_cat = 6 
		}
		else if (Line == 3) {
			vCoef = bL4_TXD
			dist = fbL4_TXD
			vTXR = oL4_TXR
			BCR_cat = 3  
		}
		else if (Line >= 4) {
			vCoef = bLX_TXD
			dist = fbLX_TXD
			vTXR = J(1, 0, .)
			BCR_cat = 3  
		}
		
		nTXR = cols(vTXR)
		vRN = runiform(rows(mMOR), 1)
		vOut = J(rows(mMOR), 1, .)
		
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
		vCons_e = vCons[idxEligible]
		
		// Build patient matrix
		mPat = (vAge_e, vAge2_e, vMale_e, 
				vECOG0_e, vECOG1_e, vECOG2_e, 
				vRISS1_e, vRISS2_e, vRISS3_e)
		
		// Add treatment regimen dummies
		currentTX = mTXR[idxEligible, Line]
		if (nTXR >= 1) {
			vTXR1_e = (currentTX :== vTXR[1,1])
			mPat = (mPat, vTXR1_e)
		}
		if (nTXR >= 2) {
			vTXR2_e = (currentTX :== vTXR[1,2])
			mPat = (mPat, vTXR2_e)
		}
		if (nTXR >= 3) {
			vTXR3_e = (currentTX :== vTXR[1,3])
			mPat = (mPat, vTXR3_e)
		}
		
		// Add previous BCR
		if (BCR_cat == 6) {
			prevBCR = mBCR[idxEligible, Line]
			vBCR_CR_e = (prevBCR :== 1)
			vBCR_VG_e = (prevBCR :== 2)
			vBCR_PR_e = (prevBCR :== 3)
			vBCR_MR_e = (prevBCR :== 4)
			vBCR_SD_e = (prevBCR :== 5)
			vBCR_PD_e = (prevBCR :== 6)
			mPat = (mPat, vBCR_CR_e, vBCR_VG_e, vBCR_PR_e, vBCR_MR_e, vBCR_SD_e, vBCR_PD_e)
		}
		else if (BCR_cat == 3) {
			prevBCR = mBCR[idxEligible, Line]
			vBCR_CR_e = (prevBCR :== 1)
			vBCR_PR_e = (prevBCR :== 3)
			vBCR_SD_e = (prevBCR :== 5)
			mPat = (mPat, vBCR_CR_e, vBCR_PR_e, vBCR_SD_e)
		}
		
		// Add constant
		mPat = (mPat, vCons_e)
		
		// Extract coefficients
		aux = vCoef[1, cols(vCoef)]
		nPredictors = cols(mPat)
		vCoef = vCoef[1, 1..nPredictors]'
		
		// Calculate XB and duration
		vXB = mPat * vCoef
		vRN_e = vRN[idxEligible]
		vOut[idxEligible] = calcSurvTime(vXB, vRN_e, dist, aux)
		
		// Store outcomes
		mTXD[idxEligible, LX+1] = round(vOut[idxEligible], 0.1)
		mTNE[idxEligible, OMC] = round(vOut[idxEligible], 0.1)
		mTSD[idxEligible, OMC+1] = mTSD[idxEligible, OMC-1] :+ mTNE[idxEligible, OMC]
	}
}
