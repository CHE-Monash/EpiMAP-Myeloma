**********
*SIM TFI (L2+) 
*
* Purpose: Treatment-free Interval at Line 2+ End (time from LXE to L(X+1)S)
* Method: Parametric survival analysis
* Outcome: Continous time (months)
**********

mata {
	// Filters
	vEligible = (mMOR[.,OMC-1] :== 0) :& (mState[.,1] :<= OMC+1)
	idx = selectindex(vEligible)
	
	if (rows(idx) > 0) {
		// Generate random numbers for all patients
		vRN = runiform(rows(mMOR), 1)
		
		// Initialize outcome vector
		vOutcome = J(rows(mMOR), 1, .)
		
		// Select coefficient matrix based on Line
		if (Line == 2) {
			vCoef = bL2_TFI
			fbCoef = fbL2_TFI
			maxTFI = maxL2_TFI
			BCR_cat = 6 
		}
		else if (Line == 3) {
			vCoef = bL3_TFI
			fbCoef = fbL3_TFI
			maxTFI = maxL3_TFI
			BCR_cat = 3 
		}
		else if (Line == 4) {
			vCoef = bL4_TFI
			fbCoef = fbL4_TFI
			maxTFI = maxL4_TFI
			BCR_cat = 3
		}
		else if (Line >= 5) {
			vCoef = bLX_TFI
			fbCoef = fbLX_TFI
			maxTFI = maxLX_TFI
			BCR_cat = 3 			
		}
		
		// Patient vectors
		vAge_e = mAge[idx, OMC]
		vAge2_e = vAge_e :^ 2
		vMale_e = vMale[idx]
		vECOG0_e = (vECOG[idx] :== 0)
		vECOG1_e = (vECOG[idx] :== 1)
		vECOG2_e = (vECOG[idx] :== 2)
		vRISS1_e = (vRISS[idx] :== 1)
		vRISS2_e = (vRISS[idx] :== 2)
		vRISS3_e = (vRISS[idx] :== 3)
		vBCR1_e = (mBCR[idx, Line] :== 1)
		vBCR2_e = (mBCR[idx, Line] :== 2)
		vBCR3_e = (mBCR[idx, Line] :== 3)
		vBCR4_e = (mBCR[idx, Line] :== 4)
		vBCR5_e = (mBCR[idx, Line] :== 5)
		vBCR6_e = (mBCR[idx, Line] :== 6)
		vCons_e = vCons[idx]
		
		// Build patient matrix
		mPat = (vAge_e, vAge2_e, vMale_e, 
				   vECOG0_e, vECOG1_e, vECOG2_e,
		           vRISS1_e, vRISS2_e, vRISS3_e)
		
		// Add previous BCR
		if (BCR_cat == 6) {
			prevBCR = mBCR[idx, Line]
			vBCR_CR_e = (prevBCR :== 1)
			vBCR_VG_e = (prevBCR :== 2)
			vBCR_PR_e = (prevBCR :== 3)
			vBCR_MR_e = (prevBCR :== 4)
			vBCR_SD_e = (prevBCR :== 5)
			vBCR_PD_e = (prevBCR :== 6)
			mPat = (mPat, vBCR_CR_e, vBCR_VG_e, vBCR_PR_e, vBCR_MR_e, vBCR_SD_e, vBCR_PD_e)
		}
		else if (BCR_cat == 3) {
			prevBCR = mBCR[idx, Line]
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
		
		// Calculate XB
		vXB = mPat * vCoef
		
		// Calculate outcome (survival time)
		vRN_e = vRN[idx]
		vOut[idx] = calcSurvTime(vXB, vRN_e, fbCoef, aux)
		
		// Curtail if beyond maximum observed
		vOut[idx] = rowmin((vOut[idx], J(rows(idx), 1, maxTFI)))
		
		// Handle prevalent patients
		prevalent = selectindex(mState[., 1] :> OMC + 1)
		if (rows(prevalent) > 0) {
			vOut[prevalent] = mTNE[prevalent, OMC]
		}
		
		// Update matrices
		mTFI[., LX+1] = round(vOut, 0.1)
		mTNE[., OMC] = round(vOut, 0.1)
		mTSD[., OMC+1] = mTSD[., OMC] + mTNE[., OMC]
	}
}
