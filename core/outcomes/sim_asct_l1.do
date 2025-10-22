**********
* SIM SCT L1 - Vectorised Implementation
* 
* Purpose: Determine ASCT eligibility at Line 1 End
* Outcome: Binary (0 = Not eligible, 1 = Eligible)
* Additional filters: CR != 7 AND BCR != 6
**********
	
mata {
	// Initialise outcome
	oSCT_L1 = J(st_nobs(), 1, .)
	
	// Filter for valid patients: Alive & State filters
	idx = selectindex((mMOR[., OMC-1] :== 0) :& (mState[., 1] :<= OMC + 1))
		
	// Calculate for valid patients
	if (rows(idx) > 0) {
	
		// Extract BCR from matrix (Line 1 BCR is in column 1 of mBCR)
		vBCR_L1 = mBCR[., 1]
		
		// Create BCR dummy variables
		vBCR2 = (vBCR_L1 :== 2)
		vBCR3 = (vBCR_L1 :== 3)
		vBCR4 = (vBCR_L1 :== 4)
		vBCR5 = (vBCR_L1 :== 5)
		
		// Assemble patient matrix
		// Note: Coefficient order matches bL1_SCT structure
		pSCT_L1 = (vAge, vAge2, vMale, vECOG1, vECOG2, vRISS2, vRISS3, vBCR2, vBCR3, vBCR4, vBCR5, vAge70, vAge75, vCMc1, vCMc2, vCMc3, vCons)

		// Extract coefficients
		nPredictors = cols(pSCT_L1)
		coefSCT_L1 = bL1_SCT[1, 1..nPredictors]'

		// Calculate XB
		xbSCT_L1 = pSCT_L1 * coefSCT_L1
		
		// Calculate probabilities
		prSCT_L1 = 1 :/ (1 :+ exp(-xbSCT_L1))
		
		// Generate random numbers
		rnSCT_L1 = runiform(rows(pSCT_L1), 1)
		
		// Extract CR from matrix (Line 1 CR is in column 1 of mTXR)
		vCR_L1 = mTXR[., 1]
			
		// Determine outcome with additional filters
		// Eligible if: probability > random number AND CR != 7 AND BCR != 6
		oSCT_L1 = (prSCT_L1 :> rnSCT_L1) :& (vCR_L1 :!= 7) :& (vBCR_L1 :!= 6)
		
		// Convert logical to numeric (0/1)
		oSCT_L1 = oSCT_L1 :* 1
	}
		
	// Update matrices 
	vSCT_L1 = oSCT_L1
	mSCT[., 2] = oSCT_L1           // Column 2 = SCT eligibility at L1E
	mCore[., cSCT] = oSCT_L1       // Update mCore for backwards compatibility
}
