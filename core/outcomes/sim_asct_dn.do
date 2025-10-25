**********
* SIM SCT DN
* 
* Purpose: Determine ASCT eligibility at diagnosis
* Method: Logistic regression
* Outcome: Binary (0 = Not eligible, 1 = Eligible)
**********
	
mata {
	// Initialize outcome
	vOC = J(st_nobs(), 1, .)
	
	// Filter for valid patients
	idx = selectindex(mState[., 1] :<= OMC + 1)
		
	// Calculate for valid patients
	if (rows(idx) > 0) {
	
		// Assemble patient matrix
		mPat = (vAge, vAge2, vMale, 
				vECOG0, vECOG1, vECOG2, 
				vRISS1, vRISS2, vRISS3, 
				vAge70, vAge75, 
				vCMc0, vCMc1, vCMc2, vCMc3,
				vCons)

		// Extract coefficients
		nPredictors = cols(mPat)
		vCoef = bDN_SCT[1, 1..nPredictors]'

		// Calculate XB
		vXB = mPat * vCoef
		
		// Calculate probabilities
		vPR = 1 :/ (1 :+ exp(-vXB))
		
		// Generate random numbers
		vRN = runiform(rows(idx), 1)
			
		// Determine outcome 
		vOC = (vPR :> vRN)
	}
		
	// Update matrices 
	vSCT_DN = vOC 
	mSCT[., 1] = vOC
} 
	