**********
* SIM SCT DN - Vectorised Implementation
* 
* Purpose: Determine ASCT eligibility at diagnosis
* Outcome: Binary (0 = Not eligible, 1 = Eligible)
**********
	
mata {
	// Initialize outcome
	oSCT_DN = J(st_nobs(), 1, .)
	
	// Filter for valid patients
	idx = selectindex(mState[., 1] :<= OMC + 1)
		
	// Calculate for valid patients
	if (rows(idx) > 0) {
	
		// Assemble patient matrix
		pSCT_DN = (vAge, vAge2, vMale, vECOG1, vECOG2, vRISS2, vRISS3, vAge70, vAge75, vCMc1, vCMc2, vCMc3, vCons)

		// Extract coefficients
		coefSCT_DN = bDN_SCT[1, (1,2,3,5,6,8,9,10,11,13,14,15,16)]'

		// Calculate XB
		xbSCT_DN = pSCT_DN * coefSCT_DN
		
		// Calculate probabilities
		prSCT_DN = 1 :/ (1 :+ exp(-xbSCT_DN))
		
		// Generate random numbers
		rnSCT_DN = runiform(rows(pSCT_DN), 1)
			
		// Determine outcome 
		oSCT_DN = (prSCT_DN :> rnSCT_DN)
	}
		
	// Update matrices 
	vSCT_DN = oSCT_DN 
	mSCT[., 1] = oSCT_DN           // Column 1 = SCT eligibility at DN
	mCore[., cSCT] = oSCT_DN       // Update mCore for use in other outcomes
} 
	