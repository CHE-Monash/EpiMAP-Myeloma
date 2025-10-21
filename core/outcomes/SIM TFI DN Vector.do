**********
* SIM TFI DN - Vectorised Implementation
* 
* Purpose: Calculate Treatment-free Interval at Diagnosis (time from diagnosis to Line 1 start)
* Outcome: Continuous survival time (in days)
**********
	
mata {
	// Initialize outcome
	oTFI_DN = J(st_nobs(), 1, .)
		
	// Filter for incident patients
	incident = selectindex(mState[., 1] :<= OMC + 1)
	
	// Calculate for incident patients
	if (rows(incident) > 0) {
		
		// Assemble patient matrix
		pTFI_DN = (vAge, vAge2, vMale, vECOG1, vECOG2, vRISS2, vRISS3, vSCT_DN, vCons)
		
		// Extract coefficients
		coefTFI_DN = bDN_CI[1, (1,2,3,5,6,8,9,10,11)]'
			
		// Calculate XB
		xbTFI_DN = pTFI_DN * coefTFI_DN
			
		// Generate random numbers
		rnTFI_DN = runiform(rows(incident), 1)
			
		// Calculate survival time using utility function
		oTFI_DN = calcSurvTime(xbTFI_DN, rnTFI_DN, fbDN_CI, bDN_CI[1, cols(bDN_CI)])
	}	
		
	// Grab prevalent patient data
	prevalent = selectindex(mState[., 1] :> OMC + 1)
	if (rows(prevalent) > 0) {
		oTFI_DN[prevalent] = mTNE[prevalent, OMC] :* 365.25  // Convert years to days
	}
		
	// Update matrices	
	mTNE[., OMC] = oTFI_DN :/ 365.25                       // Convert days to years
	mTFI[., LX+1] = oTFI_DN                                // Store in days
	mTSD[., OMC+1] = mTSD[., OMC] + mTNE[., OMC]          // Cumulative time
}
