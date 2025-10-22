**********
* SIM TFI DN - Vectorised Implementation
* 
* Purpose: Calculate Treatment-free Interval at Diagnosis (time from diagnosis to Line 1 start)
* Method: Parametric survival analysis
* Outcome: Continuous survival time (in days)
**********
	
mata {
	// Initialize outcome
	vOutcome = J(st_nobs(), 1, .)
		
	// Filter for incident patients
	incident = selectindex(mState[., 1] :<= OMC + 1)
	
	// Calculate for incident patients
	if (rows(incident) > 0) {
		
		// Assemble patient matrix
		mPatient = (vAge, vAge2, vMale, 
				   vECOG0, vECOG1, vECOG2, 
				   vRISS1, vRISS2, vRISS3, 
				   vSCT_DN,
				   vCons)
		
		// Extract coefficients
		nPredictors = cols(mPatient)
		vCoef= bDN_CI[1,1..nPredictors]'
		aux = bDN_CI[1, cols(bDN_CI)]
			
		// Calculate XB
		vXB = mPatient * vCoef
			
		// Generate random numbers
		vRN = runiform(rows(incident), 1)
			
		// Calculate survival time using utility function
		vOutcome = calcSurvTime(vXB, vRN, fbDN_CI, aux)
	}	
		
	// Grab prevalent patient data
	prevalent = selectindex(mState[., 1] :> OMC + 1)
	if (rows(prevalent) > 0) {
		vOutcome[prevalent] = mTNE[prevalent, OMC] :* 365.25  // Convert years to days
	}
		
	// Update matrices	
	mTNE[., OMC] = vOutcome :/ 365.25                       // Convert days to years
	mTFI[., LX+1] = vOutcome                                // Store in days
	mTSD[., OMC+1] = mTSD[., OMC] + mTNE[., OMC]          // Cumulative time
}
