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
		mPat = (vAge, vAge2, vMale, 
				   vECOG0, vECOG1, vECOG2, 
				   vRISS1, vRISS2, vRISS3, 
				   vSCT_DN,
				   vCons)
		
		// Extract coefficients
		nPredictors = cols(mPat)
		vCoef= bDN_TFI[1,1..nPredictors]'
		dist = fbDN_TFI
		aux = bDN_TFI[1, cols(bDN_TFI)]
			
		// Calculate XB
		vXB = mPat * vCoef
			
		// Generate random numbers
		vRN = runiform(rows(incident), 1)
			
		// Calculate outcome (survival time)
		vOut = calcSurvTime(vXB, vRN, dist, aux)
	}	
		
	// Grab prevalent patient data
	prevalent = selectindex(mState[., 1] :> OMC + 1)
	if (rows(prevalent) > 0) {
		vOut[prevalent] = mTNE[prevalent, OMC]
	}
		
	// Update matrices	
	mTFI[., LX+1] = round(vOut, 0.1)                                
	mTNE[., OMC] = round(vOut, 0.1)
	mTSD[., OMC+1] = mTSD[., OMC] + mTNE[., OMC]
}
