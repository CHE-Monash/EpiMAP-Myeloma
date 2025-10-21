**********
* SIM OS DN - Vectorised Implementation
* 
* Purpose: Calculate Overall Survival at Diagnosis
* Method: Parametric survival model (Weibull) using vectorised matrix multiplication  
* Outcome: Continuous time (in years)
**********
	
mata {
	
	// Extract current BCR from mCore for TX#BCR interaction
	// At diagnosis, BCR may be 0 or missing - these dummies handle that
	currentBCR = (Line == 0) ? J(nObs, 1, 5) : mBCR[., Line]
	bcr1 = (currentBCR :== 1)
	bcr2 = (currentBCR :== 2)
	bcr3 = (currentBCR :== 3)
	bcr4 = (currentBCR :== 4)
	bcr5 = (currentBCR :== 5)
	bcr6 = (currentBCR :== 6)
	
	// Assemble patient matrix
	// Variables: Age, Age2, Male, ECOG1, ECOG2, RISS2, RISS3, BCR1-6, Constant
	// Note: At Line 0 (diagnosis), only BCR dummies at positions 10-15 are used
	pOS_DN = (vAge, vAge2, vMale, vECOG1, vECOG2, vRISS2, vRISS3, bcr1, bcr2, bcr3, bcr4, bcr5, bcr6, vCons)
	
	// Extract coefficients
	
		coefOS_DN = bOS[1, (1,2,3,5,6,8,9,10,11,12,13,14,15,58)]'
	}
	// Additional logic for other lines would go here when vectorising those

	// Calculate XB
	xbOS_DN = pOS_DN * coefOS_DN
		
	// Generate random numbers
	rnOS_DN = runiform(rows(pOS_DN), 1)
		
	// Calculate survival time using utility function
	oOS_DN = calcSurvTime(xbOS_DN, rnOS_DN, distribution, bOS[1, cols(bOS)])
		
	// State filter (no-one dies before DN, no prevalent patients at diagnosis)
	incident = selectindex(mState[., 1] :<= OMC + 1)
	if (rows(incident) > 0) {
		// Incident patients keep their calculated values
		// (oOS_DN already calculated for all, just don't overwrite)
	}
	
	// Set any patients outside valid state to missing (shouldn't happen at DN)
	invalid = selectindex(mState[., 1] :> OMC + 1)
	if (rows(invalid) > 0) {
		oOS_DN[invalid] = .
	}
		
	// Update mOS matrix (Overall Survival times in years)
	mOS[., OMC] = oOS_DN
}
