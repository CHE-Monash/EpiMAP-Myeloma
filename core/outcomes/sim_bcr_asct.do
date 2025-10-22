**********
* SIM BCR SCT Vector

* Purpose: Determine BCR to SCT
* Outcome: Categorical 1 to 5 (1 = CR, 5 = SD), no 6 
**********

mata {
	// Initialise outcome
	oBCR_SCT = J(st_nobs(), 1, 0)
	
	// Filter for valid SCT patients
	idx = selectindex((mMOR[., OMC-1] :== 0) :& 
	                  (mState[., 1] :<= OMC + 1) :&
	                  (vSCT_L1 :== 1))
	
	// Calculate for SCT patients only
	if (rows(idx) > 0) {
		
		// Extract L1 BCR
		vBCR_L1 = mBCR[., 1]
		
		// Create BCR dummy variables
		pBCR_2 = (vBCR_L1 :== 2)
		pBCR_3 = (vBCR_L1 :== 3)
		pBCR_4 = (vBCR_L1 :== 4)
		pBCR_5 = (vBCR_L1 :== 5)
		
		// Assemble patient matrix
		pBCR_SCT = (vAge[idx], vAge2[idx], vMale[idx], 
		            vECOG1[idx], vECOG2[idx], 
		            vRISS2[idx], vRISS3[idx],
		            pBCR_2[idx], pBCR_3[idx], pBCR_4[idx], pBCR_5[idx])
		
		// Extract coefficients
		nPredictors = cols(pBCR_SCT)
		coefBCR_SCT = bSCT_BCR[1, 1..nPredictors]'
		
		// Extract cut points
		nCutPoints = 4
		cutPointIndices = (cols(bSCT_BCR) - nCutPoints + 1)..cols(bSCT_BCR)
		cutPoints = bSCT_BCR[1, cutPointIndices]
		
		// Calculate XB
		XB = pBCR_SCT * coefBCR_SCT
		
		// Calculate probabilities
		cumProbs = calcOrdLogitProbs(XB, cutPoints)
		
		// Generate random numbers
		RN = runiform(rows(idx), 1)
		
		// Assign outcomes
		categoryValues = (1, 2, 3, 4, 5)
		vOutcome = assignOrdOutcome(RN, cumProbs, categoryValues)
		
		// Update outcome vector
		oBCR_SCT[idx] = vOutcome
	}
	
	// Handle prevalent patients
	idxPrevalent = selectindex(mState[., 1] :> OMC + 1)
	if (rows(idxPrevalent) > 0) {
		oBCR_SCT[idxPrevalent] = mBCR[idxPrevalent, 10]
	}
	
	// Update matrices
	mBCR[., 10] = oBCR_SCT
}
