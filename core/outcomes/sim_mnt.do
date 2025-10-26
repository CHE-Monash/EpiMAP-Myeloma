**********
* SIM MNT Vector

* Purpose: Determine BCR to SCT
* Outcome: Categorical 1 to 5 (1 = CR, 5 = SD), no 6
**********

mata {
	// Initialise outcome
	oMNT = J(st_nobs(), 1, 0)
	
	// Filter for valid patients
	idx = selectindex((mMOR[., OMC-1] :== 0) :& (mState[., 1] :<= OMC + 1))
	
	// Calculate for valid patients
	if (rows(idx) > 0) {
		
		// Extract current treatment regimen
		currentTX = mTXR[., Line]
		
		// Determine number of treatment regimen dummies
		if (Line == 1) oVector = oL1_TXR
		else if (Line == 2) oVector = oL2_TXR
		else if (Line == 3) oVector = oL3_TXR
		else if (Line == 4) oVector = oL4_TXR
		
		nRegimens = cols(oVector)
		
		// Create treatment regimen dummies
		if (nRegimens >= 2) TXR_is_R2 = (currentTX :== oVector[1, 2])
		if (nRegimens >= 3) TXR_is_R3 = (currentTX :== oVector[1, 3])
		if (nRegimens >= 4) TXR_is_R4 = (currentTX :== oVector[1, 4])
		
		// Extract L1 BCR
		vBCR_L1 = mBCR[., 1]
		
		// Create BCR dummy variables (reference = 1)
		pBCR_2 = (vBCR_L1 :== 2)
		pBCR_3 = (vBCR_L1 :== 3)
		pBCR_4 = (vBCR_L1 :== 4)
		pBCR_5 = (vBCR_L1 :== 5)
		pBCR_6 = (vBCR_L1 :== 6)
		
		// Assemble patient matrix
		pMNT = (vAge[idx], vAge2[idx], vMale[idx], 
		        vECOG0[idx], vECOG1[idx], vECOG2[idx], 
		        vRISS1[idx], vRISS2[idx], vRISS3[idx])
		
		// Add treatment regimen dummies
		if (nRegimens >= 2) pMNT = pMNT, TXR_is_R2[idx]
		if (nRegimens >= 3) pMNT = pMNT, TXR_is_R3[idx]
		if (nRegimens >= 4) pMNT = pMNT, TXR_is_R4[idx]
		
		// Add SCT
		pMNT = pMNT, vSCT_L1[idx]
		
		// Add BCR dummies
		pMNT = pMNT, (pBCR_2[idx], pBCR_3[idx], pBCR_4[idx], 
		              pBCR_5[idx], pBCR_6[idx])
		
		// Add constant
		pMNT = pMNT, vCons[idx]
		
		// Extract coefficients
		nPredictors = cols(pMNT)
		coefMNT = bMNT[1, 1..nPredictors]'
		
		// Calculate XB
		XB = pMNT * coefMNT
		
		// Calculate probability
		prMNT = 1 :/ (1 :+ exp(-XB))
		
		// Generate random numbers
		RN = runiform(rows(idx), 1)
		
		// Determine outcome
		vOutcome = (prMNT :> RN) :* 1
		
		// Update outcome vector
		oMNT[idx] = vOutcome
	}
	
	// Update matrices
	vMNT = oMNT
}
