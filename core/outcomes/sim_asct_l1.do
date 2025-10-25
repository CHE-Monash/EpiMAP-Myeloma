**********
* SIM SCT L1
* 
* Purpose: Determine ASCT eligibility at Line 1 End
* Method: Logistic regression with outcome filters: TXR != 7 AND BCR != 6
* Outcome: Binary (0 = Not eligible, 1 = Eligible)
**********
	
mata {
	// Initialise outcome
	oSCT_L1 = J(st_nobs(), 1, .)
	
	// Filter for valid patients: Alive & State filters
	idx = selectindex((mMOR[., OMC-1] :== 0) :& (mState[., 1] :<= OMC + 1))
		
	// Calculate for valid patients
	if (rows(idx) > 0) {
		
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
		vAge70_e = vAge70[idx]
		vAge75_e = vAge75[idx]
		vCMc0_e = vCMc0[idx]
		vCMc1_e = vCMc1[idx]
		vCMc2_e = vCMc2[idx]
		vCMc3_e = vCMc3[idx]
		vCons_e = vCons[idx]
		
		// Assemble patient matrix
		mPat = (vAge_e, vAge2_e, vMale_e, 
		        vECOG0_e, vECOG1_e, vECOG2_e, 
				vRISS1_e, vRISS2_e, vRISS3_e, 
				vBCR1_e, vBCR2_e, vBCR3_e, vBCR4_e, vBCR5_e,
				vAge70_e, vAge75_e, 
				vCMc0_e, vCMc1_e, vCMc2_e, vCMc3_e, 
				vCons_e)

		// Extract coefficients
		nPredictors = cols(mPat)
		vCoef = bL1_SCT[1, 1..nPredictors]'

		// Calculate XB
		vXB = mPat * vCoef
		
		// Calculate probabilities
		vPR = 1 :/ (1 :+ exp(-vXB))
		
		// Generate random numbers
		vRN = runiform(rows(idx), 1)
		
		// Extract TXR from matrix
		vTXR_L1_e = mTXR[idx, 1]
			
		// Determine outcome: PR > RN AND TXR != 7 AND BCR != 6
		vOC = (vPR :> vRN) :& (vTXR_L1_e :!= 7) :& (vBCR1_e :!= 6)		
		
		// Convert logical to numeric (0/1)
		vOC = vOC :* 1
		
		// Update ONLY the filtered patients (0 or 1)
		oSCT_L1[idx] = vOC
	}
		
	// Update matrices 
	vSCT_L1 = oSCT_L1
	mSCT[., 2] = oSCT_L1
}
