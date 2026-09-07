**********
* Monash Myeloma Model - Sim TXD
*
* Purpose: Draw treatment duration for L2 onwards via parametric survival. Continuous time in
*          months.
**********

mata {
	// Initialise outcome
	vOC = J(Obs, 1, .)
	
	// Filter for alive and eligible
	idx = selectindex((mMOR[.,OMC-1] :== 0) :& (mState[.,1] :<= OMC))
	if (rows(idx) > 0) {
		
		// Determine model stage and BCR structure
		if (Line == 2) {
			vCoef = bL2_TXD
			dist = fbL2_TXD
			vTXR = oL2_TXR
			maxTXD = maxL2_TXD
		}
		else if (Line == 3) {
			vCoef = bL3_TXD
			dist = fbL3_TXD
			vTXR = oL3_TXR
			maxTXD = maxL3_TXD
		}
		else if (Line == 4) {
			vCoef = bL4_TXD
			dist = fbL4_TXD
			vTXR = oL4_TXR
			maxTXD = maxL4_TXD
		}
		else if (Line == 5) {
			vCoef = bL5_TXD
			dist = fbL5_TXD
			vTXR = J(1, 0, .)
			maxTXD = maxL5_TXD
		}
		else if (Line >= 6) {
			vCoef = bLX_TXD
			dist = fbLX_TXD
			vTXR = J(1, 0, .)
			maxTXD = maxLX_TXD
		}
		
		// Assemble patient matrix
		mPat = (vAge[idx], vAge2[idx], vMale[idx], 
				vECOG0[idx], vECOG1[idx], vECOG2[idx], 
				vRISS1[idx], vRISS2[idx], vRISS3[idx])
		
		// Add BCR
		vBCR = mBCR[idx, Line]
		vBCR_1 = (vBCR :== 1)
		vBCR_2 = (vBCR :== 2)
		vBCR_3 = (vBCR :== 3)
		vBCR_4 = (vBCR :== 4)
		vBCR_5 = (vBCR :== 5)
		vBCR_6 = (vBCR :== 6)
		mPat = (mPat, vBCR_1, vBCR_2, vBCR_3, vBCR_4, vBCR_5, vBCR_6)
		
		// Add TXR
		currentTX = mTXR[idx, Line]			
		if (cols(vTXR) >= 1) mPat = mPat, (currentTX :== vTXR[1, 1])
		if (cols(vTXR) >= 2) mPat = mPat, (currentTX :== vTXR[1, 2])
		if (cols(vTXR) >= 3) mPat = mPat, (currentTX :== vTXR[1, 3])
		if (cols(vTXR) >= 4) mPat = mPat, (currentTX :== vTXR[1, 4])
		
		// Add constant
		mPat = (mPat, vCons[idx])
		
		// Extract coefficients
		aux = vCoef[1, cols(vCoef)]
		nPredictors = cols(mPat)
		// Guard: design columns must equal coefficients minus the ancillary (1). Catches a
		// silent off-by-one if a factor's level count drifts (e.g. an MI-introduced category).
		if (nPredictors != cols(vCoef) - 1) {
			errprintf("sim_txd: design/coefficient mismatch at Line %g - mPat has %g columns but coefficients imply %g predictors\n", Line, nPredictors, cols(vCoef) - 1)
			exit(459)
		}
		vCoef = vCoef[1, 1..nPredictors]'
		
		// Calculate XB and OC
		vXB = mPat * vCoef
		vRN = rnDraw(idx, rn_txd(Line))
		vOC = calcSurvTime(vXB, vRN, dist, aux)
		
		// Curtail if beyond maximum observed
		vOC = rowmin((vOC, J(rows(vOC), 1, maxTXD)))

		// Update matrices
		mTXD[idx, Line] = round(vOC, 0.01)
		mTNE[idx, OMC] = round(vOC, 0.01)
		mTSD[idx, OMC+1] = mTSD[idx, OMC] :+ mTNE[idx, OMC]
	}
}
