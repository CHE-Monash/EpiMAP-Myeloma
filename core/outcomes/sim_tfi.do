**********
* Monash Myeloma Model - Sim TFI (L2+)
*
* Purpose: Draw Treatment-free Interval (time from LXE to L(X+1)S) via parametric survival.
*          Continuous time in months.
**********

mata {
    // Initialise outcome
	vOC = J(Obs, 1, .)
	
	// Filter for alive and eligible
    idx = selectindex((mMOR[.,OMC-1] :== 0) :& (mState[.,1] :<= OMC))
	if (rows(idx) > 0) {

		// Select coefficients
		if (Line == 2) {
			vCoef = bL2_TFI
			dist = fbL2_TFI
			maxTFI = maxL2_TFI
		}
		else if (Line == 3) {
			vCoef = bL3_TFI
			dist = fbL3_TFI
			maxTFI = maxL3_TFI
		}
		else if (Line == 4) {
			vCoef = bL4_TFI
			dist = fbL4_TFI
			maxTFI = maxL4_TFI
		}
		else if (Line == 5) {
			vCoef = bL5_TFI
			dist = fbL5_TFI
			maxTFI = maxL5_TFI
		}
		else if (Line >= 6) {
			vCoef = bLX_TFI
			dist = fbLX_TFI
			maxTFI = maxLX_TFI
		}
		
		// Assemble patient matrix
		mPat = (vAge[idx], vAge2[idx], vMale[idx], 
				vECOG0[idx], vECOG1[idx], vECOG2[idx],
		        vRISS1[idx], vRISS2[idx], vRISS3[idx])
		
		// Add BCR. The block WIDTH comes from the coefficient vector (9 built above + block + _cons
		// + 1 ancillary), so a line whose fit collapses the response works with no special case, and
		// the VALUE domain is clamped to match - without the clamp a simulated SD/PD patient scores 0
		// on every dummy of a narrowed block and silently falls through to the CR base level.
		// L5 is fitted on min(BCR_L5, 2): see prep/risk_equations.do and docs/refractory.md.
		nBCR = cols(vCoef) - 11
		vBCR = mBCR[idx, Line]
		if (nBCR < 6) vBCR = rowmin((vBCR, J(rows(idx), 1, nBCR)))
		for (k = 1; k <= nBCR; k++) mPat = mPat, (vBCR :== k)
	
		// Add constant
		mPat = mPat, vCons[idx]
		
		// Extract coefficients
		aux = vCoef[1, cols(vCoef)]
		nPredictors = cols(mPat)
		// Guard: design columns must equal coefficients minus the ancillary (1). Catches a
		// silent off-by-one if a factor's level count drifts (e.g. an MI-introduced category).
		if (nPredictors != cols(vCoef) - 1) {
			errprintf("sim_tfi: design/coefficient mismatch at Line %g - mPat has %g columns but coefficients imply %g predictors\n", Line, nPredictors, cols(vCoef) - 1)
			exit(459)
		}
		vCoef = vCoef[1, 1..nPredictors]'
		
		// Calculate XB and OC
		vXB = mPat * vCoef
		vRN = rnDraw(idx, rn_tfi(Line))
		vOC = calcSurvTime(vXB, vRN, dist, aux)
		
		// Curtail if beyond maximum observed
		vOC = rowmin((vOC, J(rows(vOC), 1, maxTFI)))
		
		// Update matrices
		mTFI[idx, Line+1] = round(vOC, 0.01) // Col 1 of mTFI is DN, so TFI_L3 goes in column 4
		mTNE[idx, OMC] = round(vOC, 0.01)
		mTSD[idx, OMC+1] = mTSD[idx, OMC] + mTNE[idx, OMC]
	}
}
