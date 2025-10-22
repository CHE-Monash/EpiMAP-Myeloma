**********
*SIM CR L1
*
* Purpose: Select chemotherapy regimen for Line 1 using multinomial logit
* Coefficients: bL1_CR[1, 12:22] for XB2, [1, 23:33] for XB3, [1, 34:44] for XB4
**********

mata {
	// Filters
	vEligible = (mMOR[.,OMC-1] :== 0) :& (mState[.,1] :<= OMC+1)
	idxEligible = selectindex(vEligible)
	
	if (rows(idxEligible) > 0) {
		// Get number of regimen options
		nL1 = cols(oL1_CR)
		
		// Generate random numbers
		vRN = runiform(rows(mMOR), 1)
		
		// XB1 is always 1 (reference category)
		vXB1 = J(rows(idxEligible), 1, 1)
		vXB2 = J(rows(idxEligible), 1, 0)
		vXB3 = J(rows(idxEligible), 1, 0)
		vXB4 = J(rows(idxEligible), 1, 0)
		
		// Extract patient characteristics for eligible patients
		vAge_e = mAge[idxEligible, OMC]
		vAge2_e = vAge_e :^ 2
		vMale_e = mCore[idxEligible, cMale]
		vECOG1_e = (mCore[idxEligible, cECOG] :== 1)
		vECOG2_e = (mCore[idxEligible, cECOG] :== 2)
		vRISS2_e = (mCore[idxEligible, cRISS] :== 2)
		vRISS3_e = (mCore[idxEligible, cRISS] :== 3)
		vSCT_e = mCore[idxEligible, cSCT]
		vCons_e = vCons[idxEligible]
		
		// Build patient matrix (all predictors)
		pMatrix = (vAge_e, vAge2_e, vMale_e, vECOG1_e, vECOG2_e, 
		           vRISS2_e, vRISS3_e, vSCT_e, vCons_e)
		
		if (nL1 >= 2) {
			// Extract coefficients for XB2: columns 12-22
			coef_XB2 = bL1_CR[1, (12, 13, 14, 16, 17, 19, 20, 21, 22)]'
			
			// Calculate XB2 and exponentiate
			vXB2 = exp(pMatrix * coef_XB2)
		}
		
		if (nL1 >= 3) {
			// Extract coefficients for XB3: columns 23-33
			coef_XB3 = bL1_CR[1, (23, 24, 25, 27, 28, 30, 31, 32, 33)]'
			
			// Calculate XB3 and exponentiate
			vXB3 = exp(pMatrix * coef_XB3)
		}
		
		if (nL1 >= 4) {
			// Extract coefficients for XB4: columns 34-44
			coef_XB4 = bL1_CR[1, (34, 35, 36, 38, 39, 41, 42, 43, 44)]'
			
			// Calculate XB4 and exponentiate
			vXB4 = exp(pMatrix * coef_XB4)
		}
		
		// Calculate denominator for softmax
		vDenom = vXB1 :+ vXB2 :+ vXB3 :+ vXB4
		
		// Calculate cumulative probabilities
		vPR1 = vXB1 :/ vDenom
		vPR2 = vPR1 :+ (vXB2 :/ vDenom)
		vPR3 = vPR2 :+ (vXB3 :/ vDenom)
		vPR4 = vPR3 :+ (vXB4 :/ vDenom)
		
		// Assign regimen based on random number and probabilities
		vRN_e = vRN[idxEligible]
		vOutcome = J(rows(idxEligible), 1, .)
		
		// Default to first regimen (reference)
		vOutcome = (vRN_e :< vPR1) :* oL1_CR[1,1]
		
		// Assign second regimen
		if (nL1 >= 2) {
			vOutcome = vOutcome :+ ((vRN_e :>= vPR1) :& (vRN_e :< vPR2)) :* oL1_CR[1,2]
		}
		
		// Assign third regimen
		if (nL1 >= 3) {
			vOutcome = vOutcome :+ ((vRN_e :>= vPR2) :& (vRN_e :< vPR3)) :* oL1_CR[1,3]
		}
		
		// Assign fourth regimen
		if (nL1 >= 4) {
			vOutcome = vOutcome :+ ((vRN_e :>= vPR3) :& (vRN_e :< vPR4)) :* oL1_CR[1,4]
		}
		
		// Update outcome matrices
		mCore[idxEligible, cCR] = vOutcome
		mTXR[idxEligible, OMC-1] = vOutcome  // L1 is column 1 in mCR (v2.0 has OMC-1)
	}
}
