**********
*SIM TXR

*Purpose: Select treatment regimen
*Method: Multinomial logit
*Outcome: Regimen number
**********

mata {
	vEligible = (mMOR[.,OMC-1] :== 0) :& (mState[.,1] :<= OMC+1)
	idxEligible = selectindex(vEligible)
	
	if (rows(idxEligible) > 0) {
		
		// Determine model stage
		if (Line == 0) {
			vCoef = bL1_TXR
			vTXR = oL1_TXR
		}
		else if (Line == 1) {
			vCoef = bL2_TXR
			vTXR = oL2_TXR
		}
		else if (Line == 2) {
			vCoef = bL3_TXR
			vTXR = oL3_TXR
		}
		else if (Line = 3) {
			vCoef = bL4_TXR
			vTXR = oL4_TXR
		}
	
		nRegimens = cols(vTXR)
		
		vRN = runiform(rows(mMOR), 1)
		
		vXB1 = J(rows(idxEligible), 1, 1)
		vXB2 = J(rows(idxEligible), 1, 0)
		vXB3 = J(rows(idxEligible), 1, 0)
		vXB4 = J(rows(idxEligible), 1, 0)
		
		// Patient characteristics
		vAge_e = mAge[idxEligible, OMC]
		vAge2_e = vAge_e :^ 2
		vMale_e = vMale[idxEligible]
		vECOG0_e = (vECOG[idxEligible] :== 0)
		vECOG1_e = (vECOG[idxEligible] :== 1)
		vECOG2_e = (vECOG[idxEligible] :== 2)
		vRISS1_e = (vRISS[idxEligible] :== 1)
		vRISS2_e = (vRISS[idxEligible] :== 2)
		vRISS3_e = (vRISS[idxEligible] :== 3)
		vCons_e = vCons[idxEligible]
		
		// Build patient matrix (L1)
		if (Line == 0) {
			vSCT_e = vSCT_DN[idxEligible]
			mPat = (vAge_e, vAge2_e, vMale_e, 
					   vECOG0_e, vECOG1_e, vECOG2_e, 
			           vRISS1_e, vRISS2_e, vRISS3_e, 
					   vSCT_e, 
					   vCons_e)
			
			nPredictors = cols(mPat)

			if (nRegimens >= 2) {
				coef_XB2 = vCoef[1, (nPredictors+1)..(2*nPredictors)]'
				vXB2 = exp(mPat * coef_XB2)
			}

			if (nRegimens >= 3) {
				coef_XB3 = vCoef[1, (2*nPredictors+1)..(3*nPredictors)]'
				vXB3 = exp(mPat * coef_XB3)
			}

			if (nRegimens >= 4) {
				coef_XB4 = vCoef[1, (3*nPredictors+1)..(4*nPredictors)]'
				vXB4 = exp(mPat * coef_XB4)
			}
		}	
		// Build patient matrix (L2 or L3)
		else if (Line == 1 | Line == 2) {
			prevBCR = mBCR[idxEligible, Line]
			pBCR_1 = (prevBCR :== 1)
			pBCR_2 = (prevBCR :== 2)
			pBCR_3 = (prevBCR :== 3)
			pBCR_4 = (prevBCR :== 4)
			pBCR_5 = (prevBCR :== 5)
			pBCR_6 = (prevBCR :== 6)
			
			mPat = (vAge_e, vAge2_e, vMale_e, 
					vECOG0_e, vECOG1_e, vECOG2_e, 
			        vRISS1_e, vRISS2_e, vRISS3_e, 
					pBCR_1, pBCR_2, pBCR_3, pBCR_4, pBCR_5, pBCR_6, 
					vCons_e)
			
			nPredictors = cols(mPat)

			if (nRegimens >= 2) {
				coef_XB2 = vCoef[1, (nPredictors+1)..(2*nPredictors)]'
				vXB2 = exp(mPat * coef_XB2)
			}

			if (nRegimens >= 3) {
				coef_XB3 = vCoef[1, (2*nPredictors+1)..(3*nPredictors)]'
				vXB3 = exp(mPat * coef_XB3)
			}

			if (nRegimens >= 4) {
				coef_XB4 = vCoef[1, (3*nPredictors+1)..(4*nPredictors)]'
				vXB4 = exp(mPat * coef_XB4)
			}
		}
		// Build patient matrix (L4+)
		else {
			prevBCR = mBCR[idxEligible, Line]
			pBCR_1 = (prevBCR :== 1)
			pBCR_3 = (prevBCR :== 3)
			pBCR_5 = (prevBCR :== 5)
			
			mPat = (vAge_e, vAge2_e, vMale_e, 
					vECOG0_e, vECOG1_e, vECOG2_e, 
			        vRISS1_e, vRISS2_e, vRISS3_e, 
					pBCR_1, pBCR_3, pBCR_5, 
					vCons_e)
			
			if (nRegimens >= 2) {
				coef_XB2 = vCoef[1, (nPredictors+1)..(2*nPredictors)]'
				vXB2 = exp(mPat * coef_XB2)
			}
			
			if (nRegimens >= 3) {
				coef_XB3 = vCoef[1, (2*nPredictors+1)..(3*nPredictors)]'
				vXB3 = exp(mPat * coef_XB3)
			}
			
			if (nRegimens >= 4) {
				coef_XB4 = vCoef[1, (3*nPredictors+1)..(4*nPredictors)]'
				vXB4 = exp(mPat * coef_XB4)
			}
		}
		
		// Calculate probabilities
		vPR1 = vXB1 :/ (vXB1 :+ vXB2 :+ vXB3 :+ vXB4)
		vPR2 = vPR1 :+ (vXB2 :/ (vXB1 :+ vXB2 :+ vXB3 :+ vXB4))
		vPR3 = vPR2 :+ (vXB3 :/ (vXB1 :+ vXB2 :+ vXB3 :+ vXB4))
		vPR4 = vPR3 :+ (vXB4 :/ (vXB1 :+ vXB2 :+ vXB3 :+ vXB4))
		
		vRN_e = vRN[idxEligible]
		vOut = J(rows(idxEligible), 1, .)
		
		// Determine Outcome
		vOut = (vRN_e :< vPR1) :* vTXR[1,1]

		if (nRegimens >= 2) {
			vOut = vOut :+ ((vRN_e :>= vPR1) :& (vRN_e :< vPR2)) :* vTXR[1,2]
		}
		
		if (nRegimens >= 3) {
			vOut = vOut :+ ((vRN_e :>= vPR2) :& (vRN_e :< vPR3)) :* vTXR[1,3]
		}
		
		if (nRegimens >= 4) {
			vOut = vOut :+ ((vRN_e :>= vPR3) :& (vRN_e :< vPR4)) :* vTXR[1,4]
		}
		
		// Update matrix
		mTXR[idxEligible, Line+1] = vOut
	}
}
