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
			mCoef = bL1_CR
			rOutcomes = oL1_CR
		}
		else if (Line == 1) {
			mCoef = bL2_CR
			rOutcomes = oL2_CR
		}
		else if (Line == 2) {
			mCoef = bL3_CR
			rOutcomes = oL3_CR
		}
		else if (Line == 3) {
			mCoef = bL4_CR
			rOutcomes = oL4_CR
		}
/*		else if (Line == 4) {
			mCoef = bL5_CR
			rOutcomes = oL5_CR
		}
		else if (Line == 5) {
			mCoef = bL6_CR
			rOutcomes = oL6_CR
		}
		else if (Line == 6) {
			mCoef = bL7_CR
			rOutcomes = oL7_CR
		}
		else if (Line == 7) {
			mCoef = bL8_CR
			rOutcomes = oL8_CR
		}
		else if (Line == 8) {
			mCoef = bL9_CR
			rOutcomes = oL9_CR
		}
*/		
		nRegimens = cols(rOutcomes)
		
		vRN = runiform(rows(mMOR), 1)
		
		vXB1 = J(rows(idxEligible), 1, 1)
		vXB2 = J(rows(idxEligible), 1, 0)
		vXB3 = J(rows(idxEligible), 1, 0)
		vXB4 = J(rows(idxEligible), 1, 0)
		
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

		if (Line == 0) {
			vSCT_e = vSCT_DN[idxEligible]
			pMatrix = (vAge_e, vAge2_e, vMale_e, 
					   vECOG0_e, vECOG1_e, vECOG2_e, 
			           vRISS1_e, vRISS2_e, vRISS3_e, 
					   vSCT_e, 
					   vCons_e)
			
			nPredictors = cols(pMatrix)

			if (nRegimens >= 2) {
				coef_XB2 = mCoef[1, (nPredictors+1)..(2*nPredictors)]'
				vXB2 = exp(pMatrix * coef_XB2)
			}

			if (nRegimens >= 3) {
				coef_XB3 = mCoef[1, (2*nPredictors+1)..(3*nPredictors)]'
				vXB3 = exp(pMatrix * coef_XB3)
			}

			if (nRegimens >= 4) {
				coef_XB4 = mCoef[1, (3*nPredictors+1)..(4*nPredictors)]'
				vXB4 = exp(pMatrix * coef_XB4)
			}
		}		

	else if (Line == 1) {
			prevBCR = mBCR[idxEligible, 1]
			pBCR_1 = (prevBCR :== 1)
			pBCR_2 = (prevBCR :== 2)
			pBCR_3 = (prevBCR :== 3)
			pBCR_4 = (prevBCR :== 4)
			pBCR_5 = (prevBCR :== 5)
			pBCR_6 = (prevBCR :== 6)
			
			pMatrix = (vAge_e, vAge2_e, vMale_e, 
					   vECOG0_e, vECOG1_e, vECOG2_e, 
			           vRISS1_e, vRISS2_e, vRISS3_e, 
					   pBCR_1, pBCR_2, pBCR_3, pBCR_4, pBCR_5, pBCR_6, 
					   vCons_e)
			
			nPredictors = cols(pMatrix)

			if (nRegimens >= 2) {
				coef_XB2 = mCoef[1, (nPredictors+1)..(2*nPredictors)]'
				vXB2 = exp(pMatrix * coef_XB2)
			}

			if (nRegimens >= 3) {
				coef_XB3 = mCoef[1, (2*nPredictors+1)..(3*nPredictors)]'
				vXB3 = exp(pMatrix * coef_XB3)
			}

			if (nRegimens >= 4) {
				coef_XB4 = mCoef[1, (3*nPredictors+1)..(4*nPredictors)]'
				vXB4 = exp(pMatrix * coef_XB4)
			}
		}
		else if (Line == 2 | Line == 3) {
			prevBCR = mBCR[idxEligible, Line-1]
			pBCR_1 = (prevBCR :== 1)
			pBCR_2 = (prevBCR :== 2)
			pBCR_3 = (prevBCR :== 3)
			pBCR_4 = (prevBCR :== 4)
			pBCR_5 = (prevBCR :== 5)
			pBCR_6 = (prevBCR :== 6)
			
			pMatrix = (vAge_e, vAge2_e, vMale_e, 
					   vECOG0_e, vECOG1_e, vECOG2_e, 
			           vRISS1_e, vRISS2_e, vRISS3_e, 
					   pBCR_1, pBCR_2, pBCR_3, pBCR_4, pBCR_5, pBCR_6, 
					   vCons_e)
			
			nPredictors = cols(pMatrix)

			if (nRegimens >= 2) {
				coef_XB2 = mCoef[1, (nPredictors+1)..(2*nPredictors)]'
				vXB2 = exp(pMatrix * coef_XB2)
			}

			if (nRegimens >= 3) {
				coef_XB3 = mCoef[1, (2*nPredictors+1)..(3*nPredictors)]'
				vXB3 = exp(pMatrix * coef_XB3)
			}

			if (nRegimens >= 4) {
				coef_XB4 = mCoef[1, (3*nPredictors+1)..(4*nPredictors)]'
				vXB4 = exp(pMatrix * coef_XB4)
			}
		}
		else {
			prevBCR = mBCR[idxEligible, Line-1]
			pBCR_1 = (prevBCR :== 1)
			pBCR_3 = (prevBCR :== 3)
			pBCR_5 = (prevBCR :== 5)
			
			pMatrix = (vAge_e, vAge2_e, vMale_e, 
					   vECOG0_e, vECOG1_e, vECOG2_e, 
			           vRISS1_e, vRISS2_e, vRISS3_e, 
					   pBCR_1, pBCR_3, pBCR_5, 
					   vCons_e)
			
			if (nRegimens >= 2) {
				coef_XB2 = mCoef[1, (14, 15, 16, 18, 19, 21, 22, 24, 25, 26)]'
				vXB2 = exp(pMatrix * coef_XB2)
			}
			
			if (nRegimens >= 3) {
				coef_XB3 = mCoef[1, (27, 28, 29, 31, 32, 34, 35, 37, 38, 39)]'
				vXB3 = exp(pMatrix * coef_XB3)
			}
			
			if (nRegimens >= 4) {
				coef_XB4 = mCoef[1, (40, 41, 42, 44, 45, 47, 48, 50, 51, 52)]'
				vXB4 = exp(pMatrix * coef_XB4)
			}
		}
		
		vPR1 = vXB1 :/ (vXB1 :+ vXB2 :+ vXB3 :+ vXB4)
		vPR2 = vPR1 :+ (vXB2 :/ (vXB1 :+ vXB2 :+ vXB3 :+ vXB4))
		vPR3 = vPR2 :+ (vXB3 :/ (vXB1 :+ vXB2 :+ vXB3 :+ vXB4))
		vPR4 = vPR3 :+ (vXB4 :/ (vXB1 :+ vXB2 :+ vXB3 :+ vXB4))
		
		vRN_e = vRN[idxEligible]
		vOutcome = J(rows(idxEligible), 1, .)
		
		vOutcome = (vRN_e :< vPR1) :* rOutcomes[1,1]

		if (nRegimens >= 2) {
			vOutcome = vOutcome :+ ((vRN_e :>= vPR1) :& (vRN_e :< vPR2)) :* rOutcomes[1,2]
		}
		
		if (nRegimens >= 3) {
			vOutcome = vOutcome :+ ((vRN_e :>= vPR2) :& (vRN_e :< vPR3)) :* rOutcomes[1,3]
		}
		
		if (nRegimens >= 4) {
			vOutcome = vOutcome :+ ((vRN_e :>= vPR3) :& (vRN_e :< vPR4)) :* rOutcomes[1,4]
		}
		
		mTXR[idxEligible, Line+1] = vOutcome
	}
}
