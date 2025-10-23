**********
*SIM TFI L1
*
* Purpose: Treatment-free Interval at Line 1 End (time from L1E to L2S)
* Method: Parametric survival analysis, split by ASCT status
* Outcome: Continuous time (Months)
**********

mata {
	// Filters
	vEligible = (mMOR[.,OMC-1] :== 0) :& (mState[.,1] :<= OMC+1)
	idxEligible = selectindex(vEligible)
	
	if (rows(idxEligible) > 0) {
		// Generate random numbers for all patients
		vRN = runiform(rows(mMOR), 1)
		
		// Initialize outcome vector
		vOutcome = J(rows(mMOR), 1, .)
		
		// GROUP 1: ASCT Patients
		vIsASCT = (vSCT_L1[.,1] :== 1) :& vEligible
		idxASCT = selectindex(vIsASCT)
		
		if (rows(idxASCT) > 0) {
			// Patient vectors
			vAge_ASCT = mAge[idxASCT, OMC]
			vAge2_ASCT = vAge_ASCT :^ 2
			vMale_ASCT = vMale[idxASCT]
			vECOG0_ASCT = (vECOG[idxASCT] :== 0)
			vECOG1_ASCT = (vECOG[idxASCT] :== 1)
			vECOG2_ASCT = (vECOG[idxASCT] :== 2)
			vRISS1_ASCT = (vRISS[idxASCT] :== 1)
			vRISS2_ASCT = (vRISS[idxASCT] :== 2)
			vRISS3_ASCT = (vRISS[idxASCT] :== 3)
			vMNT_ASCT = vMNT[idxASCT]
			vBCR1_ASCT = (mBCR[idxASCT, Line] :== 1)
			vBCR2_ASCT = (mBCR[idxASCT, Line] :== 2)
			vBCR3_ASCT = (mBCR[idxASCT, Line] :== 3)
			vBCR4_ASCT = (mBCR[idxASCT, Line] :== 4)
			vCons_ASCT = vCons[idxASCT]
			
			// Patient matrix (ASCT patients cannot have BCR = 5 or 6)
			mPat_ASCT = (vAge_ASCT, vAge2_ASCT, vMale_ASCT, 
						 vECOG0_ASCT, vECOG1_ASCT, vECOG2_ASCT, 
			             vRISS1_ASCT, vRISS2_ASCT, vRISS3_ASCT, 
						 vMNT_ASCT, 
			             vBCR1_ASCT, vBCR2_ASCT, vBCR3_ASCT, vBCR4_ASCT, 
						 vCons_ASCT)
			
			// Extract coefficients for ASCT
			nPredictors = cols(mPat_ASCT)
			coef_ASCT = bL1_TFI_ASCT[1, 1..nPredictors]' 
			aux_ASCT = bL1_TFI_ASCT[1, cols(bL1_TFI_ASCT)]
			
			// Calculate XB
			vXB_ASCT = mPat_ASCT * coef_ASCT
			
			// Calculate outcome (survival time)
			vRN_ASCT = vRN[idxASCT]
			vOut[idxASCT] = calcSurvTime(vXB_ASCT, vRN_ASCT, fbL1_TFI_ASCT, aux_ASCT)
			
			// Curtail if beyond maximum observed
			vOut[idxASCT] = rowmin((vOut[idxASCT], J(rows(idxASCT), 1, maxL1_TFI_ASCT)))
		}
		
		// GROUP 2: NoASCT Patients
		vIsNoASCT = (vSCT_L1[.,1] :== 0) :& vEligible
		idxNoASCT = selectindex(vIsNoASCT)
		
		if (rows(idxNoASCT) > 0) {
			// Patient vectors
			vAge_NoASCT = mAge[idxNoASCT, OMC]
			vAge2_NoASCT = vAge_NoASCT :^ 2
			vMale_NoASCT = vMale[idxNoASCT]
			vECOG0_NoASCT = (vECOG[idxNoASCT] :== 0)
			vECOG1_NoASCT = (vECOG[idxNoASCT] :== 1)
			vECOG2_NoASCT = (vECOG[idxNoASCT] :== 2)
			vRISS1_NoASCT = (vRISS[idxNoASCT] :== 1)
			vRISS2_NoASCT = (vRISS[idxNoASCT] :== 2)
			vRISS3_NoASCT = (vRISS[idxNoASCT] :== 3)
			vMNT_NoASCT = vMNT[idxNoASCT]
			vBCR1_NoASCT = (mBCR[idxNoASCT, Line] :== 1)
			vBCR2_NoASCT = (mBCR[idxNoASCT, Line] :== 2)
			vBCR3_NoASCT = (mBCR[idxNoASCT, Line] :== 3)
			vBCR4_NoASCT = (mBCR[idxNoASCT, Line] :== 4)
			vBCR5_NoASCT = (mBCR[idxNoASCT, Line] :== 5)
			vBCR6_NoASCT = (mBCR[idxNoASCT, Line] :== 6)
			vCons_NoASCT = vCons[idxNoASCT]
			
			// Patient matrix (Non-SCT includes BCR 5 and 6)
			mPat_NoASCT = (vAge_NoASCT, vAge2_NoASCT, vMale_NoASCT, 
						   vECOG0_NoASCT, vECOG1_NoASCT, vECOG2_NoASCT,
			               vRISS1_NoASCT, vRISS2_NoASCT, vRISS3_NoASCT, 
						   vMNT_NoASCT,
			               vBCR1_NoASCT, vBCR2_NoASCT, vBCR3_NoASCT, vBCR4_NoASCT, vBCR5_NoASCT, vBCR6_NoASCT, 
						   vCons_NoASCT)
			
			// Extract coefficients for NoASCT
			nPredictors = cols(mPat_NoASCT)
			vCoef_NoASCT = bL1_TFI_NoASCT[1, 1..nPredictors]'
			aux_NoASCT = bL1_TFI_NoASCT[1, cols(bL1_TFI_NoASCT)]
			
			// Calculate XB
			vXB_NoASCT = mPat_NoASCT * vCoef_NoASCT
			
			// Calculate outcome (survival time)
			vRN_NoASCT = vRN[idxNoASCT]
			vOut[idxNoASCT] = calcSurvTime(vXB_NoASCT, vRN_NoASCT, fbL1_TFI_NoASCT, aux_NoASCT)
			
			// Curtail if beyond maximum observed
			vOut[idxNoASCT] = rowmin((vOut[idxNoASCT], J(rows(idxNoASCT), 1, maxL1_TFI_NoASCT)))

		}
		
		// Handle prevalent patients
		prevalent = selectindex(mState[., 1] :> OMC + 1)
		if (rows(prevalent) > 0) {
			vOutcome[prevalent] = mTNE[prevalent, OMC]
		}
		
		// Update matrices
		mTFI[., LX+1] = round(vOut, 0.1)
		mTNE[., OMC] = round(vOut, 0.1)
		mTSD[., OMC+1] = mTSD[., OMC] + mTNE[., OMC]
	}
}
