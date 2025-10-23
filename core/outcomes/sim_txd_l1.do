**********
*SIM CD L1
*
* Purpose: Chemotherapy Duration for Line 1
* Note: No mCore references except final update for backwards compatibility
**********

mata {
	// Filters
	vEligible = (mMOR[.,OMC-1] :== 0) :& (mState[.,1] :<= OMC+1)
	idxEligible = selectindex(vEligible)
	
	if (rows(idxEligible) > 0) {
		// Generate random numbers for all eligible patients
		vRN = runiform(rows(mMOR), 1)
		
		// Initialize outcome vector
		vOut = J(rows(mMOR), 1, .)
		
		// ========================================
		// GROUP 1: Fixed + ASCT (with splines)
		// ========================================
		// CR != 7 AND SCT == 1
		vIsASCT = (mTXR[.,1] :!= 7) :& (mSCT[.,1] :== 1) :& vEligible
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
			vCons_ASCT = vCons[idxASCT]
			
			// Patient matrix
			mPat_ASCT = (vAge_ASCT, vAge2_ASCT, vMale_ASCT, 
					   vECOG0_ASCT, vECOG1_ASCT, vECOG2_ASCT, 
					   vRISS1_ASCT, vRISS2_ASCT, vRISS3_ASCT)
			nCR = cols(oL1_TXR)
			if (nCR >= 1) {
				vCR1_ASCT = (mTXR[idxASCT, 1] :== oL1_TXR[1,1])
				mPat_ASCT = (mPat_ASCT , vCR1_ASCT)
			}
			if (nCR >= 2) {
				vCR2_ASCT = (mTXR[idxASCT, 1] :== oL1_TXR[1,2])
				mPat_ASCT = (mPat_ASCT , vCR2_ASCT)
			}
			if (nCR >= 3) {
				vCR3_ASCT = (mTXR[idxASCT, 1] :== oL1_TXR[1,3])
				mPat_ASCT = (mPat_ASCT , vCR3_ASCT)	
			}
			mPat_ASCT = (mPat_ASCT, vCons_ASCT)

			// --- SPLINE 1 ---
				// Extract coefficients
				nPredictors = cols(mPat_ASCT)
				coef_S1 = bL1_TXD_ASCT_S1[1, 1..nPredictors]'
				aux_S1 = bL1_TXD_ASCT_S1[1, cols(bL1_TXD_ASCT_S1)]
				
				// Calculate XB
				vXB_S1 = mPat_ASCT * coef_S1
				
				// Calculate survival time
				vRN_ASCT = vRN[idxASCT]
				vDur_S1 = calcSurvTime(vXB_S1, vRN_ASCT, fbL1_TXD_ASCT_S1, aux_S1)
			
				// --- SPLINE 2 (for those beyond Cutoff 1) ---
					vBeyondC1 = (vDur_S1 :> L1_TXD_ASCT_C1)
					idxBeyondC1 = selectindex(vBeyondC1)

					if (rows(idxBeyondC1) > 0) {
						// Patient matrix
						mPat_ASCT_C1 = mPat_ASCT[idxBeyondC1, .]
						
						// Extract coefficients
						nPredictors = cols(mPat_ASCT_C1)
						coef_S2 = bL1_TXD_ASCT_S2[1, 1..nPredictors]'
						aux_S2 = bL1_TXD_ASCT_S2[1, cols(bL1_TXD_ASCT_S2)]
						
						// Calculate XB
						vXB_S2 = mPat_ASCT_C1 * coef_S2
						
						// Conditional RN: survive past C1
						vSurvC1 = exp(-exp(vXB_S2) :* (L1_TXD_ASCT_C1 :^ exp(aux_S2)))
						vRN_S2 = runiform(rows(idxBeyondC1), 1) :* vSurvC1
						
						// Recalculate duration with Spline 2
						vDur_S2 = calcSurvTime(vXB_S2, vRN_S2, fbL1_TXD_ASCT_S2, aux_S2)
						
						// Update durations for those beyond C1
						vDur_S1[idxBeyondC1] = vDur_S2
				
					// --- SPLINE 3 (for those beyond Cutoff 2) ---
						vBeyondC2 = (vDur_S2 :> L1_TXD_ASCT_C2)
						idxBeyondC2_local = selectindex(vBeyondC2)
						
						if (rows(idxBeyondC2_local) > 0) {
							// Patient matrix
							idxBeyondC2 = idxBeyondC1[idxBeyondC2_local]
							mPat_ASCT_C2 = mPat_ASCT[idxBeyondC2, .]
							
							// Extract coefficients
							nPredictors = cols(mPat_ASCT_C2)
							coef_S3 = bL1_TXD_ASCT_S3[1, 1..nPredictors]'
							aux_S3 = bL1_TXD_ASCT_S3[1, cols(bL1_TXD_ASCT_S3)]
							
							// Calculate XB
							vXB_S3 = mPat_ASCT_C2 * coef_S3
							
							// Conditional RN: survive past C2
							vSurvC2 = exp(-exp(vXB_S3) :* (L1_TXD_ASCT_C2 :^ exp(aux_S3)))
							vRN_S3 = runiform(rows(idxBeyondC2), 1) :* vSurvC2
							
							// Recalculate duration with Spline 3
							vDur_S3 = calcSurvTime(vXB_S3, vRN_S3, fbL1_TXD_ASCT_S3, aux_S3)
							
							// Update durations for those beyond C2
							vDur_S1[idxBeyondC2] = vDur_S3
						}
					}
			
			// Store final durations for ASCT patients
			vOut[idxASCT] = vDur_S1
		}
	
		// ========================================
		// GROUP 2: Fixed + No ASCT
		// ========================================
		vIsNoASCT = (mTXR[.,1] :!= 7) :& (mSCT[.,1] :== 0) :& vEligible
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
			vCons_NoASCT = vCons[idxNoASCT]

			// Patient matrix
			mPat_NoASCT = (vAge_NoASCT, vAge2_NoASCT, vMale_NoASCT, 
					       vECOG0_NoASCT, vECOG1_NoASCT, vECOG2_NoASCT, 
					       vRISS1_NoASCT, vRISS2_NoASCT, vRISS3_NoASCT)
			if (nCR >= 1) {
				vCR1_NoASCT = (mTXR[idxNoASCT, 1] :== oL1_TXR[1,1])
				mPat_NoASCT = (mPat_NoASCT , vCR1_NoASCT)
			}
			if (nCR >= 2) {
				vCR2_NoASCT = (mTXR[idxNoASCT, 1] :== oL1_TXR[1,2])
				mPat_NoASCT = (mPat_NoASCT , vCR2_NoASCT)
			}
			if (nCR >= 3) {
				vCR3_NoASCT = (mTXR[idxNoASCT, 1] :== oL1_TXR[1,3])
				mPat_NoASCT = (mPat_NoASCT , vCR3_NoASCT)	
			}	
			mPat_NoASCT = (mPat_NoASCT , vCons_NoASCT)
						
			// Extract coefficients 
			nPredictors = cols(mPat_NoASCT)
			coef_NoASCT = bL1_TXD_NoASCT[1, 1..nPredictors]'
			aux_NoASCT = bL1_TXD_NoASCT[1, cols(bL1_TXD_NoASCT)]

			// Calculate XB
			vXB_NoASCT = mPat_NoASCT * coef_NoASCT
			
			// Calculate survival time
			vRN_NoASCT = vRN[idxNoASCT]
			vOut[idxNoASCT] = calcSurvTime(vXB_NoASCT, vRN_NoASCT, fbL1_TXD_NoASCT, aux_NoASCT)
		}
				
		// Update outcome matrices
		mTXD[idxEligible, 1] = round(vOut[idxEligible], 0.1)
		mTNE[idxEligible, OMC] = round(vOut[idxEligible], 0.1)
		mTSD[idxEligible, OMC+1] = mTSD[idxEligible, OMC] :+ mTNE[idxEligible, OMC]
	}
}	

// ========================================
// GROUP 3: Continuous Therapy (CR == 7)
// Separate mata block to avoid parsing errors
// ========================================
capture confirm matrix bL1_TXD_Cont
if _rc == 0 {
	mata {
		vEligible = (mMOR[.,OMC-1] :== 0) :& (mState[.,1] :<= OMC+1)
		vIsCont = (mTXR[.,1] :== 7) :& vEligible
		idxCont = selectindex(vIsCont)
		
		if (rows(idxCont) > 0) {
			// Build patient matrix from source vectors (no CR dummies for continuous)
			vAge_Cont = mAge[idxCont, OMC]
			vAge2_Cont = vAge_Cont :^ 2
			vMale_Cont = vMale[idxCont]
			vECOG0_Cont = (vECOG[idxCont] :== 0)
			vECOG1_Cont = (vECOG[idxCont] :== 1)
			vECOG2_Cont = (vECOG[idxCont] :== 2)
			vRISS1_Cont = (vRISS[idxCont] :== 1)
			vRISS2_Cont = (vRISS[idxCont] :== 2)
			vRISS3_Cont = (vRISS[idxCont] :== 3)
			vCons_Cont = vCons[idxCont]
			
			mPat_Cont = (vAge_Cont, vAge2_Cont, vMale_Cont,
					     vECOG0_Cont, vECOG1_Cont, vECOG2_Cont, 
					     vRISS1_Cont, vRISS2_Cont, vRISS3_Cont, 
					     vCons_Cont)
			
			// Extract coefficients 
			nPredictors = cols(mPat_Cont)
			coef_Cont = bL1_TXD_Cont[1, 1..nPredictors]'
			dist_Cont = fbL1_TXD_Cont
			aux_Cont = bL1_TXD_Cont[1, cols(bL1_TXD_Cont)]

			// Calculate XB
			vXB_Cont = mPat_Cont * coef_Cont
			
			// Calculate survival time using mata function
			vRN_Cont = vRN_Cont[idxCont]
			vOut_Cont = calcSurvTime(vXB_Cont, vRN_Cont, dist_Cont, aux_Cont)
			
			// Curtail if beyond observed maximum
			vOut_Cont = (vOut_Cont :> maxL1_TXD_Cont) :* maxL1_TXD_Cont :+ (vOut_Cont :<= maxL1_TXD_Cont) :* vOut_Cont
			
			// Update matrices
			mTXD[idxCont, 1] = round(vOut_Cont, 0.1)
			mTNE[idxCont, OMC] = round(vOut_Cont, 0.1)
			mTSD[idxCont, OMC+1] = mTSD[idxCont, OMC] :+ mTNE[idxCont, OMC]
		}
	}

