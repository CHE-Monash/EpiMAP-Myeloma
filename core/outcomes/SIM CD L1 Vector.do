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
		vOutcome = J(rows(mMOR), 1, .)
		
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
			vECOG1_ASCT = (vECOG[idxASCT] :== 1)
			vECOG2_ASCT = (vECOG[idxASCT] :== 2)
			vRISS2_ASCT = (vRISS[idxASCT] :== 2)
			vRISS3_ASCT = (vRISS[idxASCT] :== 3)
			vCons_ASCT = vCons[idxASCT]
			
			// Patient matrix
			pMatrix = (vAge_ASCT, vAge2_ASCT, vMale_ASCT, vECOG1_ASCT, vECOG2_ASCT, vRISS2_ASCT, vRISS3_ASCT)
			nCR = cols(oL1_CR)
			if (nCR >= 1) {
				vCR1_ASCT = (mTXR[idxASCT, 1] :== oL1_CR[1,1])
				pMatrix = (pMatrix , vCR1_ASCT)
			}
			if (nCR >= 2) {
				vCR2_ASCT = (mTXR[idxASCT, 1] :== oL1_CR[1,2])
				pMatrix = (pMatrix , vCR2_ASCT)
			}
			if (nCR >= 3) {
				vCR3_ASCT = (mTXR[idxASCT, 1] :== oL1_CR[1,3])
				pMatrix = (pMatrix , vCR3_ASCT)	
			}
			pMatrix = (pMatrix, vCons_ASCT)

			// --- SPLINE 1 ---
				// Extract coefficients
				nPredictors = cols(pMatrix)
				coef_S1 = bL1F1_CD_S1[1, 1..nPredictors]'
				aux_S1 = bL1F1_CD_S1[1, cols(bL1F1_CD_S1)]
				
				// Calculate XB
				vXB_S1 = pMatrix * coef_S1
				
				// Calculate survival time
				vRN_ASCT = vRN[idxASCT]
				vDur_S1 = calcSurvTime(vXB_S1, vRN_ASCT, fbL1F1_CD_S1, aux_S1)
			
				// --- SPLINE 2 (for those beyond Cutoff 1) ---
					vBeyondC1 = (vDur_S1 :> L1F1_CD_C1)
					idxBeyondC1 = selectindex(vBeyondC1)

					if (rows(idxBeyondC1) > 0) {
						// Patient matrix
						pMatrix_C1 = pMatrix[idxBeyondC1, .]
						
						// Extract coefficients
						nPredictors = cols(pMatrix_C1)
						coef_S2 = bL1F1_CD_S2[1, 1..nPredictors]'
						aux_S2 = bL1F1_CD_S2[1, cols(bL1F1_CD_S2)]
						
						// Calculate XB
						vXB_S2 = pMatrix_C1 * coef_S2
						
						// Conditional RN: survive past C1
						vSurvC1 = exp(-exp(vXB_S2) :* (L1F1_CD_C1 :^ exp(aux_S2)))
						vRN_S2 = runiform(rows(idxBeyondC1), 1) :* vSurvC1
						
						// Recalculate duration with Spline 2
						vDur_S2 = calcSurvTime(vXB_S2, vRN_S2, fbL1F1_CD_S2, aux_S2)
						
						// Update durations for those beyond C1
						vDur_S1[idxBeyondC1] = vDur_S2
				
					// --- SPLINE 3 (for those beyond Cutoff 2) ---
						vBeyondC2 = (vDur_S2 :> L1F1_CD_C2)
						idxBeyondC2_local = selectindex(vBeyondC2)
						
						if (rows(idxBeyondC2_local) > 0) {
							// Patient matrix
							idxBeyondC2 = idxBeyondC1[idxBeyondC2_local]
							pMatrix_C2 = pMatrix[idxBeyondC2, .]
							
							// Extract coefficients
							nPredictors = cols(pMatrix_C2)
							coef_S3 = bL1F1_CD_S3[1, 1..nPredictors]'
							aux_S3 = bL1F1_CD_S3[1, cols(bL1F1_CD_S3)]
							
							// Calculate XB
							vXB_S3 = pMatrix_C2 * coef_S3
							
							// Conditional RN: survive past C2
							vSurvC2 = exp(-exp(vXB_S3) :* (L1F1_CD_C2 :^ exp(aux_S3)))
							vRN_S3 = runiform(rows(idxBeyondC2), 1) :* vSurvC2
							
							// Recalculate duration with Spline 3
							vDur_S3 = calcSurvTime(vXB_S3, vRN_S3, fbL1F1_CD_S3, aux_S3)
							
							// Update durations for those beyond C2
							vDur_S1[idxBeyondC2] = vDur_S3
						}
					}
			
			// Store final durations for ASCT patients
			vOutcome[idxASCT] = vDur_S1
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
			vECOG1_NoASCT = (vECOG[idxNoASCT] :== 1)
			vECOG2_NoASCT = (vECOG[idxNoASCT] :== 2)
			vRISS2_NoASCT = (vRISS[idxNoASCT] :== 2)
			vRISS3_NoASCT = (vRISS[idxNoASCT] :== 3)
			vCons_NoASCT = vCons[idxNoASCT]

			// Patient matrix
			pMatrix = (vAge_NoASCT, vAge2_NoASCT, vMale_NoASCT, vECOG1_NoASCT, vECOG2_NoASCT, vRISS2_NoASCT, vRISS3_NoASCT)
			if (nCR >= 1) {
				vCR1_NoASCT = (mTXR[idxNoASCT, 1] :== oL1_CR[1,1])
				pMatrix = (pMatrix , vCR1_NoASCT)
			}
			if (nCR >= 2) {
				vCR2_NoASCT = (mTXR[idxNoASCT, 1] :== oL1_CR[1,2])
				pMatrix = (pMatrix , vCR2_NoASCT)
			}
			if (nCR >= 3) {
				vCR3_NoASCT = (mTXR[idxNoASCT, 1] :== oL1_CR[1,3])
				pMatrix = (pMatrix , vCR3_NoASCT)	
			}	
			pMatrix = (pMatrix , vCons_NoASCT)
						
			// Extract coefficients 
			nPredictors = cols(pMatrix)
			coef_F0 = bL1F0_CD[1, 1..nPredictors]'
			aux_F0 = bL1F0_CD[1, cols(bL1F0_CD)]

			// Calculate XB
			vXB_F0 = pMatrix * coef_F0
			
			// Calculate survival time
			vRN_NoASCT = vRN[idxNoASCT]
			vDur_F0 = calcSurvTime(vXB_F0, vRN_NoASCT, fbL1F0_CD, aux_F0)
			
			vOutcome[idxNoASCT] = vDur_F0
		}
				
		// ========================================
		// Update outcome matrices
		// ========================================
		mTXD[idxEligible, 1] = vOutcome[idxEligible]  // L1 duration in column 1
		
		// Update mTNE and mTSD
		mTNE[idxEligible, OMC] = vOutcome[idxEligible] :/ 365.25  // Convert days to years
		mTSD[idxEligible, OMC+1] = mTSD[idxEligible, OMC] :+ mTNE[idxEligible, OMC]
		
		// Update mCore for backwards compatibility (will be removed after full vectorisation)
		mCore[idxEligible, cCD] = vOutcome[idxEligible]
			
	}
}	

// ========================================
// GROUP 3: Continuous Therapy (CR == 7)
// Separate mata block to avoid parsing errors
// ========================================
capture confirm matrix bL1C_CD
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
			vECOG1_Cont = (vECOG[idxCont] :== 1)
			vECOG2_Cont = (vECOG[idxCont] :== 2)
			vRISS2_Cont = (vRISS[idxCont] :== 2)
			vRISS3_Cont = (vRISS[idxCont] :== 3)
			vCons_Cont = vCons[idxCont]
			
			pMatrix = (vAge_Cont, vAge2_Cont, vMale_Cont, vECOG1_Cont, vECOG2_Cont, vRISS2_Cont, vRISS3_Cont, vCons_Cont)
			
			// Calculate XB with bL1C_CD coefficients
			coef_C = bL1C_CD[1, (1, 2, 3, 5, 6, 8, 9, cols(bL1C_CD)-1)]'
			vXB_C = pMatrix * coef_C
			gamma_C = bL1C_CD[1, cols(bL1C_CD)]
			vRN_Cont = runiform(rows(mMOR), 1)
			vRN_Cont = vRN_Cont[idxCont]
			
			// Calculate survival time using mata function
			vDur_C = calcSurvTime(vXB_C, vRN_Cont, fbL1C_CD, gamma_C)
			
			// Curtail if beyond observed maximum
			vDur_C = (vDur_C :> maxL1C_CD) :* maxL1C_CD :+ (vDur_C :<= maxL1C_CD) :* vDur_C
			
			// Update outcome matrices
			mTXD[idxCont, 1] = vDur_C
			mTNE[idxCont, OMC] = vDur_C :/ 365.25
			mTSD[idxCont, OMC+1] = mTSD[idxCont, OMC] :+ mTNE[idxCont, OMC]
			mCore[idxCont, cCD] = vDur_C
	}
}
		
