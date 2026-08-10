**********
* Monash Myeloma Model - Sim TFI L1
*
* Purpose: Draw Treatment-free Interval at Line 1 End (time from L1E to L2S) via parametric
*          survival, split by ASCT status. Continuous time in months.
*
* Notes:   TRUNCATED BELOW AT THE DRAWN MAINTENANCE DURATION for maintenance patients. The gap
*          cannot be shorter than the maintenance it contains, and enforcing it here is what lets
*          sim_mnd.do drop the ln(gap) covariate that was restricting its fit.
*
*          calcSurvTime maps a survivor probability U to a time and S is decreasing, so T >= L is
*          exactly U <= S(L). Drawing U' = U * S(L) satisfies it by construction - no rejection
*          sampling, no clipping, and the same random number is consumed, so CRN alignment holds.
*
*          The alternative was capping MND at the gap in process_data.do, which fired on a large
*          minority and could only ever shorten. Truncating the gap upward instead moves the
*          adjustment onto the quantity that can absorb it. That cap remains as a safety net.
*
* ORDER:   AFTER sim_mnd.do, which is a reversal of the previous order and the point of the design.
**********

mata {
	// Initialise outcome
	vOC = J(Obs, 1, .)
	
	// Filter for alive and eligible
	idx = selectindex((mMOR[.,OMC-1] :== 0) :& (mState[.,1] :<= OMC))	
	if (rows(idx) > 0) {

		// Group 1: ASCT
		idxASCT = idx[selectindex(vSCT_L1[idx] :== 1)]
		if (rows(idxASCT) > 0) {
			
			// Grab BCR to ASCT
			vBCR_1 = (mBCR[idxASCT, 10] :== 1)
			vBCR_2 = (mBCR[idxASCT, 10] :== 2)
			vBCR_3 = (mBCR[idxASCT, 10] :== 3)
			vBCR_4 = (mBCR[idxASCT, 10] :== 4)	
			
			// Assemble patient matrix (ASCT patients cannot have BCR = 5 or 6)
			mPat_ASCT = (vAge[idxASCT], vAge2[idxASCT], vMale[idxASCT], 
						 vECOG0[idxASCT], vECOG1[idxASCT], vECOG2[idxASCT], 
			             vRISS1[idxASCT], vRISS2[idxASCT], vRISS3[idxASCT], 
						 vMNT[idxASCT], 
			             vBCR_1, vBCR_2, vBCR_3, vBCR_4, 
						 vCons[idxASCT])
			
			// Extract coefficients for ASCT
			nPredictors = cols(mPat_ASCT)
			if (cols(bL1_TFI_ASCT) != nPredictors + 1) {
				errprintf("sim_tfi_l1 (ASCT): design/coefficient mismatch - mPat has %g columns so %g were expected (mean + ancillary), but bL1_TFI_ASCT has %g. A BCR level was likely empty in the fit, or the coefficient file was built by a different specification.\n",
					nPredictors, nPredictors + 1, cols(bL1_TFI_ASCT))
				exit(459)
			}
			coef_ASCT = bL1_TFI_ASCT[1, 1..nPredictors]' 
			aux_ASCT = bL1_TFI_ASCT[1, cols(bL1_TFI_ASCT)]
			
			// Calculate XB
			vXB_ASCT = mPat_ASCT * coef_ASCT
			
			// Calculate outcome (survival time)
			vRN_ASCT = rnDraw(idxASCT, rn_tfi_l1(1))

			// Truncate below at the maintenance already drawn (sim_mnd.do runs first). vMND is
			// missing for non-maintenance patients and whenever no MND model was fitted, so the
			// bound is only applied where it exists - everyone else keeps the untruncated draw.
			vLB_A = editmissing(vMND[idxASCT], 0)
			vSL_A = J(rows(idxASCT), 1, 1)
			iTr_A = selectindex(vLB_A :> 0)
			if (rows(iTr_A) > 0) {
				vSL_A[iTr_A] = calcSurvProb(vXB_ASCT[iTr_A], vLB_A[iTr_A], fbL1_TFI_ASCT, aux_ASCT)
				// Floor the survivor probability: a bound out in the far tail returns S ~ 0, and
				// U' = U * 0 inverts to +infinity for every family here.
				vSL_A[iTr_A] = rowmax((vSL_A[iTr_A], J(rows(iTr_A), 1, 1e-8)))
			}
			vRN_ASCT = vRN_ASCT :* vSL_A

			vOC[idxASCT] = calcSurvTime(vXB_ASCT, vRN_ASCT, fbL1_TFI_ASCT, aux_ASCT)
			
			// Curtail if beyond maximum observed - but NEVER below the truncation bound, or the
			// ceiling silently undoes the truncated draw above. A patient whose maintenance runs
			// past maxL1_TFI would otherwise be handed a gap shorter than the maintenance it is
			// supposed to contain, which is the bias reversing the draw order removed.
			vCap_A = rowmax((J(rows(idxASCT), 1, maxL1_TFI_ASCT), vLB_A))
			vOC[idxASCT] = rowmin((vOC[idxASCT], vCap_A))
		}
		
		// Group 2: No ASCT
		idxNoASCT = idx[selectindex(vSCT_L1[idx] :== 0)]
		if (rows(idxNoASCT) > 0) {

			//Grab BCR to L1
			vBCR_1 = (mBCR[idxNoASCT, 1] :== 1)
			vBCR_2 = (mBCR[idxNoASCT, 1] :== 2)
			vBCR_3 = (mBCR[idxNoASCT, 1] :== 3)
			vBCR_4 = (mBCR[idxNoASCT, 1] :== 4)
			vBCR_5 = (mBCR[idxNoASCT, 1] :== 5)	
			vBCR_6 = (mBCR[idxNoASCT, 1] :== 6)	
		
			// Assemble patient matrix (NoASCT includes BCR 5 and 6)
			mPat_NoASCT = (vAge[idxNoASCT], vAge2[idxNoASCT], vMale[idxNoASCT], 
						   vECOG0[idxNoASCT], vECOG1[idxNoASCT], vECOG2[idxNoASCT],
			               vRISS1[idxNoASCT], vRISS2[idxNoASCT], vRISS3[idxNoASCT], 
						   vMNT[idxNoASCT],
			               vBCR_1, vBCR_2, vBCR_3, vBCR_4, vBCR_5, vBCR_6, 
						   vCons[idxNoASCT])
			
			// Extract coefficients for NoASCT
			nPredictors = cols(mPat_NoASCT)
			if (cols(bL1_TFI_NoASCT) != nPredictors + 1) {
				errprintf("sim_tfi_l1 (NoASCT): design/coefficient mismatch - mPat has %g columns so %g were expected (mean + ancillary), but bL1_TFI_NoASCT has %g. A BCR level was likely empty in the fit, or the coefficient file was built by a different specification.\n",
					nPredictors, nPredictors + 1, cols(bL1_TFI_NoASCT))
				exit(459)
			}
			vCoef_NoASCT = bL1_TFI_NoASCT[1, 1..nPredictors]'
			aux_NoASCT = bL1_TFI_NoASCT[1, cols(bL1_TFI_NoASCT)]
			
			// Calculate XB
			vXB_NoASCT = mPat_NoASCT * vCoef_NoASCT
			
			// Calculate outcome (survival time)
			vRN_NoASCT = rnDraw(idxNoASCT, rn_tfi_l1(2))

			vLB_N = editmissing(vMND[idxNoASCT], 0)
			vSL_N = J(rows(idxNoASCT), 1, 1)
			iTr_N = selectindex(vLB_N :> 0)
			if (rows(iTr_N) > 0) {
				vSL_N[iTr_N] = calcSurvProb(vXB_NoASCT[iTr_N], vLB_N[iTr_N], fbL1_TFI_NoASCT, aux_NoASCT)
				vSL_N[iTr_N] = rowmax((vSL_N[iTr_N], J(rows(iTr_N), 1, 1e-8)))
			}
			vRN_NoASCT = vRN_NoASCT :* vSL_N

			vOC[idxNoASCT] = calcSurvTime(vXB_NoASCT, vRN_NoASCT, fbL1_TFI_NoASCT, aux_NoASCT)
			
			// Curtail if beyond maximum observed - never below the truncation bound (see the ASCT arm)
			vCap_N = rowmax((J(rows(idxNoASCT), 1, maxL1_TFI_NoASCT), vLB_N))
			vOC[idxNoASCT] = rowmin((vOC[idxNoASCT], vCap_N))
		}
		
		// Update matrices
		mTFI[idx, 2] = round(vOC[idx], 0.01)
		mTNE[idx, OMC] = round(vOC[idx], 0.01)
		mTSD[idx, OMC+1] = mTSD[idx, OMC] :+ mTNE[idx, OMC]
	}
}
