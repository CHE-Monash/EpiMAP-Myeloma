**********
* Monash Myeloma Model - Sim TXR
*
* Purpose: Select treatment regimen via multinomial logit. Outcome is the regimen number.
**********

mata {
	if (txr_model_exists(Line)) {
	// Model exists - run multinomial logit
			
	vCoef = get_txr_coef(Line)
	vTXR = get_txr_outcome(Line)
		
		// Initialise outcome
		vOC = J(Obs, 1, .)
		
		// Filter for alive and eligible
		idx = selectindex((mMOR[., OMC-1] :== 0) :& (mState[., 1] :<= OMC))
		if (rows(idx) > 0) {
		
			nRegimens = cols(vTXR)
			
			vXB1 = J(rows(idx), 1, 1)
			vXB2 = J(rows(idx), 1, 0)
			vXB3 = J(rows(idx), 1, 0)
			vXB4 = J(rows(idx), 1, 0)
			
			// Assemble patient matrix (L1) - Age + transplant status (SCT gates some regimens).
			// Must match the fit in prep/risk_equations.do: mlogit TXR_L1 Age Age2 SCT.
			if (Line == 1) {
				mPat = (vAge[idx], vAge2[idx], vSCT_DN[idx], vCons[idx])

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
			
			// Build patient matrix (L2+). L2-L4 carry the len-refractory flag: refractory status
			// drives regimen choice almost deterministically (docs/refractory.md 3.4), which Age
			// alone cannot reproduce. L5+ carry Age only - refractoriness is saturated there (5(3))
			// and only L2-L4 were fitted with it. Column order must match the fit:
			// mlogit TXR_L{2..4} Age Age2 LenRefr_any.
			else if (Line >= 2) {
				if (Line <= 4) mPat = (vAge[idx], vAge2[idx], vLenRefr_in[idx], vCons[idx])
				else           mPat = (vAge[idx], vAge2[idx], vCons[idx])

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
			
			// Calculate probabilities
			vPR1 = vXB1 :/ (vXB1 :+ vXB2 :+ vXB3 :+ vXB4)
			vPR2 = vPR1 :+ (vXB2 :/ (vXB1 :+ vXB2 :+ vXB3 :+ vXB4))
			vPR3 = vPR2 :+ (vXB3 :/ (vXB1 :+ vXB2 :+ vXB3 :+ vXB4))
			vPR4 = vPR3 :+ (vXB4 :/ (vXB1 :+ vXB2 :+ vXB3 :+ vXB4))
			
			// Determine Outcome
			vRN = rnDraw(idx, rn_txr(Line))
			vOC = (vRN :< vPR1) :* vTXR[1,1]

			if (nRegimens >= 2) {
				vOC = vOC :+ ((vRN :>= vPR1) :& (vRN :< vPR2)) :* vTXR[1,2]
			}
			
			if (nRegimens >= 3) {
				vOC = vOC :+ ((vRN :>= vPR2) :& (vRN :< vPR3)) :* vTXR[1,3]
			}
			
			if (nRegimens >= 4) {
				vOC = vOC :+ ((vRN :>= vPR3) :& (vRN :< vPR4)) :* vTXR[1,4]
			}
			
			// Update matrix
			mTXR[idx, Line] = vOC
		}	
	}
	else {
		// No model - assign pooled "Other" treatment
		idxAlive = selectindex(mMOR[., OMC-1] :== 0)
		idxDead = selectindex(mMOR[., OMC-1] :!= 0)
			
		// Assign "Other" (code 0) to alive patients
		if (rows(idxAlive) > 0) {
			mTXR[idxAlive, Line] = J(rows(idxAlive), 1, 0)
		}
			
		// Set missing for dead patients
		if (rows(idxDead) > 0) {
			mTXR[idxDead, Line] = J(rows(idxDead), 1, .)
		}
			
	}
}

// Check for override file, execute if it exists
mata: st_local("current_line", strofreal(Line))
if `current_line' == ${line} {
	local override_file "${outcomes_path}/sim_txr_override.do"
	capture confirm file "`override_file'"
	if _rc == 0 {
		qui do `override_file'
	}
}

