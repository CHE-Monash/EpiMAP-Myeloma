**********
* Monash Myeloma Model - Sim MNR
*
* Purpose: Select the L1 maintenance regimen among patients who receive maintenance (MNT == 1).
*          Lenalidomide (1) and thalidomide (5) only (docs/refractory.md 7.4), so this is
*          effectively a logit with lenalidomide the base. 'Other' regimens are not modelled -
*          the engine never produces one. The maintenance analogue of sim_txr.do.
*
* Notes:   WHICH drug codes are modelled is an analysis-level choice, declared in
*          analyses/$analysis/outcomes/txr_$coeffs.do and applied by gen_mnr in
*          prep/risk_equations.do; anything unlisted was folded into 0 = 'other' at fit
*          time. That is what lets one extraction serve both the historical window the
*          out-of-sample validation scores against (thalidomide was the majority regimen
*          until 2020) and a current-paradigm window. See docs/refractory.md 7.4.
*
*          Consumed only by cost_tx_mnt in core/process_data.do, so this runs but does
*          nothing in the $line 2 analyses, where maintenance is never costed (7.3).
*
*          Must match the fit in prep/risk_equations.do:
*              mlogit MNR_L1 Age Age2 Male i.ECOGcc i.RISS SCT, baseoutcome(1)   (len + thal only)
**********

mata {
	if (mnr_model_exists()) {
	// Model exists - run multinomial logit

		vCoef = get_mnr_coef()
		vMNRout = get_mnr_outcome()

		// Initialise outcome
		vOC = J(Obs, 1, .)

		// Filter for alive, eligible AND receiving maintenance. Unlike sim_txr this is
		// conditional on another simulated outcome (vMNT, drawn immediately above).
		idx = selectindex((mMOR[., OMC-1] :== 0) :& (mState[., 1] :<= OMC) :& (vMNT :== 1))
		if (rows(idx) > 0) {

			nRegimens = cols(vMNRout)

			vXB1 = J(rows(idx), 1, 1)
			vXB2 = J(rows(idx), 1, 0)
			vXB3 = J(rows(idx), 1, 0)
			vXB4 = J(rows(idx), 1, 0)

			// Assemble patient matrix. ECOG/RISS dummies include the base level: the
			// coefficient vector carries 0b.ECOGcc / 0b.RISS, as in sim_mnt.do.
			mPat = (vAge[idx], vAge2[idx], vMale[idx],
					vECOG0[idx], vECOG1[idx], vECOG2[idx],
					vRISS1[idx], vRISS2[idx], vRISS3[idx],
					vSCT_L1[idx], vCons[idx])

			nPredictors = cols(mPat)

			// Guard: the coefficient vector must be nRegimens blocks of nPredictors (mlogit
			// carries an all-zero block for the base outcome). A mismatch means the design
			// here has drifted from the fit in risk_equations.do.
			if (cols(vCoef) != nRegimens * nPredictors) {
				errprintf("sim_mnr: design/coefficient mismatch - mPat has %g columns and there are %g outcomes (expected %g coefficients), but the coefficient vector has %g\n",
					nPredictors, nRegimens, nRegimens * nPredictors, cols(vCoef))
				exit(459)
			}

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

			// Calculate cumulative probabilities
			vPR1 = vXB1 :/ (vXB1 :+ vXB2 :+ vXB3 :+ vXB4)
			vPR2 = vPR1 :+ (vXB2 :/ (vXB1 :+ vXB2 :+ vXB3 :+ vXB4))
			vPR3 = vPR2 :+ (vXB3 :/ (vXB1 :+ vXB2 :+ vXB3 :+ vXB4))
			vPR4 = vPR3 :+ (vXB4 :/ (vXB1 :+ vXB2 :+ vXB3 :+ vXB4))

			// Determine outcome
			vRN = rnDraw(idx, rn_mnr())
			vOC = (vRN :< vPR1) :* vMNRout[1,1]

			if (nRegimens >= 2) {
				vOC = vOC :+ ((vRN :>= vPR1) :& (vRN :< vPR2)) :* vMNRout[1,2]
			}

			if (nRegimens >= 3) {
				vOC = vOC :+ ((vRN :>= vPR2) :& (vRN :< vPR3)) :* vMNRout[1,3]
			}

			if (nRegimens >= 4) {
				vOC = vOC :+ ((vRN :>= vPR3) :& (vRN :< vPR4)) :* vMNRout[1,4]
			}

			// Update vector. Patients with MNT == 0 keep the missing they were initialised
			// with in mata_setup.do - they have no maintenance regimen.
			vMNR[idx] = vOC
		}
	}
	else {
		// No mlogit was fitted. TWO DIFFERENT CASES, and conflating them is a silent data loss:
		//
		//   (a) the analysis declared exactly ONE maintenance regimen (car_t declares $MNR_L1 "1"),
		//       so there was nothing to choose between and risk_equations.do stored oL1_MNR instead
		//       of fitting. Every maintenance patient gets THAT regimen.
		//   (b) the analysis declared NO maintenance list at all, so maintenance is not modelled and
		//       everyone falls to 'other'. This is the $line 2 case (docs/refractory.md 7.3).
		//
		// This branch used to assign 0 in both cases, which broke (a) end to end and quietly: with
		// vMNR == 0 sim_mnd draws no duration (nothing is billed), and sim_mnt_refr fires only where
		// vMNR == 1, so the whole maintenance len-refractory arm produced nobody.
		vMNRout  = get_mnr_outcome()
		idxOther = selectindex((mMOR[., OMC-1] :== 0) :& (mState[., 1] :<= OMC) :& (vMNT :== 1))
		if (rows(idxOther) > 0) {
			if (cols(vMNRout) == 1) vMNR[idxOther] = J(rows(idxOther), 1, vMNRout[1,1])
			else                    vMNR[idxOther] = J(rows(idxOther), 1, 0)
		}
	}
}

// Check for override file, execute if it exists
local override_file "${outcomes_path}/sim_mnr_override.do"
capture confirm file "`override_file'"
if _rc == 0 {
	qui do `override_file'
}
