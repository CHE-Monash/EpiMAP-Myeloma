**********
* Monash Myeloma Model - Sim LenRefr (treatment lines)
*
* Purpose: update the LATCHED lenalidomide-refractory-from-treatment state (vLenRefr_in = the
*          LenRefr_Tx_in covariate). Fires at each line's END OMC, AFTER that line's OS, so the
*          state read by this line's TXR and OS is the value from STRICTLY PRIOR lines - which is
*          how LenRefr_Tx_in was fitted (docs/refractory.md 3.5 / 4).
*
*          The only event is the 0 -> 1 flip, and it can only happen to a not-yet-refractory
*          patient on a lenalidomide regimen:
*            BCR in {5,6}  -> refractory        (definitional, no equation)
*            BCR in {1-4}  -> Bernoulli(residual logit bLENREFR_TX)
*          Already-refractory patients stay 1 (latched); non-len lines and BCR-missing rows do
*          nothing. The line enters the logit collapsed to L1 / L2 / L3+ (min(Line, 3)).
*
* Fired by: core/simulation_engine.do at OMC 3,5,7,9,11 (L1E..L5E), i.e. after sim_os at each
*           line end. Line still holds the just-completed line (it increments at the next start).
* Reads:    vLenRefr_in (state), mTXR[.,Line] (regimen), mBCR[.,Line] (response), bLENREFR_TX,
*           LENREFR_regimens (len codes), the baseline patient vectors, rn_lenrefr(Line).
* Writes:   vLenRefr_in (latched 0 -> 1).
**********

mata {
	if (lenrefr_model_exists()) {

		bLR    = get_lenrefr_coef()       // 1 x 21 logit coefficients (factor-expanded, base = 0)
		vLenReg = get_lenrefr_regimens()  // row vector of len-containing regimen codes

		// Alive and reached this line (same filter as the other line-level sims)
		idx = selectindex((mMOR[., OMC-1] :== 0) :& (mState[., 1] :<= OMC))
		if (rows(idx) > 0) {

			// Snapshot the ENTRY-to-this-line state (before the update below) for export/validation.
			// This is LenRefr_Tx_in as at line entry - refractoriness from strictly prior lines.
			mLenRefr_in[idx, Line] = vLenRefr_in[idx]

			// Regimen on THIS line, and whether it is a lenalidomide code
			vReg   = mTXR[idx, Line]
			vIsLen = J(rows(idx), 1, 0)
			for (c = 1; c <= cols(vLenReg); c++) vIsLen = vIsLen :| (vReg :== vLenReg[1, c])

			// Eligible to flip: not yet refractory AND on a len line
			vElig = (vLenRefr_in[idx] :== 0) :& vIsLen

			// This line's response, and the definitional / residual split
			vB   = mBCR[idx, Line]
			vDef = (vB :== 5) :| (vB :== 6)             // SD/PD: definitionally refractory
			vRes = (vB :>= 1) :& (vB :<= 4)             // responders: the residual logit's arm

			// Residual-arm probability from the logit. Line dummies are constant across idx
			// (the whole block is one line); BCR dummies are per-patient. Column order MUST match
			// the fit: Age Age2 Male i.ECOGcc i.RISS CM(4) i.LENREFR_line i.BCR _cons.
			lg3 = min((Line, 3))
			vLR1 = J(rows(idx), 1, lg3 == 1)
			vLR2 = J(rows(idx), 1, lg3 == 2)
			vLR3 = J(rows(idx), 1, lg3 == 3)

			mPat = (vAge[idx], vAge2[idx], vMale[idx],
			        vECOG0[idx], vECOG1[idx], vECOG2[idx],
			        vRISS1[idx], vRISS2[idx], vRISS3[idx],
			        vCKD[idx], vCRD[idx], vPLM[idx], vDBT[idx],
			        vLR1, vLR2, vLR3,
			        (vB :== 1), (vB :== 2), (vB :== 3), (vB :== 4),
			        vCons[idx])

			// Guard: design columns must equal the coefficient count (plain logit, no ancillary)
			if (cols(mPat) != cols(bLR)) {
				errprintf("sim_lenrefr: design/coefficient mismatch at Line %g - mPat has %g columns but coefficient vector has %g\n", Line, cols(mPat), cols(bLR))
				exit(459)
			}

			vXB = mPat * bLR'
			vPR = 1 :/ (1 :+ exp(-vXB))

			// One CRN draw per line; residual-arm patients flip when p > u
			vRN      = rnDraw(idx, rn_lenrefr(Line))
			vResRefr = vRes :& (vPR :> vRN)

			// Acquire refractoriness this line (definitional or residual), among the eligible
			vAcq = vElig :& (vDef :| vResRefr)

			// Latch: 1 stays 1, eligible flips where acquired
			vLenRefr_in[idx] = vLenRefr_in[idx] :| vAcq
		}
	}
}
