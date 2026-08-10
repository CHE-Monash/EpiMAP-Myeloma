**********
* Monash Myeloma Model - Sim LenRefr (treatment lines)
*
* Purpose: update the LATCHED lenalidomide-refractory-from-treatment state (vLenRefr_in = the
*          refr_len_tx_in covariate). Fires at each line's END OMC, AFTER that line's OS, so the
*          state read by this line's TXR and OS is the value from STRICTLY PRIOR lines - which is
*          how refr_len_tx_in was fitted (docs/refractory.md 3.5 / 4).
*
*          The only event is the 0 -> 1 flip, and it can only happen to a not-yet-refractory
*          patient, with probability
*            P(acquire) = P(lenalidomide exposure) x  1              if BCR in {5,6}  (definitional)
*                                                     logit p       if BCR in {1-4}  (residual)
*          Already-refractory patients stay 1 (latched). The line enters the logit collapsed to
*          L1 / L2 / L3+ (min(Line, 3)).
*
*          EXPOSURE IS PROBABILISTIC for the 'other' regimen bucket. A modelled lenalidomide code
*          gives exposure 1; TXR code 0 gives p_line, the registry's P(len | regimen not modelled).
*          Treating 'other' as non-len made this arm under-generate about four-fold (2.2% of L2
*          patients against a registry ~8.3%) and hid exactly the patients the OS-by-refractory
*          gate needs - at L1 ten years older, a third as often transplanted and 2.4x as likely to
*          be refractory. See scratch/lenrefr_other.log and scratch/refractory/_notes.md.
*
* Fired by: core/simulation_engine.do at OMC 3,5,7,9,11 (L1E..L5E), i.e. after sim_os at each
*           line end. Line still holds the just-completed line (it increments at the next start).
* Reads:    vLenRefr_in (state), mTXR[.,Line] (regimen), mBCR[.,Line] (response), bLENREFR_TX,
*           LENREFR_regimens (len codes), LENREFR_pother (p_line), the baseline patient vectors,
*           rn_lenrefr(Line).
* Writes:   vLenRefr_in (latched 0 -> 1).
**********

mata {
	if (lenrefr_model_exists()) {

		bLR    = get_lenrefr_coef()       // 1 x 21 logit coefficients (factor-expanded, base = 0)
		vLenReg = get_lenrefr_regimens()  // row vector of len-containing regimen codes
		vPOth   = get_lenrefr_pother()    // P(len | regimen not modelled), one column per line

		// SNAPSHOT and DRAW are two different populations, deliberately.
		//
		// The snapshot records the state a patient ENTERED this line with - refractoriness from
		// strictly prior lines, fixed before this line began. It is therefore well defined for a
		// patient who dies DURING the line, and must NOT be gated on survival. It used to share the
		// draw's index, which left refr_len_ll missing for everyone who died in line l - 11.4% at L2 -
		// and mata_setup.do coerces missing to 0. A line-L cohort analysis re-simulating those
		// patients from entry therefore admitted them as NON-refractory. They are the sickest
		// entrants (worse ECOG, RISS and L1 response, transplanted 23% against 30%) and so the most
		// likely to be refractory, making the loss differential rather than noise.
		//
		// Reached this line, alive or not. Same test sim_asct_dn.do uses.
		idxSnap = selectindex(mState[., 1] :<= OMC)
		if (rows(idxSnap) > 0) mLenRefr_in[idxSnap, Line] = vLenRefr_in[idxSnap]

		// The DRAW is gated on survival - only a living patient can flip (same filter as the
		// other line-level sims).
		idx = selectindex((mMOR[., OMC-1] :== 0) :& (mState[., 1] :<= OMC))
		if (rows(idx) > 0) {

			// Regimen on THIS line, and the PROBABILITY it contained lenalidomide.
			//
			// A modelled len code is certain (1). The 'other' bucket (code 0) is a MIXTURE - it holds
			// every regimen the analysis does not model, and at L1 half of it is lenalidomide, at L4
			// a third, while the L4/L5+ modelled lists contain no len regimen at all. Treating it as
			// non-len made the treatment arm under-generate about four-fold, and the patients it hid
			// are exactly the ones the OS gate needs: at L1 ten years older, a third as often
			// transplanted, 2.4x as likely to be refractory (scratch/lenrefr_other.log).
			//
			// So 'other' contributes p_line rather than 0, marginalising over the bucket's unobserved
			// drug content. TXR and TXD are untouched - this changes only what the GATE counts, which
			// is why it avoids the TXD regression that reverted the add-Rd-to-the-lists attempt.
			vReg  = mTXR[idx, Line]
			vPLen = J(rows(idx), 1, 0)
			for (c = 1; c <= cols(vLenReg); c++) vPLen = vPLen :+ (vReg :== vLenReg[1, c])
			if (cols(vPOth) >= Line) vPLen = vPLen :+ (vReg :== 0) :* vPOth[1, Line]

			// Not yet refractory - the only patients who can flip (the state is latched)
			vElig = (vLenRefr_in[idx] :== 0)

			// This line's response, and the definitional / residual split
			vB   = mBCR[idx, Line]
			vDef = (vB :== 5) :| (vB :== 6)             // SD/PD: definitionally refractory
			vRes = (vB :>= 1) :& (vB :<= 4)             // responders: the residual logit's arm

			// Residual-arm probability from the logit. Line dummies are constant across idx
			// (the whole block is one line); BCR dummies are per-patient. Column order MUST match
			// the fit: Age Age2 Male i.ECOGcc i.RISS CM(4) i.refr_len_grp i.BCR _cons.
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

			// Acquisition probability, conditional on lenalidomide exposure:
			//   definitional arm (BCR 5/6) - certain GIVEN exposure, so probability 1
			//   residual arm    (BCR 1-4)  - the logit's p
			// multiplied by P(exposure). For a modelled len regimen vPLen is 1 and this reduces
			// exactly to the previous behaviour; for 'other' it scales by p_line. The definitional
			// arm becoming probabilistic is correct rather than a side effect - a progressive-disease
			// patient on an unknown regimen is only lenalidomide-refractory if they actually had it.
			vPAcq = vPLen :* (vDef :+ vRes :* vPR)

			// One CRN draw per line (same slot as before - no rn_K change); flip when p > u
			vRN  = rnDraw(idx, rn_lenrefr(Line))
			vAcq = vElig :& (vPAcq :> vRN)

			// Latch: 1 stays 1, eligible flips where acquired
			vLenRefr_in[idx] = vLenRefr_in[idx] :| vAcq
		}
	}
}
