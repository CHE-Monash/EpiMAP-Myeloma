**********
* Monash Myeloma Model - Sim OS
*
* Purpose: Draw stage-specific parametric Overall Survival for the current pathway stage (OMC).
*          Each stage has its own fitted model (own shape, demographic, comorbidity + BCR
*          effects) clocked from that stage's own entry event, so the draw is fresh survival
*          from that clock (elapsed = 0 at a line's first stage -> UNCONDITIONAL: no
*          from-diagnosis mTSD conditioning, which removed the heavy-TFI OS lift that
*          over-predicted weak responders). Stored on the DIAGNOSIS clock (originTSD +
*          residual) so sim_mort compares like-for-like. Continuous months to mOS[.,OMC].
* Notes:   L6+ is the exception: one model clocked once from L6 start, so its later stages
*          condition on survival since L6. Covariates per stage: Age Age2 Male i.ECOGcc i.RISS
*          CM_CKD CM_CRD CM_PLM CM_DBT + BCR, FULL factor levels (base levels carry 0 coeffs);
*          the four comorbidity flags are single 0/1 columns after RISS and before the BCR
*          block, matching the streg varlist order in prep/risk_equations.do. Firing map below.
**********

* Firing map (which fitted model, clock origin and BCR source fire at each OMC):
*
*   OMC  Stage         Model            origin mTSD col   BCR mBCR col
*   ---  ------------  ---------------  ---------------  ----------------
*    1   DN            OS_DN            1  (=0)          none
*    2   L1S           OS_L1S           2  (TSD_L1S)     1  (BCR_L1)
*    3   L1E No-ASCT   OS_L1E_NoASCT    3  (TSD_L1E)     1  (BCR_L1)
*    3   L1E ASCT      OS_L1E_ASCT      3  (TSD_L1E)     10 (BCR_SCT)
*    4   L2S           OS_L2S           4  (TSD_L2S)     2  (BCR_L2)
*    5   L2E           OS_L2E           5  (TSD_L2E)     2  (BCR_L2)
*    6   L3S           OS_L3S           6  (TSD_L3S)     3  (BCR_L3)
*    7   L3E           OS_L3E           7  (TSD_L3E)     3  (BCR_L3)
*    8   L4S           OS_L4S           8  (TSD_L4S)     4  (BCR_L4)
*    9   L4E           OS_L4E           9  (TSD_L4E)     4  (BCR_L4)
*   10   L5S           OS_L5S           10 (TSD_L5S)     5  (BCR_L5)
*   11   L5E           OS_L5E           11 (TSD_L5E)     5  (BCR_L5)
* 12-19  L6S .. L9E    OS_L6plus        12 (TSD_L6S)     Line (running BCR)
*
* Pattern: origin col = OMC (each stage's own clock, elapsed = 0 -> unconditional),
* BCR col = floor(OMC/2). Exceptions: DN (no BCR), L1E ASCT (BCR_SCT col 10), and L6+
* (all clocked from L6 start col 12 -> conditional, current line's running BCR).
* BCR_SCT has 4 levels (1-4) not 6; the BCR block width is read from each coefficient
* matrix (nBCR = cols - 15), so a 4- or 6-level block is handled automatically.

mata {
	// Alive-and-eligible filter (no alive-check at DN, where mMOR[.,OMC-1] does not exist)
	baseFilter = (OMC == 1 ? (mState[.,1] :<= OMC) : ((mMOR[.,OMC-1] :== 0) :& (mState[.,1] :<= OMC)))

	// One sub-fit per stage, except L1E which splits into No-ASCT / ASCT
	nSub = (OMC == 3 ? 2 : 1)

	for (s = 1; s <= nSub; s++) {

		// ---- Clock origin, BCR source and SCT filter (computed; only the model needs a branch) ----
		// Each stage clocks from its own mTSD column (= OMC); L6+ (OMC 12-19) shares L6 start, col 12.
		originCol = (OMC < 12 ? OMC : 12)

		// BCR source column in mBCR: none at DN, running BCR at L6+, else the line's own response.
		// (ASCT overrides to BCR_SCT col 10 below.)
		if (OMC == 1) bcrCol = 0
		else if (OMC >= 12) bcrCol = Line
		else bcrCol = floor(OMC / 2)

		sctFilter = J(Obs, 1, 1)   // all patients unless L1E splits by transplant

		// Select the fitted model + functional form for this stage
		if (OMC == 1) {
			vCoef = bOS_DN
			dist  = fbOS_DN
		}
		else if (OMC == 2) {
			vCoef = bOS_L1S
			dist  = fbOS_L1S
		}
		else if (OMC == 3 & s == 1) {
			vCoef = bOS_L1E_NoASCT
			dist  = fbOS_L1E_NoASCT
			sctFilter = (vSCT_L1 :== 0)
		}
		else if (OMC == 3 & s == 2) {
			vCoef = bOS_L1E_ASCT
			dist  = fbOS_L1E_ASCT
			bcrCol = 10
			sctFilter = (vSCT_L1 :== 1)
		}
		else if (OMC == 4) {
			vCoef = bOS_L2S
			dist  = fbOS_L2S
		}
		else if (OMC == 5) {
			vCoef = bOS_L2E
			dist  = fbOS_L2E
		}
		else if (OMC == 6) {
			vCoef = bOS_L3S
			dist  = fbOS_L3S
		}
		else if (OMC == 7) {
			vCoef = bOS_L3E
			dist  = fbOS_L3E
		}
		else if (OMC == 8) {
			vCoef = bOS_L4S
			dist  = fbOS_L4S
		}
		else if (OMC == 9) {
			vCoef = bOS_L4E
			dist  = fbOS_L4E
		}
		else if (OMC == 10) {
			vCoef = bOS_L5S
			dist  = fbOS_L5S
		}
		else if (OMC == 11) {
			vCoef = bOS_L5E
			dist  = fbOS_L5E
		}
		else {
			vCoef = bOS_L6plus
			dist  = fbOS_L6plus
		}

		idx = selectindex(baseFilter :& sctFilter)
		if (rows(idx) == 0) continue

		// ---- Design: full factor levels (base levels carry 0 coefficients) ----
		// Age, Age2, Male, ECOG(0,1,2), RISS(1,2,3), then the four comorbidity flags,
		// then the BCR block, then the constant. Order matches the streg varlist.
		mPat = (vAge[idx], vAge2[idx], vMale[idx],
				vECOG0[idx], vECOG1[idx], vECOG2[idx],
				vRISS1[idx], vRISS2[idx], vRISS3[idx],
				vCKD[idx], vCRD[idx], vPLM[idx], vDBT[idx])

		// BCR block width implied by the coefficient vector:
		// L2S..L4E (OMC 4..9) carry ONE combined len-refractory covariate after the BCR block. The
		// treatment- and maintenance-dose flags carry the same conditional hazard (test Tx = Mnt
		// p = 0.97), so they collapse to one and the engine keeps a single latched state. The
		// implied BCR width therefore drops by one there.
		hasLenRefr = (OMC >= 4 & OMC <= 9)

		//   cols - 13 covariates (Age,Age2,Male,ECOGx3,RISSx3,CMx4) - cons - aux - hasLenRefr
		nBCR = cols(vCoef) - 15 - hasLenRefr
		if (nBCR > 0 & bcrCol > 0) {
			vB = mBCR[idx, bcrCol]
			for (k = 1; k <= nBCR; k++) mPat = mPat, (vB :== k)
		}
		if (hasLenRefr) mPat = mPat, (vLenRefr_in[idx])
		mPat = mPat, vCons[idx]

		// Guard: design columns must equal coefficients minus the ancillary (1). Catches a
		// silent off-by-one if a factor's level count drifts (e.g. an MI-introduced category).
		aux = vCoef[1, cols(vCoef)]
		if (cols(mPat) != cols(vCoef) - 1) {
			errprintf("sim_os: design/coefficient mismatch at OMC %g (sub %g) - mPat has %g columns but coefficients imply %g predictors\n", OMC, s, cols(mPat), cols(vCoef) - 1)
			exit(459)
		}
		vBeta = vCoef[1, 1..(cols(vCoef) - 1)]'

		// ---- Linear predictor and conditional residual survival on the line clock ----
		vXB = mPat * vBeta
		vElapsed = mTSD[idx, OMC] :- mTSD[idx, originCol]   // months since this line's start (0 at a line's first stage)
		vPR = calcSurvProb(vXB, vElapsed, dist, aux)        // S(elapsed); = 1 when elapsed = 0
		vRN = rnDraw(idx, rn_os(OMC)) :* vPR                // condition on survival to now
		vT  = calcSurvTime(vXB, vRN, dist, aux)             // total survival on the line clock

		// Store on the diagnosis clock (originTSD + residual) for sim_mort
		mOS[idx, OMC] = round(mTSD[idx, originCol] :+ vT, 0.01)
	}
}
