**********
*SIM MORT
* 
* Purpose: Determine which patients have died
* Outcome: Binary (0 = Alive, 1 = Dead)
**********
mata {
	// Filters
	vStateValid = (mState[.,1] :<= OMC+1)
	vWasAlive = (OMC == 1 ? J(rows(mMOR), 1, 1) : (mMOR[.,OMC-1] :== 0))
	vEligible = vStateValid :& vWasAlive
		
	// Mortality determination
	vDies = (mTSD[.,OMC+1] :>= mOS[.,OMC]) :& vEligible
	vLives = (mTSD[.,OMC+1] :< mOS[.,OMC]) :& vEligible
	
	// Get indices
	idxEligible = selectindex(vEligible)
	idxDies = selectindex(vDies)
	idxLives = selectindex(vLives)
		
	// Age limit: cap survival time if age at death exceeds maximum
	vExceedsLimit = ((mAge[.,1] :+ mOS[.,OMC]) :> Limit) :& vDies
	mOS[.,OMC] = vExceedsLimit :* (Limit :- mAge[.,1]) :+ (!vExceedsLimit) :* mOS[.,OMC]
		
	// Update mMOR - only update eligible patients
	if (rows(idxEligible) > 0) {
		mMOR[idxEligible, OMC] = vDies[idxEligible]
	}
		
	// Update mOC - only for deaths
	if (rows(idxDies) > 0) {
		mOC[idxDies, 1] = mOS[idxDies, OMC]
		mOC[idxDies, 2] = J(rows(idxDies), 1, 1)
	}
		
	// Update mTSD - set to missing for deaths
	if (rows(idxDies) > 0) {
		mTSD[idxDies, OMC+1] = J(rows(idxDies), 1, .)
	}
		
	// Update mTNE - only for eligible patients
	if (rows(idxEligible) > 0) {
		vTSDsafe = editmissing(mTSD[.,OMC], 0)
		mTNE[idxEligible, OMC] = vDies[idxEligible] :* (mOS[idxEligible, OMC] :- vTSDsafe[idxEligible]) :+ 
		                         vLives[idxEligible] :* mTNE[idxEligible, OMC]
	}
		
	// Update mTFI or mTXD - only if there are deaths
	if (rows(idxDies) > 0) {
		vTSDsafe = editmissing(mTSD[idxDies, OMC], 0)
		vTimeDays = (mOS[idxDies, OMC] :- vTSDsafe) :* 365.25
		
		if (mod(OMC, 2) == 0) {
			mTFI[idxDies, OMC/2] = vTimeDays
		}
		else {
			mTXD[idxDies, (OMC+1)/2] = vTimeDays
		}
	}
		
	// Update mSCT only if death before L1E (OMC <= 2)
	if (OMC <= 2 & rows(idxEligible) > 0) {
		mSCT[idxDies, 1] = J(rows(idxDies), 1, 0)
		mSCT[idxLives, 1] = mSCT[idxLives, 1]  // Keep existing values
	}
}
