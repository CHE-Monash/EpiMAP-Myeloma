**********
*SIM Age
* 
* Purpose: Update patient age and handle age limit deaths
* Outcome: Continuous age in years
**********
mata {
	// Filters
	vStateValid = (mState[.,1] :<= OMC+1)
	vWasAlive = (OMC == 1 ? J(rows(mMOR), 1, 0) : (mMOR[.,OMC-1] :== 0))  // No one alive before diagnosis
	vWasDead = (OMC == 1 ? J(rows(mMOR), 1, 0) : (mMOR[.,OMC-1] :== 1))
	vEligible = vStateValid :& vWasAlive
	
	// Get indices
	idxEligible = selectindex(vEligible)
	idxDead = selectindex(vStateValid :& vWasDead)
	
	// Update mCore age to missing for dead patients - remove after vectorisation
	if (rows(idxDead) > 0) {
		mCore[idxDead, cAge] = J(rows(idxDead), 1, .)
	}
	
	// Update age for alive patients
	if (rows(idxEligible) > 0 & OMC > 1) {
		// Age = previous age + time in previous state
		mAge[idxEligible, OMC] = mAge[idxEligible, OMC-1] :+ mTNE[idxEligible, OMC-1]
		
		// Update mCore - removed after vectorisation
		mCore[idxEligible, cAge] = mAge[idxEligible, OMC]
		
		// Check for patients exceeding age limit
		vExceedsLimit = (mAge[idxEligible, OMC] :> Limit) :& (mAge[idxEligible, OMC] :< .)
		idxExceeds = selectindex(vExceedsLimit)
		
		if (rows(idxExceeds) > 0) {
			// Map back to full index
			idxExceedsFull = idxEligible[idxExceeds]
			
			// Cap age at limit
			mAge[idxExceedsFull, OMC] = J(rows(idxExceedsFull), 1, Limit)
			
			// Mark as dead in PREVIOUS OMC (because age calculation happens AFTER previous state)
			mMOR[idxExceedsFull, OMC-1] = J(rows(idxExceedsFull), 1, 1)
			
			// Update mCore age to missing - remove after vectorisation
			mCore[idxExceedsFull, cAge] = J(rows(idxExceedsFull), 1, .)
			
			// Set outcome time (survival from diagnosis to death at age limit)
			mOC[idxExceedsFull, 1] = Limit :- mAge[idxExceedsFull, 1]
			mOC[idxExceedsFull, 2] = J(rows(idxExceedsFull), 1, 1)
			
			// Update mTFI or mTXD based on (OMC-1) parity (since death happened at previous OMC)
			vTimeDays = (mOC[idxExceedsFull, 1] :- mTSD[idxExceedsFull, OMC-1]) :* 365.25
			
			if (mod(OMC-1, 2) == 0) {
				// Even OMC-1: update CI
				mTFI[idxExceedsFull, (OMC-1)/2] = vTimeDays
			}
			else {
				// Odd OMC-1: update CD  
				mTXD[idxExceedsFull, OMC/2] = vTimeDays
			}
			
			// Set mTSD for current OMC to missing
			mTSD[idxExceedsFull, OMC] = J(rows(idxExceedsFull), 1, .)
			
			// Set mTNE for previous OMC to missing
			mTNE[idxExceedsFull, OMC-1] = J(rows(idxExceedsFull), 1, .)
		}
	}
}
