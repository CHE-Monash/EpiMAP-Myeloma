**********
*SIM Age
* 
* Purpose: Update patient age and handle age limit deaths
* Method: Comparing times
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
	
	// Update age for alive patients
	if (rows(idxEligible) > 0 & OMC > 1) {
		// Age = previous age + time in previous state
		vPrevAge = mAge[idxEligible, OMC-1]
		vPrevTNE = mTNE[idxEligible, OMC-1] / 12  // Convert months to years
		vNewAge = vPrevAge :+ vPrevTNE
		mAge[idxEligible, OMC] = round(vNewAge, 0.1)
		
		// Check for patients exceeding age limit
		vCurrentAges = mAge[idxEligible, OMC]
		vExceedsLimit = (vCurrentAges :> Limit) :& (vCurrentAges :< .)
		idxExceeds = selectindex(vExceedsLimit)
		
		if (rows(idxExceeds) > 0) {
			// Map back to full index
			idxExceedsFull = idxEligible[idxExceeds]
			
			// Cap age at limit
			mAge[idxExceedsFull, OMC] = J(rows(idxExceedsFull), 1, Limit)
			
			// Mark as dead in PREVIOUS OMC
			mMOR[idxExceedsFull, OMC-1] = J(rows(idxExceedsFull), 1, 1)
			
			// Set outcome time in months (survival from diagnosis to death at age limit)
			vDiagAge = mAge[idxExceedsFull, 1]
			mOC[idxExceedsFull, 1] = (Limit :- vDiagAge) * 12 
			mOC[idxExceedsFull, 2] = J(rows(idxExceedsFull), 1, 1)
			
			// Calculate time for TXD/TFI in months
			vPrevTSD = mTSD[idxExceedsFull, OMC-1]
			vOCTime = mOC[idxExceedsFull, 1]
			vTimeMonths = vOCTime :- vPrevTSD 
			
			if (mod(OMC-1, 2) == 0) {
				// Even OMC-1: update TXD
				lineIdx = floor((OMC-1)/2)  // Ensure integer index
				
				// Bounds check: TXD has 9 columns (L1-L9)
				if (lineIdx >= 1 & lineIdx <= 9) {
					mTXD[idxExceedsFull, lineIdx] = vTimeMonths
				}
			}
			else {
				// Odd OMC-1: update TFI
				lineIdx = floor((OMC+1)/2)  // Ensure integer index
				
				// Bounds check: TFI has 9 columns (DN, L1-L8)
				if (lineIdx >= 1 & lineIdx <= 9) {
					mTFI[idxExceedsFull, lineIdx] = vTimeMonths
				}
			}
			
			// Set mTSD for current OMC to missing
			mTSD[idxExceedsFull, OMC] = J(rows(idxExceedsFull), 1, .)
			
			// Set mTNE for previous OMC to missing
			mTNE[idxExceedsFull, OMC-1] = J(rows(idxExceedsFull), 1, .)
		}
	}
}
