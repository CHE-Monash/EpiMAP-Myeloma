**********
* SIM OS - Vectorised Implementation
* 
* Purpose: Calculate Overall Survival for any line of therapy
* Outcome: Continuous survival time (in days)
**********

mata {  
    nObs = rows(vAge)
    
	// Determine which column of mBCR to use based on OMC
	if (OMC == 1 | OMC == 2) {
		currentBCR = J(nObs, 1, 5)
		bcrCol = 0
	}
	else if (mod(OMC, 2) == 0) {
		bcrCol = floor((OMC - 2) / 2)
		currentBCR = mBCR[., bcrCol]
	}
	else {
		bcrCol = floor((OMC - 1) / 2)
		currentBCR = mBCR[., bcrCol]
	}
    
	// Create BCR dummy variables
	bcr1 = (currentBCR :== 1)
	bcr2 = (currentBCR :== 2)
	bcr3 = (currentBCR :== 3)
	bcr4 = (currentBCR :== 4)
	bcr5 = (currentBCR :== 5)
	bcr6 = (currentBCR :== 6)
	
	// Initialize outcome
    vOS = J(nObs, 1, .)
	
	// Determine coefficient segments
	// Segment maps to BCR coefficient position: BCR_start = 10 + segment*6	
	if (OMC <= 2) { // DN or L1S
		segments = 0
	}
	else if (OMC == 3) { // L1E 
		segments = (1, 2)
	}
	else {
		lineNum = floor((OMC - 2) / 2)
		segments = lineNum + 2
	}
 
    // Loop through segments
    for (s = 1; s <= cols(segments); s++) {
        segment = segments[s]
        
        // Calculate BCR coefficient start position
        bcrStart = 10 + segment * 6
        
        // Build coefficient column vector
        coefCols = (1, 2, 3, 5, 6, 8, 9,                           // Base effects
                    bcrStart, bcrStart+1, bcrStart+2,              // BCR 1-3
                    bcrStart+3, bcrStart+4, bcrStart+5,            // BCR 4-6
                    58)                                            // Constant
        
		// Determine patients for this segment
		if (segment == 0) { // DN
			if (OMC == 1) {
				idx = selectindex(mState[., 1] :<= OMC + 1)
			}
			else {
				idx = selectindex((mMOR[., OMC-1] :== 0) :& (mState[., 1] :<= OMC + 1))
			}
		}
		else if (segment == 1) { // L1 No ASCT
			idx = selectindex((mMOR[., OMC-1] :== 0) :& 
			                  (mState[., 1] :<= OMC + 1) :& 
			                  (vSCT_DN :== 0))
		}
		else if (segment == 2) { // L1 ASCT
			idx = selectindex((mMOR[., OMC-1] :== 0) :& 
			                  (mState[., 1] :<= OMC + 1) :& 
			                  (vSCT_DN :== 1))
		}
		else {
			idx = selectindex((mMOR[., OMC-1] :== 0) :& (mState[., 1] :<= OMC+1))
		}
        
        // Calculate for this segment if patients exist
        if (rows(idx) > 0) {
            // Assemble patient matrix for this segment
			pOS = (vAge[idx], vAge2[idx], vMale[idx], 
			       vECOG1[idx], vECOG2[idx], vRISS2[idx], vRISS3[idx],
			       bcr1[idx], bcr2[idx], bcr3[idx],
			       bcr4[idx], bcr5[idx], bcr6[idx],
			       vCons[idx])
            
            // Extract coefficients
            coefOS = bOS[1, coefCols]'
            
            // Calculate XB
            xbOS = pOS * coefOS
            
            // Generate random numbers
            rnOS = runiform(rows(idx), 1)
            
            // Calculate survival time
            vOS[idx] = calcSurvTime(xbOS, rnOS, fbOS, bOS[1, cols(bOS)])
        }
    }
    
    // Update matrix
    mOS[., OMC] = vOS
}
