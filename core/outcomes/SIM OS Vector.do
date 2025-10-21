**********
* SIM OS - Vectorised Implementation
* 
* Purpose: Calculate Overall Survival for any line of therapy
* Outcome: Continuous survival time (in days)
**********

mata {  
    nObs = st_nobs()
    
    // Get current BCR
    currentBCR = (Line == 0) ? J(nObs, 1, 5) : mBCR[., Line]
    
    // Create BCR dummy matrix once (nObs x 6)
    bcrDummies = ((currentBCR :== 1), (currentBCR :== 2), (currentBCR :== 3),
                  (currentBCR :== 4), (currentBCR :== 5), (currentBCR :== 6))
    
    // Initialize outcome
    oOS = J(nObs, 1, .)
    
    // Define segments to calculate
    // Segment maps to BCR coefficient position: BCR_start = 10 + segment*6
    // Segments: 0=L0, 1=L1noSCT, 2=L1SCT, 3=L2, 4=L3, 5=L4, 6=L5, 7=L6
    
    segments = J(0, 1, .)  // Will hold list of segments to calculate
    
    if (Line == 0) {
        segments = 0
    }
    else if (Line == 1) {
        segments = (1 \ 2)  // Both noSCT and SCT
    }
    else if (Line <= 6) {
        segments = Line + 1  // L2=3, L3=4, L4=5, L5=6, L6=7
    }
    
    // Loop through segments
    for (s = 1; s <= rows(segments); s++) {
        segment = segments[s]
        
        // Calculate BCR coefficient start position
        bcrStart = 10 + segment * 6
        
        // Build coefficient column vector
        coefCols = (1, 2, 3, 5, 6, 8, 9,                           // Base effects
                    bcrStart, bcrStart+1, bcrStart+2,              // BCR 1-3
                    bcrStart+3, bcrStart+4, bcrStart+5,            // BCR 4-6
                    58)                                            // Constant
        
        // Determine which patients for this segment
        if (segment == 0) {
            // Line 0 (diagnosis)
            idx = selectindex(mState[., 1] :<= OMC + 1)
        }
        else if (segment == 1) {
            // Line 1, no SCT
            idx = selectindex(mState[., 1] :<= OMC + 1 :& vSCT_DN :== 0)
        }
        else if (segment == 2) {
            // Line 1, SCT
            idx = selectindex(mState[., 1] :<= OMC + 1 :& vSCT_DN :== 1)
        }
        else {
            // Line 2+
            idx = selectindex(mState[., 1] :<= OMC + 1)
        }
        
        // Calculate for this segment if patients exist
        if (rows(idx) > 0) {
            // Assemble patient matrix for this segment
            pOS = (vAge[idx], vAge2[idx], vMale[idx], 
                       vECOG1[idx], vECOG2[idx], vRISS2[idx], vRISS3[idx],
                       bcrDummies[idx, 1], bcrDummies[idx, 2], bcrDummies[idx, 3],
                       bcrDummies[idx, 4], bcrDummies[idx, 5], bcrDummies[idx, 6],
                       vCons[idx])
            
            // Extract coefficients
            coefOS = bOS[1, coefCols]'
            
            // Calculate XB
            xbOS = pOS * coefOS
            
            // Generate random numbers
            rnOS = runiform(rows(idx), 1)
            
            // Calculate survival time
            oOS[idx] = calcSurvTime(xbOS, rnOS, fbOS, bOS[1, cols(bOS)])
        }
    }
    
    // Update matrix
    mOS[., OMC] = oOS
}
