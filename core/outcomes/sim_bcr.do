**********
*SIM BCR Vector

* Purpose: Vectorised Best Clinical Response (BCR) calculation for all lines
*          of therapy (L1-L9)
*
* Key Features:
*   - Single implementation works across all therapy lines
*   - Handles 6 categories (L1-L2): CR, nCR, VGPR, MR, PR, SD/PD
*   - Handles 3 categories (L3-L9): CR, VGPR, PR
*   - Uses vectorised matrix multiplication (no patient loops)
*   - Dynamically adjusts predictors based on line number
*
* Categories:
* 6-category (L1-L2): 1=CR, 2=VGPR, 3=PR, 4=MR, 5=SD, 6=PD
* 3-category (L3+): 1=CR/VGPR, 3=PR/MR, 5=SD/PD
*
* Author: Modernised vectorised implementation
* Date: October 2025
**********

mata:
	// Determine structure
	if (Line <= 2) {
		nCategories = 6
		categoryValues = (1, 2, 3, 4, 5, 6)
		nCutPoints = 5
	}
	else {
		nCategories = 3
		categoryValues = (1, 3, 5)
		nCutPoints = 2
	}
	
	// Extract previous BCR
	if (Line >= 2) {
		if (Line == 2) {
			// L2: BCR_L1 for all, plus BCR_SCT interactions for ASCT patients
			BCR_L1 = mBCR[., 1]
			pBCR_VGPR = (BCR_L1 :== 2)
			pBCR_PR = (BCR_L1 :== 3)
			pBCR_MR = (BCR_L1 :== 4)
			pBCR_SD = (BCR_L1 :== 5)
			pBCR_PD = (BCR_L1 :== 6)
			
			BCR_SCT = mBCR[., 10]
			pBCR_SCT_1 = (BCR_SCT :== 1)
			pBCR_SCT_2 = (BCR_SCT :== 2)
			pBCR_SCT_3 = (BCR_SCT :== 3)
			pBCR_SCT_4 = (BCR_SCT :== 4)
		}
		else {
			// L3+: Use previous line BCR
			prevLineIdx = Line - 1
			prevBCR = mBCR[., prevLineIdx]
			pBCR_PR = (prevBCR :== 3)
			pBCR_SD = (prevBCR :== 5)
		}
	}
	
	// Extract TXR dummies
	if (Line <= 4) {
		if (Line == 1) oVector = oL1_TXR
		if (Line == 2) oVector = oL2_TXR
		if (Line == 3) oVector = oL3_TXR
		if (Line == 4) oVector = oL4_TXR
		
		currentTX = mTXR[., Line]
		
		if (cols(oVector) >= 2) TXR_is_R2 = (currentTX :== oVector[1, 2])
		if (cols(oVector) >= 3) TXR_is_R3 = (currentTX :== oVector[1, 3])
		if (cols(oVector) >= 4) TXR_is_R4 = (currentTX :== oVector[1, 4])
		
		nRegimens = cols(oVector)
	}
	else {
		nRegimens = 0
	}
	
	// Build patient matrix
	pMatrix = (vAge, vAge2, vMale, vECOG0, vECOG1, vECOG2, vRISS1, vRISS2, vRISS3)
	
	// Add SCT
	if (Line == 1) pMatrix = pMatrix, vSCT_DN
	else pMatrix = pMatrix, vSCT_L1
	
	// Add TXR dummies
	if (Line <= 4) {
		if (Line == 1) oVector = oL1_TXR
		if (Line == 2) oVector = oL2_TXR
		if (Line == 3) oVector = oL3_TXR
		if (Line == 4) oVector = oL4_TXR
		
		currentTX = mTXR[., Line]
		
		if (cols(oVector) >= 1) pMatrix = pMatrix, (currentTX :== oVector[1, 1])
		if (cols(oVector) >= 2) pMatrix = pMatrix, (currentTX :== oVector[1, 2])
		if (cols(oVector) >= 3) pMatrix = pMatrix, (currentTX :== oVector[1, 3])
		if (cols(oVector) >= 4) pMatrix = pMatrix, (currentTX :== oVector[1, 4])
	}
	
	// Add previous BCR
	if (Line == 2) {
		pMatrix = pMatrix, (pBCR_VGPR, pBCR_PR, pBCR_MR, pBCR_SD, pBCR_PD)
		pMatrix = pMatrix, (pBCR_SCT_1 :* vSCT_L1, pBCR_SCT_2 :* vSCT_L1, 
		                    pBCR_SCT_3 :* vSCT_L1, pBCR_SCT_4 :* vSCT_L1)
	}
	else if (Line >= 3) {
		pMatrix = pMatrix, (pBCR_PR, pBCR_SD)
	}
	
	nPredictors = cols(pMatrix)
	
	// Extract coefficients
	if (Line == 1) bMatrix = bL1_BCR
	if (Line == 2) bMatrix = bL2_BCR
	if (Line == 3) bMatrix = bL3_BCR
	if (Line == 4) bMatrix = bL4_BCR
	if (Line >= 5) bMatrix = bLX_BCR
	
	coefVector = bMatrix[1, 1..nPredictors]'
	cutPointIndices = (cols(bMatrix) - nCutPoints + 1)..cols(bMatrix)
	cutPoints = bMatrix[1, cutPointIndices]
	
	// Calculate XB
	XB = pMatrix * coefVector
	
	// Calculate probabilities
	cumProbs = calcOrdLogitProbs(XB, cutPoints)
	
	// Assign outcomes
	RN = runiform(nObs, 1)
	vOutcome = assignOrdOutcome(RN, cumProbs, categoryValues)
	
	// Update matrices
	for (i = 1; i <= nObs; i++) {
		if (mMOR[i, OMC - 1] == 0 & mState[i, 1] <= OMC + 1) {
			finalOutcome = vOutcome[i]
		}
		else if (mState[i, 1] > OMC + 1) {
			finalOutcome = mBCR[i, LX]
		}
		else {
			finalOutcome = .
		}
			
		if (mMOR[i, OMC - 1] == 0) {
			mBCR[i, LX] = finalOutcome
		}
	}
end
