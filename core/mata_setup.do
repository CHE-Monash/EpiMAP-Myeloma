**********
* EpiMAP Myeloma - Mataq Setup
*
* Architecture:
*   - Vectors: Patient characteristics that don't change (Age, Male, ECOG, RISS)
*   - Matrices: Line-varying outcomes preserved in columns (mBCR, mCD, mOS, etc.)
*   - Dynamic assembly: Combine vectors as needed for specific calculations
*
* Version: 2.1.0
* Date: October 2025
* Author: Adam Irving + Claude
**********

capture program drop mata_setup
program define mata_setup
	
	di "Setting up matrices"
	mata: Obs = st_nobs()
	
	*mCore - Core patient characteristics
{
    mata {
        mCore = st_data(.,"ID Male ECOGcc RISS CMc")
        rmCore = J(Obs, 2, "")
        rmCore[., 2] = strofreal(1::Obs)
        cmCore = J(5, 2, "")
        cmCore[1,1] = "ID"
        cmCore[2,1] = "Male"
        cmCore[3,1] = "ECOGcc"
        cmCore[4,1] = "RISS"
		cmCore[5,1] = "CMc"
    }
}

	*mState - State matrix
{
	mata {
		mState = st_data(.,"State DateDN")
		rmState = J(Obs,2,"")
		rmState[.,2] = strofreal(mState[.,1])
		cmState = J(2,2,"")
		cmState[1,1] = "State"
		cmState[2,1] = "DateDN"
	}	
}	
	*_matrix_list(mState, rmState, cmState)

	*mCom - Comorbidity matrix
{
	mata {
		mCom = st_data(.,"Age70 Age75 CMc")
		rmCom = J(Obs,2,"")
		rmCom[.,2] = strofreal(1::Obs)
		cmCom = J(3,2,"")
		cmCom[1,1] = "Age70" 
		cmCom[2,1] = "Age75"
		cmCom[3,1] = "CMc" 
	}
}
	*mata: _matrix_list(mCom, rmCom, cmCom)
		
	*mSCT
{
	mata {
		mSCT = st_data(., "SCT SCT") // Currently using same variable twice
		rmSCT = J(Obs,2,"")
		rmSCT[.,2] = strofreal(1::Obs)
		cmSCT = J(2,2,"")
		cmSCT[1,1] = "SCT_DN"
		cmSCT[2,1] = "SCT_L1"
	}
}		
	*mata: _matrix_list(mSCT, rmSCT, cmSCT)
		
	*mMNT
{
	mata {
		mMNT = st_data(.,"MNT")
		rmMNT = J(Obs,2,"")
		rmMNT[.,2] = strofreal(1::Obs)
		cmMNT = J(1,2,"")
		cmMNT[1,1] = "MNT"
	}
}		
	*mata: _matrix_list(mMNT, rmMNT, cmMNT)	
		
	*mCons
{
	gen Cons = 1
	mata {
		mCons = st_data(.,"Cons")
		rmCons = J(Obs, 2,"")
		rmCons[.,2] = strofreal(1::Obs)
		cmCons = J(1,2,"")
		cmCons[1,1] = "Cons"
	}
}		
		*mata: _matrix_list(mCons, rmCons, cmCons)	
	
	*mAge - Age (at event)
{	
	mata {
		mAge = st_data(., "Age_DN Age_L1S Age_L1E Age_L2S Age_L2E Age_L3S Age_L3E Age_L4S Age_L4E Age_L5S Age_L5E Age_L6S Age_L6E Age_L7S Age_L7E Age_L8S Age_L8E Age_L9S Age_L9E")
		rmAge = J(Obs, 2, "")
		rmAge[., 2] = strofreal(1::Obs)
		cmAge = J(19, 2, "")
		cmAge[1,1] = "Age_DN"
		cmAge[2,1] = "Age_L1S"
		cmAge[3,1] = "Age_L1E"
		cmAge[4,1] = "Age_L2S"
		cmAge[5,1] = "Age_L2E"
		cmAge[6,1] = "Age_L3S"
		cmAge[7,1] = "Age_L3E"
		cmAge[8,1] = "Age_L4S"
		cmAge[9,1] = "Age_L4E"
		cmAge[10,1] = "Age_L5S"
		cmAge[11,1] = "Age_L5E"
		cmAge[12,1] = "Age_L6S"
		cmAge[13,1] = "Age_L6E"
		cmAge[14,1] = "Age_L7S"
		cmAge[15,1] = "Age_L7E"
		cmAge[16,1] = "Age_L8S"
		cmAge[17,1] = "Age_L8E"
		cmAge[18,1] = "Age_L9S"
		cmAge[19,1] = "Age_L9E"
	}
}		
	*mata: _matrix_list(mAge, rmAge, cmAge)		

	*mOS - Overall Survival (from DN)
{	
	mata {
		mOS = J(Obs,19,.)
		rmOS = J(Obs, 2, "")
		rmOS[., 2] = strofreal(1::Obs)
		cmOS = J(19, 2, "")
		cmOS[1,1] = "OS_DN"
		cmOS[2,1] = "OS_L1S"
		cmOS[3,1] = "OS_L1E"
		cmOS[4,1] = "OS_L2S"
		cmOS[5,1] = "OS_L2E"
		cmOS[6,1] = "OS_L3S"
		cmOS[7,1] = "OS_L3E"
		cmOS[8,1] = "OS_L4S"
		cmOS[9,1] = "OS_L4E"
		cmOS[10,1] = "OS_L5S"
		cmOS[11,1] = "OS_L5E"
		cmOS[12,1] = "OS_L6S"
		cmOS[13,1] = "OS_L6E"
		cmOS[14,1] = "OS_L7S"
		cmOS[15,1] = "OS_L7E"
		cmOS[16,1] = "OS_L8S"
		cmOS[17,1] = "OS_L8E"
		cmOS[18,1] = "OS_L9S"
		cmOS[19,1] = "OS_L9E"
	}
}			
	*mata: _matrix_list(mOS, rmOS, cmOS)		

	*mTNE - Time to Next Event
{	
	mata {
		mTNE = st_data(., "TNE_DN TNE_L1S TNE_L1E TNE_L2S TNE_L2E TNE_L3S TNE_L3E TNE_L4S TNE_L4E TNE_L5S TNE_L5E TNE_L6S TNE_L6E TNE_L7S TNE_L7E TNE_L8S TNE_L8E TNE_L9S TNE_L9E")
		rmTNE = J(Obs, 2, "")
		rmTNE[., 2] = strofreal(1::Obs)
		cmTNE = J(19, 2, "")
		cmTNE[1,1] = "TNE_DN"
		cmTNE[2,1] = "TNE_L1S"
		cmTNE[3,1] = "TNE_L1E"
		cmTNE[4,1] = "TNE_L2S"
		cmTNE[5,1] = "TNE_L2E"
		cmTNE[6,1] = "TNE_L3S"
		cmTNE[7,1] = "TNE_L3E"
		cmTNE[8,1] = "TNE_L4S"
		cmTNE[9,1] = "TNE_L4E"
		cmTNE[10,1] = "TNE_L5S"
		cmTNE[11,1] = "TNE_L5E"
		cmTNE[12,1] = "TNE_L6S"
		cmTNE[13,1] = "TNE_L6E"
		cmTNE[14,1] = "TNE_L7S"
		cmTNE[15,1] = "TNE_L7E"
		cmTNE[16,1] = "TNE_L8S"
		cmTNE[17,1] = "TNE_L8E"
		cmTNE[18,1] = "TNE_L9S"
		cmTNE[19,1] = "TNE_L9E"
	}
}			
	*mata: _matrix_list(mTNE, rmTNE, cmTNE)	
		
	*mTSD - Time Since Diagnosis
{	
	mata {
		mTSD = st_data(., "TSD_DN TSD_L1S TSD_L1E TSD_L2S TSD_L2E TSD_L3S TSD_L3E TSD_L4S TSD_L4E TSD_L5S TSD_L5E TSD_L6S TSD_L6E TSD_L7S TSD_L7E TSD_L8S TSD_L8E TSD_L9S TSD_L9E")
		mTSD[., 1] = J(Obs, 1, 0) // Set to TSD_DN to 0
		rmTSD = J(Obs, 2, "")
		rmTSD[., 2] = strofreal(1::Obs)
		cmTSD = J(19, 2, "")
		cmTSD[1,1] = "TSD_DN"
		cmTSD[2,1] = "TSD_L1S"
		cmTSD[3,1] = "TSD_L1E"
		cmTSD[4,1] = "TSD_L2S"
		cmTSD[5,1] = "TSD_L2E"
		cmTSD[6,1] = "TSD_L3S"
		cmTSD[7,1] = "TSD_L3E"
		cmTSD[8,1] = "TSD_L4S"
		cmTSD[9,1] = "TSD_L4E"
		cmTSD[10,1] = "TSD_L5S"
		cmTSD[11,1] = "TSD_L5E"
		cmTSD[12,1] = "TSD_L6S"
		cmTSD[13,1] = "TSD_L6E"
		cmTSD[14,1] = "TSD_L7S"
		cmTSD[15,1] = "TSD_L7E"
		cmTSD[16,1] = "TSD_L8S"
		cmTSD[17,1] = "TSD_L8E"
		cmTSD[18,1] = "TSD_L9S"
		cmTSD[19,1] = "TSD_L9E"
	}
}			
	*mata: _matrix_list(mTSD, rmTSD, cmTSD)		
		
	*mMOR - Mortality (=1 if patient dies before next event) 
{	
	mata {
		mMOR = st_data(., "MOR_DN MOR_L1S MOR_L1E MOR_L2S MOR_L2E MOR_L3S MOR_L3E MOR_L4S MOR_L4E MOR_L5S MOR_L5E MOR_L6S MOR_L6E MOR_L7S MOR_L7E MOR_L8S MOR_L8E MOR_L9S MOR_L9E")
		rmMOR = J(Obs, 2, "")
		rmMOR[., 2] = strofreal(1::Obs)
		cmMOR = J(19, 2, "")
		cmMOR[1,1] = "MOR_DN"
		cmMOR[2,1] = "MOR_L1S"
		cmMOR[3,1] = "MOR_L1E"
		cmMOR[4,1] = "MOR_L2S"
		cmMOR[5,1] = "MOR_L2E"
		cmMOR[6,1] = "MOR_L3S"
		cmMOR[7,1] = "MOR_L3E"
		cmMOR[8,1] = "MOR_L4S"
		cmMOR[9,1] = "MOR_L4E"
		cmMOR[10,1] = "MOR_L5S"
		cmMOR[11,1] = "MOR_L5E"
		cmMOR[12,1] = "MOR_L6S"
		cmMOR[13,1] = "MOR_L6E"
		cmMOR[14,1] = "MOR_L7S"
		cmMOR[15,1] = "MOR_L7E"
		cmMOR[16,1] = "MOR_L8S"
		cmMOR[17,1] = "MOR_L8E"
		cmMOR[18,1] = "MOR_L9S"
		cmMOR[19,1] = "MOR_L9E"
	}
}			
	*mata: _matrix_list(mMOR, rmMOR, cmMOR)

	*mOC - Outcome
{	
	mata {
		mOC = J(Obs,2,.)				
		rmOC = J(Obs, 2, "")
		rmOC[., 2] = strofreal(1::Obs)
		cmOC = J(2, 2, "")
		cmOC[1,1] = "Time"
		cmOC[2,1] = "Mort"
	}
}		
	*mata: _matrix_list(mOC, rmOC, cmOC)	
		
	*mTXR - Treatment Regimen 
{		
	mata {
		mTXR = st_data(., "CR_L1 CR_L2 CR_L3 CR_L4 CR_L5 CR_L6 CR_L7 CR_L8 CR_L9")		
		rmTXR = J(Obs, 2, "")
		rmTXR[., 2] = strofreal(1::Obs)
		cmTXR = J(9, 2, "")
		cmTXR[1,1] = "TXR_L1"
		cmTXR[2,1] = "TXR_L2"
		cmTXR[3,1] = "TXR_L3"
		cmTXR[4,1] = "TXR_L4"
		cmTXR[5,1] = "TXR_L5"
		cmTXR[6,1] = "TXR_L6"
		cmTXR[7,1] = "TXR_L7"
		cmTXR[8,1] = "TXR_L8"
		cmTXR[9,1] = "TXR_L9"
	}
}		
	*mata: _matrix_list(mTXR, rmTXR, cmTXR)
		
	*mTXD - Treatment Duration
{		
	mata {
		mTXD = st_data(., "CD_L1 CD_L2 CD_L3 CD_L4 CD_L5 CD_L6 CD_L7 CD_L8 CD_L9")		
		rmTXD = J(Obs, 2, "")
		rmTXD[., 2] = strofreal(1::Obs)
		cmTXD = J(9, 2, "")
		cmTXD[1,1] = "TXD_L1"
		cmTXD[2,1] = "TXD_L2"
		cmTXD[3,1] = "TXD_L3"
		cmTXD[4,1] = "TXD_L4"
		cmTXD[5,1] = "TXD_L5"
		cmTXD[6,1] = "TXD_L6"
		cmTXD[7,1] = "TXD_L7"
		cmTXD[8,1] = "TXD_L8"
		cmTXD[9,1] = "TXD_L9"
	}
}		
	*mata: _matrix_list(mTXD, rmTXD, cmTXD)	
		
	*mBCR - Best Clinical Response 
{	
	mata {
		mBCR = st_data(., "BCR_L1 BCR_L2 BCR_L3 BCR_L4 BCR_L5 BCR_L6 BCR_L7 BCR_L8 BCR_L9 BCR_SCT")		
		rmBCR = J(Obs, 2, "")
		rmBCR[., 2] = strofreal(1::Obs)
		cmBCR = J(10, 2, "")
		cmBCR[1,1] = "BCR_L1"
		cmBCR[2,1] = "BCR_L2"
		cmBCR[3,1] = "BCR_L3"
		cmBCR[4,1] = "BCR_L4"
		cmBCR[5,1] = "BCR_L5"
		cmBCR[6,1] = "BCR_L6"
		cmBCR[7,1] = "BCR_L7"
		cmBCR[8,1] = "BCR_L8"
		cmBCR[9,1] = "BCR_L9"
		cmBCR[10,1] = "BCR_SCT"
	}
}			
	*mata: _matrix_list(mBCR, rmBCR, cmBCR)
		
	*mTFI - Treatment-free Interval 
{		
	mata {
		mTFI = st_data(., "CI_L1 CI_L2 CI_L3 CI_L4 CI_L5 CI_L6 CI_L7 CI_L8 CI_L9")		
		rmTFI = J(Obs, 2, "")
		rmTFI[., 2] = strofreal(1::Obs)
		cmTFI = J(9, 2, "")
		cmTFI[1,1] = "TFI_DN"
		cmTFI[2,1] = "TFI_L1"
		cmTFI[3,1] = "TFI_L2"
		cmTFI[4,1] = "TFI_L3"
		cmTFI[5,1] = "TFI_L4"
		cmTFI[6,1] = "TFI_L5"
		cmTFI[7,1] = "TFI_L6"
		cmTFI[8,1] = "TFI_L7"
		cmTFI[9,1] = "TFI_L8"
	}
}		
	*mata: _matrix_list(mTFI, rmTFI, cmTFI)

	di as text _n "Setting up patient characteristic vectors..."
	
	// =========================================================================
	// CREATE DUMMY VARIABLES IN STATA
	// =========================================================================
	// Generate categorical dummy variables needed for vector creation
	
	capture gen ECOGcc0 = (ECOGcc == 0)
	capture gen ECOGcc1 = (ECOGcc == 1)
	capture gen ECOGcc2 = (ECOGcc == 2)
	capture gen RISS1 = (RISS == 1)
	capture gen RISS2 = (RISS == 2)
	capture gen RISS3 = (RISS == 3)
	capture gen CMc0 = (CMc == 0)
	capture gen CMc1 = (CMc == 1)
	capture gen CMc2 = (CMc == 2)
	capture gen CMc3 = (CMc == 3)
	capture gen Cons = 1
	
	di as text "  - Dummy variables created"
	
	mata {
		// =====================================================================
		// SYSTEM SETUP
		// =====================================================================
		// Note: Variable types are inferred by Mata, no declarations needed
		
		nObs = st_nobs()
		printf("  - Patient count: %f\n", nObs)
		
		
		// =====================================================================
		// BASELINE PATIENT CHARACTERISTICS
		// =====================================================================
		// Source: Stata variables in memory (loaded from patient dataset)
		// These vectors are created ONCE and reused across all outcome calculations
		
		displayas("text")
		printf("  - Loading demographic vectors...\n")
		
		// --- Demographics ---
		vAge = st_data(., "Age")                  // Age (continuous)
		vAge2 = vAge :^ 2                         // Age squared (for quadratic effects)
		vMale = st_data(., "Male")                // Sex (0=Female, 1=Male)
		
		
		// --- ECOG Performance Status ---
		displayas("text")
		printf("  - Loading ECOG vectors...\n")
		
		vECOG = st_data(., "ECOGcc")              // Original ECOG (0, 1, 2+)
		vECOG0 = st_data(., "ECOGcc0")            // ECOG = 0 (fully active)
		vECOG1 = st_data(., "ECOGcc1")            // ECOG = 1 (restricted strenuous activity)
		vECOG2 = st_data(., "ECOGcc2")            // ECOG = 2+ (ambulatory, self-care)
		// Note: ECOG0 typically used as reference category in regressions
		
		
		// --- Revised International Staging System (RISS) ---
		displayas("text")
		printf("  - Loading RISS vectors...\n")
		
		vRISS = st_data(., "RISS")                // Original RISS (1, 2, 3)
		vRISS1 = st_data(., "RISS1")              // RISS = I (best prognosis)
		vRISS2 = st_data(., "RISS2")              // RISS = II (intermediate)
		vRISS3 = st_data(., "RISS3")              // RISS = III (poorest prognosis)
		// Note: RISS1 typically used as reference category in regressions
		
		
		// --- Treatment Characteristics (Set Once) ---
		displayas("text")
		printf("  - Loading treatment characteristic vectors...\n")
		
		vSCT_DN = st_data(., "SCT")               // Intent for ASCT at DN (0/1)
		vSCT_L1 = st_data(., "SCT")               // Receipt of ASCT at L1 (0/1)
		vMNT = st_data(., "MNT")                  // Receipt of Maintenance therapy at L1 (0/1)
		
		// --- Utility Vector ---
		vCons = st_data(., "Cons")                // Constant vector (all ones)
		
		
		// =====================================================================
		// COMORBIDITY VECTORS
		// =====================================================================
		// Used for Intent / Receipt of ASCT 
		
		displayas("text")
		printf("  - Loading comorbidity vectors...\n")
		
		vAge70 = st_data(., "Age70")              // Age >= 70 indicator
		vAge75 = st_data(., "Age75")              // Age >= 75 indicator
		vCMc = st_data(., "CMc")                    // Comorbidity score (0, 1, 2, 3)
		vCMc0 = st_data(., "CMc0")                  // CMc = 0 (no comorbidities)
		vCMc1 = st_data(., "CMc1")                  // CMc = 1 (1 comorbidity)
		vCMc2 = st_data(., "CMc2")                  // CMc = 2 (2 comorbidities)
		vCMc3 = st_data(., "CMc3")                  // CMc = 3 (3+ comorbidities)
		// Note: CM0 typically used as reference category in regressions
		
		
		// =====================================================================
		// BEST CLINICAL RESPONSE (BCR) CATEGORY DUMMIES
		// =====================================================================
		// Pre-create BCR dummy vectors for convenient access
		// Note: These will be extracted from mBCR matrix columns as needed
		//       This section documents the BCR coding scheme
		
		// BCR Coding:
		//   1 = CR   (Complete Response)
		//   2 = VGPR (Very Good Partial Response)
		//   3 = PR   (Partial Response)
		//   4 = SD   (Stable Disease)
		//   5 = MR   (Minimal Response)
		//   6 = PD   (Progressive Disease)
		
		// These vectors can be created from mBCR columns on-the-fly:
		// Example for Line 2 BCR:
		//   vBCR2_CR   = (mBCR[., 2] :== 1)
		//   vBCR2_VGPR = (mBCR[., 2] :== 3)
		//   vBCR2_MR   = (mBCR[., 2] :== 5)
		
		
		// =====================================================================
		// VECTOR VALIDATION
		// =====================================================================
		// Ensure all vectors have correct dimensions and valid data
		
		displayas("text")
		printf("\n  - Validating vectors...\n")
		
		errorCount = 0
		
		// Check dimensions
		if (rows(vAge) != nObs) {
			displayas("error")
			printf("ERROR: vAge has incorrect dimensions (%f rows, expected %f)\n", rows(vAge), nObs)
			errorCount++
		}
		if (rows(vMale) != nObs) {
			displayas("error")
			printf("ERROR: vMale has incorrect dimensions (%f rows, expected %f)\n", rows(vMale), nObs)
			errorCount++
		}
		if (rows(vECOG) != nObs) {
			displayas("error")
			printf("ERROR: vECOG has incorrect dimensions (%f rows, expected %f)\n", rows(vECOG), nObs)
			errorCount++
		}
		if (rows(vRISS) != nObs) {
			displayas("error")
			printf("ERROR: vRISS has incorrect dimensions (%f rows, expected %f)\n", rows(vRISS), nObs)
			errorCount++
		}
		if (rows(vCons) != nObs) {
			displayas("error")
			printf("ERROR: vCons has incorrect dimensions (%f rows, expected %f)\n", rows(vCons), nObs)
			errorCount++
		}
		
		// Check for missing values in critical vectors
		if (hasmissing(vAge)) {
			displayas("error")
			printf("WARNING: vAge contains missing values\n")
		}
		if (hasmissing(vECOG)) {
			displayas("error")
			printf("WARNING: vECOG contains missing values\n")
		}
		if (hasmissing(vRISS)) {
			displayas("error")
			printf("WARNING: vRISS contains missing values\n")
		}
		
		// Check valid ranges
		if (min(vAge) < 18 | max(vAge) > 120) {
			displayas("error")
			printf("WARNING: vAge contains values outside expected range [18, 120]\n")
			printf("  Range: %f to %f\n", min(vAge), max(vAge))
		}
		if (min(vECOG) < 0 | max(vECOG) > 2) {
			displayas("error")
			printf("WARNING: vECOG contains values outside expected range [0, 2]\n")
			printf("  Range: %f to %f\n", min(vECOG), max(vECOG))
		}
		if (min(vRISS) < 1 | max(vRISS) > 3) {
			displayas("error")
			printf("WARNING: vRISS contains values outside expected range [1, 3]\n")
			printf("  Range: %f to %f\n", min(vRISS), max(vRISS))
		}
		
		// Validate dummy variables sum to 1
		ecogSum = vECOG0 + vECOG1 + vECOG2
		rissSum = vRISS1 + vRISS2 + vRISS3
		cmSum = vCMc0 + vCMc1 + vCMc2 + vCMc3
		
		if (min(ecogSum) != 1 | max(ecogSum) != 1) {
			displayas("error")
			printf("ERROR: ECOG dummy variables do not sum to 1 for all patients\n")
			errorCount++
		}
		if (min(rissSum) != 1 | max(rissSum) != 1) {
			displayas("error")
			printf("ERROR: RISS dummy variables do not sum to 1 for all patients\n")
			errorCount++
		}
		if (min(cmSum) != 1 | max(cmSum) != 1) {
			displayas("error")
			printf("ERROR: CM dummy variables do not sum to 1 for all patients\n")
			errorCount++
		}
		
		// Report validation results
		if (errorCount > 0) {
			displayas("error")
			printf("\nVector validation FAILED with %f errors\n", errorCount)
			exit(198)
		}
		else {
			displayas("result")
			printf("  - Vector validation: PASSED\n")
		}
		
		
		// =====================================================================
		// SUMMARY STATISTICS
		// =====================================================================
		
		displayas("text")
		printf("\nVector setup summary:\n")
		printf("  ----------------------------------------\n")
		
		// Calculate statistics first
		ageMean = mean(vAge)
		ageSD = sqrt(variance(vAge))
		ageMin = min(vAge)
		ageMax = max(vAge)
		malePct = 100*mean(vMale)
		ecog0Pct = 100*mean(vECOG0)
		ecog1Pct = 100*mean(vECOG1)
		ecog2Pct = 100*mean(vECOG2)
		riss1Pct = 100*mean(vRISS1)
		riss2Pct = 100*mean(vRISS2)
		riss3Pct = 100*mean(vRISS3)
		cm0Pct = 100*mean(vCMc0)
		cm1Pct = 100*mean(vCMc1)
		cm2Pct = 100*mean(vCMc2)
		cm3Pct = 100*mean(vCMc3)
		
		// Display with simple formatting
		printf("  Age:           Mean=%g, SD=%g, Range=[%g, %g]\n", 
		       ageMean, ageSD, ageMin, ageMax)
		printf("  Male:          %g%%\n", malePct)
		printf("  ECOG 0/1/2:    %g%% / %g%% / %g%%\n", 
		       ecog0Pct, ecog1Pct, ecog2Pct)
		printf("  RISS I/II/III: %g%% / %g%% / %g%%\n", 
		       riss1Pct, riss2Pct, riss3Pct)
		printf("  CM 0/1/2/3:    %g%% / %g%% / %g%% / %g%%\n", 
		       cm0Pct, cm1Pct, cm2Pct, cm3Pct)
		printf("  ----------------------------------------\n")
		
		displayas("text")
		printf("\nVector setup complete")
		displayas("result")
		printf(" ✓\n")
		
	} 
	
end

