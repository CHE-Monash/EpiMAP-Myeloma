**********
* EpiMAP Myeloma - Vector Setup
* Create reusable patient characteristic vectors
* 
* Purpose: Replaces mCore matrix approach with efficient vector-based storage
*          Patient characteristics stored ONCE and reused across all outcomes
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

capture program drop vector_setup
program define vector_setup
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
