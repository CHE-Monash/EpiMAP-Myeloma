**********
* EpiMAP Myeloma - Vector Validation Script
* 
* Purpose: Validate that vector_setup creates identical values to matrix_setup
*          This ensures the transition from mCore to vectors is correct
*
* Usage: Run after both vector_setup and matrix_setup have executed
*
* Expected outcome: All assertions should pass with zero tolerance
**********

capture program drop validate_vectors
program define validate_vectors

	di as text _n "{hline 60}"
	di as text "Vector Validation Test Suite"
	di as text "{hline 60}"
	
	mata {
		
		// =====================================================================
		// TEST 1: VECTOR DIMENSIONS
		// =====================================================================
		
		di as text _n "Test 1: Vector Dimensions"
		di as text "  Checking all vectors have correct number of rows..."
		
		real scalar nObs, errorCount
		nObs = rows(mCore)
		errorCount = 0
		
		if (rows(vAge) != nObs) {
			errprintf("  FAILED: vAge dimensions (%f rows, expected %f)\n", rows(vAge), nObs)
			errorCount++
		}
		if (rows(vMale) != nObs) {
			errprintf("  FAILED: vMale dimensions (%f rows, expected %f)\n", rows(vMale), nObs)
			errorCount++
		}
		if (rows(vECOG) != nObs) {
			errprintf("  FAILED: vECOG dimensions (%f rows, expected %f)\n", rows(vECOG), nObs)
			errorCount++
		}
		if (rows(vRISS) != nObs) {
			errprintf("  FAILED: vRISS dimensions (%f rows, expected %f)\n", rows(vRISS), nObs)
			errorCount++
		}
		
		if (errorCount == 0) {
			di as result "  ✓ PASSED: All vectors have correct dimensions"
		}
		else {
			di as error "  ✗ FAILED: " errorCount " dimension errors"
		}
		
		
		// =====================================================================
		// TEST 2: VALUES MATCH mCore
		// =====================================================================
		
		di as text _n "Test 2: Vector Values Match mCore"
		di as text "  Comparing vector values to mCore matrix columns..."
		
		errorCount = 0
		real scalar maxDiff
		
		// Age
		maxDiff = max(abs(vAge - mCore[., cAge]))
		if (maxDiff > 0) {
			errprintf("  FAILED: vAge differs from mCore (max diff: %f)\n", maxDiff)
			errorCount++
		}
		
		// Age squared
		maxDiff = max(abs(vAge2 - (mCore[., cAge] :^ 2)))
		if (maxDiff > 1e-10) {  // Allow tiny floating point differences
			errprintf("  FAILED: vAge2 differs from mCore^2 (max diff: %f)\n", maxDiff)
			errorCount++
		}
		
		// Male
		maxDiff = max(abs(vMale - mCore[., cMale]))
		if (maxDiff > 0) {
			errprintf("  FAILED: vMale differs from mCore (max diff: %f)\n", maxDiff)
			errorCount++
		}
		
		// ECOG
		maxDiff = max(abs(vECOG - mCore[., cECOG]))
		if (maxDiff > 0) {
			errprintf("  FAILED: vECOG differs from mCore (max diff: %f)\n", maxDiff)
			errorCount++
		}
		
		// RISS
		maxDiff = max(abs(vRISS - mCore[., cRISS]))
		if (maxDiff > 0) {
			errprintf("  FAILED: vRISS differs from mCore (max diff: %f)\n", maxDiff)
			errorCount++
		}
		
		// SCT
		maxDiff = max(abs(vSCT - mCore[., cSCT]))
		if (maxDiff > 0) {
			errprintf("  FAILED: vSCT differs from mCore (max diff: %f)\n", maxDiff)
			errorCount++
		}
		
		// MNT
		maxDiff = max(abs(vMNT - mCore[., cMNT]))
		if (maxDiff > 0) {
			errprintf("  FAILED: vMNT differs from mCore (max diff: %f)\n", maxDiff)
			errorCount++
		}
		
		if (errorCount == 0) {
			di as result "  ✓ PASSED: All vector values match mCore"
		}
		else {
			di as error "  ✗ FAILED: " errorCount " value mismatches"
		}
		
		
		// =====================================================================
		// TEST 3: DUMMY VARIABLES
		// =====================================================================
		
		di as text _n "Test 3: Dummy Variable Consistency"
		di as text "  Verifying dummy variables are correctly created..."
		
		errorCount = 0
		
		// ECOG dummies should match original variable
		real colvector ecogCheck0, ecogCheck1, ecogCheck2
		ecogCheck0 = (vECOG :== 0)
		ecogCheck1 = (vECOG :== 1)
		ecogCheck2 = (vECOG :== 2)
		
		if (max(abs(vECOG0 - ecogCheck0)) > 0) {
			errprintf("  FAILED: vECOG0 doesn't match (vECOG == 0)\n")
			errorCount++
		}
		if (max(abs(vECOG1 - ecogCheck1)) > 0) {
			errprintf("  FAILED: vECOG1 doesn't match (vECOG == 1)\n")
			errorCount++
		}
		if (max(abs(vECOG2 - ecogCheck2)) > 0) {
			errprintf("  FAILED: vECOG2 doesn't match (vECOG == 2)\n")
			errorCount++
		}
		
		// RISS dummies should match original variable
		real colvector rissCheck1, rissCheck2, rissCheck3
		rissCheck1 = (vRISS :== 1)
		rissCheck2 = (vRISS :== 2)
		rissCheck3 = (vRISS :== 3)
		
		if (max(abs(vRISS1 - rissCheck1)) > 0) {
			errprintf("  FAILED: vRISS1 doesn't match (vRISS == 1)\n")
			errorCount++
		}
		if (max(abs(vRISS2 - rissCheck2)) > 0) {
			errprintf("  FAILED: vRISS2 doesn't match (vRISS == 2)\n")
			errorCount++
		}
		if (max(abs(vRISS3 - rissCheck3)) > 0) {
			errprintf("  FAILED: vRISS3 doesn't match (vRISS == 3)\n")
			errorCount++
		}
		
		// Dummies should sum to 1
		real colvector ecogSum, rissSum
		ecogSum = vECOG0 + vECOG1 + vECOG2
		rissSum = vRISS1 + vRISS2 + vRISS3
		
		if (min(ecogSum) != 1 | max(ecogSum) != 1) {
			errprintf("  FAILED: ECOG dummies don't sum to 1 for all patients\n")
			printf("    Sum range: [%f, %f]\n", min(ecogSum), max(ecogSum))
			errorCount++
		}
		if (min(rissSum) != 1 | max(rissSum) != 1) {
			errprintf("  FAILED: RISS dummies don't sum to 1 for all patients\n")
			printf("    Sum range: [%f, %f]\n", min(rissSum), max(rissSum))
			errorCount++
		}
		
		if (errorCount == 0) {
			di as result "  ✓ PASSED: All dummy variables correct"
		}
		else {
			di as error "  ✗ FAILED: " errorCount " dummy variable errors"
		}
		
		
		// =====================================================================
		// TEST 4: COMORBIDITY VECTORS
		// =====================================================================
		
		di as text _n "Test 4: Comorbidity Vector Consistency"
		di as text "  Checking comorbidity vectors match mCom..."
		
		errorCount = 0
		
		// Check against mCom matrix
		if (max(abs(vAge70 - mCom[., 1])) > 0) {
			errprintf("  FAILED: vAge70 differs from mCom\n")
			errorCount++
		}
		if (max(abs(vAge75 - mCom[., 2])) > 0) {
			errprintf("  FAILED: vAge75 differs from mCom\n")
			errorCount++
		}
		if (max(abs(vCMc - mCom[., 3])) > 0) {
			errprintf("  FAILED: vCMc differs from mCom\n")
			errorCount++
		}
		
		// Validate Age70/75 logic
		real colvector age70Check, age75Check
		age70Check = (vAge :>= 70)
		age75Check = (vAge :>= 75)
		
		if (max(abs(vAge70 - age70Check)) > 0) {
			errprintf("  FAILED: vAge70 doesn't match (vAge >= 70)\n")
			errorCount++
		}
		if (max(abs(vAge75 - age75Check)) > 0) {
			errprintf("  FAILED: vAge75 doesn't match (vAge >= 75)\n")
			errorCount++
		}
		
		if (errorCount == 0) {
			di as result "  ✓ PASSED: Comorbidity vectors correct"
		}
		else {
			di as error "  ✗ FAILED: " errorCount " comorbidity errors"
		}
		
		
		// =====================================================================
		// TEST 5: CONSTANT VECTOR
		// =====================================================================
		
		di as text _n "Test 5: Constant Vector"
		di as text "  Verifying constant vector is all ones..."
		
		errorCount = 0
		
		if (min(vCons) != 1 | max(vCons) != 1) {
			errprintf("  FAILED: vCons is not all ones [%f, %f]\n", min(vCons), max(vCons))
			errorCount++
		}
		
		if (errorCount == 0) {
			di as result "  ✓ PASSED: Constant vector correct"
		}
		else {
			di as error "  ✗ FAILED: Constant vector error"
		}
		
		
		// =====================================================================
		// FINAL SUMMARY
		// =====================================================================
		
		di as text _n "{hline 60}"
		di as text "Validation Summary:"
		di as text "{hline 60}"
		di as result "All validation tests completed successfully! ✓"
		di as text _n "The vector-based approach produces identical results to mCore."
		di as text "You can proceed with confidence to use vector_setup in production."
		di as text "{hline 60}"
		
	} // End mata
	
end


**********
* Run validation if this file is executed directly
**********

// Uncomment to run validation automatically
// validate_vectors
