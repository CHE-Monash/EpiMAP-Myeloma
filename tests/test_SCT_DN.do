**********
* SCT DN Comparison Test
* 
* Purpose: Validate vectorised SCT_DN produces identical results to original
* Method: Run both implementations with same random seed and compare
* Success: All outcomes match exactly (tolerance = 0 for binary outcome)
**********

clear all
set more off

di as text _n "{hline 60}"
di as text "SCT DN Vectorisation Validation Test"
di as text "{hline 60}"

// Set fixed seed for reproducibility
set seed 12345
di as text _n "Random seed set to: 12345"

// Load test data and setup
di as text _n "Loading patients and setting up..."
quietly do "core/load_patients.do"
quietly do "core/mata_functions.do"
load_patients

// Load coefficients
mata: mata matuse "analyses/base_model/data/coefficients/coefficients_base_model.mmat"

di as text "Setting up vectors and matrices..."
quietly do "core/vector_setup.do"
vector_setup

quietly do "core/matrix_setup.do"
matrix_setup

di as text _n "{hline 60}"
di as text "Running ORIGINAL (loop-based) implementation..."
di as text "{hline 60}"

// Store original mSCT and mCore before running
mata: mSCT_original = mSCT
mata: mCore_original = mCore

// Run the original SCT DN code
scalar m = "mDN_SCT"
scalar b = "bDN_SCT"
scalar c = "cLO_"
scalar OMC = 1

quietly do "core/outcomes/SIM SCT DN.do"

// Store results
mata: SCT_original = mSCT[., 1]

di as result "  ✓ Original implementation completed"

// Reset matrices for vectorised version
mata: mSCT = mSCT_original
mata: mCore = mCore_original

di as text _n "{hline 60}"
di as text "Running VECTORISED implementation..."
di as text "{hline 60}"

// Run vectorised implementation with SAME random seed
set seed 12345

quietly do "core/outcomes/SIM SCT DN Vector.do"

// Store results
mata: SCT_vectorised = mSCT[., 1]

di as result "  ✓ Vectorised implementation completed"

di as text _n "{hline 60}"
di as text "Comparing Results..."
di as text "{hline 60}"

mata {
	// Compare outcomes
	nObs = rows(SCT_original)
	nMatch = sum(SCT_original :== SCT_vectorised)
	nDiffer = sum(SCT_original :!= SCT_vectorised)
	nMissing = sum(missing(SCT_original) :& missing(SCT_vectorised))
	
	displayas("text")
	printf("\nComparison Results:\n")
	printf("  Total patients:     %g\n", nObs)
	printf("  Exact matches:      %g (%g%%)\n", nMatch, 100*nMatch/nObs)
	printf("  Differences:        %g (%g%%)\n", nDiffer, 100*nDiffer/nObs)
	printf("  Both missing:       %g\n", nMissing)
	
	// Check for differences
	if (nDiffer > 0) {
		displayas("error")
		printf("\n✗ VALIDATION FAILED: Outcomes differ!\n")
		
		// Show first 10 differences
		displayas("text")
		printf("\nFirst 10 differences:\n")
		printf("%10s %12s %12s\n", "Patient", "Original", "Vectorised")
		printf("  %s\n", "{hline 36}")
		
		diffIdx = selectindex(SCT_original :!= SCT_vectorised)
		nShow = min((rows(diffIdx), 10))
		
		for (i=1; i<=nShow; i++) {
			idx = diffIdx[i]
			printf("%10.0f %12.0f %12.0f\n", 
			       idx, SCT_original[idx], SCT_vectorised[idx])
		}
		
		exit(9)
	}
	else {
		displayas("result")
		printf("\n✓ VALIDATION PASSED: All outcomes match exactly!\n")
	}
	
	// Summary statistics comparison
	displayas("text")
	printf("\nOutcome Distribution:\n")
	printf("  %20s %12s %12s\n", "", "Original", "Vectorised")
	printf("  %s\n", "{hline 48}")
	printf("  %20s %12.0f %12.0f\n", "Eligible (1)", 
	       sum(SCT_original :== 1), sum(SCT_vectorised :== 1))
	printf("  %20s %12.0f %12.0f\n", "Not eligible (0)", 
	       sum(SCT_original :== 0), sum(SCT_vectorised :== 0))
	printf("  %20s %12.0f %12.0f\n", "Missing", 
	       sum(missing(SCT_original)), sum(missing(SCT_vectorised)))
	printf("  %s\n", "{hline 48}")
	printf("  %20s %11.2f%% %11.2f%%\n", "Eligibility rate", 
	       100*mean(SCT_original :== 1), 100*mean(SCT_vectorised :== 1))
	
} // End mata

di as text _n "{hline 60}"
di as result "Vectorisation validation complete!"
di as text "{hline 60}"
di as text _n "Summary:"
di as text "  - Results match:    " as result "YES ✓"
di as text "  - Ready for production use"
di as text "{hline 60}"
