**********
* Test Survival Time Calculation Functions
* Tests calcSurvivalTime() against old formulas
**********

clear all

// Clear Mata memory
mata: mata clear

// Load the utility functions
quietly do "core/mata_functions.do"

//=============================================================================
// Setup test data
//=============================================================================
mata:
nObs = 1000
XB = rnormal(nObs, 1, 0, 1)
RN = runiform(nObs, 1)
aux = 0.5

printf("Generated %f observations\n", nObs)
printf("XB ~ N(0,1), RN ~ U(0,1), aux = %f\n", aux)

end

//=============================================================================
// TEST 1: Exponential Distribution
//=============================================================================
mata:
// OLD formula (from your codebase)
old_exp = ln(RN) :/ -exp(XB)

// NEW function
new_exp = calcSurvivalTime(XB, RN, "exponential")

// Compare
maxDiff_exp = max(abs(old_exp - new_exp))
meanDiff_exp = mean(abs(old_exp - new_exp))

printf("Maximum difference: %e\n", maxDiff_exp)
printf("Mean difference:    %e\n", meanDiff_exp)

if (maxDiff_exp < 1e-10) {
	printf("{result}✓ PASS: Exponential results match perfectly\n")
}
else {
	printf("{error}✗ FAIL: Exponential results differ\n")
	exit(9)
}
end

//=============================================================================
// TEST 2: Weibull Distribution
//=============================================================================
mata:
// OLD formula (from your codebase)
old_weib = ((ln(RN) :/ -exp(XB)) :^ (1 :/ exp(aux)))

// NEW function
new_weib = calcSurvivalTime(XB, RN, "weibull", aux)

// Compare
maxDiff_weib = max(abs(old_weib - new_weib))
meanDiff_weib = mean(abs(old_weib - new_weib))

printf("Maximum difference: %e\n", maxDiff_weib)
printf("Mean difference:    %e\n", meanDiff_weib)

if (maxDiff_weib < 1e-10) {
    printf("{result}✓ PASS: Weibull results match perfectly\n")
}
else {
    printf("{error}✗ FAIL: Weibull results differ\n")
    exit(9)
}
end

//=============================================================================
// TEST 3: Gompertz Distribution
//=============================================================================
mata:
// OLD formula (from your codebase)
old_gomp = (ln(1 :- ((aux :* ln(RN)) :/ exp(XB))) :/ aux)

// NEW function
new_gomp = calcSurvivalTime(XB, RN, "gompertz", aux)

// Compare
maxDiff_gomp = max(abs(old_gomp - new_gomp))
meanDiff_gomp = mean(abs(old_gomp - new_gomp))

printf("Maximum difference: %e\n", maxDiff_gomp)
printf("Mean difference:    %e\n", meanDiff_gomp)

if (maxDiff_gomp < 1e-10) {
    printf("{result}✓ PASS: Gompertz results match perfectly\n")
}
else {
    printf("{error}✗ FAIL: Gompertz results differ\n")
    exit(9)
}
end

//=============================================================================
// TEST 4: Test "ereg" alias for exponential
//=============================================================================
mata:
new_ereg = calcSurvivalTime(XB, RN, "ereg", aux)
maxDiff_ereg = max(abs(new_exp - new_ereg))

printf("Maximum difference: %e\n", maxDiff_ereg)

if (maxDiff_ereg < 1e-10) {
    printf("{result}✓ PASS: 'ereg' equals 'exponential'\n")
}
else {
    printf("{error}✗ FAIL: 'ereg' differs from 'exponential'\n")
    exit(9)
}
end

//=============================================================================
// TEST 5: Error handling for invalid distribution
//=============================================================================
mata:
bad_result = calcSurvivalTime(XB, RN, "invalid_dist", aux)
end

//=============================================================================
// TEST 6: Vectorisation check (different sized vectors)
//=============================================================================
mata:
// Test with different sized vectors
sizes = (10, 100, 1000, 5000)

for (i = 1; i <= cols(sizes); i++) {
    n = sizes[i]
    XB_test = rnormal(n, 1, 0, 1)
    RN_test = runiform(n, 1)
    
    result = calcSurvivalTime(XB_test, RN_test, "weibull", 0.5)
    
    if (rows(result) == n) {
        printf("  n=%5.0f: ✓ Correct output size\n", n)
    }
    else {
        printf("{error}  n=%5.0f: ✗ Wrong output size\n", n)
        exit(9)
    }
}

printf("{result}✓ PASS: Function works with different vector sizes\n")
end

//=============================================================================
// TEST 7: Descriptive statistics check
//=============================================================================
mata:
// Generate new test data with known parameters
nTest = 10000
XB_test = J(nTest, 1, 0)  // XB = 0 means lambda = 1
RN_test = runiform(nTest, 1)
aux_test = 1.0  // Shape = e^1 ≈ 2.718

result_test = calcSurvivalTime(XB_test, RN_test, "weibull", aux_test)

printf("\nWeibull with XB=0 (lambda=1), shape=e^1:\n")
printf("  Mean:   %8.4f\n", mean(result_test))
printf("  Min:    %8.4f\n", min(result_test))
printf("  Max:    %8.4f\n", max(result_test))
printf("  SD:     %8.4f\n", sqrt(variance(result_test)))

// All values should be positive
if (min(result_test) > 0) {
    printf("{result}✓ PASS: All survival times are positive\n")
}
else {
    printf("{error}✗ FAIL: Found non-positive survival times\n")
    exit(9)
}
end
