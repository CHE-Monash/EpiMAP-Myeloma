**********
* Test Mata Functions
**********

clear all
quietly do "core/mata_functions.do"

mata:
    // Test 1: calcSurvivalTime() matches old formula
    printf("\n{hline 80}\n")
    printf("TEST 1: Weibull survival calculation\n")
    printf("{hline 80}\n")
    
    XB = rnormal(100, 1, 0, 1)
    RN = runiform(100, 1)
    auxParam = 0.5
    
    // OLD formula
    old = ((ln(RN) :/ -exp(XB)) :^ (1 :/ exp(auxParam)))
    
    // NEW function
    new = calcSurvivalTime(XB, RN, "weibull", auxParam)
    
    // Compare
    maxDiff = max(abs(old - new))
    printf("Maximum difference: %f\n", maxDiff)
    if (maxDiff < 1e-10) {
        printf("{result}PASS: Results match perfectly\n")
    }
    else {
        printf("{error}FAIL: Results differ\n")
    }
    
    
    // Test 2: Ordered logit
    printf("\n{hline 80}\n")
    printf("TEST 2: Ordered logit assignment\n")
    printf("{hline 80}\n")
    
    XB = rnormal(100, 1, 0, 1)
    RN = runiform(100, 1)
    cutPoints = (0.5, 1.5)
    
    // Calculate using functions
    probs = calcOrderedLogitProbs(XB, cutPoints)
    outcome = assignOrderedOutcome(RN, probs, (1, 3, 5))
    
    // Validate all outcomes are valid
    validateOutcomes(outcome, (1, 3, 5), "BCR Test")
    
    // Check distribution
    printf("Outcome distribution:\n")
    printf("  CR (1):   %f%%\n", 100 * sum(outcome :== 1) / rows(outcome))
    printf("  VGPR (3): %f%%\n", 100 * sum(outcome :== 3) / rows(outcome))
    printf("  PR (5):   %f%%\n", 100 * sum(outcome :== 5) / rows(outcome))
    printf("{result}PASS: All outcomes valid\n")
    
    
    // Test 3: Curtailing
    printf("\n{hline 80}\n")
    printf("TEST 3: Value curtailing\n")
    printf("{hline 80}\n")
    
    values = (100 \ 200 \ 300 \ 400 \ 500)
    maxVal = 350
    
    curtailed = curtailValues(values, maxVal)
    
    printf("Before: "), values'
    printf("After (max=350): "), curtailed'
    
    if (max(curtailed) <= maxVal) {
        printf("{result}PASS: All values <= maximum\n")
    }
    
    printf("\n{hline 80}\n")
    printf("All tests complete!\n")
    printf("{hline 80}\n")
end
