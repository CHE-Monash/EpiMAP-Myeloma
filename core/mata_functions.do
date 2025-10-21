**********
* EpiMAP Myeloma - Mata Utility Functions
* 
* Purpose: Reusable functions for survival calculations, ordered logit,
*          filtering, and validation
*
* Usage: Load once at start of analysis
*        quietly do "core/mata_functions.do"
*
* Author: EpiMAP Team
* Date: October 2025
**********

mata:
mata clear

//=============================================================================
// SURVIVAL DISTRIBUTION FUNCTIONS
//=============================================================================

/**
 * Calculate survival time from parametric distributions
 * 
 * @param XB       Linear predictor (vector or scalar)
 * @param RN       Random numbers in (0,1) (vector or scalar)
 * @param dist     Distribution type: "exponential", "ereg", "weibull", "gompertz"
 * @param aux      Auxiliary parameter (shape for Weibull, gamma for Gompertz)
 *                 NOTE: This is always the final column in coefficient matrix
 *                 Extract as: aux = bMatrix[1, cols(bMatrix)]
 * @return         Survival times
 */
real colvector calcSurvTime(
    real colvector XB,
    real colvector RN,
    string scalar dist,
    | real scalar aux 	// Optional parameter
) {
    real colvector outcome
    
    if (dist == "exponential" | dist == "ereg") {
        // Exponential: T = -ln(U) / lambda, where lambda = exp(XB)
        outcome = ln(RN) :/ -exp(XB)
    }
    else if (dist == "weibull") {
        // Weibull: T = (-ln(U) / lambda)^(1/p), where p = exp(aux)
        outcome = ((ln(RN) :/ -exp(XB)) :^ (1 :/ exp(aux)))
    }
    else if (dist == "gompertz") {
        // Gompertz: T = (1/gamma) * ln(1 - (gamma * ln(U)) / lambda)
        outcome = (ln(1 :- ((aux :* ln(RN)) :/ exp(XB))) :/ aux)
    }
    else {
        errprintf("ERROR: Unknown distribution '%s'\n", dist)
        errprintf("Valid options: 'exponential', 'ereg', 'weibull', 'gompertz'\n")
        exit(198)
    }
    
    return(outcome)
}


/**
 * Calculate survival probability at time TSD
 * 
 * Returns S(TSD) = P(T > TSD) - the probability of surviving to time TSD
 * Used for conditional survival calculations where patient already alive at TSD
 * 
 * @param XB       Linear predictor (vector or scalar)
 * @param TSD      Time since diagnosis - calculate S(TSD)
 * @param dist     Distribution type ("exponential"/"ereg", "weibull", "gompertz")
 * @param aux      Auxiliary parameter (shape/gamma parameter)
 * @return         Survival probabilities S(TSD)
 */
real colvector calcSurvProb(
    real colvector XB,
    real colvector TSD,
    string scalar dist,
    | real scalar aux
) {
    real colvector PR
    
    // Calculate S(t) based on distribution
    if (dist == "exponential" | dist == "ereg") {
        // Exponential: S(t) = exp(-lambda * t)
        PR = exp(-exp(XB) :* TSD)
    }
    else if (dist == "weibull") {
        // Weibull: S(t) = exp(-lambda * t^shape)
        PR = exp(-exp(XB) :* (TSD :^ exp(aux)))
    }
    else if (dist == "gompertz") {
        // Gompertz: S(t) = exp(-lambda/gamma * (exp(gamma*t) - 1))
        PR = exp(-exp(XB) :* (1 :/ aux) :* (exp(aux :* TSD) :- 1))
    }
    else {
        errprintf("ERROR: Unknown distribution '%s'\n", dist)
        errprintf("Valid distributions: 'exponential', 'ereg', 'weibull', 'gompertz'\n")
        exit(198)
    }
    
    return(PR)
}


//=============================================================================
// ORDERED LOGIT FUNCTIONS
//=============================================================================

/**
 * Calculate cumulative probabilities for ordered logit model
 * 
 * For ordered logit: P(Y <= j) = 1 / (1 + exp(XB - cutpoint_j))
 * 
 * @param XB         Linear predictor (vector)
 * @param cutPoints  Row vector of cut points
 *                   NOTE: Extract from coefficient matrix as final columns
 *                   For 3-category: cutPoints = bMatrix[1, (cols(bMatrix)-1, cols(bMatrix))]
 *                   For 6-category: cutPoints = bMatrix[1, (n+1..n+5)] where n = # predictors
 * @return           Matrix of cumulative probabilities (rows = obs, cols = cutpoints)
 */
real matrix calcOrdLogitProbs(
    real colvector XB,
    real rowvector cutPoints	
) {
    real matrix probs
    real scalar i
    
    probs = J(rows(XB), cols(cutPoints), .)
    
    for (i = 1; i <= cols(cutPoints); i++) {
        probs[., i] = 1 :/ (1 :+ exp(XB :- cutPoints[i]))
    }
    
    return(probs)
}


/**
 * Assign outcomes based on ordered logit cumulative probabilities
 * 
 * Works for any number of categories (2-6)
 * 
 * @param RN         Random numbers in (0,1)
 * @param probs      Cumulative probabilities from calcOrderedLogitProbs()
 * @param values     Row vector of outcome values for each category
 *                   Example for 3-category BCR: (1, 3, 5) for CR, VGPR, PR
 *                   Example for 6-category BCR: (1, 2, 3, 4, 5, 6)
 * @return           Vector of assigned outcomes
 */
real colvector assignOrdOutcome(
    real colvector RN,
    real matrix probs,
    real rowvector values
) {
    real colvector outcome
    real scalar i
    
    // Start with highest category as default
    outcome = J(rows(RN), 1, values[cols(values)])
    
    // Work BACKWARDS from highest to lowest category
    // This ensures the lowest matching category wins
    for (i = cols(values)-1; i >= 1; i--) {
        outcome = (RN :< probs[., i]) :* values[i] + ///
                  (RN :>= probs[., i]) :* outcome
    }
    
    return(outcome)
}


//=============================================================================
// FILTERING FUNCTIONS
//=============================================================================

/**
 * Get indices of patients eligible for outcome calculation
 * 
 * Filters based on:
 * 1. Alive at previous time point (mMOR[., currentOMC-1] == 0)
 * 2. Optional: In valid state (mState[., 1] <= maxState)
 * 
 * @param mMOR       Mortality matrix
 * @param currentOMC Current outcome matrix column (OMC)
 * @param mState     Optional: State matrix
 * @param maxState   Optional: Maximum allowed state value
 * @return           Vector of row indices for valid patients
 */
real colvector getValidPatients(
    real matrix mMOR,
    real scalar currentOMC,
    | real matrix mState,
      real scalar maxState
) {
    real colvector aliveIdx, stateIdx, validIdx
    
    // Get patients alive at previous time point
    aliveIdx = selectindex(mMOR[., currentOMC-1] :== 0)
    
    // If state filter provided, apply it
    if (args() >= 3) {
        stateIdx = selectindex(mState[., 1] :<= maxState)
        // Find intersection of alive and valid state
        validIdx = select(aliveIdx, rowsum(aliveIdx :== stateIdx') :> 0)
    }
    else {
        validIdx = aliveIdx
    }
    
    return(validIdx)
}


/**
 * Check if patients are alive (returns logical vector)
 * 
 * @param mMOR       Mortality matrix
 * @param currentOMC Current outcome matrix column
 * @return           Logical vector (1 = alive, 0 = dead)
 */
real colvector isAlive(
    real matrix mMOR,
    real scalar currentOMC
) {
    return(mMOR[., currentOMC-1] :== 0)
}


//=============================================================================
// VALUE MANIPULATION FUNCTIONS
//=============================================================================

/**
 * Truncate values at maximum and/or minimum
 * 
 * Used to curtail survival times at maximum observed values
 * 
 * @param values     Vector to truncate
 * @param maxValue   Optional: Maximum value (values > max set to max)
 * @param minValue   Optional: Minimum value (values < min set to min)
 * @return           Truncated values
 */
real colvector curtailValues(
    real colvector values,
    | real scalar maxValue,
      real scalar minValue
) {
    real colvector result
    
    result = values
    
    // Apply maximum if provided and not missing
    if (args() >= 2 & maxValue < .) {
        result = result :* (result :<= maxValue) + ///
                 maxValue :* (result :> maxValue)
    }
    
    // Apply minimum if provided and not missing
    if (args() >= 3 & minValue > .) {
        result = result :* (result :>= minValue) + ///
                 minValue :* (result :< minValue)
    }
    
    return(result)
}


/**
 * Truncate values at maximum/minimum, preserving missing values
 * 
 * Safer version that doesn't modify missing values
 * 
 * @param values     Vector to truncate
 * @param maxValue   Maximum value
 * @param minValue   Optional: Minimum value
 * @return           Truncated values with missing preserved
 */
real colvector curtailValuesSafe(
    real colvector values,
    real scalar maxValue,
    | real scalar minValue
) {
    real colvector result, nonMissing
    
    result = values
    nonMissing = (values :< .)
    
    // Apply maximum to non-missing values only
    if (maxValue < .) {
        result = result :* (nonMissing :* (result :<= maxValue) + (1 :- nonMissing)) + ///
                 maxValue :* (nonMissing :* (result :> maxValue))
    }
    
    // Apply minimum to non-missing values only
    if (args() >= 3 & minValue > .) {
        result = result :* (nonMissing :* (result :>= minValue) + (1 :- nonMissing)) + ///
                 minValue :* (nonMissing :* (result :< minValue))
    }
    
    return(result)
}


//=============================================================================
// VALIDATION FUNCTIONS (Use during development, disable in production)
//=============================================================================

/**
 * Validate patient matrix dimensions match coefficient vector
 * 
 * Critical check before matrix multiplication: pMatrix * coefVector
 * If dimensions don't match, matrix multiplication will fail
 * 
 * @param pMatrix    Patient matrix (n × k)
 * @param coefs      Coefficient row vector (1 × k)
 * @param matrixName Descriptive name for error messages
 */
void validateDimensions(
    real matrix pMatrix,
    real rowvector coefs,
    string scalar matrixName
) {
    if (cols(pMatrix) != cols(coefs)) {
        errprintf("\n========================================\n")
        errprintf("DIMENSION MISMATCH ERROR: %s\n", matrixName)
        errprintf("========================================\n")
        errprintf("Patient matrix columns: %f\n", cols(pMatrix))
        errprintf("Coefficient vector columns: %f\n", cols(coefs))
        errprintf("These must match for matrix multiplication\n")
        errprintf("\nCheck your patient matrix construction:\n")
        errprintf("  - Are you including all required variables?\n")
        errprintf("  - Are reference categories handled correctly?\n")
        errprintf("  - Does coefficient extraction match patient matrix?\n")
        errprintf("========================================\n\n")
        exit(503)
    }
}


/**
 * Validate outcomes are in expected range
 * 
 * Checks that calculated outcomes match expected values
 * Useful for categorical outcomes (BCR, CR, etc.)
 * 
 * @param outcomes   Vector of outcome values
 * @param validValues Row vector of valid outcome values
 * @param outcomeName Descriptive name for warning messages
 */
void validateOutcomes(
    real colvector outcomes,
    real rowvector validValues,
    string scalar outcomeName
) {
    real colvector invalid
    real scalar nInvalid, i
    
    // Check each value against list of valid values
    invalid = J(rows(outcomes), 1, 1)
    for (i = 1; i <= cols(validValues); i++) {
        invalid = invalid :* (outcomes :!= validValues[i])
    }
    // Exclude missing values from check
    invalid = invalid :* (outcomes :< .)
    
    nInvalid = sum(invalid)
    
    if (nInvalid > 0) {
        printf("\n")
        printf("========================================\n")
        printf("WARNING: %f invalid %s values found\n", nInvalid, outcomeName)
        printf("========================================\n")
        printf("Valid values: ")
        validValues
        printf("\n")
        printf("Check your outcome calculation logic\n")
        printf("========================================\n\n")
    }
}


/**
 * Validate random numbers are in interval (0,1)
 * 
 * Random numbers outside (0,1) indicate programming error
 * 
 * @param RN         Vector of random numbers
 * @param varName    Descriptive name for error messages
 */
void validateRandomNumbers(
    real colvector RN,
    string scalar varName
) {
    real scalar nInvalid
    
    // Count non-missing values outside (0,1)
    nInvalid = sum((RN :<= 0) :| (RN :>= 1) :& (RN :< .))
    
    if (nInvalid > 0) {
        errprintf("\n========================================\n")
        errprintf("INVALID RANDOM NUMBERS: %s\n", varName)
        errprintf("========================================\n")
        errprintf("Found %f values outside interval (0,1)\n", nInvalid)
        errprintf("Random numbers must be: 0 < RN < 1\n")
        errprintf("\nCheck random number generation:\n")
        errprintf("  RN = runiform(n, 1)  // generates (0,1)\n")
        errprintf("========================================\n\n")
        exit(498)
    }
}


//=============================================================================
// COEFFICIENT EXTRACTION HELPERS
//=============================================================================
// TODO: Add functions to standardize coefficient extraction from matrices
// 
// Example patterns in current code:
//   - Main effects: bMatrix[1, 1..n]
//   - Auxiliary parameter (Weibull/Gompertz): bMatrix[1, cols(bMatrix)]
//   - Cut points (ordered logit): bMatrix[1, (n+1..n+k)]
// 
// Potential helper functions:
//   - extractMainCoefs(bMatrix, nCoefs)
//   - extractAuxParam(bMatrix) 
//   - extractCutPoints(bMatrix, nMainCoefs, nCutPoints)
//
// This would standardize the pattern and reduce errors
//=============================================================================


//=============================================================================
// EXAMPLE USAGE PATTERNS (for reference)
//=============================================================================
//
// *** SURVIVAL OUTCOME (CI, CD, TNI, OS) ***
//
// mata:
//     // Build patient matrix
//     pCI_L2 = (vAge, vAge2, vMale, vECOG1, vECOG2, vRISS2, vRISS3, 
//               (prevBCR:==3), (prevBCR:==5), vCons)
//     
//     // Extract coefficients
//     coefVector = bCI_L2[1, 1..10]'
//     aux = bCI_L2[1, cols(bCI_L2)]  // Shape/gamma parameter (final column)
//     
//     // Calculate XB
//     XB = pCI_L2 * coefVector
//     
//     // Get valid patients
//     validIdx = getValidPatients(mMOR, OMC, mState, OMC+1)
//     
//     // Calculate survival times
//     RN = runiform(nObs, 1)
//     outcome = J(nObs, 1, .)
//     outcome[validIdx] = calcSurvivalTime(XB[validIdx], RN[validIdx], 
//                                          "weibull", aux)
//     
//     // Curtail at maximum observed
//     outcome = curtailValuesSafe(outcome, maxCI_L2)
// end
//
//
// *** ORDERED LOGIT OUTCOME (BCR) ***
//
// mata:
//     // Build patient matrix
//     pBCR_L3 = (vAge, vAge2, vMale, vECOG1, vECOG2, vRISS2, vRISS3,
//                (prevBCR:==3), (prevBCR:==5))
//     
//     // Extract coefficients and cut points
//     coefVector = bBCR_L3[1, 1..9]'
//     cutPoints = bBCR_L3[1, (10, 11)]  // Final 2 columns are cut points
//     
//     // Calculate XB
//     XB = pBCR_L3 * coefVector
//     
//     // Calculate probabilities and assign outcomes
//     RN = runiform(nObs, 1)
//     probs = calcOrderedLogitProbs(XB, cutPoints)
//     outcome = assignOrderedOutcome(RN, probs, (1, 3, 5))
//     
//     // Get valid patients and update
//     validIdx = getValidPatients(mMOR, OMC, mState, OMC+1)
//     mBCR[validIdx, 3] = outcome[validIdx]
// end
//
//
// *** CONDITIONAL SURVIVAL (OS at later time points) ***
//
// mata:
//     // Build patient matrix
//     pOS_L2S = (vAge, vAge2, vMale, vECOG1, vECOG2, vRISS2, vRISS3, vCons)
//     
//     // Extract coefficients
//     coefVector = bOS[1, 1..8]'
//     aux = bOS[1, cols(bOS)]
//     
//     // Calculate XB
//     XB = pOS_L2S * coefVector
//     
//     // Calculate conditional survival (given already survived to mTSD)
//     outcome = calcConditionalSurvivalTime(XB, mTSD[., OMC], "weibull", aux)
//     
//     // Update OS matrix
//     validIdx = getValidPatients(mMOR, OMC)
//     mOS[validIdx, OMC] = outcome[validIdx]
// end
//
//=============================================================================

end

di as text "{hline 80}"
di as result "Mata utility functions loaded successfully"
di as text "{hline 80}"
di as text "Available functions:"
di as text "  Survival: calcSurvivalTime(), calcConditionalSurvivalTime()"
di as text "  Ordered logit: calcOrderedLogitProbs(), assignOrderedOutcome()"
di as text "  Filtering: getValidPatients(), isAlive()"
di as text "  Utilities: curtailValues(), curtailValuesSafe()"
di as text "  Validation: validateDimensions(), validateOutcomes(), validateRandomNumbers()"
di as text "{hline 80}"
