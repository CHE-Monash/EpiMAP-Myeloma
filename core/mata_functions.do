**********
* Monash Myeloma Model - Mata Utility Functions
*
* Purpose: Reusable Mata functions for survival calculations, ordered logit, patient
*          filtering, value truncation and validation. Used across the outcome modules.
* Usage:   Load once per (cleared) Mata state, e.g. run "core/mata_functions.do".
**********

*Functions for sim_txr.do
mata:

// Helper function: Get TXR coefficients for a given line
real matrix get_txr_coef(real scalar line) {
	external bL1_TXR, bL2_TXR, bL3_TXR, bL4_TXR, bL5_TXR, bL6_TXR, bL7_TXR, bL8_TXR, bL9_TXR
    if (line == 1 & rows(bL1_TXR) > 0) return(bL1_TXR)
    if (line == 2 & rows(bL2_TXR) > 0) return(bL2_TXR)
    if (line == 3 & rows(bL3_TXR) > 0) return(bL3_TXR)
    if (line == 4 & rows(bL4_TXR) > 0) return(bL4_TXR)
    if (line == 5 & rows(bL5_TXR) > 0) return(bL5_TXR)
    if (line == 6 & rows(bL6_TXR) > 0) return(bL6_TXR)
    if (line == 7 & rows(bL7_TXR) > 0) return(bL7_TXR)
    if (line == 8 & rows(bL8_TXR) > 0) return(bL8_TXR)
    if (line == 9 & rows(bL9_TXR) > 0) return(bL9_TXR)
    return(J(0, 0, .))  // Empty matrix if doesn't exist
}

// Helper function: Get TXR outcome codes
real rowvector get_txr_outcome(real scalar line) {
	external oL1_TXR, oL2_TXR, oL3_TXR, oL4_TXR, oL5_TXR, oL6_TXR, oL7_TXR, oL8_TXR, oL9_TXR
    if (line == 1) return(oL1_TXR)
    if (line == 2) return(oL2_TXR)
    if (line == 3) return(oL3_TXR)
    if (line == 4) return(oL4_TXR)
    if (line == 5) return(oL5_TXR)
    if (line == 6) return(oL6_TXR)
    if (line == 7) return(oL7_TXR)
    if (line == 8) return(oL8_TXR)
    if (line == 9) return(oL9_TXR)
	return(J(1, 0, .))  // Empty rowvector if not found
}

// Helper function: Check if TXR model exists
real scalar txr_model_exists(real scalar line) {
    return(rows(get_txr_coef(line)) > 0)
}

end

*Functions for sim_mnr.do / sim_mnd.do
mata:

// Same `external` trick as the TXR helpers above, and for the same reason: when an analysis
// declares no maintenance regimens, risk_equations.do stores oL1_MNR = 0 and never creates
// bL1_MNR, so it is absent from the coefficient file. Referencing it directly would fail on
// an undefined symbol; declaring it external creates it empty, so rows() can be tested.

// Helper function: Get MNR coefficients (empty matrix if the analysis modelled no regimens)
real matrix get_mnr_coef() {
	external bL1_MNR
	if (rows(bL1_MNR) > 0) return(bL1_MNR)
	return(J(0, 0, .))
}

// Helper function: Get MNR outcome codes (the drug codes the analysis modelled)
real rowvector get_mnr_outcome() {
	external oL1_MNR
	if (cols(oL1_MNR) > 0) return(oL1_MNR)
	return(J(1, 0, .))
}

// Helper function: Check if MNR model exists
real scalar mnr_model_exists() {
	return(rows(get_mnr_coef()) > 0)
}

// Helper function: Get MND (maintenance duration) coefficients, LENALIDOMIDE
// L1_MND is fitted per REGIMEN and pooled across transplant, not the other way round. The two drugs
// are different processes - thalidomide is a fixed course (9% censored), lenalidomide runs on (47%)
// - while the two transplant arms have near-identical Kaplan-Meier curves. SCT enters as a covariate
// instead. See prep/risk_equations.do.
real matrix get_mnd_coef_len() {
	external bL1_MND_LEN
	if (rows(bL1_MND_LEN) > 0) return(bL1_MND_LEN)
	return(J(0, 0, .))
}

// Helper function: Get MND coefficients, THALIDOMIDE (censored at 18 months in the fit)
real matrix get_mnd_coef_thal() {
	external bL1_MND_THAL
	if (rows(bL1_MND_THAL) > 0) return(bL1_MND_THAL)
	return(J(0, 0, .))
}

// Helper function: Get LENREFR_TX residual-arm logit coefficients (empty if not fitted)
// `external' creates the symbol empty when an analysis did not fit the block, so rows() can be
// tested and sim_lenrefr.do becomes a no-op rather than failing on an undefined symbol.
real matrix get_lenrefr_coef() {
	external bLENREFR_TX
	if (rows(bLENREFR_TX) > 0) return(bLENREFR_TX)
	return(J(0, 0, .))
}

// Helper function: Get LENREFR_MNT (L1 maintenance) logit coefficients (empty if not fitted)
real matrix get_mntrefr_coef() {
	external bLENREFR_MNT
	if (rows(bLENREFR_MNT) > 0) return(bLENREFR_MNT)
	return(J(0, 0, .))
}

// Helper function: Get the lenalidomide-containing regimen codes the analysis declared
real rowvector get_lenrefr_regimens() {
	external LENREFR_regimens
	if (cols(LENREFR_regimens) > 0) return(LENREFR_regimens)
	return(J(1, 0, .))
}

// Helper function: Get P(lenalidomide | regimen not modelled) per line - the 'other'-bucket gate
// probability, one column per line. Empty when an analysis did not fit it, which makes sim_lenrefr
// fall back to the strict modelled-regimen gate rather than failing.
real rowvector get_lenrefr_pother() {
	external LENREFR_pother
	if (cols(LENREFR_pother) > 0) return(LENREFR_pother)
	return(J(1, 0, .))
}

// Helper function: LenRefr treatment arm runs only if BOTH the logit and a len-regimen list exist.
// The maintenance arm has its own getter (get_mntrefr_coef) - it is a fitted logit too, since the
// tail rule it used to apply was replaced.
real scalar lenrefr_model_exists() {
	return(rows(get_lenrefr_coef()) > 0 & cols(get_lenrefr_regimens()) > 0)
}

// Helper function: MND model exists if EITHER transplant arm was fitted
real scalar mnd_model_exists() {
	return(cols(get_mnd_coef_len()) > 0 | cols(get_mnd_coef_thal()) > 0)
}

end

mata:
//=============================================================================
// SURVIVAL DISTRIBUTION FUNCTIONS
//=============================================================================

/**
 * Calculate survival time from parametric distributions
 * 
 * @param XB       Linear predictor (vector or scalar)
 * @param RN       Random numbers in (0,1) (vector or scalar)
 * @param dist     Distribution type: "exponential", "ereg", "weibull", "gompertz", "llogistic", "lnormal"
 * @param aux      Auxiliary parameter (ln_p for Weibull, gamma for Gompertz, ln_gam for llogistic, ln_sig for lnormal).
 *                 Usually the final column of the coefficient matrix, extracted as
 *                 aux = bMatrix[1, cols(bMatrix)] (a single shared shape). May instead
 *                 be a per-observation colvector if a model varies the shape by covariate.
 *                 All uses of aux below are elementwise, so a scalar or colvector both work.
 * @return         Survival times
 */
real colvector calcSurvTime(
    real colvector XB,
    real colvector RN,
    string scalar dist,
    | real colvector aux 	// Optional; scalar or per-observation colvector
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
    else if (dist == "llogistic") {
        // Log-logistic (AFT): T = exp(XB) * (1/U - 1)^gamma, gamma = exp(aux) [aux = ln_gam].
        // exp(XB) is the median/scale; heavier tail than Weibull. Verified against streg to machine
        // precision. Finite mean but infinite variance (gamma<1) -> very heavy simulated tail.
        outcome = exp(XB) :* ((1 :/ RN :- 1) :^ exp(aux))
    }
    else if (dist == "lnormal") {
        // Log-normal (AFT): ln(T) ~ N(XB, sigma^2), sigma = exp(aux) [aux = ln_sig]. Invert S(T)=U:
        // T = exp(XB + sigma*invnormal(1-U)). Heavier-than-Weibull tail but FINITE variance (lighter
        // extreme tail than log-logistic -- TFI uses this). Verified against streg to machine precision.
        outcome = exp(XB :+ exp(aux) :* invnormal(1 :- RN))
    }
    else {
        errprintf("ERROR: Unknown distribution '%s'\n", dist)
        errprintf("Valid options: 'exponential', 'ereg', 'weibull', 'gompertz', 'llogistic', 'lnormal'\n")
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
 * @param dist     Distribution type ("exponential"/"ereg", "weibull", "gompertz", "llogistic", "lnormal")
 * @param aux      Auxiliary parameter (shape/gamma parameter). Scalar, or a per-observation
 *                 colvector when the shape varies by covariate (OS: per-BCR Weibull shape).
 * @return         Survival probabilities S(TSD)
 */
real colvector calcSurvProb(
    real colvector XB,
    real colvector TSD,
    string scalar dist,
    | real colvector aux
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
    else if (dist == "llogistic") {
        // Log-logistic (AFT): S(t) = 1 / (1 + (t/exp(XB))^(1/gamma)), gamma = exp(aux) [aux = ln_gam].
        PR = 1 :/ (1 :+ ((TSD :/ exp(XB)) :^ (1 :/ exp(aux))))
    }
    else if (dist == "lnormal") {
        // Log-normal (AFT): S(t) = 1 - Phi((ln t - XB)/sigma) = Phi((XB - ln t)/sigma), sigma = exp(aux)
        PR = normal((XB :- ln(TSD)) :/ exp(aux))
    }
    else {
        errprintf("ERROR: Unknown distribution '%s'\n", dist)
        errprintf("Valid distributions: 'exponential', 'ereg', 'weibull', 'gompertz', 'llogistic', 'lnormal'\n")
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

end
