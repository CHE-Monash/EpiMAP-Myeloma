**********
* Monash Myeloma Model - Mata Setup
*
* Purpose: Load the cohort into Mata and allocate the working structures for the
*          simulation pass: fixed patient characteristics as vectors (Age, Male, ECOG,
*          RISS, ...) and line-varying outcomes as matrices with one column per pathway
*          point (mOS, mBCR, mTXR, ...). Also builds the common-random-number matrix mRN.
* Notes:   CRN seed = $crn_seed_base (+ bootstrap iteration $b); requires canonical
*          ID == _n row order (asserted) for cross-arm alignment. See core/rng_slots.do.
**********

cap program drop mata_setup
program define mata_setup
	
	di "Setting up matrices"
	mata: Obs = st_nobs()

	**********
	* Common Random Numbers - build the pre-drawn uniform matrix mRN (Obs x K)
	*   Keyed by patient row x event slot (see core/rng_slots.do). Identical across
	*   arms within a replication (same seed, same cohort order) and independent
	*   across replications / bootstrap iterations (seed = base + b). CRN is
	*   unconditional: every migrated event reads its column from mRN.
	**********
	// Ensure the slot registry is loaded (run_pipeline loads it; guard direct callers)
	capture mata: st_numscalar("__rn_K", rn_K())
	if _rc {
		run "core/rng_slots.do"
		mata: st_numscalar("__rn_K", rn_K())
	}

	// Seed base (override by setting global crn_seed_base before the run)
	if ("$crn_seed_base" == "") global crn_seed_base 20260615

	// Replication offset: the bootstrap iteration when bootstrapping, else 0.
	//   Same value for both arms within a replication; differs across replications.
	local _b = 0
	if ("$boot" == "1" & "$b" != "") local _b = $b
	local _crn_seed = $crn_seed_base + `_b'

	// Cross-arm alignment requires canonical row order (load_patients resets ID=_n).
	//   A missing reset or a re-sorted cohort would silently break CRN - fail loud.
	capture assert ID == _n
	if _rc {
		di as error "mata_setup: cohort is not in canonical ID==_n order - CRN alignment unsafe."
		error 459
	}

	set seed `_crn_seed'
	mata: mRN = runiform(Obs, rn_K())
	di as text "CRN: built mRN " _N " x " scalar(__rn_K) " (seed `_crn_seed')"


	*mState - State matrix
	mata {
		mState = st_data(.,"State DateDN")
		rmState = J(Obs,2,"")
		rmState[.,2] = strofreal(mState[.,1])
		cmState = J(2,2,"")
		cmState[1,1] = "State"
		cmState[2,1] = "DateDN"
	}		
	*_matrix_list(mState, rmState, cmState)
		
	*mAge - Age (at event)
	mata {
		mAge = st_data(., "Age_DN Age_L1S Age_L1E Age_L2S Age_L2E Age_L3S Age_L3E Age_L4S Age_L4E Age_L5S Age_L5E Age_L6S Age_L6E Age_L7S Age_L7E Age_L8S Age_L8E Age_L9S Age_L9E")
		mAge[., 1] = round(mAge[., 1], 0.1)  // Round Age_DN to fix floating-point precision
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
	*mata: _matrix_list(mAge, rmAge, cmAge)		

	*mOS - Overall Survival (from DN)
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
	*mata: _matrix_list(mOS, rmOS, cmOS)		

	*mTNE - Time to Next Event	
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
	*mata: _matrix_list(mTNE, rmTNE, cmTNE)	
		
	*mTSD - Time Since Diagnosis
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
	*mata: _matrix_list(mTSD, rmTSD, cmTSD)		
		
	*mMOR - Mortality (=1 if patient dies before next event) 
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
	*mata: _matrix_list(mMOR, rmMOR, cmMOR)

	*mOC - Outcome
	mata {
		mOC = J(Obs,2,.)				
		rmOC = J(Obs, 2, "")
		rmOC[., 2] = strofreal(1::Obs)
		cmOC = J(2, 2, "")
		cmOC[1,1] = "Time"
		cmOC[2,1] = "Mort"
	}
	*mata: _matrix_list(mOC, rmOC, cmOC)	
		
	*mTXR - Treatment Regimen 	
	mata {
		mTXR = st_data(., "TXR_L1 TXR_L2 TXR_L3 TXR_L4 TXR_L5 TXR_L6 TXR_L7 TXR_L8 TXR_L9")		
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
	*mata: _matrix_list(mTXR, rmTXR, cmTXR)
		
	*mTXD - Treatment Duration
	mata {
		mTXD = st_data(., "TXD_L1 TXD_L2 TXD_L3 TXD_L4 TXD_L5 TXD_L6 TXD_L7 TXD_L8 TXD_L9")		
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
	*mata: _matrix_list(mTXD, rmTXD, cmTXD)	
		
	*mBCR - Best Clinical Response 
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
	*mata: _matrix_list(mBCR, rmBCR, cmBCR)
		
	*mTFI - Treatment-free Interval 
	mata {
		mTFI = st_data(., "TFI_DN TFI_L1 TFI_L2 TFI_L3 TFI_L4 TFI_L5 TFI_L6 TFI_L7 TFI_L8")		
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
	*mata: _matrix_list(mTFI, rmTFI, cmTFI)

	di "Setting up vectors"
	
	*Generate categorical dummy variables 
		cap gen ECOGcc0 = (ECOGcc == 0)
		cap gen ECOGcc1 = (ECOGcc == 1)
		cap gen ECOGcc2 = (ECOGcc == 2)
		cap gen RISS1 = (RISS == 1)
		cap gen RISS2 = (RISS == 2)
		cap gen RISS3 = (RISS == 3)
		cap gen ISS1 = (ISS == 1)
		cap gen ISS2 = (ISS == 2)
		cap gen ISS3 = (ISS == 3)
		cap gen Cons = 1
	
	*Vectors
	mata {
		vID = st_data(., "ID")					  // Patient ID
		vAge = st_data(., "Age_DN")               // Age at diagnosis (continuous)
		vAge2 = vAge :^ 2                         // Age squared (for quadratic effects)
		vMale = st_data(., "Male")                // Sex (0=Female, 1=Male)	
		vECOG = st_data(., "ECOGcc")              // Original ECOG (0, 1, 2+)
		vECOG0 = st_data(., "ECOGcc0")            // ECOG = 0 (fully active)
		vECOG1 = st_data(., "ECOGcc1")            // ECOG = 1 (restricted strenuous activity)
		vECOG2 = st_data(., "ECOGcc2")            // ECOG = 2+ (ambulatory, self-care)
		vRISS = st_data(., "RISS")                // Original RISS (1, 2, 3)
		vRISS1 = st_data(., "RISS1")              // RISS = I (best prognosis)
		vRISS2 = st_data(., "RISS2")              // RISS = II (intermediate)
		vRISS3 = st_data(., "RISS3")              // RISS = III (poorest prognosis)
		
		vISS = st_data(., "ISS")
		vISS1 = st_data(., "ISS1")              
		vISS2 = st_data(., "ISS2")              
		vISS3 = st_data(., "ISS3") 	

		vAge70 = st_data(., "Age70")              // Age >= 70 indicator
		vAge75 = st_data(., "Age75")              // Age >= 75 indicator
		vCKD = st_data(., "CM_CKD")				  // Chronic Kidney Disease
		vCRD = st_data(., "CM_CRD")				  // Moderate-to-severe cardiac disease
		vPLM = st_data(., "CM_PLM")				  // Moderate-to-severe chronic lung disease
		vDBT = st_data(., "CM_DBT")				  // Diabetes (insulin / oral hypoglycaemics)

		vCons = st_data(., "Cons")                // Constant vector (all ones)

		vSCT_DN = st_data(., "SCT_DN")            // Intent for ASCT at DN (0/1)
		vSCT_L1 = st_data(., "SCT_L1")            // Receipt of ASCT at L1 (0/1)
		vMNT = st_data(., "MNT")                  // Receipt of Maintenance therapy at L1 (0/1)

		// Maintenance regimen and duration share at L1. Pure OUTPUTS - unlike the vectors
		// above there is no cohort column to read, so they are initialised missing and
		// filled by sim_mnr.do / sim_mnd.do among MNT == 1. Missing is the correct resting
		// value: a patient with MNT == 0 has no maintenance regimen and no duration.
		// vMND is the maintenance DURATION in months (the raw survival draw); process_data.do
		// caps it at the realised TFI_L1 at billing time. See docs/refractory.md 4.4.
		vMNR = J(Obs, 1, .)                       // L1 maintenance regimen (1 = len, 5 = thal)
		vMND = J(Obs, 1, .)

		// Lenalidomide-refractory, the LATCHED at-line-entry state. 0 = not yet refractory, which is
		// the correct resting value and is 0 at L1 entry by construction. sim_lenrefr.do flips it at
		// a line's END (after that line's OS), so every equation reads the state from STRICTLY PRIOR
		// lines. See docs/refractory.md 3.5 / 4.
		//
		// ONE vector, not two. The treatment-dose and maintenance-dose flags were kept separate
		// while it was open whether they carried different prognostic weight; os_lenrefr_check.do
		// settled it at p = 0.97, so they collapse and the union is what the OS equations read.
		//
		// READ IT IN for a line-entry analysis. At $line 1 the state is 0 by construction (no prior
		// lines), but a $line >= 2 analysis starts its patients mid-pathway, where refractoriness
		// acquired earlier is real and already recorded: process_data.do exports LenRefr_L1..L9 and
		// cohort_pool.do carries them through, so the value exists in the pool. Initialising to 0
		// regardless would hand every patient the NON-refractory branch of the L2-L4 OS equations,
		// which now carry a LenRefr_any coefficient - over-predicting survival relative to the
		// population those equations were fitted on. mBCR is read from the cohort for exactly this
		// reason; this vector was the one piece of line-entry state that was not.
		//
		// TWO CASES, distinguished by $data_type (set by load_patients.do), not by $line alone.
		// A pool build sets $line to the decision line but SIMULATES FROM DIAGNOSIS on a synthetic
		// incidence cohort, so the state must start at 0 and be generated - warning there would be
		// noise on every pool build. Only a cohort that genuinely ENTERS at the line needs it read in.
		lrLine = strtoreal(st_global("line"))
		if (lrLine < 1 | lrLine >= .) lrLine = 1          // $line 0 = full pathway from diagnosis
		lrName = sprintf("LenRefr_L%g", lrLine)
		lrFromDx = (st_global("data_type") == "synthetic")

		if (lrFromDx | lrLine == 1) {
			// Simulating the pathway from diagnosis: 0 is correct and the engine fills it in.
			vLenRefr_in = J(Obs, 1, 0)
		}
		else if (_st_varindex(lrName) < .) {
			vLenRefr_in = editmissing(st_data(., lrName), 0)
			printf("  mata_setup: len-refractory state read in from %s (mean %5.3f)\n",
			       lrName, mean(vLenRefr_in))
		}
		else {
			vLenRefr_in = J(Obs, 1, 0)
			errprintf("mata_setup: %s is not in this line-entry cohort, so every patient enters\n", lrName)
			errprintf("            NON-refractory. The L2-L4 OS equations carry a LenRefr_any term,\n")
			errprintf("            so this over-predicts survival. Rebuild the cohort pool with a\n")
			errprintf("            model version that exports LenRefr_L*.\n")
		}

		// Per-line SNAPSHOT for export and validation: mLenRefr_in[.,l] is the entry-to-Ll value,
		// recorded at each line end before the state updates. Missing for a line never reached,
		// matching mBCR. process_data.do writes it out as LenRefr_L1..L9.
		mLenRefr_in = J(Obs, 9, .)                       // L1 maintenance duration, months (capped at TFI when billed)

	}
end
