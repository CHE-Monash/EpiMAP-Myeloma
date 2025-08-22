**********
	*EpiMAP VRd Line 1 Post-Market Analysis
**********

**********
*Analysis Configuration
	local analysis_name "vrd_l1_post"
	local coefficient_file "coefficients_vrd_l1_post"

**********
*Load Coefficients
	if("$Boot" == "0") {
		mata: mata matuse "$coefficients_path/`coefficient_file'"
	}
	else {
		forvalues b = $MinBS / $MaxBS {
			mata: mata matuse "$coefficients_path/`coefficient_file'_B`b'"
		}
	}

**********
*Load Patient Data
	// Determine data source based on $Data parameter
	if ("$Data" == "Population") {
		// Use one of the base population datasets
		local pop_number = 1  // Default to population 1, could be parameterized
		use "data/2025-2030/population_`pop_number'.dta", clear
		di as text "Loaded base population `pop_number'"
	}
	else if ("$Data" == "Predicted") {
		// Use analysis-specific predicted data
		use "$patients_path/patients_vrd_l1_post.dta", clear
		di as text "Loaded VRd L1 specific patient data"
	}
	else if ("$Data" == "Cohort10") {
		// Use a specific cohort from population data
		use "data/2025-2030/population_1.dta", clear
		keep if _n <= 1000  // Take first 1000 patients
		di as text "Loaded cohort of 1000 patients from population 1"
	}
	else {
		di as error "Unknown data type: $Data"
		exit 198
	}

**********
*VRd-specific data modifications
	// Apply VRd-specific patient selection criteria or modifications
	if ("$Int" == "VRd") {
		// Example: VRd might have specific eligibility criteria
		// keep if age >= 18 & age <= 80  // Age restrictions
		// replace some_variable = new_value if condition
		
		di as text "Applied VRd intervention modifications"
	}

**********
*Standard data preparation
	*State filter
		keep if State <= ($Line * 2) + 1
		replace ID = _n

	*ID filter
		keep if ID >= $MinID & ID <= $MaxID
		quietly sum ID
		scalar Obs = r(N)
		di as text "Final sample size: `=Obs' patients"

	*Age limit
		scalar Limit = 100
		
**********
*Matrix Setup
	quietly do "core/matrix_setup.do"
	matrix_setup

**********
*Execute simulation
	quietly do "core/simulation_engine.do"
	simulation
	
**********
*Process results
	quietly do "core/process_data.do"
	process_data

**********
*Save results
	save "$results_path/$Int $Line $Data $MinID $MaxID.dta", replace
	di as text "Results saved to: $results_path/$Int $Line $Data $MinID $MaxID.dta"
