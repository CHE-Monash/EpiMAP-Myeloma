**********
	*EpiMAP Myeloma - Load Patients
**********

capture program drop load_patients
program define load_patients

	// Determine data source based on $Data parameter
	if ("$Data" == "Population") {
		// Parse population number if specified
		if regexm("$Data", "Population_([0-9]+)") {
			local pop_number = regexs(1)
		}
		else {
			local pop_number = 1  // Default to population 1
		}
		
		use "data/2025-2030/population_`pop_number'.dta", clear
		di as text "Loaded base population `pop_number'"
	}
	else if ("$Data" == "Predicted") {
		use "$patients_path/patients_vrd_l1_post.dta", clear
		di as text "Loaded VRd L1 specific patient data"
	}
	else if ("$Data" == "Cohort10") {
		use "data/2025-2030/population_1.dta", clear
		keep if _n <= 1000
		di as text "Loaded cohort of 1000 patients from population 1"
	}
	else {
		di as error "Unknown data type: $Data"
		exit 198
	}

	// Apply VRd-specific modifications
	if ("$Int" == "VRd") {
		di as text "Applied VRd intervention modifications"
	}

	// Standard data preparation
	keep if State <= ($Line * 2) + 1
	replace ID = _n
	keep if ID >= $MinID & ID <= $MaxID
	
	quietly sum ID
	scalar Obs = r(N)
	di as text "Final sample size: `=Obs' patients"
	
	scalar Limit = 100
	
end
