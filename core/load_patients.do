**********
*EpiMAP Myeloma - Load Patients
**********

capture program drop load_patients
program define load_patients

	// Determine data source based on $DataType (not $Data)
	if ("$DataType" == "population") {
		// Use population number extracted by main dispatcher
		use "patients/population/2025-2030/population_${PopulationNumber}.dta", clear
		di as text "Loaded population ${PopulationNumber}"
	}
	else if ("$DataType" == "predicted") {
		use "$patients_path/patients_${Analysis}.dta", clear
		di as text "Loaded ${Analysis} specific patient data"
	}
	else {
		di as error "Unknown data type: $DataType"
		di as error "Valid options: Population_#, Predicted"
		exit 198
	}

	// Apply intervention-specific modifications if needed
	if ("$Int" == "VRd") {
		di as text "Applied VRd intervention modifications"
	}

	// Standard data preparation and filtering
	keep if State <= ($Line * 2) + 2
	replace ID = _n
	keep if ID >= $MinID & ID <= $MaxID
	
	quietly sum ID
	scalar Obs = r(N)
	di as text "Final sample size: `=Obs' patients (IDs $MinID to $MaxID)"
	
	scalar Limit = 100 // Can be removed once model is fully matarised
	mata: Limit = 100
	
end
