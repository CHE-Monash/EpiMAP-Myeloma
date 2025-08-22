**********
	*EpiMAP Myeloma v2.0 - Main Dispatcher
**********
	
	clear
	clear mata
	clear matrix
	set more off
	
**********
*Capture arguments 
	global Analysis `1'      // Analysis name (base_model, vrd_l1_post, etc.)
	global Int `2'           // Intervention
	global Line `3'          // Line of therapy
	global Coeffs `4'        // Coefficient set
	global Data `5'          // Data type (Population, Predicted, Cohort10, Population_3, etc.)
	global MinID `6'         // Min patient ID
	global MaxID `7'         // Max patient ID
	global Boot `8'          // Bootstrap flag
	global MinBS `9'         // Min bootstrap sample
	global MaxBS `10'        // Max bootstrap sample

**********
*Parse Data parameter for population-specific requests
	if regexm("$Data", "Population_([0-9]+)") {
		global PopulationNumber = regexs(1)
		global DataType = "Population"
		di "Using specific population: ${PopulationNumber}"
	}
	else {
		global DataType = "$Data"
		global PopulationNumber = 1  // Default population
	}

**********
*Validate Analysis exists
	capture confirm file "analyses/$Analysis/${Analysis}.do"
	if _rc != 0 {
		di as error "Error: Analysis '$Analysis' not found."
		di as error "Available analyses:"
		
		local analyses : dir "analyses" dirs "*"
		foreach analysis of local analyses {
			di "  - `analysis'"
		}
		exit 601
	}

**********
*Validate Population data if needed
	if ("$DataType" == "Population") {
		capture confirm file "data/2025-2030/population_${PopulationNumber}.dta"
		if _rc != 0 {
			di as error "Error: Population ${PopulationNumber} not found."
			di as error "Available populations:"
			local populations : dir "data/populations" files "EpiMAP_Population_*.dta"
			foreach pop of local populations {
				di "  - `pop'"
			}
			exit 601
		}
	}

**********
*Set paths based on analysis
	global analysis_path "analyses/$Analysis"
	global coefficients_path "$analysis_path/data/coefficients"
	global patients_path "$analysis_path/data/patients" 
	global simulated_path "$analysis_path/data/simulated"
	global populations_path "data/populations"
	
	// Create results directory if it doesn't exist
	capture mkdir "$results_path"

**********
*Run analysis
	di as text "Running EpiMAP Myeloma Analysis: $Analysis"
	di as text "Intervention: $Int, Line: $Line, Data: $DataType"
	if ("$DataType" == "Population") di as text "Population: ${PopulationNumber}"
	di as text "Patient range: $MinID to $MaxID"
	
	// Run the specific analysis
	do "$analysis_path/${Analysis}.do"
	
	di "Analysis completed successfully."
