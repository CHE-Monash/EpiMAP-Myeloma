**********
*EpiMAP Myeloma - Base Model Analysis
**********

**********
*Analysis Configuration
	local analysis_name "base_model"
	local coefficient_file "coefficients_base_model"

**********
*Load Programs
	quietly do "core/mata_functions.do"
	quietly do "core/load_patients.do"
	quietly do "core/mata_setup.do"
	quietly do "core/simulation_engine.do"
	quietly do "core/process_data.do"
	
**********
*Determine processing approach
	if("$Boot" == "0") {
		// No Bootstrapping
		mata: mata matuse "$coefficients_path/`coefficient_file'"
				
		load_patients
		mata_setup
		simulation
		process_data
		
		save "$simulated_path/$Int $Line $Data $MinID $MaxID.dta", replace
		di as text "Results saved to: $simulated_path/$Int $Line $Data $MinID $MaxID.dta"
	}
	else {
		// Bootstrapping
		di as text "Running bootstrap analysis with $MaxBS iterations..."
		
		forvalues b = $MinBS / $MaxBS {
			di as text _newline "========================================" 
			di as text "Processing bootstrap iteration `b' of $MaxBS..."
			di as text "========================================"
			
			mata: mata clear
			mata: mata matuse "$coefficients_path/bootstrap/`coefficient_file'_B`b'"
			
			load_patients
			mata_setup
			simulation
			process_data
			
			save "$simulated_path/bootstrap/$Int $Line $Data $MinID $MaxID Bootstrap_B`b'.dta", replace
			di as text "Bootstrap iteration `b' completed: Bootstrap_B`b'.dta"
		}
		
		di as result _newline "All bootstrap iterations completed successfully!"
		di as text "Bootstrap results saved in: $simulated_path/bootstrap/"
		di as text "Files saved as: Bootstrap_B[$MinBS-$MaxBS].dta"
	}
