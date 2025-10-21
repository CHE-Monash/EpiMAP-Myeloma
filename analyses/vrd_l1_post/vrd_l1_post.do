**********
*EpiMAP Myeloma - VRd Line 1 Post-Market Analysis
**********

**********
*Analysis Configuration
	local analysis_name "vrd_l1_post"
	local coefficient_file "coefficients_vrd_l1_post"

**********
*Load Programs
	quietly do "core/mata_functions.do"
	quietly do "core/load_patients.do"
	quietly do "core/matrix_setup.do"	
	quietly do "core/simulation_engine.do"
	quietly do "core/process_data.do"
	
**********
*Determine processing approach
	if("$Boot" == "0") {
		// No Bootstrapping
		mata: mata matuse "$coefficients_path/`coefficient_file'"
				
		mata_functions
		load_patients		
		matrix_setup
		simulation
		process_data
		
		save "$simulated_path/$Int $Line $Data $MinID $MaxID.dta", replace
		di "Analysis completed: $simu;ated_path/$Int $Line $Data $MinID $MaxID.dta"
	}
	else {
		// Bootstrapping
		forvalues b = $MinBS / $MaxBS {
			di as text "Processing bootstrap iteration `b' of $MaxBS..."
			
			mata: mata clear
			mata: mata matuse "$coefficients_path/bootstrap/`coefficient_file'_B`b'"
			
			mata_functions
			load_patients		
			matrix_setup
			simulation
			process_data
			
			save "$simulated_path/bootstrap/$Int $Line $Data $MinID $MaxID Bootstrap_B`b'.dta", replace
			di "Bootstrap iteration `b' completed: Bootstrap_B`b'.dta"
		}
		
		di "All bootstrap iterations completed. Files saved as Bootstrap_B[MinBS-MaxBS].dta"
	}
