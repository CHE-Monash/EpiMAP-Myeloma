**********
	*EpiMAP Myeloma Analysis
**********
	
	clear
	set more off
		
	*Set Analysis Directory
		 cd "~/Documents/Monash/Research/Blood Disorders/Myeloma/EpiMAP/Admin/Github/v2.0"
		
	*Globals
		global core ID Event0 Date0 Event1 Date1
		global pt Age Male ECOGc ISS
		global tx SCT MNT BCR Regimen
		global surv _st _d _origin _t _t0
		
**********			
*Simulation

	*Arg 1 - Analysis
	*Arg 2 - Int // If Int != All, SIM CR $Analysis
	*Arg 3 - Line
	*Arg 4 - Coeffs	
	*Arg 5 - Data // Population or Predicted
	*Arg 6 - MinID
	*Arg 7 - MaxID
	*Arg 8 - Bootstrap
	*Arg 9 - MinBS
	*Arg 10 - MaxBS		
	
	
	*VRd-Post				
		*No Bootstrap 
			do "EpiMAP_Myeloma_v2.0.do" VRd-Post VRd 1 VRd Predicted 1 4884 0
			do "EpiMAP_Myeloma_v2.0.do" VRd-Post NoVRd 1 VRd Predicted 1 4884 0
			
			do "Simulation.do" VRd-Post NoVRd 1 VRd Cohort10 1 1000 0	// 48741	
				
		*Bootstrap
			do "Simulation.do" VRd-Post VRd 1 VRd Predicted 1 100 1 1 5
			do "Simulation.do" VRd-Post NoVRd 1 VRd Predicted 1 100 1 1 5
