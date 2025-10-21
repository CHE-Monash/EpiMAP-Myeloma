**********
	*SIM CR L1
**********

	*Patient matrix	
		mata: `=m' = mCore , mML
		mata: r`=m' = rmCore
		mata: c`=m' = cmCore \ cmML
		mata: `=m'[.,`=c'RN] = runiform(`=Obs',1)	
		mata: _matrix_list(`=m', r`=m', c`=m')
		
	*Create nL1
		mata: st_matrix("stL1_CR", oL1_CR)
		scalar nL1 = colsof(stL1_CR)
		
	*Determine outcome 
		forvalues i = 1/`=Obs' {
			mata {
				if (mMOR[`i',`=OMC'-1] == 0 & mState[`i',1] <= `=OMC'+1) { // Alive & State filters
				
					`=m'[`i',`=c'XB1] = 1 // Set e(XB1) to 1
					
					if (`=nL1' >= 2) { // Calculate e(XB2)
						*Age
							`=m'[`i',`=c'XB2] = `=m'[`i',`=c'XB2] + (`=m'[`i',cAge] * `=b'[1,12])
						*Age2
							`=m'[`i',`=c'XB2] = `=m'[`i',`=c'XB2] + (`=m'[`i',cAge]^2 * `=b'[1,13])							
						*Male
							`=m'[`i',`=c'XB2] = `=m'[`i',`=c'XB2] + (`=m'[`i',cMale] * `=b'[1,14])
						*ECOG 
							if (`=m'[`i',cECOG] == 1) `=m'[`i',`=c'XB2] = `=m'[`i',`=c'XB2] + `=b'[1,16] 
							if (`=m'[`i',cECOG] == 2) `=m'[`i',`=c'XB2] = `=m'[`i',`=c'XB2] + `=b'[1,17] 
						*RISS 
							if (`=m'[`i',cRISS] == 2) `=m'[`i',`=c'XB2] = `=m'[`i',`=c'XB2] + `=b'[1,19] 
							if (`=m'[`i',cRISS] == 3) `=m'[`i',`=c'XB2] = `=m'[`i',`=c'XB2] + `=b'[1,20]
						*SCT
							`=m'[`i',`=c'XB2] = `=m'[`i',`=c'XB2] + (`=m'[`i',cSCT] * `=b'[1,21])
						*cons
							`=m'[`i',`=c'XB2] = `=m'[`i',`=c'XB2] + `=b'[1,22]
						*Exponent
							`=m'[`i',`=c'XB2] = exp(`=m'[`i',`=c'XB2])
					}
					
					if (`=nL1' >= 3) { // Calculate e(XB3)
						*Age
							`=m'[`i',`=c'XB3] = `=m'[`i',`=c'XB3] + (`=m'[`i',cAge] * `=b'[1,23])
						*Age2
							`=m'[`i',`=c'XB3] = `=m'[`i',`=c'XB3] + (`=m'[`i',cAge]^2 * `=b'[1,24])								
						*Male
							`=m'[`i',`=c'XB3] = `=m'[`i',`=c'XB3] + (`=m'[`i',cMale] * `=b'[1,25])
						*ECOG 
							if (`=m'[`i',cECOG] == 1) `=m'[`i',`=c'XB3] = `=m'[`i',`=c'XB3] + `=b'[1,27] 
							if (`=m'[`i',cECOG] == 2) `=m'[`i',`=c'XB3] = `=m'[`i',`=c'XB3] + `=b'[1,28] 
						*RISS 
							if (`=m'[`i',cRISS] == 2) `=m'[`i',`=c'XB3] = `=m'[`i',`=c'XB3] + `=b'[1,30] 
							if (`=m'[`i',cRISS] == 3) `=m'[`i',`=c'XB3] = `=m'[`i',`=c'XB3] + `=b'[1,31]
						*SCT
							`=m'[`i',`=c'XB3] = `=m'[`i',`=c'XB3] + (`=m'[`i',cSCT] * `=b'[1,32])	
						*cons
							`=m'[`i',`=c'XB3] = `=m'[`i',`=c'XB3] + `=b'[1,33]
						*Exponent
							`=m'[`i',`=c'XB3] = exp(`=m'[`i',`=c'XB3])
					}
					
					if (`=nL1' >= 4) {	// Calculate e(XB4)
						*Age
							`=m'[`i',`=c'XB4] = `=m'[`i',`=c'XB4] + (`=m'[`i',cAge] * `=b'[1,34])
						*Age2
							`=m'[`i',`=c'XB4] = `=m'[`i',`=c'XB4] + (`=m'[`i',cAge]^2 * `=b'[1,35])								
						*Male
							`=m'[`i',`=c'XB4] = `=m'[`i',`=c'XB4] + (`=m'[`i',cMale] * `=b'[1,36])
						*ECOG 
							if (`=m'[`i',cECOG] == 1) `=m'[`i',`=c'XB4] = `=m'[`i',`=c'XB4] + `=b'[1,38] 
							if (`=m'[`i',cECOG] == 2) `=m'[`i',`=c'XB4] = `=m'[`i',`=c'XB4] + `=b'[1,39] 
						*RISS 
							if (`=m'[`i',cRISS] == 2) `=m'[`i',`=c'XB4] = `=m'[`i',`=c'XB4] + `=b'[1,41] 
							if (`=m'[`i',cRISS] == 3) `=m'[`i',`=c'XB4] = `=m'[`i',`=c'XB4] + `=b'[1,42]
						*SCT
							`=m'[`i',`=c'XB4] = `=m'[`i',`=c'XB4] + (`=m'[`i',cSCT] * `=b'[1,43])	
						*cons
							`=m'[`i',`=c'XB4] = `=m'[`i',`=c'XB4] + `=b'[1,44]
						*Exponent
							`=m'[`i',`=c'XB4] = exp(`=m'[`i',`=c'XB4])
					}
										
					*Calculate probabilities
						`=m'[`i',`=c'PR1] = `=m'[`i',`=c'XB1]/(`=m'[`i',`=c'XB1] + `=m'[`i',`=c'XB2] + `=m'[`i',`=c'XB3] + `=m'[`i',`=c'XB4])
						`=m'[`i',`=c'PR2] = `=m'[`i',`=c'PR1] + `=m'[`i',`=c'XB2]/(`=m'[`i',`=c'XB1] + `=m'[`i',`=c'XB2] + `=m'[`i',`=c'XB3] + `=m'[`i',`=c'XB4])
						`=m'[`i',`=c'PR3] = `=m'[`i',`=c'PR2] + `=m'[`i',`=c'XB3]/(`=m'[`i',`=c'XB1] + `=m'[`i',`=c'XB2] + `=m'[`i',`=c'XB3] + `=m'[`i',`=c'XB4])
						`=m'[`i',`=c'PR4] = `=m'[`i',`=c'PR3] + `=m'[`i',`=c'XB4]/(`=m'[`i',`=c'XB1] + `=m'[`i',`=c'XB2] + `=m'[`i',`=c'XB3] + `=m'[`i',`=c'XB4])
								
					*Compare to RN
						if (`=m'[`i',`=c'PR1] != . & `=m'[`i',`=c'RN] < `=m'[`i',`=c'PR1]) `=m'[`i',`=c'OC] = `=o'[1,1]
						if (`=m'[`i',`=c'PR2] != . & `=m'[`i',`=c'RN] > `=m'[`i',`=c'PR1]  & `=m'[`i',`=c'RN] < `=m'[`i',`=c'PR2]) `=m'[`i',`=c'OC] = `=o'[1,2]
						if (`=m'[`i',`=c'PR3] != . & `=m'[`i',`=c'RN] > `=m'[`i',`=c'PR2]  & `=m'[`i',`=c'RN] < `=m'[`i',`=c'PR3]) `=m'[`i',`=c'OC] = `=o'[1,3]
						if (`=m'[`i',`=c'PR4] != . & `=m'[`i',`=c'RN] > `=m'[`i',`=c'PR3]  & `=m'[`i',`=c'RN] < `=m'[`i',`=c'PR4]) `=m'[`i',`=c'OC] = `=o'[1,4]
						
					*Update outcome matrices
						mCore[`i',cCR] = `=m'[`i',`=c'OC]
						mCR[`i',`=LX'+1] = `=m'[`i',`=c'OC]
				}
			}
		}

	*VRd 
{
		if ("$Line" == "1" & "$Int" == "VRd") {	
			forvalues i = 1/`=Obs' { 
				mata {
					if (mMOR[`i',`=OMC'-1] == 0) {
						`=m'[`i',`=c'OC] = 31
						mCore[`i',cCR] = 31
						mCR[`i',`=LX'+1] = 31
					}	
				}
			}
		}
}
