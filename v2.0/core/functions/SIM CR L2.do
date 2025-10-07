**********
	*SIM CR L2
**********

	*Patient matrix	
		mata: `=m' = mCore , mML
		mata: r`=m' = rmCore
		mata: c`=m' = cmCore \ cmML
		mata: `=m'[.,`=c'RN] = runiform(`=Obs',1)	
		mata: _matrix_list(`=m', r`=m', c`=m')
		
	*Create nL2
		mata: st_matrix("stL2_CR", oL2_CR)
		scalar nL2 = colsof(stL2_CR)
		
	*Determine outcome
		forvalues i = 1/`=Obs' {
			mata {
				if (mMOR[`i',`=OMC'-1] == 0 & mState[`i',1] <= `=OMC'+1) { // Alive & State filters
				
					`=m'[`i',`=c'XB1] = 1 // Set e(XB1) to 1
					
					if (`=nL2' >= 2) { // Calculate e(XB2)
						*Age
							`=m'[`i',`=c'XB2] = `=m'[`i',`=c'XB2] + (`=m'[`i',cAge] * `=b'[1,17])
						*Age2
							`=m'[`i',`=c'XB2] = `=m'[`i',`=c'XB2] + (`=m'[`i',cAge]^2 * `=b'[1,18])							
						*Male
							`=m'[`i',`=c'XB2] = `=m'[`i',`=c'XB2] + (`=m'[`i',cMale] * `=b'[1,19])
						*ECOG
							if (`=m'[`i',cECOG] == 1) `=m'[`i',`=c'XB2] = `=m'[`i',`=c'XB2] + `=b'[1,21] 
							if (`=m'[`i',cECOG] == 2) `=m'[`i',`=c'XB2] = `=m'[`i',`=c'XB2] + `=b'[1,22] 
						*RISS
							if (`=m'[`i',cRISS] == 2) `=m'[`i',`=c'XB2] = `=m'[`i',`=c'XB2] + `=b'[1,24] 
							if (`=m'[`i',cRISS] == 3) `=m'[`i',`=c'XB2] = `=m'[`i',`=c'XB2] + `=b'[1,25]
						*BCR
							if (`=m'[`i',cBCR] == 2) `=m'[`i',`=c'XB2] = `=m'[`i',`=c'XB2] + `=b'[1,27] 
							if (`=m'[`i',cBCR] == 3) `=m'[`i',`=c'XB2] = `=m'[`i',`=c'XB2] + `=b'[1,28]
							if (`=m'[`i',cBCR] == 4) `=m'[`i',`=c'XB2] = `=m'[`i',`=c'XB2] + `=b'[1,29] 
							if (`=m'[`i',cBCR] == 5) `=m'[`i',`=c'XB2] = `=m'[`i',`=c'XB2] + `=b'[1,30]
							if (`=m'[`i',cBCR] == 6) `=m'[`i',`=c'XB2] = `=m'[`i',`=c'XB2] + `=b'[1,31]
						*cons
							`=m'[`i',`=c'XB2] = `=m'[`i',`=c'XB2] + `=b'[1,32]
						*Exponent
							`=m'[`i',`=c'XB2] = exp(`=m'[`i',`=c'XB2])
					}
					
					if (`=nL2' >= 3) { //Calculate e(XB3)
						*Age
							`=m'[`i',`=c'XB3] = `=m'[`i',`=c'XB3] + (`=m'[`i',cAge] * `=b'[1,33])
						*Age2
							`=m'[`i',`=c'XB3] = `=m'[`i',`=c'XB3] + (`=m'[`i',cAge]^2 * `=b'[1,34])							
						*Male
							`=m'[`i',`=c'XB3] = `=m'[`i',`=c'XB3] + (`=m'[`i',cMale] * `=b'[1,35])
						*ECOG
							if (`=m'[`i',cECOG] == 1) `=m'[`i',`=c'XB3] = `=m'[`i',`=c'XB3] + `=b'[1,37]
							if (`=m'[`i',cECOG] == 2) `=m'[`i',`=c'XB3] = `=m'[`i',`=c'XB3] + `=b'[1,38]
						*RISS 
							if (`=m'[`i',cRISS] == 2) `=m'[`i',`=c'XB3] = `=m'[`i',`=c'XB3] + `=b'[1,40]
							if (`=m'[`i',cRISS] == 3) `=m'[`i',`=c'XB3] = `=m'[`i',`=c'XB3] + `=b'[1,41]
						*BCR
							if (`=m'[`i',cBCR] == 2) `=m'[`i',`=c'XB3] = `=m'[`i',`=c'XB3] + `=b'[1,43]
							if (`=m'[`i',cBCR] == 3) `=m'[`i',`=c'XB3] = `=m'[`i',`=c'XB3] + `=b'[1,44]
							if (`=m'[`i',cBCR] == 4) `=m'[`i',`=c'XB3] = `=m'[`i',`=c'XB3] + `=b'[1,45] 
							if (`=m'[`i',cBCR] == 5) `=m'[`i',`=c'XB3] = `=m'[`i',`=c'XB3] + `=b'[1,46]
							if (`=m'[`i',cBCR] == 6) `=m'[`i',`=c'XB3] = `=m'[`i',`=c'XB3] + `=b'[1,47]
						*cons
							`=m'[`i',`=c'XB3] = `=m'[`i',`=c'XB3] + `=b'[1,48]
						*Exponent
							`=m'[`i',`=c'XB3] = exp(`=m'[`i',`=c'XB3])
					}
					
					if (`=nL2' >= 4) {	// Calculate e(XB4)
						*Age
							`=m'[`i',`=c'XB4] = `=m'[`i',`=c'XB4] + (`=m'[`i',cAge] * `=b'[1,49])
						*Age2
							`=m'[`i',`=c'XB4] = `=m'[`i',`=c'XB4] + (`=m'[`i',cAge]^2 * `=b'[1,50])							
						*Male
							`=m'[`i',`=c'XB4] = `=m'[`i',`=c'XB4] + (`=m'[`i',cMale] * `=b'[1,51])
						*ECOG
							if (`=m'[`i',cECOG] == 1) `=m'[`i',`=c'XB4] = `=m'[`i',`=c'XB4] + `=b'[1,53] 
							if (`=m'[`i',cECOG] == 2) `=m'[`i',`=c'XB4] = `=m'[`i',`=c'XB4] + `=b'[1,54]
						*RISS
							if (`=m'[`i',cRISS] == 2) `=m'[`i',`=c'XB4] = `=m'[`i',`=c'XB4] + `=b'[1,56] 
							if (`=m'[`i',cRISS] == 3) `=m'[`i',`=c'XB4] = `=m'[`i',`=c'XB4] + `=b'[1,57]
						*BCR
							if (`=m'[`i',cBCR] == 2) `=m'[`i',`=c'XB4] = `=m'[`i',`=c'XB4] + `=b'[1,59]
							if (`=m'[`i',cBCR] == 3) `=m'[`i',`=c'XB4] = `=m'[`i',`=c'XB4] + `=b'[1,60]
							if (`=m'[`i',cBCR] == 4) `=m'[`i',`=c'XB4] = `=m'[`i',`=c'XB4] + `=b'[1,61] 
							if (`=m'[`i',cBCR] == 5) `=m'[`i',`=c'XB4] = `=m'[`i',`=c'XB4] + `=b'[1,62]
							if (`=m'[`i',cBCR] == 6) `=m'[`i',`=c'XB4] = `=m'[`i',`=c'XB4] + `=b'[1,63]
						*cons
							`=m'[`i',`=c'XB4] = `=m'[`i',`=c'XB4] + `=b'[1,64]
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
/*	
	*DVd analysis
{
		*Pre
			if ("$Analysis" == "DVd-Pre" & "$Int" == "DVd" & "$Coeff" == "Pre") { 
				forvalues i = 1/`=Obs' { 
					mata { 
						if (mMOR[`i',`=OMC'-1] == 0) {
							`=m'[`i',`=c'OC] = 80
							mCore[`i',cCR] = 80
							mCR[`i',`=LX'+2] = 80
						}
					}
				}
			}
	
		*Post - includes $Data as a filter as there is an initial prediction using Population data
			if ("$Analysis" == "DVd-Post" & "$Int" == "DVd" & "$Data" == "Predicted"  & "$Coeff" == "DVd") {
				forvalues i = 1/`=Obs' { 
					mata {
						if (mMOR[`i',`=OMC'-1] == 0) {
							`=m'[`i',`=c'OC] = 80
							mCore[`i',cCR] = 80
							mCR[`i',`=LX'+2] = 80
						}	
					}
				}
			}
}		
