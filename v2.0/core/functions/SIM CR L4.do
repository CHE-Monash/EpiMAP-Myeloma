**********
	*SIM CR L4
**********

	*Patient matrix	
		mata: `=m' = mCore , mML
		mata: r`=m' = rmCore
		mata: c`=m' = cmCore \ cmML
		mata: `=m'[.,`=c'RN] = runiform(`=Obs',1)	
		mata: _matrix_list(`=m', r`=m', c`=m')
		
	*Create nL4 and update nL
		mata: st_matrix("stL4_CR", oL4_CR)
		scalar nL4 = colsof(stL4_CR)
		
	*Determine outcome
		forvalues i = 1/`=Obs' {
			mata {
				if (mMOR[`i',`=OMC'-1] == 0 & mState[`i',2] <= `=OMC') { // Alive & State filters		
					
					`=m'[`i',`=c'XB1] = 1 // Set e(XB1) to 1
					
					if (`=nL4' >= 2) { // Calculate e(XB2)
						*Age
							`=m'[`i',`=c'XB2] = `=m'[`i',`=c'XB2] + (`=m'[`i',cAge] * `=b'[1,14])
						*Age2
							`=m'[`i',`=c'XB2] = `=m'[`i',`=c'XB2] + (`=m'[`i',cAge]^2 * `=b'[1,15])							
						*Male
							`=m'[`i',`=c'XB2] = `=m'[`i',`=c'XB2] + (`=m'[`i',cMale] * `=b'[1,16])
						*ECOG
							if (`=m'[`i',cECOG] == 1) `=m'[`i',`=c'XB2] = `=m'[`i',`=c'XB2] + `=b'[1,18] 
							if (`=m'[`i',cECOG] == 2) `=m'[`i',`=c'XB2] = `=m'[`i',`=c'XB2] + `=b'[1,19] 
						*RISS
							if (`=m'[`i',cRISS] == 2) `=m'[`i',`=c'XB2] = `=m'[`i',`=c'XB2] + `=b'[1,21] 
							if (`=m'[`i',cRISS] == 3) `=m'[`i',`=c'XB2] = `=m'[`i',`=c'XB2] + `=b'[1,22]
						*BCR
							if (`=m'[`i',cBCR] == 3) `=m'[`i',`=c'XB2] = `=m'[`i',`=c'XB2] + `=b'[1,24]
							if (`=m'[`i',cBCR] == 5) `=m'[`i',`=c'XB2] = `=m'[`i',`=c'XB2] + `=b'[1,25]
						*cons
							`=m'[`i',`=c'XB2] = `=m'[`i',`=c'XB2] + `=b'[1,26]
						*Exponent
							`=m'[`i',`=c'XB2] = exp(`=m'[`i',`=c'XB2])
					}
					
					if (`=nL4' >= 3) { //Calculate e(XB3)
						*Age
							`=m'[`i',`=c'XB3] = `=m'[`i',`=c'XB3] + (`=m'[`i',cAge] * `=b'[1,27])
						*Age2
							`=m'[`i',`=c'XB3] = `=m'[`i',`=c'XB3] + (`=m'[`i',cAge]^2 * `=b'[1,28])							
						*Male
							`=m'[`i',`=c'XB3] = `=m'[`i',`=c'XB3] + (`=m'[`i',cMale] * `=b'[1,29])
						*ECOG
							if (`=m'[`i',cECOG] == 1) `=m'[`i',`=c'XB3] = `=m'[`i',`=c'XB3] + `=b'[1,31]
							if (`=m'[`i',cECOG] == 2) `=m'[`i',`=c'XB3] = `=m'[`i',`=c'XB3] + `=b'[1,32]
						*RISS 
							if (`=m'[`i',cRISS] == 2) `=m'[`i',`=c'XB3] = `=m'[`i',`=c'XB3] + `=b'[1,34]
							if (`=m'[`i',cRISS] == 3) `=m'[`i',`=c'XB3] = `=m'[`i',`=c'XB3] + `=b'[1,35]
						*BCR
							if (`=m'[`i',cBCR] == 3) `=m'[`i',`=c'XB3] = `=m'[`i',`=c'XB3] + `=b'[1,37]
							if (`=m'[`i',cBCR] == 5) `=m'[`i',`=c'XB3] = `=m'[`i',`=c'XB3] + `=b'[1,38]
						*cons
							`=m'[`i',`=c'XB3] = `=m'[`i',`=c'XB3] + `=b'[1,39]
						*Exponent
							`=m'[`i',`=c'XB3] = exp(`=m'[`i',`=c'XB3])
					}
					
					if (`=nL4' >= 4) {	// Calculate e(XB4)
						*Age
							`=m'[`i',`=c'XB4] = `=m'[`i',`=c'XB4] + (`=m'[`i',cAge] * `=b'[1,40])
						*Age2
							`=m'[`i',`=c'XB4] = `=m'[`i',`=c'XB4] + (`=m'[`i',cAge]^2 * `=b'[1,41])							
						*Male
							`=m'[`i',`=c'XB4] = `=m'[`i',`=c'XB4] + (`=m'[`i',cMale] * `=b'[1,42])
						*ECOG
							if (`=m'[`i',cECOG] == 1) `=m'[`i',`=c'XB4] = `=m'[`i',`=c'XB4] + `=b'[1,44] 
							if (`=m'[`i',cECOG] == 2) `=m'[`i',`=c'XB4] = `=m'[`i',`=c'XB4] + `=b'[1,45]
						*RISS
							if (`=m'[`i',cRISS] == 2) `=m'[`i',`=c'XB4] = `=m'[`i',`=c'XB4] + `=b'[1,47] 
							if (`=m'[`i',cRISS] == 3) `=m'[`i',`=c'XB4] = `=m'[`i',`=c'XB4] + `=b'[1,48]
						*BCR
							if (`=m'[`i',cBCR] == 3) `=m'[`i',`=c'XB4] = `=m'[`i',`=c'XB4] + `=b'[1,50]
							if (`=m'[`i',cBCR] == 5) `=m'[`i',`=c'XB4] = `=m'[`i',`=c'XB4] + `=b'[1,51]
						*cons
							`=m'[`i',`=c'XB4] = `=m'[`i',`=c'XB4] + `=b'[1,52]
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
						mCR[`i',`=LX'+2] = `=m'[`i',`=c'OC]
				}
			}
		}
