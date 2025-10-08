**********
	*SIM BCR L2
**********

	*Patient matrix
		mata: `=m' = mCore , mOL
		mata: r`=m' = rmCore
		mata: c`=m' = cmCore \ cmOL
		mata: `=m'[.,`=c'RN] = runiform(`=Obs',1)					
		mata: _matrix_list(`=m', r`=m', c`=m')
			
	*Determine outcome
		forvalues i = 1/`=Obs' {
			mata {
				if (mMOR[`i',`=OMC'-1] == 0 & mState[`i',1] <= `=OMC'+1) { // Alive & State filters
				
				*Calculate XB
					*Age
						`=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + (`=m'[`i',cAge] * `=b'[1,1])
					*Age2
						`=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + (`=m'[`i',cAge]^2 * `=b'[1,2])						
					*Male 
						`=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + (`=m'[`i',cMale] * `=b'[1,3])
					*ECOG
						if (`=m'[`i',cECOG] == 1) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,5] 
						if (`=m'[`i',cECOG] == 2) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,6]  
					*RISS 
						if (`=m'[`i',cRISS] == 2) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,8] 
						if (`=m'[`i',cRISS] == 3) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,9]
					*CR
						if (cols(`=o') >= 2) {
							if (`=m'[`i',cCR] == `=o'[1,2]) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,11]
						}
						if (cols(`=o') >= 3) {
							if (`=m'[`i',cCR] == `=o'[1,3]) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,12]
						}
						if (cols(`=o') >= 4) {
							if (`=m'[`i',cCR] == `=o'[1,4]) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,13]
						}
					*BCR_L1 - mBCR column 2
						if (mBCR[`i',1] == 2) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,11 + cols(`=o')] 
						if (mBCR[`i',1] == 3) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,12 + cols(`=o')] 
						if (mBCR[`i',1] == 4) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,13 + cols(`=o')]
						if (mBCR[`i',1] == 5) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,14 + cols(`=o')] 
						if (mBCR[`i',1] == 6) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,15 + cols(`=o')] 
					*BCR_SCT - mBCR column 10
						if (`=m'[`i',cSCT] == 1) {
							if (mBCR[`i',10] == 1) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,17 + cols(`=o')] 
							if (mBCR[`i',10] == 2) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,18 + cols(`=o')] 
							if (mBCR[`i',10] == 3) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,19 + cols(`=o')] 
							if (mBCR[`i',10] == 4) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,20 + cols(`=o')]
						}
				
				*Extract cut points	
					cutPoints = `=b'[1, (cols(`=b')-4, cols(`=b')-3, cols(`=b')-2, cols(`=b')-1, cols(`=b'))]
				
				*Calculate probabilities
					probMatrix = calcOrdLogitProbs(`=m'[`i', `=c'XB], cutPoints)
					
				*Assign BCR outcome 
					`=m'[`i',`=c'OC] = assignOrdOutcome(`=m'[`i',`=c'RN], probMatrix, (1, 2, 3, 4, 5, 6))[1,1]
				}	
			}			
		}
/*					
		*Pre-market analysis - sort by XB, overwrite RN
			mata {
				if ("$Analysis" == "DVd-Pre"){
					_sort(`=m', `=c'XB) 
					N = colnonmissing(`=m')[1,`=c'OC] 
					for (n = 1; n <= N; n++) {
					`=m'[n,`=c'RN] = n / N
					}
				}
			}

			*DVd
				forvalues i = 1/`=Obs' {
					mata {	
						if ("$Analysis" == "DVd-Pre" & "$Int" == "DVd") {
							if (`=m'[`i',`=c'RN] < 0.186) `=m'[`i',`=c'OC] = 1
							if ((`=m'[`i',`=c'RN] > 0.186) & (`=m'[`i',`=c'RN] < 0.530)) `=m'[`i',`=c'OC] = 2
							if ((`=m'[`i',`=c'RN] > 0.530) & (`=m'[`i',`=c'RN] < 0.821)) `=m'[`i',`=c'OC] = 3
							if ((`=m'[`i',`=c'RN] > 0.821) & (`=m'[`i',`=c'RN] < 0.866)) `=m'[`i',`=c'OC] = 4
							if ((`=m'[`i',`=c'RN] > 0.866) & (`=m'[`i',`=c'RN] < 0.970)) `=m'[`i',`=c'OC] = 5
							if (`=m'[`i',`=c'RN] > 0.970 & `=m'[`i',`=c'RN] != .) `=m'[`i',`=c'OC] = 6
						}
					}
				}
			
			mata: _sort(`=m', 1) // sort by ID
*/
			forvalues i = 1/`=Obs' {
				mata {	
				*Grab prevalent patient data	
					if (mState[`i',1] > `=OMC'+1) `=m'[`i',`=c'OC] = mBCR[`i',`=LX']
				*Update outcome matrices	
					if (mMOR[`i',`=OMC'-1] == 0) mCore[`i',cBCR] = `=m'[`i',`=c'OC]					
					if (mMOR[`i',`=OMC'-1] == 0) mBCR[`i',`=LX'] = `=m'[`i',`=c'OC]
				}
			}
