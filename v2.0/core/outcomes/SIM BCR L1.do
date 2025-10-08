**********
	*SIM BCR L1
**********

	*Patient matrix
		mata: `=m' = mCore , mOL
		mata: r`=m' = rmCore
		mata: c`=m' = cmCore \ cmOL
		mata: `=m'[.,`=c'RN] = runiform(`=Obs',1)					
		mata: _matrix_list(`=m', r`=m', c`=m')
		
	*VRd-Post - Add VRd to the list of regimens
		if "$Coeffs" == "VRd" {
			mata: `=o' = `=o' , 31
		}
			
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
					*SCT
						`=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + (`=m'[`i',cSCT] * `=b'[1,10])
					*CR
						if (cols(`=o') >= 2) {
							if	(`=m'[`i',cCR] == `=o'[1,2]) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,12]
						}
						if (cols(`=o') >= 3) {
							if	(`=m'[`i',cCR] == `=o'[1,3]) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,13]
						}
						if (cols(`=o') >= 4) {
							if	(`=m'[`i',cCR] == `=o'[1,4]) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,14]
						}
				
				*Extract cut points	
					cutPoints = `=b'[1, (cols(`=b')-4, cols(`=b')-3, cols(`=b')-2, cols(`=b')-1, cols(`=b'))]
				
				*Calculate probabilities
					probMatrix = calcOrdLogitProbs(`=m'[`i', `=c'XB], cutPoints)
					
				*Assign BCR outcome 
					`=m'[`i',`=c'OC] = assignOrdOutcome(`=m'[`i',`=c'RN], probMatrix, (1, 2, 3, 4, 5, 6))[1,1]
				}	
				
			*Grab prevalent patient data
				if (mState[`i',1] > `=OMC'+1) `=m'[`i',`=c'OC] = mBCR[`i',`=LX'] 
					
			*Update matrices
				mCore[`i',cBCR] = `=m'[`i',`=c'OC]					
				mBCR[`i',`=LX'] = `=m'[`i',`=c'OC]
			}			
		}
							
