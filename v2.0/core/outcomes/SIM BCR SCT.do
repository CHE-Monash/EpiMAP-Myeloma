**********
	*SIM BCR SCT
**********
		
	*Patient Matrix
		mata: `=m' = mCore , mOL
		mata: r`=m' = rmCore
		mata: c`=m' = cmCore \ cmOL
		mata: `=m'[.,`=c'RN] = runiform(`=Obs',1)			
		mata: _matrix_list(`=m', r`=m', c`=m')
		
	*Determine outcome
		forvalues i = 1/`=Obs' {
			mata {
				if (mMOR[`i',`=OMC'-1] == 0 & mState[`i',1] <= `=OMC'+1) { // Alive & State filters
					if (`=m'[`i',cSCT] == 1) { // SCT patients only
					
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
						*pBCR
							if (`=m'[`i',cBCR] == 2) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,11] 
							if (`=m'[`i',cBCR] == 3) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,12] 
							if (`=m'[`i',cBCR] == 4) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,13]
							if (`=m'[`i',cBCR] == 5) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,14] 
						
					*Probabilities	
						`=m'[`i',`=c'PR1] = 1/(1+exp(`=m'[`i',`=c'XB] - `=b'[1,cols(`=b')-2]))
						`=m'[`i',`=c'PR2] = 1/(1+exp(`=m'[`i',`=c'XB] - `=b'[1,cols(`=b')-1]))	
						`=m'[`i',`=c'PR3] = 1/(1+exp(`=m'[`i',`=c'XB] - `=b'[1,cols(`=b')]))
						`=m'[`i',`=c'PR4] = 1
					
					*BCR outcome
						if (`=m'[`i',`=c'RN] < `=m'[`i',`=c'PR1]) 											`=m'[`i',`=c'OC] = 1
						if (`=m'[`i',`=c'RN] > `=m'[`i',`=c'PR1] & `=m'[`i',`=c'RN] < `=m'[`i',`=c'PR2]) 	`=m'[`i',`=c'OC] = 2
						if (`=m'[`i',`=c'RN] > `=m'[`i',`=c'PR2] & `=m'[`i',`=c'RN] < `=m'[`i',`=c'PR3])	`=m'[`i',`=c'OC] = 3
						if (`=m'[`i',`=c'RN] > `=m'[`i',`=c'PR3])											`=m'[`i',`=c'OC] = 4
					}
				}
				
			*Grab prevalent patient data	
				else if (mState[`i',1] > `=OMC'+1) `=m'[`i',`=c'OC] = mBCR[`i',10] // BCR_SCT is column 10 of mBCR
				
			*Update outcome matrices		
				if (`=m'[`i',cSCT] == 0) `=m'[`i',`=c'OC] = 0 // Set BCR to 0 for No SCT
				if (mMOR[`i',`=OMC'-1] == 0) mCore[`i',cBCR] = `=m'[`i',`=c'OC]					
				if (mMOR[`i',`=OMC'-1] == 0) mBCR[`i',10] = `=m'[`i',`=c'OC] // BCR_SCT is column 10 of mBCR
			}
		}
