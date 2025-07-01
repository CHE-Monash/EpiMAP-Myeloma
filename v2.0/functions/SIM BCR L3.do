**********
	*SIM BCR L3
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
				if (mMOR[`i',`=OMC'-1] == 0 & mState[`i',2] <= `=OMC') { // Alive & State filters
				
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
					*pBCR
						if (`=m'[`i',cBCR] == 2) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,11 + cols(`=o')] 
						if (`=m'[`i',cBCR] == 3) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,12 + cols(`=o')] 
						if (`=m'[`i',cBCR] == 4) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,13 + cols(`=o')]
						if (`=m'[`i',cBCR] == 5) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,14 + cols(`=o')] 
						if (`=m'[`i',cBCR] == 6) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,15 + cols(`=o')] 
				
				*Probabilities
					`=m'[`i',`=c'PR1] = 1/(1+exp(`=m'[`i',`=c'XB] - `=b'[1,cols(`=b')]-1))
					`=m'[`i',`=c'PR2] = 1/(1+exp(`=m'[`i',`=c'XB] - `=b'[1,cols(`=b')]))
					if (`=m'[`i',`=c'PR2] != .) `=m'[`i',`=c'PR3] = 1
				
				*BCRc outcome
					if ((`=m'[`i',`=c'RN] < `=m'[`i',`=c'PR1]) & (`=m'[`i',`=c'PR3] == 1)) `=m'[`i',`=c'OC] = 1
					if ((`=m'[`i',`=c'RN] > `=m'[`i',`=c'PR1]) & (`=m'[`i',`=c'RN] < `=m'[`i',`=c'PR2])) `=m'[`i',`=c'OC] = 3
					if (`=m'[`i',`=c'RN] > `=m'[`i',`=c'PR2]) `=m'[`i',`=c'OC] = 5
				}
	
			*Update outcome matrices	
				if (mState[`i',2] > `=OMC') `=m'[`i',`=c'OC] = mBCR[`i',`=LX'+1] // Grab prevalent patient data
				mCore[`i',cBCR] = `=m'[`i',`=c'OC]					
				mBCR[`i',`=LX'+1] = `=m'[`i',`=c'OC]
			}
		}
