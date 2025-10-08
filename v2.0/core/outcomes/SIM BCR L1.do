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
				
				*Probabilities
					`=m'[`i',`=c'PR1] = 1/(1+exp(`=m'[`i',`=c'XB] - `=b'[1,cols(`=b')-4]))
					`=m'[`i',`=c'PR2] = 1/(1+exp(`=m'[`i',`=c'XB] - `=b'[1,cols(`=b')-3]))
					`=m'[`i',`=c'PR3] = 1/(1+exp(`=m'[`i',`=c'XB] - `=b'[1,cols(`=b')-2]))
					`=m'[`i',`=c'PR4] = 1/(1+exp(`=m'[`i',`=c'XB] - `=b'[1,cols(`=b')-1]))
					`=m'[`i',`=c'PR5] = 1/(1+exp(`=m'[`i',`=c'XB] - `=b'[1,cols(`=b')]))
					if (`=m'[`i',`=c'PR5] != .) `=m'[`i',`=c'PR6] = 1
				
				*BCR outcome
					if ((`=m'[`i',`=c'RN] < `=m'[`i',`=c'PR1]) & (`=m'[`i',`=c'PR6] == 1)) `=m'[`i',`=c'OC] = 1
					if ((`=m'[`i',`=c'RN] > `=m'[`i',`=c'PR1]) & (`=m'[`i',`=c'RN] < `=m'[`i',`=c'PR2])) `=m'[`i',`=c'OC] = 2
					if ((`=m'[`i',`=c'RN] > `=m'[`i',`=c'PR2]) & (`=m'[`i',`=c'RN] < `=m'[`i',`=c'PR3])) `=m'[`i',`=c'OC] = 3
					if ((`=m'[`i',`=c'RN] > `=m'[`i',`=c'PR3]) & (`=m'[`i',`=c'RN] < `=m'[`i',`=c'PR4])) `=m'[`i',`=c'OC] = 4
					if ((`=m'[`i',`=c'RN] > `=m'[`i',`=c'PR4]) & (`=m'[`i',`=c'RN] < `=m'[`i',`=c'PR5])) `=m'[`i',`=c'OC] = 5
					if (`=m'[`i',`=c'RN] > `=m'[`i',`=c'PR5]) `=m'[`i',`=c'OC] = 6
				}	
				
			*Grab prevalent patient data
				if (mState[`i',1] > `=OMC'+1) `=m'[`i',`=c'OC] = mBCR[`i',`=LX'] 
					
			*Update matrices
				mCore[`i',cBCR] = `=m'[`i',`=c'OC]					
				mBCR[`i',`=LX'] = `=m'[`i',`=c'OC]
			}			
		}
							
