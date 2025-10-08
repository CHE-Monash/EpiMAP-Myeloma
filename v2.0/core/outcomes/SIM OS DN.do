**********
	*SIM OS DN
**********

	*Patient matrix	
		mata: `=m' = mCore , mSU
		mata: r`=m' = rmCore
		mata: c`=m' = cmCore \ cmSU
		mata: `=m'[.,`=c'RN] = runiform(`=Obs',1)
		mata: _matrix_list(`=m', r`=m', c`=m')
	
	*Determine outcome
		forvalues i = 1/`=Obs' {
			mata {
				if (mState[`i',1] <= `=OMC'+1) { // State filter only (no-one dies before DN)	
				
				*Calculate XB
					*Age
						`=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + (`=m'[`i',cAge] * `=b'[1,1])
					*Age^2
						`=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + (`=m'[`i',cAge]^2  * `=b'[1,2])		
					*Male 
						`=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + (`=m'[`i',cMale] * `=b'[1,3])
					*ECOG
						if (`=m'[`i',cECOG] == 1) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,5] 
						if (`=m'[`i',cECOG] == 2) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,6]
					*RISS
						if (`=m'[`i',cRISS] == 2) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,8] 
						if (`=m'[`i',cRISS] == 3) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,9]					
					*TX#BCR
						if (`=Line' == 0) 						`=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,(10 + `=m'[`i',cBCR] - 1)]
						if (`=Line' == 1 & `=m'[`i',cSCT] == 0)	`=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,(16 + `=m'[`i',cBCR] - 1)]
						if (`=Line' == 1 & `=m'[`i',cSCT] == 1)	`=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,(22 + `=m'[`i',cBCR] - 1)]
						if (`=Line' == 2)						`=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,(28 + `=m'[`i',cBCR] - 1)]
						if (`=Line' == 3)						`=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,(34 + `=m'[`i',cBCR] - 1)]
						if (`=Line' == 4)						`=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,(40 + `=m'[`i',cBCR] - 1)]
						if (`=Line' == 5)						`=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,(46 + `=m'[`i',cBCR] - 1)]
						if (`=Line' == 6)						`=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,(52 + `=m'[`i',cBCR] - 1)]		
					*cons
						`=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,58]
					
				*Calculate survival time
					`=m'[`i',`=c'OC] = calcSurvivalTime(`=m'[`i',`=c'XB], `=m'[`i',`=c'RN], f`=b', `=b'[1,cols(`=b')])
			
				*Update mOS
					mOS[`i',`=OMC'] = `=m'[`i',`=c'OC]
				}
			}
		}
