**********
	*SIM SCT DN
**********

	*Patient matrix	
		mata: `=m' = mCore , mLO , mCom
		mata: r`=m' = rmCore
		mata: c`=m' = cmCore \ cmLO \ cmCom 
		mata: `=m'[.,`=c'RN] = runiform(`=Obs',1)	
		mata: _matrix_list(`=m', r`=m', c`=m')
		
	*Determine outcome
		forvalues i = 1/`=Obs' {
			mata {
				if (mState[`i',1] <= `=OMC'+1) { // State filter only (no-one dies before DN)	
				
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
					*Age70
						`=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + (`=m'[`i',15] * `=b'[1,10])
					*Age75
						`=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + (`=m'[`i',16] * `=b'[1,11])
					*CM
						if (`=m'[`i',17] == 1) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,13]
						if (`=m'[`i',17] == 2) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,14]
						if (`=m'[`i',17] == 3) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,15]
					*cons
						`=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,16]
		
				*Calculate probability
					`=m'[`i',`=c'PR] = 1 / (1 + exp(-`=m'[`i',`=c'XB]))	

				*Compare with RN
					if (`=m'[`i',`=c'PR] > `=m'[`i',`=c'RN]) `=m'[`i',`=c'OC] = 1
					else `=m'[`i',`=c'OC] = 0	
						
				*Update matrices
					mSCT[`i',1] = `=m'[`i',`=c'OC]
					mCore[`i',cSCT] = `=m'[`i',`=c'OC]
				}	
			}		
		}
