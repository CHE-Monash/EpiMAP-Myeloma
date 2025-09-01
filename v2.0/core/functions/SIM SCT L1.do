**********
	*SIM SCT L1
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
					*BCR
						if (`=m'[`i',cBCR] == 2) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,11]
						if (`=m'[`i',cBCR] == 3) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,12] 
						if (`=m'[`i',cBCR] == 4) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,13] 
						if (`=m'[`i',cBCR] == 5) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,14]  
					*Age70
						`=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + (`=m'[`i',15] * `=b'[1,15])
					*Age75
						`=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + (`=m'[`i',16] * `=b'[1,16])
					*CM
						if (`=m'[`i',17] == 1) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,18]
						if (`=m'[`i',17] == 2) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,19]
						if (`=m'[`i',17] == 3) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,20]
					*cons
						`=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,21]
		
				*Calculate probability
					`=m'[`i',`=c'PR] = 1 / (1 + exp(-`=m'[`i',`=c'XB]))	

				*Compare with RN
					if (`=m'[`i',`=c'PR] > `=m'[`i',`=c'RN] & `=m'[`i',cCR] != 7 & `=m'[`i',cBCR] != 6) `=m'[`i',`=c'OC] = 1
					else `=m'[`i',`=c'OC] = 0					
					
				*Update matrices
					mSCT[`i',2] = `=m'[`i',`=c'OC]
					mCore[`i',cSCT] = `=m'[`i',`=c'OC]
				}	
			}
		}
