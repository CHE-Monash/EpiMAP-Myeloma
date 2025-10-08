**********
	*SIM CI DN
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
					*cons
						`=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,11]
						
				*Calculate survival time
					`=m'[`i',`=c'OC] = calcSurvivalTime(`=m'[`i',`=c'XB], `=m'[`i',`=c'RN], f`=b', `=b'[1,cols(`=b')])
				}

			*Grab prevalent patient data	
				else if (mState[`i',1] > `=OMC'+1) `=m'[`i',`=c'OC] = mTNE[`i',`=OMC'] * 365.25
				
			*Update outcome matrices
				mTNE[`i',`=OMC'] = `=m'[`i',`=c'OC] / 365.25
				mTSD[`i',`=OMC'+1] = mTSD[`i',`=OMC'] + mTNE[`i',`=OMC']
				mCI[`i',`=LX'+1] = `=m'[`i',`=c'OC]
			}
		}
