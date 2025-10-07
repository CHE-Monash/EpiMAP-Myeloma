**********
	*SIM CD LX (5 - 9)
**********

	*Patient Matrix
		mata: `=m' = mCore , mSU
		mata: r`=m' = rmCore
		mata: c`=m' = cmCore \ cmSU
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
					*BCRc
						if (`=m'[`i',cBCR] == 3) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,11] 
						if (`=m'[`i',cBCR] == 5) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,12] 
					*cons
						`=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,13] 
		
				*Calculate survival time
					if (f`=b' == "ereg") 		`=m'[`i',`=c'OC] = (ln(`=m'[`i',`=c'RN]))/-exp(`=m'[`i',`=c'XB])
					if (f`=b' == "weibull")  	`=m'[`i',`=c'OC] = ((ln(`=m'[`i',`=c'RN])):/-exp(`=m'[`i',`=c'XB])):^(1:/exp(`=b'[1,cols(`=b')]))
					if (f`=b' == "gompertz")	`=m'[`i',`=c'OC] = (ln(1-((`=b'[1,cols(`=b')]:*(ln(`=m'[`i',`=c'RN])))/exp(`=m'[`i',`=c'XB]))))/`=b'[1,cols(`=b')]
					
				*Curtail if outcome beyond last observed in the data
					if (`=m'[`i',`=c'OC] != . & `=m'[`i',`=c'OC] > maxLX_CD)	`=m'[`i',`=c'OC] = maxLX_CD
				}
					
			*Grab prevalent patient data
				else if (mState[`i',1] > `=OMC'+1) `=m'[`i',`=c'OC] = mTNE[`i',`=OMC'] * 365.25
				
			*Update outcome matrices
				mTNE[`i',`=OMC'] = `=m'[`i',`=c'OC] / 365.25
				mTSD[`i',`=OMC'+1] = mTSD[`i',`=OMC'] + mTNE[`i',`=OMC']
				mCore[`i',cCD] = `=m'[`i',`=c'OC]
				mCD[`i',`=LX'+1] = `=m'[`i',`=c'OC]	
			}
		}		
