**********
	*SIM CI LXX
**********

	*Patient Matrix	
		mata: `=m' = mCore , mSU
		mata: r`=m' = rmCore
		mata: c`=m' = cmCore \ cmSU
		mata: `=m'[.,`=c'RN] = runiform(`=Obs', 1, 0, 1)
		mata: _matrix_list(`=m', r`=m', c`=m')

	*Determine outcome
		forvalues i = 1/`=Obs'{
			*Calculate xb
				*Age
					mata: `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + (`=m'[`i',cAge] * `=b'[1,1])
				*Male 
					mata: `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + (`=m'[`i',cMale] * `=b'[1,2])
				*ECOGcc
					mata {
						if 		(`=m'[`i',cECOGc] == 1) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,4] 
						else if (`=m'[`i',cECOGc] == 2) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,5] 
					}
				*ISS 
					mata {
						if 		(`=m'[`i',cISS] == 2) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,7] 
						else if (`=m'[`i',cISS] == 3) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,8]
					}
				*BCRc
					mata {
						if		(`=m'[`i',cBCR] == 3) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,10] 
						else if	(`=m'[`i',cBCR] == 5) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,11] 
					}
				*cons
					mata: `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,12] 			
			
			*Calculate survival time
				mata {
					if 		(f`=b' == "exponential") `=m'[`i',`=c'OC] = (ln(`=m'[`i',`=c'RN])):/-exp(`=m'[`i',`=c'XB])
					else if (f`=b' == "weibull")  	 `=m'[`i',`=c'OC] = ((ln(`=m'[`i',`=c'RN])):/-exp(`=m'[`i',`=c'XB])):^(1:/exp(`=b'[1,cols(`=b')]))
					else if (f`=b' == "gompertz")	 `=m'[`i',`=c'OC] = (ln(1:-((`=b'[1,cols(`=b')]:*(ln(`=m'[`i',`=c'RN]))):/exp(`=m'[`i',`=c'XB])))):/`=b'[1,cols(`=b')]
				}
							
			*Curtail if beyond last observed in the data
				mata {
					if	(`=m'[`i',`=c'OC] > maxLX_CI)	`=m'[`i',`=c'OC] = maxLX_CI
				}
				
			*Update outcome matrices (for MOR == 0 only)
				mata {
					if 	(mMOR[`i',`=OMC'-1] == 0) mTNE[`i',`=OMC'] = `=m'[`i',`=c'OC]/365.25
					if 	(mMOR[`i',`=OMC'-1] == 0) mTSD[`i',`=OMC'+1] = mTSD[`i',`=OMC'] + mTNE[`i',`=OMC']
					if 	(mMOR[`i',`=OMC'-1] == 0) mCI[`i',`=LX'+1] = `=m'[`i',`=c'OC]
				}
		}
		
