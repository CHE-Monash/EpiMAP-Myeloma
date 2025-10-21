**********
	*SIM CI DN
**********

	*Patient matrix	
		mata: `=m' = mCore , mSU
		mata: r`=m' = rmCore
		mata: c`=m' = cmCore \ cmSU
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
				*SCT
					mata: `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,9]	
				*cons
					mata: `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,10]
					
			*Calculate survival time
				mata: `=m'[`i',`=c'RN] = runiform(1, 1, 0, 1)	
				mata {
					if 		(f`=b' == "exponential")	`=m'[`i',`=c'OC] = (ln(`=m'[`i',`=c'RN])):/-exp(`=m'[`i',`=c'XB])
					else if (f`=b' == "weibull") 		`=m'[`i',`=c'OC] = ((ln(`=m'[`i',`=c'RN])):/-exp(`=m'[`i',`=c'XB])):^(1:/exp(`=b'[1,cols(`=b')]))
					else if (f`=b' == "gompertz")		`=m'[`i',`=c'OC] = (ln(1:-((`=b'[1,cols(`=b')]:*(ln(`=m'[`i',`=c'RN]))):/exp(`=m'[`i',`=c'XB])))):/`=b'[1,cols(`=b')]
				}
		
			*Update mTNE & mTSD
				mata: mTNE[`i',`=OMC'] = `=m'[`i',`=c'OC]/365.25
				mata: mTSD[`i',`=OMC'+1] = mTNE[`i',`=OMC']
		}
