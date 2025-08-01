**********
	*SIM OS DN
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
				*Age^2
					mata: `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + (`=m'[`i',cAge]^2  * `=b'[1,2])		
				*Male 
					mata: `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + (`=m'[`i',cMale] * `=b'[1,3])
				*ECOGcc 
					mata {
						if 		(`=m'[`i',cECOGc] == 1) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,5] 
						else if (`=m'[`i',cECOGc] == 2) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,6]
					}
				*ISS
					mata {
						if 		(`=m'[`i',cISS] == 2) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,8] 
						else if (`=m'[`i',cISS] == 3) `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,9]
					}					
				*TX#BCR
					mata {
						if		(`=Line' == 0) 							`=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,(10 + `=m'[`i',cBCR] - 1)]
						else if (`=Line' == 1 & `=m'[`i',cSCT] == 0)	`=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,(16 + `=m'[`i',cBCR] - 1)]
						else if (`=Line' == 1 & `=m'[`i',cSCT] == 1)	`=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,(22 + `=m'[`i',cBCR] - 1)]
						else if (`=Line' == 2)							`=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,(28 + `=m'[`i',cBCR] - 1)]
						else if (`=Line' == 3)							`=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,(34 + `=m'[`i',cBCR] - 1)]
						else if (`=Line' == 4)							`=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,(40 + `=m'[`i',cBCR] - 1)]
						else if (`=Line' == 5)							`=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,(46 + `=m'[`i',cBCR] - 1)]
						else if (`=Line' == 6)							`=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,(52 + `=m'[`i',cBCR] - 1)]
					}		
				*cons
					mata: `=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,58]
			
			*Draw RN
				mata: `=m'[.,`=c'RN] = runiform(`=Obs',1)
				
			*Calculate survival time
				mata {
					if 		(f`=b' == "ereg") 		`=m'[`i',`=c'OC] = (ln(`=m'[`i',`=c'RN])):/-exp(`=m'[`i',`=c'XB])
					else if (f`=b' == "weibull")  	`=m'[`i',`=c'OC] = ((ln(`=m'[`i',`=c'RN])):/-exp(`=m'[`i',`=c'XB])):^(1:/exp(`=b'[1,cols(`=b')]))
					else if (f`=b' == "gompertz")	`=m'[`i',`=c'OC] = (ln(1:-((`=b'[1,cols(`=b')]:*(ln(`=m'[`i',`=c'RN]))):/exp(`=m'[`i',`=c'XB])))):/`=b'[1,cols(`=b')]
				}
		
			*Update mOS
				mata: mOS[.,`=OMC'] = `=m'[.,`=c'OC]
		}
