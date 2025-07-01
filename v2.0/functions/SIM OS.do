**********
	*SIM OS
**********

	*Patient matrix	
		mata: `=m' = mCore , mLO
		mata: r`=m' = rmCore
		mata: c`=m' = cmCore \ cmLO
		mata: _matrix_list(`=m', r`=m', c`=m')
		
	*Determine outcome
		forvalues i = 1/`=Obs' {
			mata {
				if (mMOR[`i',`=OMC'-1] == 0 & mState[`i',2] <= `=OMC') { // Alive & State filters
				
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
					
				*Calculate probability of survival to mTSD
					if (f`=b' == "ereg") 		`=m'[`i',`=c'PR] = exp(-(exp(`=m'[`i',`=c'XB]))*mTSD[`i',`=OMC'])
					if (f`=b' == "weibull")  	`=m'[`i',`=c'PR] = exp(-(exp(`=m'[`i',`=c'XB]))*(mTSD[`i',`=OMC']:^exp(`=b'[1,cols(`=b')])))
					if (f`=b' == "gompertz")	`=m'[`i',`=c'PR] = exp(-(exp(`=m'[`i',`=c'XB]))*(1:/`=b'[1,cols(`=b')])*(exp(`=b'[1,cols(`=b')]*mTSD[`i',`=OMC']):-1))
							
				*Draw RN, conditional on survival to mTSD
					`=m'[`i',`=c'RN] = runiform(1, 1, 0, `=m'[`i',`=c'PR])
					
				*Calculate survival time
					if (f`=b' == "ereg") 		`=m'[`i',`=c'OC] = (ln(`=m'[`i',`=c'RN])):/-exp(`=m'[`i',`=c'XB])
					if (f`=b' == "weibull")  	`=m'[`i',`=c'OC] = ((ln(`=m'[`i',`=c'RN])):/-exp(`=m'[`i',`=c'XB])):^(1:/exp(`=b'[1,cols(`=b')]))
					if (f`=b' == "gompertz")	`=m'[`i',`=c'OC] = (ln(1:-((`=b'[1,cols(`=b')]:*(ln(`=m'[`i',`=c'RN]))):/exp(`=m'[`i',`=c'XB])))):/`=b'[1,cols(`=b')]
				
				*Update outcome matrix
					mOS[`i',`=OMC'] = `=m'[`i',`=c'OC]
				}
			}
		}
