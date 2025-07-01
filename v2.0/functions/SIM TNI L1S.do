**********
	*SIM TNI L1S
**********

	*Patient matrix	
		
		
	*OLD STYLE OF CODE LOOPING OVER EACH PATIENT	
		
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
					*Chemo
						`=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,5] 
					*cons
						`=m'[`i',`=c'XB] = `=m'[`i',`=c'XB] + `=b'[1,6]
					
				*Calculate probability of survival to mTSD
					if (f`=b' == "weibull")  	`=m'[`i',`=c'PR] = exp(-(exp(`=m'[`i',`=c'XB]))*(mTSD[`i',`=OMC']:^exp(`=b'[1,cols(`=b')])))
							
				*Draw RN, conditional on survival to mTSD
					*`=m'[`i',`=c'RN] = runiform(1, 1, 0, `=m'[`i',`=c'PR])
					
				*Calculate survival time
					if (f`=b' == "weibull")  	`=m'[`i',`=c'OC] = ((ln(`=m'[`i',`=c'RN])):/-exp(`=m'[`i',`=c'XB])):^(1:/exp(`=b'[1,cols(`=b')]))
				
				*Update outcome matrix
					mOS[`i',`=OMC'] = `=m'[`i',`=c'OC]
				}
			}
		}

		
		
**********
	*Matrix multiplication style
**********

	*Patient matrix
		mata: p`=m' = mCore2 , mCons // 10 columns
		mata: rp`=m' = rmCore2
		mata: cp`=m' = cmCore2 \ cmCons
		mata: _matrix_list(p`=m', rp`=m', cp`=m')
		
	*Adjust coefficient matrix // 59 columns
		mata: aux = `=b'[.,59]
		mata: adj = `=b'[.,1..9] , `=b'[.,58] // 10 columns
		
	*Outcome matrix
		mata: o`=m' = mSU
		mata: ro`=m' = rmCore2
		mata: co`=m' = cmSU
		
		mata: o`=m'[.,1] = p`=m' * adj' // XB
		mata: o`=m'[.,2] = runiform(`=Obs',1) // RN
	
		*Loop to determine outcome
			forvalues i = 1/`=Obs' {
				mata {
					if (mState[`i',2] <= `=OMC') { // State filter
						
					*Calculate survival time
						if (f`=b' == "ereg") 		o`=m'[`i',3] = (ln(o`=m'[`i',2])):/-exp(o`=m'[`i',1])
						if (f`=b' == "weibull")  	o`=m'[`i',3] = ((ln(o`=m'[`i',2])):/-exp(o`=m'[`i',1])):^(1:/exp(aux))
						if (f`=b' == "gompertz")	o`=m'[`i',3] = (ln(1:-((aux:*(ln(o`=m'[`i',2]))):/exp(o`=m'[`i',1])))):/aux
				
					*Update mOS
						mOS[`i',`=OMC'] = o`=m'[`i',3]
					}
				}	
			}
		
		mata: _matrix_list(o`=m', ro`=m', co`=m')

		