**********
	*SIM OS DN // OMC = 2, DN = 2, L1S = 3, L1E = 4
**********

	*REPLACE CORE2 WITH CORE ONCE UPDATE COMPLETE

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

