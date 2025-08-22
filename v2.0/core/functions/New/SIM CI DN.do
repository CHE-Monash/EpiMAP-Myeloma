**********
	*SIM CI DN // OMC = 2, DN = 2, L1S = 3, L1E = 4
**********

	*REPLACE CORE2 WITH CORE ONCE UPDATE COMPLETE

	*Patient matrix
		mata: p`=m' = mCore2 , mSCT[.,1] , mCons // 11 columns
		mata: rp`=m' = rmCore2
		mata: cp`=m' = cmCore2 \ cmSCT[1.,] \ cmCons
		mata: _matrix_list(p`=m', rp`=m', cp`=m')
		
	*Adjust coefficient matrix
		mata: aux = `=b'[.,12]
		mata: adj = `=b'[.,1..11]
		
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
						if (f`=b' == "exponential")	o`=m'[`i',3] = (ln(o`=m'[`i',2])):/-exp(o`=m'[`i',1])
						if (f`=b' == "weibull") 	o`=m'[`i',3] = ((ln(o`=m'[`i',2])):/-exp(o`=m'[`i',1])):^(1:/exp(aux))
						if (f`=b' == "gompertz")	o`=m'[`i',3] = (ln(1:-((aux:*(ln(o`=m'[`i',2]))):/exp(o`=m'[`i',1])))):/aux
				
					*Update outcome matrices
						mTNE[`i',`=OMC'] = o`=m'[`i',3]/365.25
						mTSD[`i',`=OMC'+1] = mTNE[`i',`=OMC']
					}
				
				*Grab Prevalent patient data
					if (mState[`i',2] > `=OMC') o`=m'[`i',3] = mTNE[`i',`=OMC']*365.25
				}
			}

		mata: _matrix_list(o`=m', ro`=m', co`=m')
